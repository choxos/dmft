#' Fit a Bayesian hierarchical model for DMFT/dmft
#'
#' Implements a GBD-style spatial-temporal-age model with BYM2 spatial
#' random effects, RW2 temporal trends, RW1 age effects, and IID sex
#' effects. Uses Gaussian meta-regression (identity link) by default.
#'
#' @param data A data frame (one of the dentition subsets from [dmft_clean()]).
#' @param adjacency Adjacency object from [dmft_adjacency()].
#' @param dentition `"deciduous"` or `"permanent"`.
#' @param config A [dmft_config] object.
#' @param family Likelihood family: `"gaussian"` (default, meta-regression)
#'   or `"nbinomial"`.
#' @param include_interactions Logical; include space-time, space-age, and
#'   time-age IID interaction random effects.
#'
#' @returns A fitted model object (INLA or wrapped cmdstanr fit).
#' @export
dmft_fit <- function(data,
                      adjacency,
                      dentition = c("deciduous", "permanent"),
                      config,
                      family = "gaussian",
                      include_interactions = TRUE) {

  dentition <- match.arg(dentition)
  stopifnot(inherits(config, "dmft_config"))

  # Prepare model data
  mdata <- prepare_model_data(data, adjacency, dentition, config)

  if (config$backend == "inla") {
    check_suggests("INLA", "for the INLA backend")
    fit_inla(mdata, config, family, include_interactions)
  } else {
    check_suggests("cmdstanr", "for the Stan backend")
    fit_stan(mdata, config, family, include_interactions, dentition)
  }
}


# -- data preparation ---------------------------------------------------------

#' Prepare model data
#' @keywords internal
prepare_model_data <- function(data, adjacency, dentition, config) {
  mean_col <- if (dentition == "deciduous") "mean_dmft" else "mean_DMFT"
  se_col   <- "se_imputed"
  age_groups <- if (dentition == "deciduous") config$age_groups_deciduous
                else config$age_groups_permanent

  # Filter valid records
  data <- data[!is.na(data[[mean_col]]) & !is.na(data$region_std) &
                 data$region_std %in% config$regions &
                 !is.na(data$year) &
                 data$year >= config$year_start &
                 data$year <= config$year_end, ]

  if (nrow(data) == 0) cli_abort("No valid records after filtering.")

  # Indices
  data$province_idx   <- match(data$region_std, config$regions)
  data$year_idx       <- data$year - config$year_start + 1L
  data$age_group_idx  <- match(data$age_group, age_groups)
  data$sex_idx        <- ifelse(is.na(data$sex_std) | data$sex_std == "Both",
                                1L, match(data$sex_std, c("Male", "Female")))
  data$sex_idx[is.na(data$sex_idx)] <- 1L

  # Interaction indices
  data$st_idx <- as.integer(factor(paste(data$province_idx, data$year_idx, sep = "_")))
  data$sa_idx <- as.integer(factor(paste(data$province_idx, data$age_group_idx, sep = "_")))
  data$ta_idx <- as.integer(factor(paste(data$year_idx, data$age_group_idx, sep = "_")))

  # Response
  data$y <- data[[mean_col]]

  # Precision weights
  if (se_col %in% names(data)) {
    se_vals <- data[[se_col]]
    data$precision <- ifelse(!is.na(se_vals) & se_vals > 0, 1 / se_vals^2, 1)
  } else {
    data$precision <- 1
  }
  data$precision <- data$precision / max(data$precision, na.rm = TRUE)

  cli_alert_info("Model data: {nrow(data)} obs, {length(unique(data$province_idx))} regions, years {min(data$year)}-{max(data$year)}")

  list(
    data       = data,
    n_regions  = config$n_regions,
    n_years    = config$n_years,
    n_age      = length(age_groups),
    age_groups = age_groups,
    adjacency  = adjacency,
    dentition  = dentition,
    mean_col   = mean_col
  )
}


# -- INLA backend --------------------------------------------------------------

#' @keywords internal
fit_inla <- function(mdata, config, family, include_interactions) {
  g <- INLA::inla.read.graph(mdata$adjacency$graph_path)

  # Build formula
  parts <- c("y ~ 1")
  parts <- c(parts, sprintf(
    "f(province_idx, model='bym2', graph=g, scale.model=TRUE, constr=TRUE, hyper=list(phi=list(prior='%s', param=c(%s)), prec=list(prior='%s', param=c(%s))))",
    config$priors$bym2_phi$prior,
    paste(config$priors$bym2_phi$param, collapse = ","),
    config$priors$bym2_prec$prior,
    paste(config$priors$bym2_prec$param, collapse = ",")
  ))
  parts <- c(parts, sprintf(
    "f(year_idx, model='rw2', constr=TRUE, hyper=list(prec=list(prior='%s', param=c(%s))))",
    config$priors$rw2_prec$prior,
    paste(config$priors$rw2_prec$param, collapse = ",")
  ))
  parts <- c(parts, sprintf(
    "f(age_group_idx, model='rw1', constr=TRUE, hyper=list(prec=list(prior='%s', param=c(%s))))",
    config$priors$rw1_prec$prior,
    paste(config$priors$rw1_prec$param, collapse = ",")
  ))
  parts <- c(parts, sprintf(
    "f(sex_idx, model='iid', hyper=list(prec=list(prior='%s', param=c(%s))))",
    config$priors$iid_prec$prior,
    paste(config$priors$iid_prec$param, collapse = ",")
  ))

  if (include_interactions) {
    for (idx in c("st_idx", "sa_idx", "ta_idx")) {
      parts <- c(parts, sprintf(
        "f(%s, model='iid', hyper=list(prec=list(prior='%s', param=c(%s))))",
        idx, config$priors$iid_prec$prior,
        paste(config$priors$iid_prec$param, collapse = ",")
      ))
    }
  }

  formula <- as.formula(paste(parts, collapse = " + "))

  # Family setup
  if (family == "gaussian") {
    control_family <- list(
      hyper = list(prec = list(initial = 0, fixed = TRUE))
    )
    scale_vals <- mdata$data$precision
  } else {
    control_family <- list(
      hyper = list(theta = list(prior = "loggamma", param = c(1, 0.1)))
    )
    scale_vals <- NULL
  }

  cli_alert_info("Fitting INLA model ({mdata$dentition})...")
  t0 <- Sys.time()

  # Assign graph to calling environment so INLA formula can find it
  env <- environment()
  assign("g", g, envir = env)

  fit <- INLA::inla(
    formula,
    family  = family,
    data    = mdata$data,
    scale   = scale_vals,
    control.family    = control_family,
    control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
    control.predictor = list(compute = TRUE, link = 1),
    control.inla      = list(strategy = "adaptive", int.strategy = "auto"),
    num.threads       = config$n_cores,
    verbose           = FALSE
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cli_alert_success("INLA fit complete in {round(elapsed, 1)} min (DIC={round(fit$dic$dic, 1)}, WAIC={round(fit$waic$waic, 1)})")

  # Attach metadata
  attr(fit, "dmft_mdata")    <- mdata
  attr(fit, "dmft_config")   <- config
  attr(fit, "dmft_backend")  <- "inla"
  attr(fit, "dmft_dentition") <- mdata$dentition
  fit
}


# -- Stan backend ---------------------------------------------------------------

#' @keywords internal
fit_stan <- function(mdata, config, family, include_interactions, dentition) {
  stan_file <- system.file("stan", "dmft_model.stan", package = "dmft")
  if (stan_file == "") cli_abort("Stan model file not found in package.")

  model <- cmdstanr::cmdstan_model(stan_file)

  adj <- mdata$adjacency
  scaling_factor <- compute_bym2_scaling_factor(adj$node1, adj$node2, mdata$n_regions)

  d <- mdata$data
  se_col <- "se_imputed"
  se_vals <- if (se_col %in% names(d)) d[[se_col]] else rep(1, nrow(d))
  se_vals <- ifelse(is.na(se_vals) | se_vals <= 0, 1, se_vals)

  stan_data <- list(
    N = nrow(d), N_prov = mdata$n_regions, N_year = mdata$n_years,
    N_age = mdata$n_age, N_edges = length(adj$node1),
    use_negbin = as.integer(family == "nbinomial"),
    y_continuous = as.numeric(d[[mdata$mean_col]]),
    y_count = as.integer(pmax(0, round(d$y))),
    sample_size = pmax(1, ifelse(is.na(d$n), 100, d$n)),
    se = pmax(0.01, se_vals),
    prov_idx = as.integer(d$province_idx),
    year_idx = as.integer(d$year_idx),
    age_idx  = as.integer(pmin(mdata$n_age, pmax(1L, d$age_group_idx))),
    sex_idx  = as.integer(d$sex_idx),
    node1 = adj$node1, node2 = adj$node2,
    K = 0L, X = matrix(0, nrow(d), 0),
    scaling_factor = scaling_factor
  )
  stan_data$age_idx[is.na(stan_data$age_idx)] <- 1L

  cli_alert_info("Sampling with Stan ({dentition})...")
  ss <- config$stan_settings
  fit <- model$sample(
    data = stan_data, seed = config$seed,
    chains = ss$chains,
    parallel_chains = min(ss$parallel_chains, config$n_cores),
    iter_warmup = ss$iter_warmup, iter_sampling = ss$iter_sampling,
    adapt_delta = ss$adapt_delta, max_treedepth = ss$max_treedepth,
    refresh = 100
  )

  wrapped <- list(
    stan_fit = fit, backend = "cmdstanr",
    summary.fixed = data.frame(
      mean = fit$summary("alpha")$mean,
      sd   = fit$summary("alpha")$sd,
      row.names = "(Intercept)"
    ),
    summary.random = list(
      province_idx  = fit$summary("province_effect"),
      year_idx      = fit$summary("year_effect"),
      age_group_idx = fit$summary("age_effect")
    ),
    dic  = list(dic = NA, p.eff = NA),
    waic = list(waic = tryCatch(fit$loo()$estimates["looic", "Estimate"],
                                 error = function(e) NA))
  )
  class(wrapped) <- c("dmft_stan_fit", "list")
  attr(wrapped, "dmft_mdata")    <- mdata
  attr(wrapped, "dmft_config")   <- config
  attr(wrapped, "dmft_backend")  <- "cmdstanr"
  attr(wrapped, "dmft_dentition") <- dentition
  cli_alert_success("Stan fit complete")
  wrapped
}
