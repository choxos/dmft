#' Fit a Bayesian random intercept model for DMFT/dmft via cmdstanr
#'
#' Experimental Bayesian alternative to [dmft_fit()]. Uses the same
#' random intercept structure (`y ~ 1 + covariates + (1|region) + (1|year)`)
#' but estimates parameters via Hamiltonian Monte Carlo (Stan) rather than
#' REML. AST smoothing (Stage 2) is applied downstream via [dmft_predict_bayes()].
#'
#' @param data A data frame (one of the dentition subsets from [dmft_clean()]).
#' @param adjacency Adjacency object from [dmft_adjacency()].
#' @param dentition `"deciduous"` or `"permanent"`.
#' @param config A [dmft_config] object.
#' @param stan_settings Named list overriding default Stan sampling settings:
#'   `chains`, `iter_warmup`, `iter_sampling`, `adapt_delta`, `max_treedepth`.
#'   If `NULL`, uses defaults from config or built-in defaults.
#'
#' @returns A list of class `"dmft_fit_bayes"` with elements `fit` (CmdStanMCMC),
#'   `mdata` (prepared model data), `config`, and `dentition`.
#' @export
dmft_fit_bayes <- function(data,
                            adjacency,
                            dentition = c("deciduous", "permanent"),
                            config,
                            stan_settings = NULL) {

  check_suggests("cmdstanr")
  dentition <- match.arg(dentition)
  stopifnot(inherits(config, "dmft_config"))

  # Reuse the same data preparation as frequentist

  mdata <- prepare_lmer_data(data, adjacency, dentition, config)

  # Build Stan data list
  stan_data <- prepare_stan_data(mdata, config)

  # Merge settings: explicit > config > defaults
  defaults <- list(
    chains = 4L,
    iter_warmup = 1000L,
    iter_sampling = 1000L,
    adapt_delta = 0.95,
    max_treedepth = 12L
  )
  cfg_settings <- config$stan_settings %||% list()
  settings <- utils::modifyList(defaults, cfg_settings)
  if (!is.null(stan_settings)) {
    settings <- utils::modifyList(settings, stan_settings)
  }

  # Compile model
  stan_file <- system.file("stan", "dmft_ast.stan", package = "dmft")
  if (!nzchar(stan_file)) {
    cli_abort("Stan model file not found. Ensure dmft is installed with inst/stan/.")
  }

  cli_alert_info("Compiling Stan model...")
  model <- cmdstanr::cmdstan_model(stan_file)

  cli_alert_info(
    "Sampling ({dentition}): {settings$chains} chains, {settings$iter_sampling} draws each"
  )
  t0 <- Sys.time()

  fit <- model$sample(
    data            = stan_data,
    chains          = settings$chains,
    parallel_chains = min(settings$chains, parallel::detectCores(logical = FALSE)),
    iter_warmup     = settings$iter_warmup,
    iter_sampling   = settings$iter_sampling,
    adapt_delta     = settings$adapt_delta,
    max_treedepth   = settings$max_treedepth,
    seed            = config$seed,
    refresh         = max(1L, settings$iter_sampling %/% 5L)
  )

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cli_alert_success("Stan fit complete in {round(elapsed, 1)}s ({dentition})")

  result <- list(
    fit       = fit,
    mdata     = mdata,
    config    = config,
    dentition = dentition
  )
  class(result) <- "dmft_fit_bayes"
  result
}


# -- helpers -------------------------------------------------------------------

#' Prepare Stan data list from model data
#' @keywords internal
prepare_stan_data <- function(mdata, config) {
  d <- mdata$data

  # Identify covariates present in the data
  covariates <- config$covariates
  K <- 0L
  X <- matrix(0, nrow = nrow(d), ncol = 0)

  if (!is.null(covariates)) {
    available <- intersect(covariates, names(d))
    if (length(available) > 0) {
      X <- as.matrix(d[, available, drop = FALSE])
      K <- ncol(X)
      cli_alert_info("Bayesian model includes {K} covariate(s): {paste(available, collapse=', ')}")
    }
  }

  # SE: use se_imputed if available, otherwise use a default
  se_vals <- if ("se_imputed" %in% names(d)) d$se_imputed else rep(1, nrow(d))
  se_vals[is.na(se_vals) | se_vals <= 0] <- stats::median(se_vals[se_vals > 0], na.rm = TRUE)

  list(
    N          = nrow(d),
    N_regions  = mdata$n_regions,
    N_years    = mdata$n_years,
    K          = K,
    y          = d$y,
    se         = se_vals,
    region_idx = d$province_idx,
    year_idx   = d$year_idx,
    X          = X
  )
}
