#' Generate predictions from a fitted DMFT model
#'
#' Produces predictions for all region-year-age combinations in the
#' historical period. Uses posterior sampling when available, falling
#' back to a delta-method approximation.
#'
#' @param fit A fitted model from [dmft_fit()].
#' @param config A [dmft_config] object.
#' @param n_posterior Number of posterior samples for uncertainty.
#'   Set to 0 to use the delta method.
#'
#' @returns A tibble with columns: `region`, `year`, `age_group`,
#'   `predicted`, `lower`, `upper`.
#' @export
dmft_predict <- function(fit, config, n_posterior = 100L) {
  stopifnot(inherits(config, "dmft_config"))

  mdata    <- attr(fit, "dmft_mdata")
  dentition <- attr(fit, "dmft_dentition")
  backend  <- attr(fit, "dmft_backend") %||% "inla"
  age_groups <- mdata$age_groups
  dmft_max <- if (dentition == "deciduous") config$dmft_max_deciduous
              else config$dmft_max_permanent

  # Prediction grid
  grid <- expand.grid(
    region    = config$regions,
    year      = config$historical_years,
    age_group = age_groups,
    stringsAsFactors = FALSE
  )
  grid$province_idx   <- match(grid$region, config$regions)
  grid$year_idx       <- grid$year - config$year_start + 1L
  grid$age_group_idx  <- match(grid$age_group, age_groups)
  grid <- dplyr::as_tibble(grid)

  if (backend == "inla" && n_posterior > 0) {
    preds <- predict_posterior_inla(fit, grid, n_posterior, dmft_max)
  } else if (backend == "inla") {
    preds <- predict_delta_inla(fit, grid, dmft_max)
  } else {
    preds <- predict_stan(fit, grid, dmft_max)
  }

  preds
}


# -- Internal prediction functions ---------------------------------------------

#' @keywords internal
predict_posterior_inla <- function(fit, grid, n_posterior, dmft_max) {
  samples <- tryCatch(
    INLA::inla.posterior.sample(n_posterior, fit),
    error = function(e) NULL
  )
  if (is.null(samples)) {
    cli_alert_warning("Posterior sampling failed; using delta method")
    return(predict_delta_inla(fit, grid, dmft_max))
  }

  n_pred <- nrow(grid)
  pred_mat <- matrix(NA_real_, n_pred, n_posterior)

  for (s in seq_len(n_posterior)) {
    lat <- samples[[s]]$latent
    rn  <- rownames(lat)
    intercept_s <- lat[grep("^\\(Intercept\\)", rn), 1]
    prov_s <- lat[grep("^province_idx:", rn), 1]
    year_s <- lat[grep("^year_idx:", rn), 1]
    age_s  <- lat[grep("^age_group_idx:", rn), 1]

    for (i in seq_len(n_pred)) {
      eta <- intercept_s
      pi <- grid$province_idx[i]; yi <- grid$year_idx[i]; ai <- grid$age_group_idx[i]
      if (!is.na(pi) && pi <= length(prov_s)) eta <- eta + prov_s[pi]
      if (!is.na(yi) && yi <= length(year_s)) eta <- eta + year_s[yi]
      if (!is.na(ai) && ai <= length(age_s))  eta <- eta + age_s[ai]
      pred_mat[i, s] <- eta
    }
  }

  grid$predicted <- pmax(0, rowMeans(pred_mat, na.rm = TRUE))
  grid$lower     <- pmax(0, apply(pred_mat, 1, quantile, 0.025, na.rm = TRUE))
  grid$upper     <- pmin(dmft_max, apply(pred_mat, 1, quantile, 0.975, na.rm = TRUE))
  grid
}


#' @keywords internal
predict_delta_inla <- function(fit, grid, dmft_max) {
  spatial  <- fit$summary.random$province_idx
  temporal <- fit$summary.random$year_idx
  age_eff  <- fit$summary.random$age_group_idx
  intercept    <- fit$summary.fixed["(Intercept)", "mean"]
  intercept_sd <- fit$summary.fixed["(Intercept)", "sd"]

  eta <- rep(intercept, nrow(grid))
  eta_var <- rep(intercept_sd^2, nrow(grid))

  for (i in seq_len(nrow(grid))) {
    pi <- grid$province_idx[i]; yi <- grid$year_idx[i]; ai <- grid$age_group_idx[i]
    if (!is.na(pi) && pi <= nrow(spatial)) {
      eta[i]     <- eta[i] + spatial$mean[pi]
      eta_var[i] <- eta_var[i] + spatial$sd[pi]^2
    }
    if (!is.na(yi) && yi <= nrow(temporal)) {
      eta[i]     <- eta[i] + temporal$mean[yi]
      eta_var[i] <- eta_var[i] + temporal$sd[yi]^2
    }
    if (!is.na(ai) && ai <= nrow(age_eff)) {
      eta[i]     <- eta[i] + age_eff$mean[ai]
      eta_var[i] <- eta_var[i] + age_eff$sd[ai]^2
    }
  }

  eta_sd <- sqrt(eta_var)
  grid$predicted <- pmax(0, eta)
  grid$lower     <- pmax(0, eta - 1.96 * eta_sd)
  grid$upper     <- pmin(dmft_max, eta + 1.96 * eta_sd)
  grid
}


#' @keywords internal
predict_stan <- function(fit, grid, dmft_max) {
  draws <- fit$stan_fit$draws(format = "matrix")
  alpha   <- draws[, "alpha"]
  prov_cols <- grep("^province_effect\\[", colnames(draws))
  year_cols <- grep("^year_effect\\[", colnames(draws))
  age_cols  <- grep("^age_effect\\[", colnames(draws))
  prov_d <- draws[, prov_cols, drop = FALSE]
  year_d <- draws[, year_cols, drop = FALSE]
  age_d  <- draws[, age_cols, drop = FALSE]

  n_draws <- length(alpha)
  n_pred  <- nrow(grid)
  pred_mat <- matrix(NA_real_, n_pred, n_draws)

  for (i in seq_len(n_pred)) {
    eta <- alpha
    pi <- grid$province_idx[i]; yi <- grid$year_idx[i]; ai <- grid$age_group_idx[i]
    if (!is.na(pi) && pi <= ncol(prov_d)) eta <- eta + prov_d[, pi]
    if (!is.na(yi) && yi <= ncol(year_d)) eta <- eta + year_d[, yi]
    if (!is.na(ai) && ai <= ncol(age_d))  eta <- eta + age_d[, ai]
    pred_mat[i, ] <- eta
  }

  grid$predicted <- pmax(0, rowMeans(pred_mat))
  grid$lower     <- pmax(0, apply(pred_mat, 1, quantile, 0.025))
  grid$upper     <- pmin(dmft_max, apply(pred_mat, 1, quantile, 0.975))
  grid
}
