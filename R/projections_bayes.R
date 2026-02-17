#' Project DMFT trends using Bayesian posterior draws
#'
#' Experimental Bayesian alternative to [dmft_project()]. Extracts year
#' random effect posterior draws, estimates a temporal trend, and generates
#' scenario projections with uncertainty propagated from the posterior.
#'
#' @param fit_bayes A fitted Bayesian model from [dmft_fit_bayes()].
#' @param estimates Predictions from [dmft_predict_bayes()].
#' @param config A [dmft_config] object with projection range set.
#' @param damping_factor Trend damping per year (default 0.95).
#'
#' @returns A tibble with columns: `year`, `scenario`, `mean_proj`,
#'   `lower`, `upper`.
#' @export
dmft_project_bayes <- function(fit_bayes, estimates, config,
                                damping_factor = 0.95) {
  check_suggests("cmdstanr")
  check_suggests("posterior")
  stopifnot(inherits(fit_bayes, "dmft_fit_bayes"))
  stopifnot(inherits(config, "dmft_config"))

  if (is.null(config$projection_start)) {
    cli_abort("Projection range not set in config. Use `projection_range` in dmft_config().")
  }

  dentition <- fit_bayes$dentition
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent
  n_forecast <- config$projection_end - config$projection_start + 1L

  # Extract posterior draws of year random effects
  stan_fit <- fit_bayes$fit
  draws <- stan_fit$draws(format = "draws_df")

  re_year_draws <- as.matrix(draws[, grep("^re_year\\[", names(draws))])
  alpha_draws   <- draws[["alpha"]]
  sigma_draws   <- draws[["sigma"]]
  n_draws       <- length(alpha_draws)

  # Posterior mean of year effects for trend estimation
  re_year_means <- colMeans(re_year_draws)
  n_years <- length(re_year_means)

  # Use last 10 years for trend estimation
  n_recent <- min(10, n_years)
  recent_idx <- (n_years - n_recent + 1):n_years
  recent_effects <- re_year_means[recent_idx]

  lm_fit <- stats::lm(recent_effects ~ seq_along(recent_effects))
  annual_trend <- stats::coef(lm_fit)[2]
  trend_se <- summary(lm_fit)$coefficients[2, 2]

  # Last observed value (posterior mean)
  last_intercept <- mean(alpha_draws)
  last_national <- last_intercept + re_year_means[n_years]
  sigma_resid <- mean(sigma_draws)

  trend <- list(
    last_value   = last_national,
    annual_trend = annual_trend,
    trend_se     = trend_se,
    sigma_resid  = sigma_resid
  )

  # Generate scenarios using existing project_linear helper
  scenarios <- config$scenarios
  proj_list <- list()

  for (sname in names(scenarios)) {
    adj <- scenarios[[sname]]$adjustment
    proj <- project_linear(
      trend, n_forecast, adj, damping_factor, dmft_max, config$seed
    )

    proj_years <- seq.int(config$projection_start, config$projection_end)
    proj_list[[sname]] <- dplyr::tibble(
      year      = proj_years,
      scenario  = scenarios[[sname]]$name,
      mean_proj = proj$mean,
      lower     = proj$lower,
      upper     = proj$upper
    )
  }

  result <- dplyr::bind_rows(proj_list)
  cli_alert_success(
    "Bayesian projections: {n_forecast} years across {length(scenarios)} scenarios"
  )
  result
}
