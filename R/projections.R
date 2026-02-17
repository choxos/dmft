#' Project DMFT trends into the future
#'
#' Extracts the temporal trend from a fitted model and generates three
#' scenarios (reference, optimistic, pessimistic) using linear trend
#' extrapolation with optional damping and bounded uncertainty.
#'
#' @param fit A fitted model from [dmft_fit()].
#' @param estimates Predictions from [dmft_predict()].
#' @param config A [dmft_config] object. Must have `projection_start` and
#'   `projection_end` set.
#' @param damping_factor Trend damping per year (default 0.95).
#'
#' @returns A tibble with columns: `year`, `scenario`, `mean_proj`,
#'   `lower`, `upper`.
#' @export
dmft_project <- function(fit, estimates, config, damping_factor = 0.95) {
  stopifnot(inherits(config, "dmft_config"))
  if (is.null(config$projection_start)) {
    cli_abort("Projection range not set in config. Use `projection_range` in dmft_config().")
  }

  dentition <- attr(fit, "dmft_dentition")
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent
  n_forecast <- config$projection_end - config$projection_start + 1L

  # Extract temporal trend from year random effects
  trend <- extract_year_trend(fit, config)

  # Generate scenarios
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
  cli_alert_success("Projected {n_forecast} years across {length(scenarios)} scenarios")
  result
}


# -- Internal helpers ----------------------------------------------------------

#' Extract temporal trend from year random effects
#' @keywords internal
extract_year_trend <- function(fit, config) {
  re_year <- lme4::ranef(fit)$year_factor
  year_effects <- re_year[, 1]
  year_labels <- as.integer(rownames(re_year))

  # Use last 10 years for trend estimation
  n_years <- length(year_effects)
  n_recent <- min(10, n_years)
  recent_idx <- (n_years - n_recent + 1):n_years
  recent_effects <- year_effects[recent_idx]

  # Linear trend from recent period
  lm_fit <- stats::lm(recent_effects ~ seq_along(recent_effects))
  annual_trend <- stats::coef(lm_fit)[2]
  trend_se <- summary(lm_fit)$coefficients[2, 2]

  # Last observed values
  intercept <- lme4::fixef(fit)[["(Intercept)"]]
  last_national <- intercept + year_effects[n_years]

  # Residual SD for uncertainty
  sigma_resid <- stats::sigma(fit)

  list(
    last_value   = last_national,
    annual_trend = annual_trend,
    trend_se     = trend_se,
    sigma_resid  = sigma_resid
  )
}


#' Project using linear extrapolation with damping
#' @keywords internal
project_linear <- function(trend, n_forecast, adjustment,
                            damping_factor, dmft_max, seed) {
  set.seed(seed)

  mean_proj  <- numeric(n_forecast)
  lower_proj <- numeric(n_forecast)
  upper_proj <- numeric(n_forecast)

  for (t in seq_len(n_forecast)) {
    # Damped linear extrapolation
    damped_trend <- trend$annual_trend * damping_factor^t
    expected <- trend$last_value + t * damped_trend + t * adjustment

    # Growing uncertainty with time
    uncertainty_sd <- sqrt(trend$sigma_resid^2 + (t * trend$trend_se)^2)

    mean_proj[t]  <- pmax(0, pmin(dmft_max, expected))
    lower_proj[t] <- pmax(0, expected - qnorm(0.975) * uncertainty_sd)
    upper_proj[t] <- pmin(dmft_max, expected + qnorm(0.975) * uncertainty_sd)
  }

  list(mean = mean_proj, lower = lower_proj, upper = upper_proj)
}
