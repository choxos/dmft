#' Project DMFT trends into the future
#'
#' Extracts the temporal trend from a fitted model and generates three
#' scenarios (reference, optimistic, pessimistic) using damped RW2
#' extrapolation with bounded uncertainty growth.
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

  backend   <- attr(fit, "dmft_backend") %||% "inla"
  dentition <- attr(fit, "dmft_dentition")
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent
  n_forecast <- config$projection_end - config$projection_start + 1L

  # Extract temporal trend
  trend <- extract_trend(fit, backend)

  # Generate scenarios
  scenarios <- config$scenarios
  proj_list <- list()

  for (sname in names(scenarios)) {
    adj <- scenarios[[sname]]$adjustment
    proj_mat <- project_scenario(
      trend, n_forecast, config$n_draws, adj, damping_factor, dmft_max,
      config$seed + match(sname, names(scenarios)) - 1L
    )

    proj_years <- seq.int(config$projection_start, config$projection_end)
    proj_list[[sname]] <- dplyr::tibble(
      year      = proj_years,
      scenario  = scenarios[[sname]]$name,
      mean_proj = rowMeans(proj_mat),
      lower     = apply(proj_mat, 1, quantile, 0.025),
      upper     = apply(proj_mat, 1, quantile, 0.975)
    )
  }

  result <- dplyr::bind_rows(proj_list)
  cli_alert_success("Projected {n_forecast} years across {length(scenarios)} scenarios")
  result
}


# -- Internal helpers ----------------------------------------------------------

#' @keywords internal
extract_trend <- function(fit, backend) {
  if (backend == "inla") {
    te <- fit$summary.random$year_idx
    n_years <- nrow(te)
    recent <- max(1, n_years - 9):n_years
    recent_eff <- te$mean[recent]
    lm_fit <- lm(recent_eff ~ seq_along(recent_eff))

    hyper_names <- rownames(fit$summary.hyperpar)
    yr_idx <- grep("year", hyper_names, ignore.case = TRUE)
    rw2_sd <- if (length(yr_idx) > 0) 1 / sqrt(fit$summary.hyperpar$mean[yr_idx[1]])
              else sd(diff(diff(te$mean)))

    list(
      last_effect = te$mean[n_years],
      last_sd     = te$sd[n_years],
      prev_effect = te$mean[n_years - 1],
      prev_sd     = te$sd[n_years - 1],
      rw2_sd      = rw2_sd
    )
  } else {
    # Stan: use year_effect draws
    draws <- fit$stan_fit$draws(format = "matrix")
    ye_cols <- grep("^year_effect\\[", colnames(draws))
    ye <- draws[, ye_cols, drop = FALSE]
    n_years <- ncol(ye)

    list(
      last_effect = mean(ye[, n_years]),
      last_sd     = sd(ye[, n_years]),
      prev_effect = mean(ye[, n_years - 1]),
      prev_sd     = sd(ye[, n_years - 1]),
      rw2_sd      = sd(diff(diff(colMeans(ye))))
    )
  }
}


#' @keywords internal
project_scenario <- function(trend, n_forecast, n_samples, adjustment,
                              damping_factor, dmft_max, seed) {
  set.seed(seed)
  proj <- matrix(NA_real_, n_forecast, n_samples)

  for (s in seq_len(n_samples)) {
    current <- trend$last_effect + rnorm(1, 0, trend$last_sd)
    prev    <- trend$prev_effect + rnorm(1, 0, trend$prev_sd)

    for (t in seq_len(n_forecast)) {
      tr <- current - prev
      expected <- current + damping_factor^t * tr + adjustment
      noise_sd <- trend$rw2_sd * sqrt(t) / sqrt(1 + t / n_forecast)
      projected <- expected + rnorm(1, 0, noise_sd)
      projected <- max(0, min(dmft_max, projected))
      proj[t, s] <- projected
      prev <- current
      current <- projected
    }
  }

  proj
}
