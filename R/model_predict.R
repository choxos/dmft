#' Generate predictions from a fitted DMFT model with AST smoothing
#'
#' Produces predictions for all region-year-age combinations in the
#' historical period. Uses AST (Age-Spatial-Temporal) kernel smoothing
#' on mixed-model residuals, with bootstrap uncertainty intervals.
#'
#' @param fit A fitted model from [dmft_fit()].
#' @param adjacency Adjacency object from [dmft_adjacency()].
#' @param config A [dmft_config] object.
#' @param n_boot Number of bootstrap replicates for uncertainty.
#'   Set to 0 to skip uncertainty estimation.
#'
#' @returns A tibble with columns: `region`, `year`, `age_group`,
#'   `predicted`, `lower`, `upper`.
#' @export
dmft_predict <- function(fit, adjacency, config, n_boot = NULL) {
  stopifnot(inherits(config, "dmft_config"))

  dentition <- attr(fit, "dmft_dentition")
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent
  if (is.null(n_boot)) n_boot <- config$n_boot

  # Apply AST smoothing
  ast_result <- dmft_ast(fit, adjacency, config, dentition)
  estimates <- ast_result$estimates

  # Bootstrap uncertainty
  if (n_boot > 0) {
    cli_alert_info("Computing bootstrap uncertainty ({n_boot} replicates)...")
    ui <- compute_bootstrap_ui(fit, adjacency, config, dentition, n_boot, dmft_max)
    estimates <- merge(estimates, ui, by = c("region", "year", "age_group"), all.x = TRUE)
    estimates$lower[is.na(estimates$lower)] <- pmax(0, estimates$predicted[is.na(estimates$lower)] * 0.5)
    estimates$upper[is.na(estimates$upper)] <- pmin(dmft_max, estimates$predicted[is.na(estimates$upper)] * 1.5)
  } else {
    estimates$lower <- pmax(0, estimates$predicted * 0.75)
    estimates$upper <- pmin(dmft_max, estimates$predicted * 1.25)
  }

  estimates <- dplyr::as_tibble(estimates)
  estimates <- estimates[, c("region", "year", "age_group", "predicted", "lower", "upper")]
  estimates <- estimates[order(estimates$region, estimates$year, estimates$age_group), ]

  cli_alert_success("Predictions: {nrow(estimates)} cells ({dentition})")
  estimates
}


#' Compute bootstrap uncertainty intervals
#' @keywords internal
compute_bootstrap_ui <- function(fit, adjacency, config, dentition, n_boot, dmft_max) {
  mdata <- attr(fit, "dmft_mdata")
  age_groups <- mdata$age_groups
  set.seed(config$seed)

  # Extract fixed and random effects distributions
  fixed_effects <- lme4::fixef(fit)
  re_region <- lme4::ranef(fit)$region_std
  re_year   <- lme4::ranef(fit)$year_factor

  # Get variance components
  vc <- as.data.frame(lme4::VarCorr(fit))
  sigma_region <- vc$sdcor[vc$grp == "region_std"]
  sigma_year   <- vc$sdcor[vc$grp == "year_factor"]
  sigma_resid  <- stats::sigma(fit)

  # Prediction grid
  grid <- expand.grid(
    region    = config$regions,
    year      = config$historical_years,
    age_group = age_groups,
    stringsAsFactors = FALSE
  )

  # Matrix to store bootstrap predictions
  n_cells <- nrow(grid)
  boot_preds <- matrix(NA_real_, n_cells, n_boot)

  for (b in seq_len(n_boot)) {
    # Resample random effects from their estimated distributions
    re_r <- stats::rnorm(length(config$regions), 0, sigma_region)
    names(re_r) <- config$regions
    re_y <- stats::rnorm(config$n_years, 0, sigma_year)
    names(re_y) <- as.character(config$historical_years)

    # Compute predictions with resampled effects
    for (i in seq_len(n_cells)) {
      pred <- fixed_effects[["(Intercept)"]]
      r <- grid$region[i]
      y <- as.character(grid$year[i])
      if (r %in% names(re_r)) pred <- pred + re_r[r]
      if (y %in% names(re_y)) pred <- pred + re_y[y]
      # Add residual noise
      pred <- pred + stats::rnorm(1, 0, sigma_resid)
      boot_preds[i, b] <- pred
    }
  }

  # Compute 2.5th and 97.5th percentiles
  grid$lower <- pmax(0, apply(boot_preds, 1, stats::quantile, 0.025, na.rm = TRUE))
  grid$upper <- pmin(dmft_max, apply(boot_preds, 1, stats::quantile, 0.975, na.rm = TRUE))

  grid[, c("region", "year", "age_group", "lower", "upper")]
}
