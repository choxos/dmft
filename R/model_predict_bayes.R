#' Generate predictions from a Bayesian DMFT model with AST smoothing
#'
#' Experimental Bayesian alternative to [dmft_predict()]. Extracts posterior
#' draws of fixed and random effects, computes predictions for the full
#' region-year-age grid, applies AST smoothing on posterior mean residuals,
#' and returns credible intervals from the posterior.
#'
#' @param fit_bayes A fitted Bayesian model from [dmft_fit_bayes()].
#' @param adjacency Adjacency object from [dmft_adjacency()].
#' @param config A [dmft_config] object.
#'
#' @returns A tibble with columns: `region`, `year`, `age_group`,
#'   `predicted`, `lower`, `upper`.
#' @export
dmft_predict_bayes <- function(fit_bayes, adjacency, config) {
  check_suggests("cmdstanr")
  check_suggests("posterior")
  stopifnot(inherits(fit_bayes, "dmft_fit_bayes"))

  dentition <- fit_bayes$dentition
  mdata     <- fit_bayes$mdata
  stan_fit  <- fit_bayes$fit
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent
  age_groups <- mdata$age_groups

  # Apply AST smoothing (now handles dmft_fit_bayes objects)
  ast_result <- dmft_ast(
    fit       = fit_bayes,
    adjacency = adjacency,
    config    = config,
    dentition = dentition
  )
  estimates <- ast_result$estimates

  # Compute posterior credible intervals
  cli_alert_info("Computing posterior credible intervals...")
  draws <- stan_fit$draws(format = "draws_df")
  alpha_draws     <- draws[["alpha"]]
  n_draws         <- length(alpha_draws)
  re_region_draws <- as.matrix(draws[, grep("^re_region\\[", names(draws))])
  re_year_draws   <- as.matrix(draws[, grep("^re_year\\[",   names(draws))])

  # Prediction grid from estimates
  grid <- estimates
  grid$province_idx <- match(grid$region, config$regions)
  grid$year_idx     <- grid$year - config$year_start + 1L
  n_cells <- nrow(grid)

  # Compute posterior predictions for CI
  pred_matrix <- matrix(NA_real_, n_cells, n_draws)
  for (d_idx in seq_len(n_draws)) {
    a <- alpha_draws[d_idx]
    re_r <- re_region_draws[d_idx, ]
    re_y <- re_year_draws[d_idx, ]
    for (i in seq_len(n_cells)) {
      r_idx <- grid$province_idx[i]
      y_idx <- grid$year_idx[i]
      re_r_val <- if (!is.na(r_idx) && r_idx <= ncol(re_region_draws)) re_r[r_idx] else 0
      re_y_val <- if (!is.na(y_idx) && y_idx <= ncol(re_year_draws)) re_y[y_idx] else 0
      pred_matrix[i, d_idx] <- a + re_r_val + re_y_val
    }
  }

  # Add AST residual to each posterior draw for proper uncertainty propagation
  ast_resid <- grid$residual_ast
  ast_resid[is.na(ast_resid)] <- 0
  pred_matrix <- pred_matrix + ast_resid

  # Posterior summaries
  grid$lower <- pmax(0, apply(pred_matrix, 1, stats::quantile, 0.025))
  grid$upper <- pmin(dmft_max, apply(pred_matrix, 1, stats::quantile, 0.975))

  # Clip predicted
  grid$predicted <- pmax(0, pmin(dmft_max, grid$predicted))

  # Clean up
  grid <- dplyr::as_tibble(grid[, c("region", "year", "age_group",
                                      "predicted", "lower", "upper")])
  grid <- grid[order(grid$region, grid$year, grid$age_group), ]

  cli_alert_success("Bayesian predictions: {nrow(grid)} cells ({dentition})")
  grid
}
