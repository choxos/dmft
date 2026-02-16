#' Run model diagnostics
#'
#' Computes fit statistics (DIC, WAIC), residual diagnostics (RMSE, MAE),
#' spatial autocorrelation of residuals (Moran's I), and prediction validity
#' checks.
#'
#' @param fit A fitted model from [dmft_fit()].
#' @param config A [dmft_config] object.
#'
#' @returns A list with elements: `fit_stats`, `residuals`, `spatial`,
#'   `validity`.
#' @export
dmft_diagnose <- function(fit, config) {
  stopifnot(inherits(config, "dmft_config"))

  mdata     <- attr(fit, "dmft_mdata")
  dentition <- attr(fit, "dmft_dentition")
  backend   <- attr(fit, "dmft_backend") %||% "inla"
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent

  results <- list()

  # -- Fit statistics --
  if (backend == "inla") {
    results$fit_stats <- data.frame(
      metric = c("DIC", "WAIC", "p_eff", "log_mlik"),
      value  = c(fit$dic$dic, fit$waic$waic, fit$dic$p.eff, fit$mlik[1, 1])
    )
  } else {
    loo_val <- tryCatch(fit$stan_fit$loo()$estimates["looic", "Estimate"],
                        error = function(e) NA)
    results$fit_stats <- data.frame(metric = "LOO-IC", value = loo_val)
  }

  # -- Residuals --
  if (backend == "inla" && !is.null(fit$summary.fitted.values)) {
    d <- mdata$data
    fitted_vals <- fit$summary.fitted.values$mean[seq_len(nrow(d))]
    observed    <- d$y
    resid <- observed - fitted_vals
    std_resid <- resid / sd(resid, na.rm = TRUE)

    results$residuals <- list(
      rmse = sqrt(mean(resid^2, na.rm = TRUE)),
      mae  = mean(abs(resid), na.rm = TRUE),
      bias = mean(resid, na.rm = TRUE),
      normality_p = tryCatch(shapiro.test(std_resid[!is.na(std_resid)])$p.value,
                             error = function(e) NA)
    )
  } else {
    results$residuals <- list(rmse = NA, mae = NA, bias = NA, normality_p = NA)
  }

  # -- Spatial autocorrelation --
  results$spatial <- tryCatch({
    if (backend != "inla" || is.null(fit$summary.fitted.values)) {
      list(statistic = NA, p.value = NA)
    } else {
      d <- mdata$data
      fitted_vals <- fit$summary.fitted.values$mean[seq_len(nrow(d))]
      d$residual <- d$y - fitted_vals
      prov_resid <- dplyr::summarize(
        dplyr::group_by(d, province_idx),
        mean_resid = mean(residual, na.rm = TRUE),
        .groups = "drop"
      )
      prov_resid <- dplyr::arrange(prov_resid, province_idx)

      nb  <- mdata$adjacency$nb
      lw  <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
      if (nrow(prov_resid) == length(lw$neighbours)) {
        mt <- spdep::moran.test(prov_resid$mean_resid, lw,
                                 alternative = "two.sided")
        list(statistic = as.numeric(mt$statistic), p.value = mt$p.value)
      } else {
        list(statistic = NA, p.value = NA)
      }
    }
  }, error = function(e) list(statistic = NA, p.value = NA))

  # -- Prediction validity --
  if (backend == "inla" && !is.null(fit$summary.fitted.values)) {
    fitted <- fit$summary.fitted.values$mean[seq_len(nrow(mdata$data))]
    results$validity <- list(
      n_negative  = sum(fitted < 0, na.rm = TRUE),
      n_above_max = sum(fitted > dmft_max, na.rm = TRUE),
      range       = range(fitted, na.rm = TRUE)
    )
  } else {
    results$validity <- list(n_negative = NA, n_above_max = NA, range = c(NA, NA))
  }

  # -- Print summary --
  cli_alert_info("Diagnostics ({toupper(dentition)})")
  if (!is.na(results$residuals$rmse)) {
    cat(sprintf("  RMSE=%.3f  MAE=%.3f  Bias=%.4f\n",
                results$residuals$rmse, results$residuals$mae,
                results$residuals$bias))
  }
  if (!is.na(results$spatial$p.value)) {
    cat(sprintf("  Moran's I: stat=%.4f, p=%.4f %s\n",
                results$spatial$statistic, results$spatial$p.value,
                ifelse(results$spatial$p.value < 0.05,
                       "(spatial autocorrelation detected)", "(OK)")))
  }
  if (!is.na(results$validity$n_negative)) {
    cat(sprintf("  Predictions: range [%.2f, %.2f], %d<0, %d>%d\n",
                results$validity$range[1], results$validity$range[2],
                results$validity$n_negative, results$validity$n_above_max,
                dmft_max))
  }

  results
}
