#' Run model diagnostics
#'
#' Computes fit statistics (AIC, BIC), residual diagnostics (RMSE, MAE),
#' spatial autocorrelation of residuals (Moran's I), and variance components.
#'
#' @param fit A fitted model from [dmft_fit()].
#' @param config A [dmft_config] object.
#'
#' @returns A list with elements: `fit_stats`, `residuals`, `spatial`,
#'   `variance_components`, `validity`.
#' @export
dmft_diagnose <- function(fit, config) {
  stopifnot(inherits(config, "dmft_config"))

  mdata     <- attr(fit, "dmft_mdata")
  dentition <- attr(fit, "dmft_dentition")
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent

  results <- list()

  # -- Fit statistics --
  results$fit_stats <- data.frame(
    metric = c("AIC", "BIC", "logLik", "deviance"),
    value  = c(stats::AIC(fit), stats::BIC(fit),
               as.numeric(stats::logLik(fit)),
               stats::deviance(fit))
  )

  # -- Variance components --
  vc <- as.data.frame(lme4::VarCorr(fit))
  results$variance_components <- data.frame(
    group    = vc$grp,
    variance = vc$vcov,
    sd       = vc$sdcor
  )

  # -- Residuals --
  d <- mdata$data
  fitted_vals <- stats::fitted(fit)
  observed    <- d$y
  resid <- observed - fitted_vals
  std_resid <- resid / stats::sd(resid, na.rm = TRUE)

  results$residuals <- list(
    rmse = sqrt(mean(resid^2, na.rm = TRUE)),
    mae  = mean(abs(resid), na.rm = TRUE),
    bias = mean(resid, na.rm = TRUE),
    normality_p = tryCatch(
      stats::shapiro.test(utils::head(std_resid[!is.na(std_resid)], 5000))$p.value,
      error = function(e) NA
    )
  )

  # -- Spatial autocorrelation --
  results$spatial <- tryCatch({
    d$residual <- resid
    prov_resid <- dplyr::summarize(
      dplyr::group_by(d, province_idx),
      mean_resid = mean(.data$residual, na.rm = TRUE),
      .groups = "drop"
    )
    prov_resid <- dplyr::arrange(prov_resid, .data$province_idx)

    nb  <- mdata$adjacency$nb
    lw  <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
    if (nrow(prov_resid) == length(lw$neighbours)) {
      mt <- spdep::moran.test(prov_resid$mean_resid, lw,
                               alternative = "two.sided")
      list(statistic = as.numeric(mt$statistic), p.value = mt$p.value)
    } else {
      list(statistic = NA, p.value = NA)
    }
  }, error = function(e) list(statistic = NA, p.value = NA))

  # -- Prediction validity --
  results$validity <- list(
    n_negative  = sum(fitted_vals < 0, na.rm = TRUE),
    n_above_max = sum(fitted_vals > dmft_max, na.rm = TRUE),
    range       = range(fitted_vals, na.rm = TRUE)
  )

  # -- Print summary --
  cli_alert_info("Diagnostics ({toupper(dentition)})")
  cat(sprintf("  AIC=%.1f  BIC=%.1f\n",
              results$fit_stats$value[1], results$fit_stats$value[2]))
  cat(sprintf("  RMSE=%.3f  MAE=%.3f  Bias=%.4f\n",
              results$residuals$rmse, results$residuals$mae,
              results$residuals$bias))
  for (i in seq_len(nrow(results$variance_components))) {
    cat(sprintf("  Var(%s): %.4f (SD=%.4f)\n",
                results$variance_components$group[i],
                results$variance_components$variance[i],
                results$variance_components$sd[i]))
  }
  if (!is.na(results$spatial$p.value)) {
    cat(sprintf("  Moran's I: stat=%.4f, p=%.4f %s\n",
                results$spatial$statistic, results$spatial$p.value,
                ifelse(results$spatial$p.value < 0.05,
                       "(spatial autocorrelation detected)", "(OK)")))
  }
  cat(sprintf("  Predictions: range [%.2f, %.2f], %d<0, %d>%d\n",
              results$validity$range[1], results$validity$range[2],
              results$validity$n_negative, results$validity$n_above_max,
              dmft_max))

  results
}
