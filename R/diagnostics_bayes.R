#' Run Bayesian model diagnostics
#'
#' Computes MCMC diagnostics (Rhat, ESS, divergences), model comparison
#' statistics (LOO-CV via PSIS), posterior predictive checks, variance
#' component summaries, and spatial autocorrelation of posterior mean residuals.
#'
#' @param fit_bayes A fitted Bayesian model from [dmft_fit_bayes()].
#' @param config A [dmft_config] object.
#'
#' @returns A list with elements: `mcmc`, `fit_stats`, `variance_components`,
#'   `residuals`, `spatial`, `validity`.
#' @export
dmft_diagnose_bayes <- function(fit_bayes, config) {
  check_suggests("cmdstanr")
  check_suggests("posterior")
  stopifnot(inherits(fit_bayes, "dmft_fit_bayes"))

  stan_fit  <- fit_bayes$fit
  mdata     <- fit_bayes$mdata
  dentition <- fit_bayes$dentition
  dmft_max  <- if (dentition == "deciduous") config$dmft_max_deciduous
               else config$dmft_max_permanent

  results <- list()

  # -- MCMC diagnostics --
  summ <- stan_fit$summary(
    variables = c("alpha", "sigma_region", "sigma_year", "sigma"),
    posterior::default_summary_measures()
  )

  diag_df <- stan_fit$diagnostic_summary()
  n_divergent <- sum(diag_df$num_divergent)
  max_treedepth_hits <- sum(diag_df$num_max_treedepth)

  results$mcmc <- list(
    summary = summ,
    n_divergent = n_divergent,
    max_treedepth_hits = max_treedepth_hits,
    all_rhat_ok = all(summ$rhat < 1.01, na.rm = TRUE),
    min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(summ$ess_tail, na.rm = TRUE)
  )

  # -- LOO-CV --
  results$fit_stats <- tryCatch({
    check_suggests("loo")
    log_lik <- stan_fit$draws("log_lik", format = "matrix")
    loo_result <- loo::loo(log_lik)
    list(
      loo     = loo_result,
      looic   = loo_result$estimates["looic", "Estimate"],
      p_loo   = loo_result$estimates["p_loo", "Estimate"],
      elpd    = loo_result$estimates["elpd_loo", "Estimate"]
    )
  }, error = function(e) {
    cli_alert_warning("LOO-CV failed: {e$message}")
    list(loo = NULL, looic = NA, p_loo = NA, elpd = NA)
  })

  # -- Variance components --
  draws <- stan_fit$draws(format = "draws_df")
  results$variance_components <- data.frame(
    group    = c("region_std", "year_factor", "Residual"),
    mean_sd  = c(
      mean(draws$sigma_region),
      mean(draws$sigma_year),
      mean(draws$sigma)
    ),
    sd_sd    = c(
      stats::sd(draws$sigma_region),
      stats::sd(draws$sigma_year),
      stats::sd(draws$sigma)
    ),
    mean_var = c(
      mean(draws$sigma_region^2),
      mean(draws$sigma_year^2),
      mean(draws$sigma^2)
    )
  )

  # -- Residuals --
  d <- mdata$data
  y_rep <- stan_fit$draws("y_rep", format = "matrix")
  y_rep_mean <- colMeans(y_rep)
  observed <- d$y
  resid <- observed - y_rep_mean
  std_resid <- resid / stats::sd(resid, na.rm = TRUE)

  results$residuals <- list(
    rmse = sqrt(mean(resid^2, na.rm = TRUE)),
    mae  = mean(abs(resid), na.rm = TRUE),
    bias = mean(resid, na.rm = TRUE),
    normality_p = tryCatch(
      stats::shapiro.test(utils::head(std_resid[!is.na(std_resid)], 5000))$p.value,
      error = function(e) NA
    ),
    pp_check = list(
      mean_obs  = mean(observed, na.rm = TRUE),
      mean_pred = mean(y_rep_mean, na.rm = TRUE),
      sd_obs    = stats::sd(observed, na.rm = TRUE),
      sd_pred   = stats::sd(y_rep_mean, na.rm = TRUE)
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
    n_negative  = sum(y_rep_mean < 0, na.rm = TRUE),
    n_above_max = sum(y_rep_mean > dmft_max, na.rm = TRUE),
    range       = range(y_rep_mean, na.rm = TRUE)
  )

  # -- Print summary --
  cli_alert_info("Bayesian diagnostics ({toupper(dentition)})")
  cat(sprintf("  Divergences: %d | Max treedepth hits: %d\n",
              n_divergent, max_treedepth_hits))
  cat(sprintf("  Rhat OK (<1.01): %s | Min ESS bulk: %.0f | Min ESS tail: %.0f\n",
              results$mcmc$all_rhat_ok, results$mcmc$min_ess_bulk,
              results$mcmc$min_ess_tail))
  if (!is.na(results$fit_stats$looic)) {
    cat(sprintf("  LOOIC: %.1f | p_loo: %.1f | ELPD: %.1f\n",
                results$fit_stats$looic, results$fit_stats$p_loo,
                results$fit_stats$elpd))
  }
  cat(sprintf("  RMSE=%.3f  MAE=%.3f  Bias=%.4f\n",
              results$residuals$rmse, results$residuals$mae,
              results$residuals$bias))
  for (i in seq_len(nrow(results$variance_components))) {
    cat(sprintf("  SD(%s): %.4f (posterior SD=%.4f)\n",
                results$variance_components$group[i],
                results$variance_components$mean_sd[i],
                results$variance_components$sd_sd[i]))
  }
  if (!is.na(results$spatial$p.value)) {
    cat(sprintf("  Moran's I: stat=%.4f, p=%.4f %s\n",
                results$spatial$statistic, results$spatial$p.value,
                ifelse(results$spatial$p.value < 0.05,
                       "(spatial autocorrelation detected)", "(OK)")))
  }

  results
}
