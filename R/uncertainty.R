#' Impute uncertainty for DMFT data
#'
#' Applies a hierarchical imputation strategy:
#' 1. Use SE if reported
#' 2. Compute SE from SD and n
#' 3. Compute SE from confidence intervals
#' 4. Estimate SD from coefficient of variation, then compute SE
#'
#' @param data A data frame.
#' @param mean_col Name of the mean column.
#' @param sd_col Name of the SD column (or `NULL`).
#' @param se_col Name of the SE column (or `NULL`).
#' @param n_col Name of the sample-size column (default `"n"`).
#' @param lci_col Name of the lower CI column (or `NULL`).
#' @param uci_col Name of the upper CI column (or `NULL`).
#' @param cv_default Default coefficient of variation when SD is missing.
#'
#' @returns The input data frame with additional columns:
#'   `se_imputed`, `sd_imputed`, `se_source`, `n_effective`.
#' @export
impute_uncertainty <- function(data,
                                mean_col,
                                sd_col = NULL,
                                se_col = NULL,
                                n_col = "n",
                                lci_col = NULL,
                                uci_col = NULL,
                                cv_default = 0.30) {

  data$se_imputed <- NA_real_
  data$sd_imputed <- NA_real_
  data$se_source  <- NA_character_
  data$n_effective <- NA_real_

  mean_vals <- data[[mean_col]]
  n_vals <- if (!is.null(n_col) && n_col %in% names(data)) data[[n_col]] else rep(NA_real_, nrow(data))
  sd_vals <- if (!is.null(sd_col) && sd_col %in% names(data)) data[[sd_col]] else rep(NA_real_, nrow(data))
  se_vals <- if (!is.null(se_col) && se_col %in% names(data)) data[[se_col]] else rep(NA_real_, nrow(data))
  lci_vals <- if (!is.null(lci_col) && lci_col %in% names(data)) data[[lci_col]] else rep(NA_real_, nrow(data))
  uci_vals <- if (!is.null(uci_col) && uci_col %in% names(data)) data[[uci_col]] else rep(NA_real_, nrow(data))

  # Median sample size from reported values (for fallback)
  reported_n <- n_vals[!is.na(n_vals) & n_vals > 0]
  assumed_n <- if (length(reported_n) > 0) median(reported_n) else 100

  for (i in seq_len(nrow(data))) {
    # Priority 1: SE reported
    if (!is.na(se_vals[i]) && se_vals[i] > 0) {
      data$se_imputed[i] <- se_vals[i]
      data$se_source[i]  <- "reported"
      if (!is.na(n_vals[i]) && n_vals[i] > 0) {
        data$sd_imputed[i] <- se_vals[i] * sqrt(n_vals[i])
      }
      next
    }
    # Priority 2: SD + n
    if (!is.na(sd_vals[i]) && sd_vals[i] > 0 &&
        !is.na(n_vals[i]) && n_vals[i] > 0) {
      data$se_imputed[i] <- sd_vals[i] / sqrt(n_vals[i])
      data$sd_imputed[i] <- sd_vals[i]
      data$se_source[i]  <- "from_sd_n"
      next
    }
    # Priority 3: CI
    if (!is.na(lci_vals[i]) && !is.na(uci_vals[i])) {
      z <- qnorm(0.975)
      data$se_imputed[i] <- (uci_vals[i] - lci_vals[i]) / (2 * z)
      data$se_source[i]  <- "from_ci"
      next
    }
    # Priority 4: CV
    if (!is.na(mean_vals[i]) && mean_vals[i] > 0) {
      imp_sd <- mean_vals[i] * cv_default
      data$sd_imputed[i] <- imp_sd
      ni <- if (!is.na(n_vals[i]) && n_vals[i] > 0) n_vals[i] else assumed_n
      data$se_imputed[i] <- imp_sd / sqrt(ni)
      data$se_source[i]  <- if (!is.na(n_vals[i]) && n_vals[i] > 0) "from_cv_n"
                            else "from_cv_assumed_n"
      data$n_effective[i] <- ni
    }
  }

  # Effective n where missing
  has_se <- !is.na(data$se_imputed) & data$se_imputed > 0
  data$n_effective <- ifelse(
    is.na(data$n_effective) & !is.na(mean_vals) & has_se,
    (mean_vals / data$se_imputed)^2,
    data$n_effective
  )

  data
}
