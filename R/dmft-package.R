#' @title dmft: Age-Spatial-Temporal Modeling of Dental Caries
#'
#' @description Implements Age-Spatial-Temporal (AST) models for estimating
#'   dental caries burden (DMFT/dmft indices) at subnational level. Uses a
#'   two-stage approach: (1) random intercept mixed-effects model via lme4
#'   (or Stan for the experimental Bayesian backend), and (2) AST
#'   kernel-smoothing of residuals across space, time, and age dimensions.
#'
#' @details
#' The main workflow functions are:
#'
#' - [dmft_config()]: Set up analysis parameters
#' - [dmft_load()]: Load study-level data
#' - [dmft_clean()]: Clean and standardize data
#' - [dmft_adjacency()]: Build spatial adjacency from shapefile
#' - [dmft_fit()]: Fit mixed-effects model (frequentist)
#' - [dmft_predict()]: AST smoothing + bootstrap uncertainty
#' - [dmft_diagnose()]: Model diagnostics
#' - [dmft_project()]: Future trend projections
#' - [dmft_run()]: Full pipeline in one call
#'
#' Bayesian alternatives (experimental): [dmft_fit_bayes()],
#' [dmft_predict_bayes()], [dmft_diagnose_bayes()], [dmft_project_bayes()].
#'
#' @references
#' Shoaee S, et al. (2022). Subnational estimation of dental caries burden.
#' *BMC Oral Health*, 22:634.
#'
#' Shoaee S, et al. (2025). Subnational estimation using AST models.
#' *BMC Oral Health*, 25:1490.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats as.formula coef fitted lm median predict qnorm quantile
#' @importFrom stats residuals rnorm runif sd setNames shapiro.test sigma var
#' @importFrom stats deviance
#' @importFrom utils head
#' @importFrom rlang .data .env := sym `%||%`
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning cli_abort
NULL
