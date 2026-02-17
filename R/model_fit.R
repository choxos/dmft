#' Fit a random intercept mixed-effects model for DMFT/dmft
#'
#' Stage 1 of the AST methodology: fits a random intercept mixed-effects model
#' via [lme4::lmer()] with optional fixed covariates and random province/year
#' intercepts. Inverse-variance weighting is applied when standard errors are
#' available.
#'
#' @param data A data frame (one of the dentition subsets from [dmft_clean()]).
#' @param adjacency Adjacency object from [dmft_adjacency()].
#' @param dentition `"deciduous"` or `"permanent"`.
#' @param config A [dmft_config] object.
#'
#' @returns A fitted lme4 model object with metadata attributes.
#' @export
dmft_fit <- function(data,
                      adjacency,
                      dentition = c("deciduous", "permanent"),
                      config) {

  dentition <- match.arg(dentition)
  stopifnot(inherits(config, "dmft_config"))

  # Prepare model data
  mdata <- prepare_lmer_data(data, adjacency, dentition, config)

  # Build formula
  mean_col <- mdata$mean_col
  formula_parts <- paste0("y ~ 1")

  # Add fixed covariates if specified
  if (!is.null(config$covariates)) {
    available_covs <- intersect(config$covariates, names(mdata$data))
    if (length(available_covs) > 0) {
      formula_parts <- paste(formula_parts, "+",
                              paste(available_covs, collapse = " + "))
      cli_alert_info("Including covariates: {paste(available_covs, collapse=', ')}")
    }
    missing_covs <- setdiff(config$covariates, available_covs)
    if (length(missing_covs) > 0) {
      cli_alert_warning("Covariates not found in data: {paste(missing_covs, collapse=', ')}")
    }
  }

  # Add random intercepts for region and year
  formula_parts <- paste(formula_parts, "+ (1|region_std) + (1|year_factor)")
  model_formula <- stats::as.formula(formula_parts)

  cli_alert_info("Fitting lme4 model ({dentition}): {deparse(model_formula)}")
  t0 <- Sys.time()

  # Fit with inverse-variance weights if available
  if ("precision" %in% names(mdata$data) && any(mdata$data$precision != 1)) {
    fit <- lme4::lmer(
      model_formula,
      data    = mdata$data,
      weights = mdata$data$precision,
      REML    = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 100000))
    )
  } else {
    fit <- lme4::lmer(
      model_formula,
      data    = mdata$data,
      REML    = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa",
                                   optCtrl = list(maxfun = 100000))
    )
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  aic_val <- stats::AIC(fit)
  bic_val <- stats::BIC(fit)
  cli_alert_success("lme4 fit complete in {round(elapsed, 1)}s (AIC={round(aic_val, 1)}, BIC={round(bic_val, 1)})")

  # Attach metadata
  attr(fit, "dmft_mdata")     <- mdata
  attr(fit, "dmft_config")    <- config
  attr(fit, "dmft_dentition") <- dentition
  fit
}


# -- data preparation ---------------------------------------------------------

#' Prepare data for lme4 mixed model
#' @keywords internal
prepare_lmer_data <- function(data, adjacency, dentition, config) {
  mean_col <- if (dentition == "deciduous") "mean_dmft" else "mean_DMFT"
  se_col   <- "se_imputed"
  age_groups <- if (dentition == "deciduous") config$age_groups_deciduous
                else config$age_groups_permanent

  # Filter valid records
  data <- data[!is.na(data[[mean_col]]) & !is.na(data$region_std) &
                 data$region_std %in% config$regions &
                 !is.na(data$year) &
                 data$year >= config$year_start &
                 data$year <= config$year_end, ]

  if (nrow(data) == 0) cli_abort("No valid records after filtering.")

  # Response
  data$y <- data[[mean_col]]

  # Indices
  data$province_idx   <- match(data$region_std, config$regions)
  data$year_idx       <- data$year - config$year_start + 1L
  data$age_group_idx  <- match(data$age_group, age_groups)

  # Year as factor for random intercept
  data$year_factor <- factor(data$year)

  # Precision weights (inverse-variance)
  if (se_col %in% names(data)) {
    se_vals <- data[[se_col]]
    data$precision <- ifelse(!is.na(se_vals) & se_vals > 0, 1 / se_vals^2, 1)
  } else {
    data$precision <- 1
  }
  data$precision <- data$precision / max(data$precision, na.rm = TRUE)

  cli_alert_info("Model data: {nrow(data)} obs, {length(unique(data$province_idx))} regions, years {min(data$year)}-{max(data$year)}")

  list(
    data       = data,
    n_regions  = config$n_regions,
    n_years    = config$n_years,
    n_age      = length(age_groups),
    age_groups = age_groups,
    adjacency  = adjacency,
    dentition  = dentition,
    mean_col   = mean_col
  )
}
