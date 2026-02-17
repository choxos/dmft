#' Run the full DMFT analysis pipeline
#'
#' Orchestrates the complete workflow: data loading, cleaning, spatial
#' structure, mixed-effects model fitting, AST smoothing, bootstrap
#' uncertainty, diagnostics, and optional projections.
#'
#' @param data_path Path to study-level CSV/XLSX data.
#' @param shapefile_path Path to regional shapefile.
#' @param config A [dmft_config] object.
#' @param dentition Which dentition(s) to model: `"both"`, `"deciduous"`,
#'   or `"permanent"`.
#' @param backend Modeling backend: `"frequentist"` (default, lme4) or
#'   `"bayesian"` (experimental, cmdstanr). The Bayesian backend requires
#'   cmdstanr, posterior, and loo to be installed.
#' @param region_name_col Column in the shapefile containing region names.
#' @param skip_projections Skip the projection step.
#'
#' @returns A list with elements `config`, `clean_data`, `adjacency`,
#'   `fits`, `estimates`, `diagnostics`, and optionally `projections`.
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- dmft_config(
#'   regions = c("Region_A", "Region_B", "Region_C"),
#'   region_col = "province",
#'   year_range = c(2000, 2020),
#'   projection_range = c(2021, 2030)
#' )
#' results <- dmft_run("data.csv", "regions.shp", config = cfg)
#' }
dmft_run <- function(data_path,
                      shapefile_path,
                      config,
                      dentition = c("both", "deciduous", "permanent"),
                      backend = c("frequentist", "bayesian"),
                      region_name_col = NULL,
                      skip_projections = FALSE) {

  dentition <- match.arg(dentition)
  backend   <- match.arg(backend)
  stopifnot(inherits(config, "dmft_config"))

  if (backend == "bayesian") {
    check_suggests("cmdstanr")
    check_suggests("posterior")
  }

  cli_alert_info("Starting DMFT analysis pipeline (AST method, backend={backend})")

  # 1. Load data
  raw_data <- dmft_load(data_path, config)

  # 2. Clean
  clean <- dmft_clean(raw_data, config)

  # 3. Adjacency
  adj <- dmft_adjacency(
    shapefile_path = shapefile_path,
    region_name_col = region_name_col,
    config = config
  )

  results <- list(
    config     = config,
    clean_data = clean,
    adjacency  = adj,
    fits       = list(),
    estimates  = list(),
    diagnostics = list(),
    projections = list()
  )

  # 4. Model fitting + AST smoothing + predictions
  dents <- switch(dentition,
    both      = c("deciduous", "permanent"),
    deciduous = "deciduous",
    permanent = "permanent"
  )

  for (dent in dents) {
    d <- clean[[dent]]
    if (is.null(d) || nrow(d) == 0) {
      cli_alert_warning("No {dent} data; skipping")
      next
    }

    if (backend == "frequentist") {
      # Stage 1: Fit mixed-effects model (lme4)
      fit <- dmft_fit(d, adj, dentition = dent, config = config)

      # Stage 2: AST smoothing + bootstrap uncertainty
      est <- dmft_predict(fit, adj, config)

      # Diagnostics
      diag <- dmft_diagnose(fit, config)

      results$fits[[dent]]       <- fit
      results$estimates[[dent]]  <- est
      results$diagnostics[[dent]] <- diag

      # Projections
      if (!skip_projections && !is.null(config$projection_start)) {
        proj <- dmft_project(fit, est, config)
        results$projections[[dent]] <- proj
      }
    } else {
      # Bayesian backend (experimental)
      fit <- dmft_fit_bayes(d, adj, dentition = dent, config = config)

      # AST smoothing + posterior credible intervals
      est <- dmft_predict_bayes(fit, adj, config)

      # Bayesian diagnostics
      diag <- dmft_diagnose_bayes(fit, config)

      results$fits[[dent]]       <- fit
      results$estimates[[dent]]  <- est
      results$diagnostics[[dent]] <- diag

      # Projections
      if (!skip_projections && !is.null(config$projection_start)) {
        proj <- dmft_project_bayes(fit, est, config)
        results$projections[[dent]] <- proj
      }
    }
  }

  cli_alert_success("Pipeline complete")
  results
}
