#' Create a DMFT analysis configuration
#'
#' Defines all parameters for a DMFT/dmft analysis: regions, time period,
#' age groups, AST smoothing parameters, and projection settings.
#'
#' @param regions Character vector of region names (e.g., provinces).
#' @param region_col Column name in the data containing region names.
#' @param year_range Length-2 integer vector `c(start, end)` for historical period.
#' @param projection_range Optional length-2 vector `c(start, end)` for projections.
#'   If `NULL`, no projections are configured.
#' @param age_groups_deciduous Character vector of deciduous age group labels.
#' @param age_groups_permanent Character vector of permanent age group labels.
#' @param dmft_max_deciduous Maximum biologically plausible dmft (default 20).
#' @param dmft_max_permanent Maximum biologically plausible DMFT (default 28).
#' @param cv_default Default coefficient of variation for uncertainty imputation.
#' @param covariates Character vector of covariate column names to include as
#'   fixed effects in the mixed model. If `NULL`, only the intercept is used.
#' @param ast_params Named list of AST smoothing parameters:
#'   \describe{
#'     \item{par_space}{Spatial correlation parameter (default 0.9).}
#'     \item{par_time}{Temporal smoothing parameter, lambda (default 2).}
#'     \item{par_age}{Age smoothing parameter, omega (default 1).}
#'     \item{weight_coverage}{Weight for high-coverage data sources (default 0.9).}
#'   }
#' @param n_boot Number of bootstrap replicates for uncertainty intervals.
#' @param seed Random seed for reproducibility.
#' @param region_mapping Optional named list mapping alternative region names to
#'   canonical names in `regions`.
#'
#' @returns An object of class `"dmft_config"`.
#' @export
#'
#' @examples
#' cfg <- dmft_config(
#'   regions = c("Region_A", "Region_B", "Region_C"),
#'   region_col = "province",
#'   year_range = c(2000, 2020)
#' )
#' cfg
dmft_config <- function(
    regions,
    region_col = "province",
    year_range = c(1990, 2025),
    projection_range = NULL,
    age_groups_deciduous = c("0-4", "5-9", "10-14"),
    age_groups_permanent = c("5-9", "10-14", "15-19", "20-24", "25-29",
                             "30-34", "35-39", "40-44", "45-49", "50-54",
                             "55-59", "60+"),
    dmft_max_deciduous = 20L,
    dmft_max_permanent = 28L,
    cv_default = 0.30,
    covariates = NULL,
    ast_params = list(
      par_space = 0.9,
      par_time = 2,
      par_age = 1,
      weight_coverage = 0.9
    ),
    n_boot = 1000L,
    seed = 12345L,
    region_mapping = NULL) {

  stopifnot(
    is.character(regions), length(regions) >= 2,
    is.character(region_col), length(region_col) == 1,
    is.numeric(year_range), length(year_range) == 2,
    year_range[1] < year_range[2]
  )

  if (!is.null(projection_range)) {
    stopifnot(
      is.numeric(projection_range), length(projection_range) == 2,
      projection_range[1] > year_range[2]
    )
  }

  if (!is.null(covariates)) {
    stopifnot(is.character(covariates))
  }

  # Fill in AST defaults for any missing params
  ast_defaults <- list(par_space = 0.9, par_time = 2, par_age = 1,
                        weight_coverage = 0.9)
  for (nm in names(ast_defaults)) {
    if (is.null(ast_params[[nm]])) ast_params[[nm]] <- ast_defaults[[nm]]
  }

  # Parse age group boundaries
  parse_age_bounds <- function(ag) {
    lower <- upper <- integer(length(ag))
    for (i in seq_along(ag)) {
      if (grepl("\\+$", ag[i])) {
        lower[i] <- as.integer(sub("\\+$", "", ag[i]))
        upper[i] <- 100L
      } else {
        parts <- as.integer(strsplit(ag[i], "-")[[1]])
        lower[i] <- parts[1]
        upper[i] <- parts[2]
      }
    }
    list(lower = lower, upper = upper)
  }

  dec_bounds <- parse_age_bounds(age_groups_deciduous)
  perm_bounds <- parse_age_bounds(age_groups_permanent)

  cfg <- list(
    regions = regions,
    n_regions = length(regions),
    region_col = region_col,
    region_mapping = region_mapping,

    year_start = as.integer(year_range[1]),
    year_end = as.integer(year_range[2]),
    historical_years = seq.int(year_range[1], year_range[2]),
    n_years = as.integer(year_range[2] - year_range[1] + 1L),

    projection_start = if (!is.null(projection_range)) as.integer(projection_range[1]) else NULL,
    projection_end   = if (!is.null(projection_range)) as.integer(projection_range[2]) else NULL,

    age_groups_deciduous = age_groups_deciduous,
    age_groups_deciduous_lower = dec_bounds$lower,
    age_groups_deciduous_upper = dec_bounds$upper,
    n_age_deciduous = length(age_groups_deciduous),

    age_groups_permanent = age_groups_permanent,
    age_groups_permanent_lower = perm_bounds$lower,
    age_groups_permanent_upper = perm_bounds$upper,
    n_age_permanent = length(age_groups_permanent),

    dmft_max_deciduous = as.integer(dmft_max_deciduous),
    dmft_max_permanent = as.integer(dmft_max_permanent),
    cv_default = cv_default,

    covariates = covariates,
    ast_params = ast_params,
    n_boot = as.integer(n_boot),
    seed = as.integer(seed),

    # Projection scenario adjustments
    scenarios = list(
      reference   = list(name = "Reference",   adjustment = 0),
      optimistic  = list(name = "Optimistic",  adjustment = -0.02),
      pessimistic = list(name = "Pessimistic", adjustment = 0.02)
    )
  )

  class(cfg) <- "dmft_config"
  cfg
}


#' @export
print.dmft_config <- function(x, ...) {
  cli_alert_info("DMFT analysis configuration")
  cat(sprintf("  Regions:      %d (%s, ...)\n", x$n_regions,
              paste(head(x$regions, 3), collapse = ", ")))
  cat(sprintf("  Historical:   %d-%d (%d years)\n",
              x$year_start, x$year_end, x$n_years))
  if (!is.null(x$projection_start)) {
    cat(sprintf("  Projections:  %d-%d\n", x$projection_start, x$projection_end))
  }
  cat(sprintf("  Age groups:   deciduous=%d, permanent=%d\n",
              x$n_age_deciduous, x$n_age_permanent))
  cat(sprintf("  AST params:   space=%.1f, time=%g, age=%g\n",
              x$ast_params$par_space, x$ast_params$par_time, x$ast_params$par_age))
  cat(sprintf("  Covariates:   %s\n",
              if (is.null(x$covariates)) "none" else paste(x$covariates, collapse = ", ")))
  cat(sprintf("  Bootstrap:    %d replicates\n", x$n_boot))
  cat(sprintf("  Seed:         %d\n", x$seed))
  invisible(x)
}
