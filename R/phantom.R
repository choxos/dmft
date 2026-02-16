#' Generate synthetic (phantom) DMFT data for testing
#'
#' Creates realistic test data following known age-by-region patterns,
#' useful for validating the analysis pipeline before using real data.
#'
#' @param config A [dmft_config] object.
#' @param dentition `"deciduous"`, `"permanent"`, or `"both"`.
#' @param n_per_cell Number of synthetic studies per region-year-age cell.
#' @param temporal_trend Annual proportional change (negative = improvement).
#' @param missing_prob Probability of a cell being missing (sparse data).
#'
#' @returns A tibble of synthetic study-level data.
#' @export
#'
#' @examples
#' cfg <- dmft_config(
#'   regions = c("A", "B", "C"),
#'   region_col = "province",
#'   year_range = c(2000, 2020)
#' )
#' phantom <- dmft_phantom(cfg, dentition = "permanent")
#' head(phantom)
dmft_phantom <- function(config,
                          dentition = c("both", "deciduous", "permanent"),
                          n_per_cell = 1L,
                          temporal_trend = -0.015,
                          missing_prob = 0.15) {

  dentition <- match.arg(dentition)
  stopifnot(inherits(config, "dmft_config"))
  set.seed(config$seed)

  results <- list()

  if (dentition %in% c("both", "deciduous")) {
    results$deciduous <- generate_phantom(
      config, "deciduous", n_per_cell, temporal_trend, missing_prob
    )
  }
  if (dentition %in% c("both", "permanent")) {
    results$permanent <- generate_phantom(
      config, "permanent", n_per_cell, temporal_trend, missing_prob
    )
  }

  out <- dplyr::bind_rows(results)
  cli_alert_success("Generated {nrow(out)} phantom records")
  out
}


#' @keywords internal
generate_phantom <- function(config, dent, n_per_cell, trend, missing_prob) {
  if (dent == "deciduous") {
    age_groups <- config$age_groups_deciduous
    baselines  <- c(2.5, 4.5, 6.0)[seq_along(age_groups)]
    sds        <- c(2.0, 2.5, 3.0)[seq_along(age_groups)]
    mean_col   <- "mean_dmft"
  } else {
    age_groups <- config$age_groups_permanent
    n_ag <- length(age_groups)
    baselines  <- seq(0.5, 22, length.out = n_ag)
    sds        <- seq(1.0, 7.0, length.out = n_ag)
    mean_col   <- "mean_DMFT"
  }

  # Region effects (random, centered at 0)
  region_effects <- rnorm(config$n_regions, 0, 0.15)
  names(region_effects) <- config$regions

  rows <- list()
  ref_year <- config$year_start

  for (reg in config$regions) {
    for (yr in config$historical_years) {
      for (a_idx in seq_along(age_groups)) {
        if (runif(1) < missing_prob) next

        for (rep in seq_len(n_per_cell)) {
          years_since <- yr - ref_year
          base <- baselines[a_idx] *
            exp(trend * years_since + region_effects[reg])
          n_sample <- sample(50:500, 1)
          obs_mean <- max(0, rnorm(1, base, sds[a_idx] / sqrt(n_sample) * 5))
          obs_sd   <- abs(rnorm(1, sds[a_idx], sds[a_idx] * 0.2))

          # Parse age bounds
          ag <- age_groups[a_idx]
          if (grepl("\\+$", ag)) {
            a_start <- as.integer(sub("\\+$", "", ag))
            a_end   <- a_start + 20L
          } else {
            parts   <- as.integer(strsplit(ag, "-")[[1]])
            a_start <- parts[1]
            a_end   <- parts[2]
          }

          row <- data.frame(
            study     = sprintf("phantom_%s_%d_%s", reg, yr, ag),
            province  = reg,
            year      = yr,
            age_start = a_start,
            age_end   = a_end,
            sex       = sample(c("Male", "Female", "Both"), 1),
            n         = n_sample,
            stringsAsFactors = FALSE
          )
          row[[mean_col]] <- round(obs_mean, 2)
          row[[sub("mean", "sd", mean_col)]] <- round(obs_sd, 2)

          rows[[length(rows) + 1]] <- row
        }
      }
    }
  }

  dplyr::bind_rows(rows)
}
