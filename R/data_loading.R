#' Load DMFT study-level data
#'
#' Reads a CSV or Excel file containing study-level DMFT/dmft data.
#' The file should have at minimum: a region column, `year`,
#' `age_start`, `age_end`, and at least one of `mean_dmft` or `mean_DMFT`.
#'
#' @param file_path Path to CSV or XLSX file.
#' @param config A [dmft_config] object (used to identify the region column).
#'
#' @returns A tibble of raw data.
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- dmft_config(regions = c("A", "B"), region_col = "province",
#'                    year_range = c(2000, 2020))
#' dat <- dmft_load("my_data.csv", config = cfg)
#' }
dmft_load <- function(file_path, config) {
  stopifnot(inherits(config, "dmft_config"))
  stopifnot(file.exists(file_path))

  ext <- tolower(tools::file_ext(file_path))

  if (ext == "xlsx" || ext == "xls") {
    data <- readxl::read_excel(file_path, sheet = 1)
  } else {
    data <- readr::read_csv(file_path, show_col_types = FALSE)
  }

  cli_alert_success("Loaded {nrow(data)} rows from {basename(file_path)}")
  data
}


#' Load a shapefile for spatial analysis
#'
#' @param shapefile_path Path to `.shp` file (or any format [sf::st_read()]
#'   supports).
#'
#' @returns An `sf` object.
#' @export
dmft_load_shapefile <- function(shapefile_path) {
  stopifnot(file.exists(shapefile_path))
  sf_obj <- sf::st_read(shapefile_path, quiet = TRUE)
  cli_alert_success("Loaded shapefile with {nrow(sf_obj)} features")
  sf_obj
}
