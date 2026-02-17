#' Create spatial adjacency structure from a shapefile
#'
#' Builds the neighbourhood graph needed for the AST spatial weighting.
#' Uses Queen contiguity (shared boundary or vertex). Isolated regions
#' (islands) are connected to their nearest neighbour.
#'
#' @param shapefile_path Path to a shapefile (`.shp`).
#' @param sf_obj Alternatively, an `sf` object already in memory.
#'   Exactly one of `shapefile_path` or `sf_obj` must be provided.
#' @param region_name_col Column in the shapefile containing region names.
#'   If `NULL`, the function tries to detect it automatically.
#' @param config A [dmft_config] object. Region ordering in the adjacency
#'   matrix will match `config$regions`.
#'
#' @returns A list with elements:
#'   \describe{
#'     \item{adj_matrix}{Binary adjacency matrix (n x n).}
#'     \item{nb}{`spdep` neighbourhood object.}
#'     \item{sf}{The `sf` object (reordered to match config).}
#'     \item{location_ids}{Named integer vector mapping region names to numeric IDs.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' adj <- dmft_adjacency(shapefile_path = "provinces.shp", config = cfg)
#' }
dmft_adjacency <- function(shapefile_path = NULL,
                            sf_obj = NULL,
                            region_name_col = NULL,
                            config) {

  stopifnot(inherits(config, "dmft_config"))
  stopifnot(!is.null(shapefile_path) || !is.null(sf_obj))

  # Load shapefile if needed
  if (is.null(sf_obj)) {
    sf_obj <- sf::st_read(shapefile_path, quiet = TRUE)
  }
  sf_obj <- sf::st_make_valid(sf_obj)

  # Disable s2 for poly2nb
  s2_was_on <- sf::sf_use_s2()
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(s2_was_on), add = TRUE)

  # Detect region name column
  if (is.null(region_name_col)) {
    region_name_col <- detect_region_col(sf_obj, config)
  }

  # Reorder sf to match config$regions
  sf_names <- sf_obj[[region_name_col]]
  match_idx <- match(config$regions, sf_names)
  if (any(is.na(match_idx))) {
    missing <- config$regions[is.na(match_idx)]
    cli_alert_warning("Regions not found in shapefile: {paste(missing, collapse=', ')}")
    match_idx <- match_idx[!is.na(match_idx)]
  }
  sf_ordered <- sf_obj[match_idx, ]

  # Build neighbourhood
  nb <- spdep::poly2nb(sf_ordered, queen = TRUE)

  # Handle islands
  islands <- which(spdep::card(nb) == 0)
  if (length(islands) > 0) {
    cli_alert_info("Connecting {length(islands)} isolated region(s) to nearest neighbour")
    centroids <- sf::st_centroid(sf::st_geometry(sf_ordered))
    for (isl in islands) {
      dists <- as.numeric(sf::st_distance(centroids[isl], centroids))
      dists[isl] <- Inf
      nearest <- which.min(dists)
      nb[[isl]] <- as.integer(nearest)
      if (!isl %in% nb[[nearest]]) {
        nb[[nearest]] <- sort(as.integer(c(nb[[nearest]], isl)))
      }
    }
  }

  # Adjacency matrix (binary: 0/1)
  adj_matrix <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  n <- nrow(adj_matrix)
  matched_names <- config$regions[config$regions %in% sf_names]
  rownames(adj_matrix) <- matched_names
  colnames(adj_matrix) <- matched_names

  # Numeric location IDs for AST (1, 2, ..., n)
  location_ids <- stats::setNames(seq_len(n), matched_names)

  cli_alert_success("Adjacency: {n} regions, {sum(adj_matrix) / 2} edges")

  list(
    adj_matrix   = adj_matrix,
    nb           = nb,
    sf           = sf_ordered,
    location_ids = location_ids
  )
}


#' Detect the region name column in an sf object
#' @keywords internal
detect_region_col <- function(sf_obj, config) {
  nms <- names(sf_obj)
  # Try common column names
  for (candidate in c("PRENAME", "PRNAME", "NAME", "name", "NAME_1",
                       "region", "Region", "province", "Province",
                       "state", "State", "admin1")) {
    if (candidate %in% nms) {
      vals <- sf_obj[[candidate]]
      if (any(config$regions %in% vals)) return(candidate)
    }
  }
  # Try every character column
  for (col in nms) {
    if (is.character(sf_obj[[col]]) || is.factor(sf_obj[[col]])) {
      vals <- as.character(sf_obj[[col]])
      if (sum(config$regions %in% vals) >= 2) return(col)
    }
  }
  cli_abort("Could not detect region name column in shapefile.
             Set `region_name_col` explicitly.")
}
