#' Create spatial adjacency structure from a shapefile
#'
#' Builds the neighbourhood graph needed for BYM2 spatial random effects.
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
#'     \item{graph_path}{Path to the INLA graph file (temp file).}
#'     \item{node1, node2}{Edge list for Stan / BYM2 scaling.}
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

  # Adjacency matrix
  adj_matrix <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  n <- nrow(adj_matrix)
  matched_names <- config$regions[config$regions %in% sf_names]
  rownames(adj_matrix) <- matched_names
  colnames(adj_matrix) <- matched_names

  # INLA graph file (temp)
  graph_path <- tempfile(fileext = ".graph")
  spdep::nb2INLA(graph_path, nb)

  # Edge list for Stan
  node1 <- node2 <- integer(0)
  for (i in seq_len(n)) {
    for (j in nb[[i]]) {
      if (j > i) {
        node1 <- c(node1, i)
        node2 <- c(node2, j)
      }
    }
  }

  cli_alert_success("Adjacency: {n} regions, {length(node1)} edges")


  list(
    adj_matrix = adj_matrix,
    nb         = nb,
    sf         = sf_ordered,
    graph_path = graph_path,
    node1      = node1,
    node2      = node2
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
