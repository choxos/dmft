# Create spatial adjacency structure from a shapefile

Builds the neighbourhood graph needed for the AST spatial weighting.
Uses Queen contiguity (shared boundary or vertex). Isolated regions
(islands) are connected to their nearest neighbour.

## Usage

``` r
dmft_adjacency(
  shapefile_path = NULL,
  sf_obj = NULL,
  region_name_col = NULL,
  config
)
```

## Arguments

- shapefile_path:

  Path to a shapefile (`.shp`).

- sf_obj:

  Alternatively, an `sf` object already in memory. Exactly one of
  `shapefile_path` or `sf_obj` must be provided.

- region_name_col:

  Column in the shapefile containing region names. If `NULL`, the
  function tries to detect it automatically.

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object. Region ordering in the adjacency matrix will match
  `config$regions`.

## Value

A list with elements:

- adj_matrix:

  Binary adjacency matrix (n x n).

- nb:

  `spdep` neighbourhood object.

- sf:

  The `sf` object (reordered to match config).

- location_ids:

  Named integer vector mapping region names to numeric IDs.

## Examples

``` r
if (FALSE) { # \dontrun{
adj <- dmft_adjacency(shapefile_path = "provinces.shp", config = cfg)
} # }
```
