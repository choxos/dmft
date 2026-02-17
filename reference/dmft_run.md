# Run the full DMFT analysis pipeline

Orchestrates the complete workflow: data loading, cleaning, spatial
structure, mixed-effects model fitting, AST smoothing, bootstrap
uncertainty, diagnostics, and optional projections.

## Usage

``` r
dmft_run(
  data_path,
  shapefile_path,
  config,
  dentition = c("both", "deciduous", "permanent"),
  backend = c("frequentist", "bayesian"),
  region_name_col = NULL,
  skip_projections = FALSE
)
```

## Arguments

- data_path:

  Path to study-level CSV/XLSX data.

- shapefile_path:

  Path to regional shapefile.

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

- dentition:

  Which dentition(s) to model: `"both"`, `"deciduous"`, or
  `"permanent"`.

- backend:

  Modeling backend: `"frequentist"` (default, lme4) or `"bayesian"`
  (experimental, cmdstanr). The Bayesian backend requires cmdstanr,
  posterior, and loo to be installed.

- region_name_col:

  Column in the shapefile containing region names.

- skip_projections:

  Skip the projection step.

## Value

A list with elements `config`, `clean_data`, `adjacency`, `fits`,
`estimates`, `diagnostics`, and optionally `projections`.

## Examples

``` r
if (FALSE) { # \dontrun{
cfg <- dmft_config(
  regions = c("Region_A", "Region_B", "Region_C"),
  region_col = "province",
  year_range = c(2000, 2020),
  projection_range = c(2021, 2030)
)
results <- dmft_run("data.csv", "regions.shp", config = cfg)
} # }
```
