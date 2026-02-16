# Load DMFT study-level data

Reads a CSV or Excel file containing study-level DMFT/dmft data. The
file should have at minimum: a region column, `year`, `age_start`,
`age_end`, and at least one of `mean_dmft` or `mean_DMFT`.

## Usage

``` r
dmft_load(file_path, config)
```

## Arguments

- file_path:

  Path to CSV or XLSX file.

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object (used to identify the region column).

## Value

A tibble of raw data.

## Examples

``` r
if (FALSE) { # \dontrun{
cfg <- dmft_config(regions = c("A", "B"), region_col = "province",
                   year_range = c(2000, 2020))
dat <- dmft_load("my_data.csv", config = cfg)
} # }
```
