# Generate synthetic (phantom) DMFT data for testing

Creates realistic test data following known age-by-region patterns,
useful for validating the analysis pipeline before using real data.

## Usage

``` r
dmft_phantom(
  config,
  dentition = c("both", "deciduous", "permanent"),
  n_per_cell = 1L,
  temporal_trend = -0.015,
  missing_prob = 0.15
)
```

## Arguments

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

- dentition:

  `"deciduous"`, `"permanent"`, or `"both"`.

- n_per_cell:

  Number of synthetic studies per region-year-age cell.

- temporal_trend:

  Annual proportional change (negative = improvement).

- missing_prob:

  Probability of a cell being missing (sparse data).

## Value

A tibble of synthetic study-level data.

## Examples

``` r
cfg <- dmft_config(
  regions = c("A", "B", "C"),
  region_col = "province",
  year_range = c(2000, 2020)
)
phantom <- dmft_phantom(cfg, dentition = "permanent")
#> ✔ Generated 645 phantom records
head(phantom)
#>                  study province year age_start age_end    sex   n mean_DMFT
#> 1   phantom_A_2000_5-9        A 2000         5       9 Female 142      0.80
#> 2 phantom_A_2000_15-19        A 2000        15      19   Both 152      4.57
#> 3 phantom_A_2000_20-24        A 2000        20      24 Female 335      7.21
#> 4 phantom_A_2000_25-29        A 2000        25      29 Female 308      6.80
#> 5 phantom_A_2000_30-34        A 2000        30      34   Both  65      8.66
#> 6 phantom_A_2000_35-39        A 2000        35      39   Both 410     11.54
#>   sd_DMFT
#> 1    0.64
#> 2    1.71
#> 3    2.91
#> 4    3.48
#> 5    3.59
#> 6    4.55
```
