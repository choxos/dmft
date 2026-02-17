# Project DMFT trends into the future

Extracts the temporal trend from a fitted model and generates three
scenarios (reference, optimistic, pessimistic) using linear trend
extrapolation with optional damping and bounded uncertainty.

## Usage

``` r
dmft_project(fit, estimates, config, damping_factor = 0.95)
```

## Arguments

- fit:

  A fitted model from
  [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md).

- estimates:

  Predictions from
  [`dmft_predict()`](https://choxos.github.io/dmft/reference/dmft_predict.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object. Must have `projection_start` and `projection_end` set.

- damping_factor:

  Trend damping per year (default 0.95).

## Value

A tibble with columns: `year`, `scenario`, `mean_proj`, `lower`,
`upper`.
