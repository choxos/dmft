# Run model diagnostics

Computes fit statistics (DIC, WAIC), residual diagnostics (RMSE, MAE),
spatial autocorrelation of residuals (Moran's I), and prediction
validity checks.

## Usage

``` r
dmft_diagnose(fit, config)
```

## Arguments

- fit:

  A fitted model from
  [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

## Value

A list with elements: `fit_stats`, `residuals`, `spatial`, `validity`.
