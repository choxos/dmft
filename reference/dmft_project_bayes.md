# Project DMFT trends using Bayesian posterior draws

Experimental Bayesian alternative to
[`dmft_project()`](https://choxos.github.io/dmft/reference/dmft_project.md).
Extracts year random effect posterior draws, estimates a temporal trend,
and generates scenario projections with uncertainty propagated from the
posterior.

## Usage

``` r
dmft_project_bayes(fit_bayes, estimates, config, damping_factor = 0.95)
```

## Arguments

- fit_bayes:

  A fitted Bayesian model from
  [`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md).

- estimates:

  Predictions from
  [`dmft_predict_bayes()`](https://choxos.github.io/dmft/reference/dmft_predict_bayes.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object with projection range set.

- damping_factor:

  Trend damping per year (default 0.95).

## Value

A tibble with columns: `year`, `scenario`, `mean_proj`, `lower`,
`upper`.
