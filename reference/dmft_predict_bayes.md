# Generate predictions from a Bayesian DMFT model with AST smoothing

Experimental Bayesian alternative to
[`dmft_predict()`](https://choxos.github.io/dmft/reference/dmft_predict.md).
Extracts posterior draws of fixed and random effects, computes
predictions for the full region-year-age grid, applies AST smoothing on
posterior mean residuals, and returns credible intervals from the
posterior.

## Usage

``` r
dmft_predict_bayes(fit_bayes, adjacency, config)
```

## Arguments

- fit_bayes:

  A fitted Bayesian model from
  [`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md).

- adjacency:

  Adjacency object from
  [`dmft_adjacency()`](https://choxos.github.io/dmft/reference/dmft_adjacency.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

## Value

A tibble with columns: `region`, `year`, `age_group`, `predicted`,
`lower`, `upper`.
