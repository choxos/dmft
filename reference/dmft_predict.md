# Generate predictions from a fitted DMFT model

Produces predictions for all region-year-age combinations in the
historical period. Uses posterior sampling when available, falling back
to a delta-method approximation.

## Usage

``` r
dmft_predict(fit, config, n_posterior = 100L)
```

## Arguments

- fit:

  A fitted model from
  [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

- n_posterior:

  Number of posterior samples for uncertainty. Set to 0 to use the delta
  method.

## Value

A tibble with columns: `region`, `year`, `age_group`, `predicted`,
`lower`, `upper`.
