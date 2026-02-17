# Generate predictions from a fitted DMFT model with AST smoothing

Produces predictions for all region-year-age combinations in the
historical period. Uses AST (Age-Spatial-Temporal) kernel smoothing on
mixed-model residuals, with bootstrap uncertainty intervals.

## Usage

``` r
dmft_predict(fit, adjacency, config, n_boot = NULL)
```

## Arguments

- fit:

  A fitted model from
  [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md).

- adjacency:

  Adjacency object from
  [`dmft_adjacency()`](https://choxos.github.io/dmft/reference/dmft_adjacency.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

- n_boot:

  Number of bootstrap replicates for uncertainty. Set to 0 to skip
  uncertainty estimation.

## Value

A tibble with columns: `region`, `year`, `age_group`, `predicted`,
`lower`, `upper`.
