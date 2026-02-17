# Apply AST (Age-Spatial-Temporal) smoothing to model residuals

Implements the AST kernel-smoothing algorithm following Foreman et al.
(2012) and the CRAN AST package methodology. Residuals from a
mixed-effects model are smoothed across three dimensions (space, time,
age) using kernel weights, then added back to predictions.

## Usage

``` r
dmft_ast(fit, adjacency, config, dentition = c("deciduous", "permanent"))
```

## Arguments

- fit:

  A fitted model from
  [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md)
  (lme4) or
  [`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md)
  (Bayesian). Both are supported transparently.

- adjacency:

  Adjacency object from
  [`dmft_adjacency()`](https://choxos.github.io/dmft/reference/dmft_adjacency.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

- dentition:

  `"deciduous"` or `"permanent"`.

## Value

A list with elements:

- estimates:

  Tibble with columns: region, year, age_group, predicted, residual_ast.

- space_matrix:

  Spatial weight matrix.

- time_matrix:

  Temporal weight matrix.

- age_matrix:

  Age weight matrix.
