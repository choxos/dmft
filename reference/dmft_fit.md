# Fit a Bayesian hierarchical model for DMFT/dmft

Implements a GBD-style spatial-temporal-age model with BYM2 spatial
random effects, RW2 temporal trends, RW1 age effects, and IID sex
effects. Uses Gaussian meta-regression (identity link) by default.

## Usage

``` r
dmft_fit(
  data,
  adjacency,
  dentition = c("deciduous", "permanent"),
  config,
  family = "gaussian",
  include_interactions = TRUE
)
```

## Arguments

- data:

  A data frame (one of the dentition subsets from
  [`dmft_clean()`](https://choxos.github.io/dmft/reference/dmft_clean.md)).

- adjacency:

  Adjacency object from
  [`dmft_adjacency()`](https://choxos.github.io/dmft/reference/dmft_adjacency.md).

- dentition:

  `"deciduous"` or `"permanent"`.

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

- family:

  Likelihood family: `"gaussian"` (default, meta-regression) or
  `"nbinomial"`.

- include_interactions:

  Logical; include space-time, space-age, and time-age IID interaction
  random effects.

## Value

A fitted model object (INLA or wrapped cmdstanr fit).
