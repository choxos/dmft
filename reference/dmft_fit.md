# Fit a random intercept mixed-effects model for DMFT/dmft

Stage 1 of the AST methodology: fits a random intercept mixed-effects
model via [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) with
optional fixed covariates and random province/year intercepts.
Inverse-variance weighting is applied when standard errors are
available.

## Usage

``` r
dmft_fit(data, adjacency, dentition = c("deciduous", "permanent"), config)
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

## Value

A fitted lme4 model object with metadata attributes.
