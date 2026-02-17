# Fit a Bayesian random intercept model for DMFT/dmft via cmdstanr

Experimental Bayesian alternative to
[`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md).
Uses the same random intercept structure
(`y ~ 1 + covariates + (1|region) + (1|year)`) but estimates parameters
via Hamiltonian Monte Carlo (Stan) rather than REML. AST smoothing
(Stage 2) is applied downstream via
[`dmft_predict_bayes()`](https://choxos.github.io/dmft/reference/dmft_predict_bayes.md).

## Usage

``` r
dmft_fit_bayes(
  data,
  adjacency,
  dentition = c("deciduous", "permanent"),
  config,
  stan_settings = NULL
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

- stan_settings:

  Named list overriding default Stan sampling settings: `chains`,
  `iter_warmup`, `iter_sampling`, `adapt_delta`, `max_treedepth`. If
  `NULL`, uses defaults from config or built-in defaults.

## Value

A list of class `"dmft_fit_bayes"` with elements `fit` (CmdStanMCMC),
`mdata` (prepared model data), `config`, and `dentition`.
