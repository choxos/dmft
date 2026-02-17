# Run Bayesian model diagnostics

Computes MCMC diagnostics (Rhat, ESS, divergences), model comparison
statistics (LOO-CV via PSIS), posterior predictive checks, variance
component summaries, and spatial autocorrelation of posterior mean
residuals.

## Usage

``` r
dmft_diagnose_bayes(fit_bayes, config)
```

## Arguments

- fit_bayes:

  A fitted Bayesian model from
  [`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

## Value

A list with elements: `mcmc`, `fit_stats`, `variance_components`,
`residuals`, `spatial`, `validity`.
