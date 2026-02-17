# dmft: Age-Spatial-Temporal Modeling of Dental Caries

Implements Age-Spatial-Temporal (AST) models for estimating dental
caries burden (DMFT/dmft indices) at subnational level. Uses a two-stage
approach: (1) random intercept mixed-effects model via lme4 (or Stan for
the experimental Bayesian backend), and (2) AST kernel-smoothing of
residuals across space, time, and age dimensions.

## Details

The main workflow functions are:

- [`dmft_config()`](https://choxos.github.io/dmft/reference/dmft_config.md):
  Set up analysis parameters

- [`dmft_load()`](https://choxos.github.io/dmft/reference/dmft_load.md):
  Load study-level data

- [`dmft_clean()`](https://choxos.github.io/dmft/reference/dmft_clean.md):
  Clean and standardize data

- [`dmft_adjacency()`](https://choxos.github.io/dmft/reference/dmft_adjacency.md):
  Build spatial adjacency from shapefile

- [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md):
  Fit mixed-effects model (frequentist)

- [`dmft_predict()`](https://choxos.github.io/dmft/reference/dmft_predict.md):
  AST smoothing + bootstrap uncertainty

- [`dmft_diagnose()`](https://choxos.github.io/dmft/reference/dmft_diagnose.md):
  Model diagnostics

- [`dmft_project()`](https://choxos.github.io/dmft/reference/dmft_project.md):
  Future trend projections

- [`dmft_run()`](https://choxos.github.io/dmft/reference/dmft_run.md):
  Full pipeline in one call

Bayesian alternatives (experimental):
[`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md),
[`dmft_predict_bayes()`](https://choxos.github.io/dmft/reference/dmft_predict_bayes.md),
[`dmft_diagnose_bayes()`](https://choxos.github.io/dmft/reference/dmft_diagnose_bayes.md),
[`dmft_project_bayes()`](https://choxos.github.io/dmft/reference/dmft_project_bayes.md).

## References

Shoaee S, et al. (2022). Subnational estimation of dental caries burden.
*BMC Oral Health*, 22:634.

Shoaee S, et al. (2025). Subnational estimation using AST models. *BMC
Oral Health*, 25:1490.

## See also

Useful links:

- <https://github.com/choxos/dmft>

- Report bugs at <https://github.com/choxos/dmft/issues>

## Author

**Maintainer**: Ahmad Sofi-Mahmudi <ahmad.sofimahmudi@gmail.com>
([ORCID](https://orcid.org/0000-0002-2460-2394)) \[copyright holder\]
