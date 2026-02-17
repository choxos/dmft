# Package index

## Configuration

Set up your analysis

- [`dmft_config()`](https://choxos.github.io/dmft/reference/dmft_config.md)
  : Create a DMFT analysis configuration

## Data Pipeline

Load, clean, and prepare data

- [`dmft_load()`](https://choxos.github.io/dmft/reference/dmft_load.md)
  : Load DMFT study-level data
- [`dmft_load_shapefile()`](https://choxos.github.io/dmft/reference/dmft_load_shapefile.md)
  : Load a shapefile for spatial analysis
- [`dmft_clean()`](https://choxos.github.io/dmft/reference/dmft_clean.md)
  : Clean and standardize DMFT data
- [`impute_uncertainty()`](https://choxos.github.io/dmft/reference/impute_uncertainty.md)
  : Impute uncertainty for DMFT data

## Spatial Structure

Build adjacency for spatial modeling

- [`dmft_adjacency()`](https://choxos.github.io/dmft/reference/dmft_adjacency.md)
  : Create spatial adjacency structure from a shapefile

## Model Fitting (Frequentist)

lme4 mixed-effects models (default)

- [`dmft_fit()`](https://choxos.github.io/dmft/reference/dmft_fit.md) :
  Fit a random intercept mixed-effects model for DMFT/dmft
- [`dmft_predict()`](https://choxos.github.io/dmft/reference/dmft_predict.md)
  : Generate predictions from a fitted DMFT model with AST smoothing
- [`dmft_diagnose()`](https://choxos.github.io/dmft/reference/dmft_diagnose.md)
  : Run model diagnostics
- [`dmft_ast()`](https://choxos.github.io/dmft/reference/dmft_ast.md) :
  Apply AST (Age-Spatial-Temporal) smoothing to model residuals

## Model Fitting (Bayesian)

Experimental Stan/cmdstanr backend

- [`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md)
  : Fit a Bayesian random intercept model for DMFT/dmft via cmdstanr
- [`dmft_predict_bayes()`](https://choxos.github.io/dmft/reference/dmft_predict_bayes.md)
  : Generate predictions from a Bayesian DMFT model with AST smoothing
- [`dmft_diagnose_bayes()`](https://choxos.github.io/dmft/reference/dmft_diagnose_bayes.md)
  : Run Bayesian model diagnostics

## Projections

Future trend scenarios

- [`dmft_project()`](https://choxos.github.io/dmft/reference/dmft_project.md)
  : Project DMFT trends into the future
- [`dmft_project_bayes()`](https://choxos.github.io/dmft/reference/dmft_project_bayes.md)
  : Project DMFT trends using Bayesian posterior draws

## Visualization

Maps and trend plots

- [`dmft_plot_map()`](https://choxos.github.io/dmft/reference/dmft_plot_map.md)
  : Plot a choropleth map of DMFT estimates
- [`dmft_plot_trends()`](https://choxos.github.io/dmft/reference/dmft_plot_trends.md)
  : Plot temporal trends by region
- [`theme_dmft()`](https://choxos.github.io/dmft/reference/theme_dmft.md)
  : Default ggplot2 theme for DMFT plots

## Pipeline

Full analysis workflow

- [`dmft_run()`](https://choxos.github.io/dmft/reference/dmft_run.md) :
  Run the full DMFT analysis pipeline
- [`dmft_phantom()`](https://choxos.github.io/dmft/reference/dmft_phantom.md)
  : Generate synthetic (phantom) DMFT data for testing
