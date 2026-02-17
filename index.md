# dmft

**dmft** implements the Age-Spatial-Temporal (AST) model for estimating
and projecting dental caries burden (DMFT/dmft indices) at subnational
level for any country.

**Documentation & Tutorial**: <https://choxos.github.io/dmft/>

## Features

- **Two-stage AST modeling**: mixed-effects model + kernel-smoothed
  residuals across age, space, and time
- **Two backends**: frequentist (lme4, default) and Bayesian (cmdstanr,
  experimental)
- **Spatial smoothing** via adjacency-based weight matrices from
  shapefiles
- **Temporal smoothing** using LOESS-inspired cubic power weights
- **Age smoothing** with exponential decay kernels
- **Uncertainty quantification**: bootstrap (frequentist) or posterior
  credible intervals (Bayesian)
- **Projection framework** with reference, optimistic, and pessimistic
  scenarios
- **Publication-ready visualizations**: choropleth maps, temporal
  trends, uncertainty ribbons
- **Country-agnostic**: works with any shapefile and study-level data

## Installation

``` r
# Install from GitHub
pak::pak("choxos/dmft")

# Or with remotes
remotes::install_github("choxos/dmft")
```

For the experimental Bayesian backend, also install cmdstanr:

``` r
install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", getOption("repos")))
cmdstanr::install_cmdstan()
```

## Quick Start

``` r
library(dmft)

# 1. Configure the analysis for your country
cfg <- dmft_config(
  regions = c("Region_A", "Region_B", "Region_C", "Region_D"),
  region_col = "province",
  year_range = c(2000, 2020),
  projection_range = c(2021, 2030)
)

# 2. Run the full pipeline
results <- dmft_run(
  data_path = "my_data.csv",
  shapefile_path = "my_regions.shp",
  config = cfg
)

# 3. Visualize
dmft_plot_map(results$estimates$permanent, results$adjacency$sf, year = 2020)
dmft_plot_trends(results$estimates$permanent, results$projections$permanent)
```

See the [full
tutorial](https://choxos.github.io/dmft/articles/tutorial.html) for a
complete walkthrough.

### Step-by-Step Usage

``` r
# Load and clean data
raw <- dmft_load("data.csv", config = cfg)
clean <- dmft_clean(raw, config = cfg)

# Create spatial structure
adj <- dmft_adjacency(shapefile_path = "regions.shp", config = cfg)

# Fit mixed-effects model (frequentist)
fit <- dmft_fit(clean$permanent, adj, dentition = "permanent", config = cfg)

# AST smoothing + bootstrap uncertainty
estimates <- dmft_predict(fit, adj, config = cfg)

# Diagnostics
diagnostics <- dmft_diagnose(fit, config = cfg)

# Project future trends
projections <- dmft_project(fit, estimates, config = cfg)
```

### Bayesian Backend (Experimental)

``` r
# Fit with Stan instead of lme4
fit_b <- dmft_fit_bayes(clean$permanent, adj, dentition = "permanent", config = cfg)

# Predictions with posterior credible intervals
estimates_b <- dmft_predict_bayes(fit_b, adj, config = cfg)

# Bayesian diagnostics (Rhat, ESS, LOO-CV)
diag_b <- dmft_diagnose_bayes(fit_b, config = cfg)

# Or run the full pipeline with backend = "bayesian"
results <- dmft_run("data.csv", "regions.shp", config = cfg, backend = "bayesian")
```

## Input Data Format

The input CSV/XLSX should contain study-level summary statistics with at
minimum:

| Column                            | Description                                               |
|-----------------------------------|-----------------------------------------------------------|
| `province` (or your `region_col`) | Region name                                               |
| `year`                            | Year of data collection                                   |
| `age_start`                       | Start of age range                                        |
| `age_end`                         | End of age range                                          |
| `mean_dmft` and/or `mean_DMFT`    | Mean dmft (deciduous) / DMFT (permanent)                  |
| `n`                               | Sample size (optional, used for weighting)                |
| `sd_dmft` / `se_dmft`             | Standard deviation / error (optional, imputed if missing) |

## Model

The package uses a two-stage AST (Age-Spatial-Temporal) approach
following Shoaee et al. (2022, 2025):

**Stage 1 — Mixed-effects model**

    y_i ~ Normal(mu_i, sigma^2 + se_i^2)

    mu_i = beta_0 + [covariates] + b_region + b_year

    b_region ~ Normal(0, sigma_region^2)   [random intercepts]
    b_year   ~ Normal(0, sigma_year^2)     [random intercepts]

The frequentist backend estimates parameters via REML (lme4). The
Bayesian backend uses HMC/NUTS sampling (Stan) with weakly informative
priors.

**Stage 2 — AST kernel smoothing of residuals**

Residuals from Stage 1 are smoothed using weighted averaging across
three dimensions:

- **Spatial**: adjacency-based weights (par_space = 0.9)
- **Temporal**: LOESS cubic power decay (par_time = 2)
- **Age**: exponential decay (par_age = 1)

The combined weight matrix is the Kronecker product of these three
weight matrices, row-normalized to sum to 1.

## References

- Shoaee S, et al. (2022). Subnational estimation of dental caries
  burden. *BMC Oral Health*, 22:634.
- Shoaee S, et al. (2025). Subnational estimation using AST models. *BMC
  Oral Health*, 25:1490.

## License

MIT
