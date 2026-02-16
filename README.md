# dmft <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/choxos/dmft/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/choxos/dmft/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**dmft** implements GBD-style Bayesian hierarchical models for estimating and projecting dental caries burden (DMFT/dmft indices) at subnational level for any country.

## Features

- **Bayesian spatial-temporal-age modeling** with BYM2 spatial effects, RW2 temporal trends, and RW1 age effects
- **Two inference backends**: INLA (fast approximate Bayesian) and Stan (full MCMC via cmdstanr)
- **Gaussian meta-regression** with inverse-variance weighting for study-level summary data
- **Projection framework** with reference, optimistic, and pessimistic scenarios
- **Publication-ready visualizations**: choropleth maps, temporal trends, uncertainty ribbons
- **Country-agnostic**: works with any shapefile and study-level data

## Installation

``` r
# Install from GitHub
pak::pak("choxos/dmft")

# Or with remotes
remotes::install_github("choxos/dmft")
```

The INLA backend (default) requires the INLA package:

``` r
install.packages("INLA",
  repos = c(CRAN = "https://cloud.r-project.org",
            INLA = "https://inla.r-inla-download.org/R/stable"))
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

### Step-by-step usage

``` r
# Load and clean data
raw <- dmft_load("data.csv", config = cfg)
clean <- dmft_clean(raw, config = cfg)

# Create spatial structure
adj <- dmft_adjacency(shapefile_path = "regions.shp", config = cfg)

# Fit model
fit <- dmft_fit(clean$permanent, adj, dentition = "permanent", config = cfg)

# Predict and diagnose
estimates <- dmft_predict(fit, config = cfg)
diagnostics <- dmft_diagnose(fit, config = cfg)

# Project future trends
projections <- dmft_project(fit, estimates, config = cfg)
```

## Input Data Format

The input CSV/XLSX should contain study-level summary statistics with at minimum:

| Column | Description |
|--------|-------------|
| `province` (or your `region_col`) | Region name |
| `year` | Year of data collection |
| `age_start` | Start of age range |
| `age_end` | End of age range |
| `mean_dmft` and/or `mean_DMFT` | Mean dmft (deciduous) / DMFT (permanent) |
| `n` | Sample size (optional, used for weighting) |
| `sd_dmft` / `se_dmft` | Standard deviation / error (optional, imputed if missing) |

## Model

The package fits a Bayesian hierarchical model:

```
y_i ~ Normal(mu_i, se_i^2)    [Gaussian meta-regression]

mu_i = beta_0 + S_region + T_year + A_age + G_sex + interactions

S ~ BYM2(rho, sigma_spatial)   [spatial: ICAR + IID mixture]
T ~ RW2(sigma_temporal)        [temporal: second-order random walk]
A ~ RW1(sigma_age)             [age: first-order random walk]
G ~ IID(sigma_sex)             [sex: exchangeable]
```

## License

MIT
