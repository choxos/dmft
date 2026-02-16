# Create a DMFT analysis configuration

Defines all parameters for a DMFT/dmft analysis: regions, time period,
age groups, model priors, and projection settings.

## Usage

``` r
dmft_config(
  regions,
  region_col = "province",
  year_range = c(1990, 2025),
  projection_range = NULL,
  age_groups_deciduous = c("0-4", "5-9", "10-14"),
  age_groups_permanent = c("5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39",
    "40-44", "45-49", "50-54", "55-59", "60+"),
  dmft_max_deciduous = 20L,
  dmft_max_permanent = 28L,
  cv_default = 0.3,
  backend = c("inla", "cmdstanr"),
  seed = 12345L,
  n_draws = 250L,
  n_cores = max(1L, parallel::detectCores() - 1L),
  region_mapping = NULL,
  stan_settings = list(chains = 4, parallel_chains = 4, iter_warmup = 1000, iter_sampling
    = 1000, adapt_delta = 0.95, max_treedepth = 12)
)
```

## Arguments

- regions:

  Character vector of region names (e.g., provinces).

- region_col:

  Column name in the data containing region names.

- year_range:

  Length-2 integer vector `c(start, end)` for historical period.

- projection_range:

  Optional length-2 vector `c(start, end)` for projections. If `NULL`,
  no projections are configured.

- age_groups_deciduous:

  Character vector of deciduous age group labels.

- age_groups_permanent:

  Character vector of permanent age group labels.

- dmft_max_deciduous:

  Maximum biologically plausible dmft (default 20).

- dmft_max_permanent:

  Maximum biologically plausible DMFT (default 28).

- cv_default:

  Default coefficient of variation for uncertainty imputation.

- backend:

  Inference backend: `"inla"` or `"cmdstanr"`.

- seed:

  Random seed for reproducibility.

- n_draws:

  Number of posterior draws.

- n_cores:

  Number of CPU cores.

- region_mapping:

  Optional named list mapping alternative region names to canonical
  names in `regions`.

- stan_settings:

  Named list of Stan MCMC settings. Only used when
  `backend = "cmdstanr"`.

## Value

An object of class `"dmft_config"`.

## Examples

``` r
cfg <- dmft_config(
  regions = c("Region_A", "Region_B", "Region_C"),
  region_col = "province",
  year_range = c(2000, 2020)
)
cfg
#> ℹ DMFT analysis configuration
#>   Regions:      3 (Region_A, Region_B, Region_C, ...)
#>   Historical:   2000-2020 (21 years)
#>   Age groups:   deciduous=3, permanent=12
#>   Backend:      inla
#>   Seed:         12345
```
