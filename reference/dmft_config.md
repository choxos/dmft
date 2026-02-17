# Create a DMFT analysis configuration

Defines all parameters for a DMFT/dmft analysis: regions, time period,
age groups, AST smoothing parameters, and projection settings.

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
  covariates = NULL,
  ast_params = list(par_space = 0.9, par_time = 2, par_age = 1, weight_coverage = 0.9),
  n_boot = 1000L,
  stan_settings = NULL,
  seed = 12345L,
  region_mapping = NULL
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

- covariates:

  Character vector of covariate column names to include as fixed effects
  in the mixed model. If `NULL`, only the intercept is used.

- ast_params:

  Named list of AST smoothing parameters:

  par_space

  :   Spatial correlation parameter (default 0.9).

  par_time

  :   Temporal smoothing parameter, lambda (default 2).

  par_age

  :   Age smoothing parameter, omega (default 1).

  weight_coverage

  :   Weight for high-coverage data sources (default 0.9).

- n_boot:

  Number of bootstrap replicates for uncertainty intervals.

- stan_settings:

  Optional named list of Stan sampling settings for the experimental
  Bayesian backend (see
  [`dmft_fit_bayes()`](https://choxos.github.io/dmft/reference/dmft_fit_bayes.md)).
  Elements: `chains`, `iter_warmup`, `iter_sampling`, `adapt_delta`,
  `max_treedepth`. If `NULL`, sensible defaults are used.

- seed:

  Random seed for reproducibility.

- region_mapping:

  Optional named list mapping alternative region names to canonical
  names in `regions`.

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
#>   AST params:   space=0.9, time=2, age=1
#>   Covariates:   none
#>   Bootstrap:    1000 replicates
#>   Seed:         12345
```
