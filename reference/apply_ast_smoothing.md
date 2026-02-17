# Apply AST smoothing to residuals

Core algorithm: combines spatial, temporal, and age weight matrices via
Kronecker products, applies coverage weighting, normalizes, and computes
weighted average of residuals for each prediction cell.

## Usage

``` r
apply_ast_smoothing(
  resid_df,
  space_mat,
  time_mat,
  age_mat,
  n_age,
  min_year,
  max_year,
  weight_coverage = 0.9
)
```

## Arguments

- resid_df:

  Data frame with columns: age, year, location, residual, and optionally
  coverage.

- space_mat:

  Spatial weight matrix.

- time_mat:

  Temporal weight matrix.

- age_mat:

  Age weight matrix.

- n_age:

  Number of age groups.

- min_year:

  Start year.

- max_year:

  End year.

- weight_coverage:

  Coverage weight parameter (default 0.9).

## Value

Data frame with columns: year, age, location, residual_AST.
