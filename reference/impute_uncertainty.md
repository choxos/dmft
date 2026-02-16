# Impute uncertainty for DMFT data

Applies a hierarchical imputation strategy:

1.  Use SE if reported

2.  Compute SE from SD and n

3.  Compute SE from confidence intervals

4.  Estimate SD from coefficient of variation, then compute SE

## Usage

``` r
impute_uncertainty(
  data,
  mean_col,
  sd_col = NULL,
  se_col = NULL,
  n_col = "n",
  lci_col = NULL,
  uci_col = NULL,
  cv_default = 0.3
)
```

## Arguments

- data:

  A data frame.

- mean_col:

  Name of the mean column.

- sd_col:

  Name of the SD column (or `NULL`).

- se_col:

  Name of the SE column (or `NULL`).

- n_col:

  Name of the sample-size column (default `"n"`).

- lci_col:

  Name of the lower CI column (or `NULL`).

- uci_col:

  Name of the upper CI column (or `NULL`).

- cv_default:

  Default coefficient of variation when SD is missing.

## Value

The input data frame with additional columns: `se_imputed`,
`sd_imputed`, `se_source`, `n_effective`.
