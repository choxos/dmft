# Clean and standardize DMFT data

Performs column standardization, region name mapping, sex
standardization, age group assignment, outlier flagging, and value
validation. Returns separate deciduous and permanent datasets ready for
modeling.

## Usage

``` r
dmft_clean(data, config)
```

## Arguments

- data:

  A data frame from
  [`dmft_load()`](https://choxos.github.io/dmft/reference/dmft_load.md).

- config:

  A
  [dmft_config](https://choxos.github.io/dmft/reference/dmft_config.md)
  object.

## Value

A list with elements `deciduous` and `permanent`, each a tibble with
standardized columns and imputed uncertainty.
