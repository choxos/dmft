# Calculate temporal weight matrix for AST

LOESS-style cubic power weighting: `W_T(i,j) = (1 - (|i-j|/T)^lambda)^3`

## Usage

``` r
calc_time_matrix(min_year, max_year, par_time = 2)
```

## Arguments

- min_year:

  Start year.

- max_year:

  End year.

- par_time:

  Temporal smoothing parameter (lambda, default 2).

## Value

Square temporal weight matrix.
