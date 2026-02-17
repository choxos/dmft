# Calculate age weight matrix for AST

Exponential decay: `W_A(i,j) = 1 / exp(omega * |i-j|)`

## Usage

``` r
calc_age_matrix(n_age, par_age = 1)
```

## Arguments

- n_age:

  Number of age groups.

- par_age:

  Age smoothing parameter (omega, default 1).

## Value

Square age weight matrix.
