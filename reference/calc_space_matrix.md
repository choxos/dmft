# Calculate spatial weight matrix for AST

Transforms binary adjacency matrix into spatial weights. Diagonal =
par_space, adjacent cells = par_space \* (1 - par_space), non-adjacent =
0.

## Usage

``` r
calc_space_matrix(adj_matrix, par_space = 0.9)
```

## Arguments

- adj_matrix:

  Binary adjacency matrix (0/1).

- par_space:

  Spatial correlation parameter (default 0.9).

## Value

Square spatial weight matrix.
