# Compute BYM2 scaling factor from ICAR precision matrix

The scaling factor is the geometric mean of the non-zero eigenvalues of
the ICAR precision matrix Q (Riebler et al. 2016). This ensures proper
variance interpretation of the BYM2 spatial random effect.

## Usage

``` r
compute_bym2_scaling_factor(node1, node2, n_nodes)
```

## Arguments

- node1:

  Integer vector of first nodes in each edge.

- node2:

  Integer vector of second nodes in each edge.

- n_nodes:

  Total number of nodes (regions).

## Value

Scaling factor (numeric scalar).

## Examples

``` r
# Simple chain graph: 1-2-3-4
compute_bym2_scaling_factor(c(1, 2, 3), c(2, 3, 4), 4)
#> [1] 1.587401
```
