# dmft: Bayesian Hierarchical Modeling of Dental Caries (DMFT/dmft)

Implements GBD-style Bayesian hierarchical models for estimating and
projecting dental caries burden (DMFT/dmft indices) at subnational
level. Combines BYM2 spatial random effects, RW2 temporal trends, and
RW1 age effects with Gaussian meta-regression. Supports any country
given study-level data and a regional shapefile. Provides a complete
pipeline from data cleaning through model fitting, diagnostics,
projections, and publication-ready visualizations. Two inference
backends are supported: INLA (fast approximate Bayesian) and Stan (full
MCMC).

## See also

Useful links:

- <https://github.com/choxos/dmft>

- Report bugs at <https://github.com/choxos/dmft/issues>

## Author

**Maintainer**: Ahmad Sofi-Mahmudi <ahmad.sofimahmudi@gmail.com>
([ORCID](https://orcid.org/0000-0002-2460-2394)) \[copyright holder\]
