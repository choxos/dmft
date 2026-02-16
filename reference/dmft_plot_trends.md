# Plot temporal trends by region

Plot temporal trends by region

## Usage

``` r
dmft_plot_trends(
  estimates,
  projections = NULL,
  regions = NULL,
  title = "DMFT Temporal Trends"
)
```

## Arguments

- estimates:

  Predictions from
  [`dmft_predict()`](https://choxos.github.io/dmft/reference/dmft_predict.md).

- projections:

  Optional projections from
  [`dmft_project()`](https://choxos.github.io/dmft/reference/dmft_project.md).

- regions:

  Character vector of regions to plot (default: all).

- title:

  Optional title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
