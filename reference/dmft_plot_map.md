# Plot a choropleth map of DMFT estimates

Plot a choropleth map of DMFT estimates

## Usage

``` r
dmft_plot_map(
  estimates,
  sf_obj,
  year,
  region_name_col = NULL,
  value_col = "predicted",
  title = NULL
)
```

## Arguments

- estimates:

  Predictions from
  [`dmft_predict()`](https://choxos.github.io/dmft/reference/dmft_predict.md).

- sf_obj:

  An `sf` object for region boundaries.

- year:

  Year to display (estimates are averaged across age groups).

- region_name_col:

  Column in `sf_obj` matching `estimates$region`.

- value_col:

  Column to map (default `"predicted"`).

- title:

  Optional plot title.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
