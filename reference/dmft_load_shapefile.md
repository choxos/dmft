# Load a shapefile for spatial analysis

Load a shapefile for spatial analysis

## Usage

``` r
dmft_load_shapefile(shapefile_path)
```

## Arguments

- shapefile_path:

  Path to `.shp` file (or any format
  [`sf::st_read()`](https://r-spatial.github.io/sf/reference/st_read.html)
  supports).

## Value

An `sf` object.
