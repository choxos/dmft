#' Plot a choropleth map of DMFT estimates
#'
#' @param estimates Predictions from [dmft_predict()].
#' @param sf_obj An `sf` object for region boundaries.
#' @param year Year to display (estimates are averaged across age groups).
#' @param region_name_col Column in `sf_obj` matching `estimates$region`.
#' @param value_col Column to map (default `"predicted"`).
#' @param title Optional plot title.
#'
#' @returns A [ggplot2::ggplot] object.
#' @export
dmft_plot_map <- function(estimates,
                           sf_obj,
                           year,
                           region_name_col = NULL,
                           value_col = "predicted",
                           title = NULL) {

  est_year <- estimates[estimates$year == year, ]
  est_agg <- dplyr::summarize(
    dplyr::group_by(est_year, region),
    value = mean(.data[[value_col]], na.rm = TRUE),
    .groups = "drop"
  )

  if (is.null(region_name_col)) {
    # Try to find
    for (col in names(sf_obj)) {
      if (is.character(sf_obj[[col]]) || is.factor(sf_obj[[col]])) {
        if (any(est_agg$region %in% as.character(sf_obj[[col]]))) {
          region_name_col <- col
          break
        }
      }
    }
    if (is.null(region_name_col)) {
      cli_abort("Cannot detect region name column in sf_obj. Set `region_name_col`.")
    }
  }

  sf_data <- dplyr::left_join(sf_obj, est_agg,
                               by = stats::setNames("region", region_name_col))

  p <- ggplot2::ggplot(sf_data) +
    ggplot2::geom_sf(ggplot2::aes(fill = value),
                     color = "grey30", linewidth = 0.3) +
    ggplot2::scale_fill_viridis_c(na.value = "grey80", name = value_col) +
    theme_dmft() +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.5, "cm")
    )

  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  else p <- p + ggplot2::labs(title = sprintf("DMFT Estimates - %d", year))
  p
}


#' Plot temporal trends by region
#'
#' @param estimates Predictions from [dmft_predict()].
#' @param projections Optional projections from [dmft_project()].
#' @param regions Character vector of regions to plot (default: all).
#' @param title Optional title.
#'
#' @returns A [ggplot2::ggplot] object.
#' @export
dmft_plot_trends <- function(estimates,
                              projections = NULL,
                              regions = NULL,
                              title = "DMFT Temporal Trends") {

  # Aggregate across age groups
  est_agg <- dplyr::summarize(
    dplyr::group_by(estimates, region, year),
    predicted = mean(predicted, na.rm = TRUE),
    lower = mean(lower, na.rm = TRUE),
    upper = mean(upper, na.rm = TRUE),
    .groups = "drop"
  )

  if (!is.null(regions)) est_agg <- est_agg[est_agg$region %in% regions, ]

  p <- ggplot2::ggplot(est_agg, ggplot2::aes(x = year, y = predicted)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = 0.2, fill = "steelblue") +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::facet_wrap(~ region, scales = "free_y") +
    ggplot2::labs(title = title, x = "Year", y = "Predicted DMFT") +
    theme_dmft()

  # Add projections if available
  if (!is.null(projections)) {
    p <- p +
      ggplot2::geom_ribbon(
        data = projections,
        ggplot2::aes(x = year, ymin = lower, ymax = upper, fill = scenario),
        alpha = 0.15, inherit.aes = FALSE
      ) +
      ggplot2::geom_line(
        data = projections,
        ggplot2::aes(x = year, y = mean_proj, color = scenario),
        linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE
      )
  }

  p
}
