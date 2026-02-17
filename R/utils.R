#' Default ggplot2 theme for DMFT plots
#'
#' @returns A [ggplot2::theme()] object.
#' @export
theme_dmft <- function() {
  ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 10),
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
}


#' Check that a package is available
#'
#' @param pkg Package name.
#' @param reason Why it is needed.
#' @returns `TRUE` invisibly, or aborts.
#' @keywords internal
check_suggests <- function(pkg, reason = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    msg <- sprintf("Package '%s' is required", pkg)
    if (!is.null(reason)) msg <- paste0(msg, " ", reason)
    cli_abort(c(msg, i = "Install it with: install.packages('{pkg}')"))
  }
  invisible(TRUE)
}
