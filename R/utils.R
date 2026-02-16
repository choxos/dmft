#' Compute BYM2 scaling factor from ICAR precision matrix
#'
#' The scaling factor is the geometric mean of the non-zero eigenvalues
#' of the ICAR precision matrix Q (Riebler et al. 2016). This ensures
#' proper variance interpretation of the BYM2 spatial random effect.
#'
#' @param node1 Integer vector of first nodes in each edge.
#' @param node2 Integer vector of second nodes in each edge.
#' @param n_nodes Total number of nodes (regions).
#'
#' @returns Scaling factor (numeric scalar).
#' @export
#'
#' @examples
#' # Simple chain graph: 1-2-3-4
#' compute_bym2_scaling_factor(c(1, 2, 3), c(2, 3, 4), 4)
compute_bym2_scaling_factor <- function(node1, node2, n_nodes) {
  Q <- matrix(0, n_nodes, n_nodes)
  for (e in seq_along(node1)) {
    i <- node1[e]
    j <- node2[e]
    Q[i, j] <- -1
    Q[j, i] <- -1
    Q[i, i] <- Q[i, i] + 1
    Q[j, j] <- Q[j, j] + 1
  }

  eigenvalues <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  nonzero_eig <- eigenvalues[eigenvalues > 1e-10]

  if (length(nonzero_eig) == 0) {
    warning("No non-zero eigenvalues; using fallback scaling factor = 1")
    return(1.0)
  }

  exp(mean(log(nonzero_eig)))
}


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
