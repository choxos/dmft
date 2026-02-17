#' Apply AST (Age-Spatial-Temporal) smoothing to model residuals
#'
#' Implements the AST kernel-smoothing algorithm following Foreman et al. (2012)
#' and the CRAN AST package methodology. Residuals from a mixed-effects model
#' are smoothed across three dimensions (space, time, age) using kernel weights,
#' then added back to predictions.
#'
#' @param fit A fitted lme4 model from [dmft_fit()].
#' @param adjacency Adjacency object from [dmft_adjacency()].
#' @param config A [dmft_config] object.
#' @param dentition `"deciduous"` or `"permanent"`.
#'
#' @returns A list with elements:
#'   \describe{
#'     \item{estimates}{Tibble with columns: region, year, age_group, predicted, residual_ast.}
#'     \item{space_matrix}{Spatial weight matrix.}
#'     \item{time_matrix}{Temporal weight matrix.}
#'     \item{age_matrix}{Age weight matrix.}
#'   }
#' @export
dmft_ast <- function(fit, adjacency, config, dentition = c("deciduous", "permanent")) {
  dentition <- match.arg(dentition)
  stopifnot(inherits(config, "dmft_config"))

  mdata <- attr(fit, "dmft_mdata")
  ap <- config$ast_params

  age_groups <- if (dentition == "deciduous") config$age_groups_deciduous
                else config$age_groups_permanent
  n_age <- length(age_groups)

  # Build weight matrices
  space_mat <- calc_space_matrix(adjacency$adj_matrix, ap$par_space)
  time_mat  <- calc_time_matrix(config$year_start, config$year_end, ap$par_time)
  age_mat   <- calc_age_matrix(n_age, ap$par_age)

  # Extract residuals from observed data
  obs_data <- mdata$data
  obs_data$residual <- stats::residuals(fit)
  obs_data$location <- match(obs_data$region_std, config$regions)
  obs_data$age_num  <- obs_data$age_group_idx

  # Prepare residual data frame for AST (columns: age, year, location, residual)
  resid_df <- data.frame(
    age      = obs_data$age_num,
    year     = obs_data$year,
    location = obs_data$location,
    residual = obs_data$residual,
    stringsAsFactors = FALSE
  )

  # Add coverage column if available
  if ("coverage" %in% names(obs_data)) {
    resid_df$coverage <- obs_data$coverage
  }

  # Apply AST smoothing
  ast_result <- apply_ast_smoothing(
    resid_df, space_mat, time_mat, age_mat,
    n_age, config$year_start, config$year_end,
    ap$weight_coverage
  )

  # Create full prediction grid with AST-smoothed residuals
  grid <- expand.grid(
    region    = config$regions,
    year      = config$historical_years,
    age_group = age_groups,
    stringsAsFactors = FALSE
  )
  grid$location <- match(grid$region, config$regions)
  grid$age_num  <- match(grid$age_group, age_groups)

  # Merge AST residuals into grid
  grid$id <- paste(grid$year, grid$age_num, grid$location, sep = "-")
  ast_result$id <- paste(ast_result$year, ast_result$age, ast_result$location, sep = "-")
  grid <- merge(grid, ast_result[, c("id", "residual_AST")], by = "id", all.x = TRUE)
  grid$residual_AST[is.na(grid$residual_AST)] <- 0

  # Get mixed-model predictions for the grid
  grid$predicted_base <- predict_lmer_grid(fit, grid, config, dentition)

  # Final AST estimate
  dmft_max <- if (dentition == "deciduous") config$dmft_max_deciduous
              else config$dmft_max_permanent
  grid$predicted <- pmax(0, pmin(dmft_max, grid$predicted_base + grid$residual_AST))

  # Clean up and return as tibble
  result <- dplyr::tibble(
    region       = grid$region,
    year         = grid$year,
    age_group    = grid$age_group,
    predicted    = grid$predicted,
    residual_ast = grid$residual_AST
  )
  result <- result[order(result$region, result$year, result$age_group), ]

  list(
    estimates    = result,
    space_matrix = space_mat,
    time_matrix  = time_mat,
    age_matrix   = age_mat
  )
}


# -- AST weight matrix functions -----------------------------------------------

#' Calculate spatial weight matrix for AST
#'
#' Transforms binary adjacency matrix into spatial weights.
#' Diagonal = par_space, adjacent cells = par_space * (1 - par_space),
#' non-adjacent = 0.
#'
#' @param adj_matrix Binary adjacency matrix (0/1).
#' @param par_space Spatial correlation parameter (default 0.9).
#' @returns Square spatial weight matrix.
#' @keywords internal
calc_space_matrix <- function(adj_matrix, par_space = 0.9) {
  mat <- adj_matrix
  diag(mat) <- par_space
  mat[mat == 1] <- par_space * (1 - par_space)
  # Ensure numeric rownames for AST compatibility
  rownames(mat) <- as.character(seq_len(nrow(mat)))
  colnames(mat) <- as.character(seq_len(ncol(mat)))
  mat
}


#' Calculate temporal weight matrix for AST
#'
#' LOESS-style cubic power weighting: `W_T(i,j) = (1 - (|i-j|/T)^lambda)^3`
#'
#' @param min_year Start year.
#' @param max_year End year.
#' @param par_time Temporal smoothing parameter (lambda, default 2).
#' @returns Square temporal weight matrix.
#' @keywords internal
calc_time_matrix <- function(min_year, max_year, par_time = 2) {
  years <- min_year:max_year
  n <- length(years)
  mat <- matrix(NA_real_, n, n)
  rownames(mat) <- as.character(years)
  colnames(mat) <- as.character(years)

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mat[i, j] <- (1 - (abs(years[i] - years[j]) / n)^par_time)^3
    }
  }
  mat
}


#' Calculate age weight matrix for AST
#'
#' Exponential decay: `W_A(i,j) = 1 / exp(omega * |i-j|)`
#'
#' @param n_age Number of age groups.
#' @param par_age Age smoothing parameter (omega, default 1).
#' @returns Square age weight matrix.
#' @keywords internal
calc_age_matrix <- function(n_age, par_age = 1) {
  mat <- matrix(NA_real_, n_age, n_age)
  rownames(mat) <- as.character(seq_len(n_age))
  colnames(mat) <- as.character(seq_len(n_age))

  for (i in seq_len(n_age)) {
    for (j in seq_len(n_age)) {
      mat[i, j] <- 1 / exp(par_age * abs(i - j))
    }
  }
  mat
}


#' Apply AST smoothing to residuals
#'
#' Core algorithm: combines spatial, temporal, and age weight matrices via
#' Kronecker products, applies coverage weighting, normalizes, and computes
#' weighted average of residuals for each prediction cell.
#'
#' @param resid_df Data frame with columns: age, year, location, residual,
#'   and optionally coverage.
#' @param space_mat Spatial weight matrix.
#' @param time_mat Temporal weight matrix.
#' @param age_mat Age weight matrix.
#' @param n_age Number of age groups.
#' @param min_year Start year.
#' @param max_year End year.
#' @param weight_coverage Coverage weight parameter (default 0.9).
#' @returns Data frame with columns: year, age, location, residual_AST.
#' @keywords internal
apply_ast_smoothing <- function(resid_df, space_mat, time_mat, age_mat,
                                 n_age, min_year, max_year,
                                 weight_coverage = 0.9) {
  n_loc <- nrow(space_mat)
  n_year <- max_year - min_year + 1L
  loc_ids <- as.integer(rownames(space_mat))

  # Create full prediction grid
  pred_grid <- expand.grid(
    year     = min_year:max_year,
    age      = seq_len(n_age),
    location = loc_ids
  )

  # Handle coverage
  if (!"coverage" %in% names(resid_df)) {
    resid_df$coverage <- 1L
    weight_coverage <- 1
  }

  # Remove NAs and filter to valid ranges
  resid_df <- resid_df[!is.na(resid_df$year) & !is.na(resid_df$age) &
                          !is.na(resid_df$location) & !is.na(resid_df$residual), ]
  resid_df <- resid_df[resid_df$age %in% seq_len(n_age) &
                          resid_df$location %in% loc_ids &
                          resid_df$year %in% (min_year:max_year), ]

  if (nrow(resid_df) == 0) {
    pred_grid$residual_AST <- 0
    return(pred_grid)
  }

  # Step 1: Combine age and time matrices via Kronecker product
  # The AST package uses: zy = kronecker(age_mat, time_mat)
  zy <- kronecker(age_mat, time_mat)

  # Normalize rows of zy
  row_sums <- rowSums(zy, na.rm = TRUE)
  row_sums[row_sums == 0] <- 1
  zy <- zy / row_sums

  # Step 2: Full weight matrix = kronecker(space_mat, normalized_zy)
  W <- kronecker(space_mat, zy)

  # Step 3: Match columns to observed data
  # Row/column ordering: for each location, for each age, for each year
  n_cells <- n_loc * n_age * n_year
  cell_ids <- character(n_cells)
  idx <- 0L
  for (loc in loc_ids) {
    for (ag in seq_len(n_age)) {
      for (yr in min_year:max_year) {
        idx <- idx + 1L
        cell_ids[idx] <- paste(yr, ag, loc, sep = "-")
      }
    }
  }
  colnames(W) <- cell_ids
  rownames(W) <- cell_ids

  # Observed data IDs
  resid_df$id <- paste(resid_df$year, resid_df$age, resid_df$location, sep = "-")
  resid_df <- resid_df[order(resid_df$location, resid_df$age, resid_df$year), ]

  # Subset W columns to only observed data
  obs_ids <- resid_df$id
  valid_ids <- obs_ids[obs_ids %in% cell_ids]

  if (length(valid_ids) == 0) {
    pred_grid$residual_AST <- 0
    return(pred_grid)
  }

  W_obs <- W[, valid_ids, drop = FALSE]
  resid_matched <- resid_df[resid_df$id %in% valid_ids, ]

  # Step 4: Apply coverage weighting
  cov_vals <- resid_matched$coverage
  for (col_idx in seq_len(ncol(W_obs))) {
    w <- if (cov_vals[col_idx] == 1) weight_coverage else (1 - weight_coverage)
    W_obs[, col_idx] <- W_obs[, col_idx] * w
  }

  # Step 5: Re-normalize rows
  row_sums <- rowSums(W_obs, na.rm = TRUE)
  row_sums[row_sums == 0] <- 1
  W_final <- W_obs / row_sums

  # Step 6: Compute weighted average of residuals
  resid_vec <- resid_matched$residual
  ast_residuals <- as.numeric(W_final %*% resid_vec)

  # Map back to prediction grid
  pred_grid$id <- paste(pred_grid$year, pred_grid$age, pred_grid$location, sep = "-")
  pred_grid$residual_AST <- ast_residuals[match(pred_grid$id, cell_ids)]
  pred_grid$residual_AST[is.na(pred_grid$residual_AST)] <- 0

  pred_grid[, c("year", "age", "location", "residual_AST")]
}


#' Predict from lme4 model for a full grid
#' @keywords internal
predict_lmer_grid <- function(fit, grid, config, dentition) {
  # Build new data for prediction
  newdata <- data.frame(
    region_std  = grid$region,
    year_factor = factor(grid$year, levels = levels(attr(fit, "dmft_mdata")$data$year_factor)),
    age_group   = grid$age_group,
    stringsAsFactors = FALSE
  )

  # Add covariates (use means for missing values)
  if (!is.null(config$covariates)) {
    orig_data <- attr(fit, "dmft_mdata")$data
    for (cv in config$covariates) {
      if (cv %in% names(orig_data)) {
        newdata[[cv]] <- mean(orig_data[[cv]], na.rm = TRUE)
      } else {
        newdata[[cv]] <- 0
      }
    }
  }

  # Predict with random effects (allow new levels)
  preds <- tryCatch(
    stats::predict(fit, newdata = newdata, allow.new.levels = TRUE),
    error = function(e) {
      # Fallback: use fixed effects only
      stats::predict(fit, newdata = newdata, re.form = NA)
    }
  )

  as.numeric(preds)
}
