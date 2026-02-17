test_that("calc_space_matrix produces correct dimensions and values", {
  # 3 regions: A-B adjacent, B-C adjacent, A-C not adjacent
  adj <- matrix(c(0, 1, 0,
                  1, 0, 1,
                  0, 1, 0), nrow = 3)
  rownames(adj) <- colnames(adj) <- c("A", "B", "C")

  sp <- calc_space_matrix(adj, par_space = 0.9)

  expect_equal(nrow(sp), 3)
  expect_equal(ncol(sp), 3)
  # Diagonal = par_space
  expect_equal(sp[1, 1], 0.9)
  expect_equal(sp[2, 2], 0.9)
  # Adjacent = par_space * (1 - par_space)
  expect_equal(sp[1, 2], 0.9 * 0.1)
  expect_equal(sp[2, 1], 0.9 * 0.1)
  # Non-adjacent = 0
  expect_equal(sp[1, 3], 0)
})

test_that("calc_time_matrix produces symmetric positive matrix", {
  tm <- calc_time_matrix(2000, 2005, par_time = 2)

  expect_equal(nrow(tm), 6)
  expect_equal(ncol(tm), 6)
  # Symmetric
  expect_equal(tm, t(tm))
  # Diagonal = 1
  expect_equal(unname(diag(tm)), rep(1, 6))
  # All values positive
  expect_true(all(tm > 0))
  # Off-diagonal < 1
  expect_true(all(tm[row(tm) != col(tm)] < 1))
})

test_that("calc_age_matrix produces correct exponential decay", {
  am <- calc_age_matrix(4, par_age = 1)

  expect_equal(nrow(am), 4)
  # Diagonal = 1
  expect_equal(unname(diag(am)), rep(1, 4))
  # Adjacent = 1/exp(1)
  expect_equal(am[1, 2], 1 / exp(1))
  # Distance 2 = 1/exp(2)
  expect_equal(am[1, 3], 1 / exp(2))
  # Symmetric
  expect_equal(am, t(am))
})

test_that("apply_ast_smoothing returns correct structure", {
  # Simple case: 2 regions, 2 years, 2 age groups
  adj <- matrix(c(0, 1, 1, 0), nrow = 2)
  rownames(adj) <- colnames(adj) <- c("R1", "R2")

  space_mat <- calc_space_matrix(adj, 0.9)
  time_mat  <- calc_time_matrix(2000, 2001, 2)
  age_mat   <- calc_age_matrix(2, 1)

  resid_df <- data.frame(
    region_std  = rep(c("R1", "R2"), each = 4),
    year        = rep(rep(c(2000, 2001), each = 2), 2),
    age_group   = rep(c("0-4", "5-9"), 4),
    residual    = rnorm(8)
  )

  result <- apply_ast_smoothing(resid_df, space_mat, time_mat, age_mat,
                                 n_age = 2, min_year = 2000, max_year = 2001)

  expect_true("residual_AST" %in% names(result))
  expect_equal(nrow(result), 8)
  # AST-smoothed residuals should be finite
  expect_true(all(is.finite(result$residual_AST)))
})
