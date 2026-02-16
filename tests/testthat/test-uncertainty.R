test_that("impute_uncertainty uses SE when reported", {
  d <- data.frame(mean_val = 5, se_val = 0.5, sd_val = NA, n = 100)
  result <- impute_uncertainty(d, "mean_val", se_col = "se_val", sd_col = "sd_val")
  expect_equal(result$se_imputed, 0.5)
  expect_equal(result$se_source, "reported")
})

test_that("impute_uncertainty computes SE from SD and n", {
  d <- data.frame(mean_val = 5, se_val = NA, sd_val = 2, n = 100)
  result <- impute_uncertainty(d, "mean_val", se_col = "se_val", sd_col = "sd_val")
  expect_equal(result$se_imputed, 2 / sqrt(100))
  expect_equal(result$se_source, "from_sd_n")
})

test_that("impute_uncertainty uses CV fallback", {
  d <- data.frame(mean_val = 5, n = 100)
  result <- impute_uncertainty(d, "mean_val", cv_default = 0.30)
  expect_equal(result$se_imputed, 5 * 0.30 / sqrt(100))
  expect_equal(result$se_source, "from_cv_n")
})

test_that("impute_uncertainty uses median n when n is missing", {
  d <- data.frame(mean_val = c(5, 10), n = c(200, NA))
  result <- impute_uncertainty(d, "mean_val", cv_default = 0.30)
  # Second row should use median of reported n (200)
  expect_equal(result$se_source[2], "from_cv_assumed_n")
  expect_equal(result$n_effective[2], 200)
})

test_that("impute_uncertainty handles all-NA input", {
  d <- data.frame(mean_val = NA_real_, n = NA_real_)
  result <- impute_uncertainty(d, "mean_val")
  expect_true(is.na(result$se_imputed))
})
