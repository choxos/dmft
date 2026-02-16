test_that("dmft_clean produces deciduous and permanent splits", {
  cfg <- dmft_config(
    regions = c("A", "B"),
    region_col = "province",
    year_range = c(2000, 2020)
  )

  d <- data.frame(
    province  = c("A", "A", "B", "B"),
    year      = c(2010, 2010, 2010, 2010),
    age_start = c(2, 6, 6, 10),
    age_end   = c(4, 9, 9, 14),
    sex       = "Both",
    mean_dmft = c(2.5, 4.0, NA, NA),
    mean_DMFT = c(NA, NA, 1.0, 2.0),
    stringsAsFactors = FALSE
  )

  result <- dmft_clean(d, cfg)
  expect_type(result, "list")
  expect_named(result, c("deciduous", "permanent"))
  expect_equal(nrow(result$deciduous), 2)
  expect_equal(nrow(result$permanent), 2)
})

test_that("standardize_sex works", {
  expect_equal(
    standardize_sex(c("M", "female", "Both", "unknown")),
    c("Male", "Female", "Both", NA)
  )
})

test_that("dmft_clean errors on missing region column", {
  cfg <- dmft_config(regions = c("A", "B"), region_col = "state",
                     year_range = c(2000, 2020))
  d <- data.frame(province = "A", year = 2010, age_start = 5,
                  age_end = 9, mean_dmft = 1.5)
  expect_error(dmft_clean(d, cfg), "state")
})
