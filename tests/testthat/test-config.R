test_that("dmft_config creates valid config object", {
  cfg <- dmft_config(
    regions = c("A", "B", "C"),
    region_col = "province",
    year_range = c(2000, 2020)
  )

  expect_s3_class(cfg, "dmft_config")
  expect_equal(cfg$n_regions, 3)
  expect_equal(cfg$year_start, 2000L)
  expect_equal(cfg$year_end, 2020L)
  expect_equal(cfg$n_years, 21L)
  expect_equal(cfg$backend, "inla")
  expect_equal(length(cfg$historical_years), 21)
})

test_that("dmft_config validates inputs", {
  expect_error(dmft_config(regions = "A", year_range = c(2000, 2020)))
  expect_error(dmft_config(regions = c("A", "B"), year_range = c(2020, 2000)))
})

test_that("dmft_config handles projections", {
  cfg <- dmft_config(
    regions = c("A", "B", "C"),
    year_range = c(2000, 2020),
    projection_range = c(2021, 2030)
  )

  expect_equal(cfg$projection_start, 2021L)
  expect_equal(cfg$projection_end, 2030L)
})

test_that("dmft_config parses age groups", {
  cfg <- dmft_config(
    regions = c("A", "B"),
    year_range = c(2000, 2020),
    age_groups_deciduous = c("0-4", "5-9"),
    age_groups_permanent = c("5-9", "10-14", "60+")
  )

  expect_equal(cfg$age_groups_deciduous_lower, c(0L, 5L))
  expect_equal(cfg$age_groups_deciduous_upper, c(4L, 9L))
  expect_equal(cfg$age_groups_permanent_lower, c(5L, 10L, 60L))
  expect_equal(cfg$age_groups_permanent_upper, c(9L, 14L, 100L))
})

test_that("print.dmft_config works", {
  cfg <- dmft_config(regions = c("A", "B"), year_range = c(2000, 2020))
  # cli_alert_info writes to stderr; capture both stdout and stderr
  out <- capture.output(print(cfg), type = "message")
  expect_true(any(grepl("DMFT analysis configuration", out)))
})
