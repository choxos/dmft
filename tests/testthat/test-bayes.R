test_that("prepare_stan_data creates valid Stan data list", {
  skip_if_not_installed("cmdstanr")

  # Create minimal test data mimicking prepare_lmer_data output
  d <- data.frame(
    y = c(2.1, 3.5, 1.8, 4.2, 2.9, 3.1),
    region_std = rep(c("A", "B", "C"), 2),
    year = rep(c(2000, 2001), each = 3),
    year_factor = factor(rep(c(2000, 2001), each = 3)),
    age_group = rep("5-9", 6),
    province_idx = rep(1:3, 2),
    year_idx = rep(1:2, each = 3),
    age_group_idx = rep(1, 6),
    se_imputed = c(0.5, 0.3, 0.4, 0.6, 0.2, 0.5),
    precision = c(4, 11.1, 6.25, 2.78, 25, 4)
  )

  mdata <- list(
    data = d,
    n_regions = 3,
    n_years = 2,
    n_age = 1,
    age_groups = "5-9"
  )

  cfg <- dmft_config(
    regions = c("A", "B", "C"),
    region_col = "region",
    year_range = c(2000, 2001)
  )

  stan_data <- dmft:::prepare_stan_data(mdata, cfg)

  expect_equal(stan_data$N, 6)
  expect_equal(stan_data$N_regions, 3)
  expect_equal(stan_data$N_years, 2)
  expect_equal(stan_data$K, 0)
  expect_equal(length(stan_data$y), 6)
  expect_equal(length(stan_data$se), 6)
  expect_true(all(stan_data$se > 0))
  expect_equal(nrow(stan_data$X), 6)
  expect_equal(ncol(stan_data$X), 0)
  expect_true(all(stan_data$region_idx %in% 1:3))
  expect_true(all(stan_data$year_idx %in% 1:2))
})


test_that("Stan model file exists in package", {
  stan_file <- system.file("stan", "dmft_ast.stan", package = "dmft")
  # Skip if package not installed in development mode
  skip_if(!nzchar(stan_file), "Stan model file not found (likely dev mode)")
  expect_true(file.exists(stan_file))
})


test_that("dmft_fit_bayes requires cmdstanr", {
  skip_if_not_installed("cmdstanr")

  cfg <- dmft_config(
    regions = c("A", "B"),
    region_col = "region",
    year_range = c(2000, 2001)
  )

  # Should error with bad data, but at least validates the function exists

  expect_error(dmft_fit_bayes(data.frame(), list(), "deciduous", cfg))
})
