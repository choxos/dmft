test_that("dmft_phantom generates data for all regions", {
  cfg <- dmft_config(
    regions = c("A", "B", "C"),
    region_col = "province",
    year_range = c(2010, 2015)
  )

  phantom <- dmft_phantom(cfg, dentition = "permanent", missing_prob = 0)

  expect_true(nrow(phantom) > 0)
  expect_true(all(c("province", "year", "age_start", "age_end",
                     "mean_DMFT") %in% names(phantom)))
  expect_true(all(unique(phantom$province) %in% cfg$regions))
  expect_true(all(phantom$year >= 2010 & phantom$year <= 2015))
})

test_that("dmft_phantom generates deciduous data", {
  cfg <- dmft_config(
    regions = c("X", "Y"),
    region_col = "province",
    year_range = c(2000, 2005)
  )

  phantom <- dmft_phantom(cfg, dentition = "deciduous", missing_prob = 0)
  expect_true("mean_dmft" %in% names(phantom))
  expect_true(all(phantom$mean_dmft >= 0))
})

test_that("dmft_phantom respects seed for reproducibility", {
  cfg <- dmft_config(regions = c("A", "B"), year_range = c(2000, 2005),
                     seed = 42)
  p1 <- dmft_phantom(cfg, dentition = "permanent")
  p2 <- dmft_phantom(cfg, dentition = "permanent")
  expect_equal(p1, p2)
})
