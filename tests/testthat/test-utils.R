test_that("compute_bym2_scaling_factor works for chain graph", {
  # 1-2-3-4 chain
  sf <- compute_bym2_scaling_factor(c(1, 2, 3), c(2, 3, 4), 4)
  expect_true(is.numeric(sf))
  expect_true(sf > 0)
})

test_that("compute_bym2_scaling_factor works for complete graph", {
  # Triangle: 1-2, 2-3, 1-3
  sf <- compute_bym2_scaling_factor(c(1, 2, 1), c(2, 3, 3), 3)
  expect_true(sf > 0)
})

test_that("compute_bym2_scaling_factor warns on empty graph", {
  # No edges: all isolated nodes
  expect_warning(sf <- compute_bym2_scaling_factor(integer(0), integer(0), 3))
  expect_equal(sf, 1.0)
})

test_that("theme_dmft returns a ggplot2 theme", {
  th <- theme_dmft()
  expect_s3_class(th, "theme")
})
