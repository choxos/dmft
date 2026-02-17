test_that("theme_dmft returns a ggplot2 theme", {
  th <- theme_dmft()
  expect_s3_class(th, "theme")
})

test_that("check_suggests errors for missing package", {
  expect_error(check_suggests("nonexistent_pkg_xyz123"))
})
