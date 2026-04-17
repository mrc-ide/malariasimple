test_that("get_peak_cc is working as expected", {
  g0 <- 0.28
  g <- c(-0.3, -0.03, 0.17)
  h <- c(-0.35, 0.33, -0.08)
  expect_equal(get_peak_cc(g0,g,h), 233)
})
