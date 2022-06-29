test_that("dmvnorm", {
  expect_equal(
    dmvnorm(1:5, 0:4, diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:5, 0:4, sqrt(100:104), log = TRUE))
  )
})
