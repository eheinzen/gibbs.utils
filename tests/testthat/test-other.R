test_that("dmvnorm", {
  expect_equal(
    dmvnorm(1:5, 0:4, diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:5, 0:4, sqrt(100:104), log = TRUE))
  )
})

test_that("dmatnorm", {
  expect_equal(
    dmatnorm(matrix(1:10, nrow = 5), matrix(0:9, nrow = 5), V = diag(1, 2), U = diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:10, 0:9, sqrt(c(100:104, 100:104)), log = TRUE))
  )
  expect_equal(
    dmatnorm(matrix(1:10, nrow = 5), matrix(0:9, nrow = 5), V = diag(c(1, 1/2)), U = diag(1/(100:104)), log = TRUE),
    sum(dnorm(1:10, 0:9, sqrt(c(100:104, 2*(100:104))), log = TRUE))
  )
})
