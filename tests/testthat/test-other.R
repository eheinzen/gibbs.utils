test_that("dmvnorm", {
  expect_equal(
    dmvnorm(1:5, 0:4, diag(1/(100:104)), use_trace = FALSE),
    dmvnorm(1:5, 0:4, diag(1/(100:104)), use_trace = TRUE)
  )
})
