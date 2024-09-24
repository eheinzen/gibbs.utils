test_that("bss `from:to` works", {
  expect_equal(bss(1:7, n_basis = 6L)[1:3, ], bss(1:3, n_basis = 6L, from = 1L, to = 7L))
})
