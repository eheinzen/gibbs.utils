


set.seed(20211203)
p <- runif(24*8*138)
n <- rpois(24*8*138, 50)
k <- rbinom(24*8*138, n, p)
dim(n) <- dim(k) <- dim(p) <- c(24*8, 138)
prec <- diag(1, 138)
n2 <- n3 <- n
k2 <- k3 <- k
k2[1:100, ] <- n2[1:100, ] <- 0
k3[] <- n3[] <- 0
microbenchmark::microbenchmark(
  none = {set.seed(99); ss_binom_reg(p, k, n, mean = 2, precision = prec)},
  none2 = {set.seed(99); ss_binom_reg2(p, k, n, mean = 2, precision = prec)},
  check = "equal",
  times = 10
)
microbenchmark::microbenchmark(
  some = {set.seed(99); ss_binom_reg(p, k2, n2, mean = 2, precision = prec)},
  some2 = {set.seed(99); ss_binom_reg2(p, k2, n2, mean = 2, precision = prec)},
  check = "equal",
  times = 10
)
microbenchmark::microbenchmark(
  all = {set.seed(99); ss_binom_reg(p, k3, n3, mean = 2, precision = prec)},
  all2 = {set.seed(99); ss_binom_reg2(p, k3, n3, mean = 2, precision = prec)},
  check = "equal",
  times = 10
)

microbenchmark::microbenchmark(
  none = {set.seed(99); mh_binom_reg(p, k, n, mean = 2, precision = prec)},
  none2 = {set.seed(99); mh_binom_reg2(p, k, n, mean = 2, precision = prec)},
  check = "equal",
  times = 10
)
