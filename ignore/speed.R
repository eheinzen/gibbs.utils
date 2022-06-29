


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
microbenchmark::microbenchmark(
  mh_norm = {set.seed(99); mh_binom_reg(p, k, n, mean = 2, precision = prec)},
  mh_unif = {set.seed(99); ss_binom_reg(p, k, n, mean = 2, precision = prec, proposal = "uniform")},
  mh_qt = {set.seed(99); ss_binom_reg(p, k, n, mean = 2, precision = prec, proposal = "quadratic taylor")},
  ss = {set.seed(99); ss_binom_reg(p, k, n, mean = 2, precision = prec)},
  times = 100
)







microbenchmark::microbenchmark(
  none = {set.seed(99); ss_pois_reg(p, n, mean = 2, precision = prec)},
  none2 = {set.seed(99); ss_pois_reg2(p, n, mean = 2, precision = prec)},
  check = "equal",
  times = 10
)












tmp <- cbind(0, matrix(c(0, 0, 0.7, 0.7), nrow = 2))
p <- array(rep(tmp, each = 400), dim = c(400, 2, 3))
mean <- round(p)
set.seed(99)
n <- rpois(400*2, 1200)
n[1:100] <- 0
k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
dim(k) <- c(400, 2, 3)
k[, 2, 2:3] <- 0

z <- matrix(1, nrow = 2, ncol = 3)
z[2, 2:3] <- 0
Q <- array(0, c(2, 2, 3))
for(j in 1:3) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
mean <- mean[, , -1]
Q <- Q[, , -1]
microbenchmark::microbenchmark(
  old = {set.seed(99); ss_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q)},
  new = {set.seed(99); ss_multinom_reg2(p = p, z = z, k = k, mean = mean, precision = Q)},
  # check = "equal",
  times = 100
)

set.seed(99)
Nr <- 10
Nc <- 10
mu <- matrix(0, Nr, Nc)
y <- matrix(rnorm(Nr*Nc), Nr, Nc)
Q <- rWishart(1, 200, diag(1000, Nc))[, , 1]
microbenchmark::microbenchmark(
  dmvnorm(y, mu, Q, log = TRUE),
  dmvnorm(y, mu, Q, log = TRUE, use_trace = TRUE),
  dmvnorm2(y, mu, Q, log = TRUE),
  times = 1000,
  check = "equal"
)




set.seed(99)
Nr <- 100
Nc <- 180
mu <- matrix(0, Nr, Nc)
y <- matrix(rnorm(Nr*Nc), Nr, Nc)
V <- rWishart(1, 200, diag(1000, Nc))[, , 1]
U <- diag(1, Nr)
microbenchmark::microbenchmark(
  dmatnorm(y, mu, V = V, U = U, log = TRUE),
  dmatnorm2(y, mu, V = V, U = NULL, log = TRUE),
  dmatnorm2(y, mu, V = V, U = U, log = TRUE),
  times = 100,
  check = "equal"
)



