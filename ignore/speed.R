


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





set.seed(99)
Nr <- 100
Nc <- 22
y <- matrix(rnorm(Nr*Nc), Nr, Nc)
Q <- diag(1, Nc)
Q0 <- diag(1, Nc)
microbenchmark::microbenchmark(
  conj_mvnorm_mu2(y = y, Q = Q, mu0 = NULL, Q0 = Q0, use = "chol"),
  conj_mvnorm_mu2(y = y, Q = Q, mu0 = NULL, Q0 = Q0, use = "spam"),
  conj_mvnorm_mu2(y = y, Q = Q, mu0 = NULL, Q0 = Q0, use = "MASS"),
  times = 1000
)



set.seed(99)
n <- 1000
p <- 154
one <- matrix(rnorm(n*p), n, p)
Q0 <- diag(1, p)
ten <- diag(10, p)
X <- matrix(sample(0:1, n*p, replace = TRUE), n, p)
zero <- diag(1, p)
old <- options(gu_ignore_diag_case = TRUE)
microbenchmark::microbenchmark(
  newer = conj_matlm_beta(y = one, X = X, V = ten, U = NULL, mu0 = NULL, Q0 = Q0, zero = zero, params.only = TRUE),
  older = conj_diagmatlm_beta(y = one, X = X, V = ten, U = NULL, mu0 = NULL, Q0 = Q0, params.only = TRUE),
  check = "equal", times = 10
)
options(old)





set.seed(88)
Q <- rWishart(1, 360, diag(1, 300))[, , 1]

x <- matrix(rnorm(1000*300), 1000, 300)
y <- matrix(rnorm(1000*300), 1000, 300)
mu <- matrix(rnorm(1000*300), 1000, 300)
expect_equal(
  dmvnorm(x, mu, Q, log = TRUE) - dmvnorm(y, mu, Q, log = TRUE),
  dmvnorm_diff(x, y, mu, Q, log = TRUE)
)
microbenchmark::microbenchmark(
  dmvnorm(x, mu, Q, log = TRUE) - dmvnorm(y, mu, Q, log = TRUE),
  dmvnorm_diff(x, y, mu, Q, log = TRUE)
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
  new3 = {set.seed(99); ss_multinom_reg3(p = p, z = z, k = k, mean = mean, precision = Q)},
  new4 = {set.seed(99); ss_multinom_reg3(p = p, z = array(rep(z, each = 400), dim = c(400, 2, 3)), k = k, mean = mean, precision = Q)},
  check = "equal",
  times = 100
)








set.seed(20210119)
L <- matrix(rnorm(50*151), nrow = 50, ncol = 151)
k <- matrix(rpois(50*151, exp(L)), nrow = 50, ncol = 151)
prec <- diag(rgamma(151, shape = 0.1, rate = 0.1))

microbenchmark::microbenchmark(
  five = {set.seed(99); sample_pois_reg(L, k, mean = 0, precision = prec, method = "quad")},
  six = {set.seed(99); sample_pois_reg2(L, k, mean = 0, precision = prec, method = "quad")},
  times = 10,
  check = "equal"
)





set.seed(20221202)
N <- 650000
M <- 5
K <- 60
Y <- matrix(rnorm(N*M), N)
mean <- matrix(rnorm(N*M), N)
sigma <- diag(1:5)
X <- matrix(rbinom(N*K, 1, prob = 0.1), N)
Xt <- t(X)
XtX <- Xt %*% X
mu0 <- matrix(0, K, M)
prec0 <- diag(10, K*M)

microbenchmark::microbenchmark(
  five = {set.seed(99); conj_matlm_beta(Y, V = sigma, mu0 = mu0, Q0 = prec0, XtU = Xt, XtUX = XtX, diag = TRUE)},
  six = {set.seed(99); conj_matlm_beta2(Y, V = sigma, mu0 = mu0, Q0 = prec0, XtU = Xt, XtUX = XtX, diag = TRUE)},
  times = 5,
  check = "equal"
)

mu0 <- matrix(0, 10*K, M)
prec0 <- diag(10, 10*K*M)
microbenchmark::microbenchmark(
  five = {set.seed(99); conj_matnorm_mu(Y[seq_len(10*K), ], V = sigma, mu0 = mu0, Q0 = prec0, diag = TRUE)},
  six = {set.seed(99); conj_matnorm_mu2(Y[seq_len(10*K), ], V = sigma, mu0 = mu0, Q0 = prec0, diag = TRUE)},
  times = 100,
  check = "equal"
)






R <- 100
I <- 22
J <- 23
mean <- p <- array(0, dim = c(R, I, J))
set.seed(99)
k <- array(rpois(R*I*J, 10), dim = c(R, I, J))
z <- array(1, dim = c(R, I, J))
Q <- array(0, c(I, I, J))
for(j in 1:J) Q[, , j] <- diag(2, I)
mean <- mean[, , -1]
Q <- Q[, , -1]

microbenchmark::microbenchmark(
  regular = sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q, diag = FALSE),
  none = sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q),
  diagtrue = sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q, diag = TRUE),
  times = 5
)





R <- 100
I <- 22
J <- 23
mean <- p <- array(0, dim = c(R, I, J))
set.seed(121)
k <- array(rpois(R*I*J, 10), dim = c(R, I, J))
z <- array(1, dim = c(R, I, J))
Q <- array(0, c(I, I, J))
for(j in 1:J) Q[, , j] <- diag(2, I)
mean <- mean[, , -1]
Q <- Q[, , -1]

microbenchmark::microbenchmark(
  regular = {set.seed(124); sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q)},
  two = {set.seed(124); sample_multinom_reg2(p = p, z = z, k = k, mean = mean, precision = Q)},
  three = {set.seed(124); sample_multinom_reg3(p = p, z = z, k = k, mean = mean, precision = Q)},
  four = {set.seed(124); sample_multinom_reg4(p = p, z = z, k = k, mean = mean, precision = Q)},
  five = {set.seed(124); sample_multinom_reg5(p = p, z = z, k = k, mean = mean, precision = Q)},
  check = "equal",
  times = 5
)






set.seed(20230124)
L <- matrix(rnorm(1000*151), nrow = 50, ncol = 151)
k <- matrix(rpois(1000*151, exp(L)), nrow = 50, ncol = 151)
prec <- diag(rgamma(151, shape = 0.1, rate = 0.1))

microbenchmark::microbenchmark(
  five = {set.seed(99); sample_pois_reg(L, k, mean = 0, precision = prec, method = "slice")},
  six = {set.seed(99); sample_pois_reg2(L, k, mean = 0, precision = prec, method = "slice")},
  times = 100,
  check = "equal"
)
microbenchmark::microbenchmark(
  five = {set.seed(99); sample_pois_reg(L, k, mean = 0, precision = prec, method = "slice")},
  six = {set.seed(99); sample_pois_reg2(L, k, mean = 0, precision = prec, method = "slice")},
  times = 10,
  check = "equal"
)

set.seed(20230124)
p <- runif(24*8*138)
n <- rpois(24*8*138, 50)
k <- rbinom(24*8*138, n, p)
dim(n) <- dim(k) <- dim(p) <- c(24*8, 138)
prec <- diag(1, 138)
microbenchmark::microbenchmark(
  five = {set.seed(99); sample_binom_reg(p, k, n, mean = 2, precision = prec, method = "quad")},
  six = {set.seed(99); sample_binom_reg2(p, k, n, mean = 2, precision = prec, method = "quad")},
  check = "equal",
  times = 10
)


R <- 100
I <- 22
J <- 23
mean <- p <- array(0, dim = c(R, I, J))
set.seed(20230124)
k <- array(rpois(R*I*J, 10), dim = c(R, I, J))
z <- array(1, dim = c(R, I, J))
Q <- array(0, c(I, I, J))
for(j in 1:J) Q[, , j] <- diag(2, I)
mean <- mean[, , -1]
Q <- Q[, , -1]

microbenchmark::microbenchmark(
  regular = {set.seed(124); sample_multinom_reg(p = p, z = z, k = k, mean = mean, precision = Q)},
  two = {set.seed(124); sample_multinom_reg2(p = p, z = z, k = k, mean = mean, precision = Q)},
  check = "equal",
  times = 5
)
