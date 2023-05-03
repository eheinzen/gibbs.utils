library(tidyverse)
y <- matrix(1)
m <- matrix(0)
Q <- diag(1, 1)
out <- out2 <- out3 <- out4 <- out5 <- out6 <- out7 <- matrix(0, 10000, 1)
for(i in seq_len(10000-1)) {
  out[i+1, ] <- sample_pois_reg(out[i, , drop = FALSE], y, m, Q, method = "slice")
  out2[i+1, ] <- sample_pois_reg(out2[i, , drop = FALSE], y, m, Q, method = "mv ind quad")
  out3[i+1, ] <- sample_pois_reg(out3[i, , drop = FALSE], y, m, Q, method = "mv trunc")
  out4[i+1, ] <- sample_pois_reg(out4[i, , drop = FALSE], y, m, Q, method = "mv quadratic")
  out5[i+1, ] <- sample_pois_reg(out5[i, , drop = FALSE], y, m, Q, method = "normal")
  out6[i+1, ] <- sample_pois_reg(out6[i, , drop = FALSE], y, m, as.vector(Q), method = "gamma")
  out7[i+1, ] <- sample_pois_reg(out7[i, , drop = FALSE], y, m, Q, method = "mv gamma")
}

abind::abind(out, out2, out4, out5, out6, out7, along = 0L) %>%
  reshape2::melt(c("method", "i", "x")) %>%
  filter(i > 1000) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ x, scales = "free") +
  geom_density()



# -------------------------------------------------------------------------


y <- matrix(3)
m <- matrix(0)
Q <- diag(1, 1)
out0 <- out <- out2 <- out3 <- out4 <- out5 <- out6 <- matrix(0, 10000, 1)
for(i in seq_len(10000-1)) {
  out0[i+1, ] <- sample_pois_reg(out0[i, , drop = FALSE], y, m, Q, method = "slice")
  out[i+1, ] <- sample_pois_reg(out[i, , drop = FALSE], y, m, Q, method = "slice", truncate = list(at = 0, allow = "above"))
  out2[i+1, ] <- sample_pois_reg(out2[i, , drop = FALSE], y, m, Q, method = "normal", truncate = list(at = 0, allow = "above"))
  out3[i+1, ] <- sample_pois_reg(out3[i, , drop = FALSE], y, m, Q, method = "unif", truncate = list(at = 0, allow = "above"))
  out4[i+1, ] <- sample_pois_reg(out4[i, , drop = FALSE], y, m, Q, method = "slice", truncate = list(at = 3, allow = "below"))
  out5[i+1, ] <- sample_pois_reg(out5[i, , drop = FALSE], y, m, Q, method = "normal", truncate = list(at = 3, allow = "below"))
  out6[i+1, ] <- sample_pois_reg(out6[i, , drop = FALSE], y, m, Q, method = "unif", truncate = list(at = 3, allow = "below"))
}

abind::abind(out0, out, out2, out3, along = 0L) %>%
  reshape2::melt(c("method", "i", "x")) %>%
  filter(i > 1000) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ x, scales = "free") +
  geom_density()


abind::abind(out0, out4, out5, out6, along = 0L) %>%
  reshape2::melt(c("method", "i", "x")) %>%
  filter(i > 1000) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ x, scales = "free") +
  geom_density()








# -------------------------------------------------------------------------



k <- matrix(1)
n <- matrix(20)
m <- matrix(-3)
Q <- diag(1, 1)
out <- out2 <- out3 <- out4 <- out5 <- out7 <- matrix(-4, 10000, 1)
for(i in seq_len(10000-1)) {
  out[i+1, ] <-  sample_binom_reg(out[i, , drop = FALSE],  k, n, m, Q, method = "slice")
  out2[i+1, ] <- sample_binom_reg(out2[i, , drop = FALSE], k, n, m, Q, method = "mv ind quad")
  out3[i+1, ] <- sample_binom_reg(out3[i, , drop = FALSE], k, n, m, Q, method = "mv trunc")
  out4[i+1, ] <- sample_binom_reg(out4[i, , drop = FALSE], k, n, m, Q, method = "mv quadratic")
  out5[i+1, ] <- sample_binom_reg(out5[i, , drop = FALSE], k, n, m, Q, method = "normal")
  out7[i+1, ] <- sample_binom_reg(out7[i, , drop = FALSE], k, n, m, Q, method = "mv beta")
}

abind::abind(out, out2, out3, out4, out5, out7, along = 0L) %>%
  reshape2::melt(c("method", "i", "x")) %>%
  filter(i > 1000) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ x, scales = "free") +
  geom_density()


# -------------------------------------------------------------------------

tmp <- cbind(0, matrix(c(0, 0, 0.7, 0.7), nrow = 2))
p <- array(rep(tmp, each = 4), dim = c(4, 2, 3))
mean <- round(p)
set.seed(99)
n <- rpois(4*2, 1200)
k <- t(sapply(n, function(x) rmultinom(1, size = x, prob = c(0.25, 0.25, 0.5))))
dim(k) <- c(4, 2, 3)
stopifnot(rowSums(k, dims = 2) == n)

z <- matrix(1, nrow = 2, ncol = 3)
Q <- array(0, c(2, 2, 3))
for(j in 1:3) Q[, , j] <- matrix(c(2, 1, 1, 2), 2)
mean <- mean[, , -1]
Q <- Q[, , -1]

k[4, , ] <- 0

out <- out2 <- out3 <- out4 <- out5 <- array(0, c(10000, dim(p)))
for(i in seq_len(10000-1)) {
  out[i+1, , , ] <-  sample_multinom_reg(p = out[i, , , ], z = z, k = k, mean = mean, precision = Q, method = "slice")
  out2[i+1, , , ] <-  sample_multinom_reg(p = out2[i, , , ], z = z, k = k, mean = mean, precision = Q, method = "unif")
  out3[i+1, , , ] <-  sample_multinom_reg(p = out3[i, , , ], z = z, k = k, mean = mean, precision = Q, method = "norm")
  out4[i+1, , , ] <-  sample_multinom_reg(p = out4[i, , , ], z = z, k = k, mean = mean, precision = Q, method = "quad")
  out5[i+1, , , ] <-  sample_multinom_reg(p = out5[i, , , ], z = z, k = k, mean = mean, precision = Q, method = "mv trunc", width = 0.01)
}

stopifnot(out[, , , 1] == 0, out2[, , , 1] == 0, out3[, , , 1] == 0, out4[, , , 1] == 0, out5[, , , 1] == 0)

library(tidyverse)
abind::abind(out, out2, out3, out4, out5, along = 0L) %>%
  reshape2::melt(c("method", "iter", "r", "i", "j")) %>%
  filter(iter > 9000, j != 1) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ r + i + j, scales = "free") +
  geom_density()













# -------------------------------------------------------------------------

set.seed(99)
beta <- matrix(c(1, 0, 1, 0, 3, 3), nrow = 2)
V <- 10*chol_inv(matrix(c(1, 0.95, 0.1, 0.95, 1, 0.25, 0.1, 0.25, 1), nrow = 3))
stopifnot(isSymmetric(V))
X <- scale(matrix(rnorm(2*100000), ncol = 2))
z <- +(beta != 0)
err <- spam::rmvnorm.prec(100000, Q = V)
err <- err - matrix(colMeans(err), nrow = 100000, ncol = 3, byrow = TRUE)
y <- X %*% beta + err
Q0 <- diag(0.1, sum(z))
Q1 <- diag(ifelse(as.vector(z) == 1, 0.1, Inf), 6)


out <- out2 <- out3 <- out4 <- array(0, c(100000, prod(dim(beta))))
for(i in seq_len(100000-1)) {
  out[i+1, as.vector(z == 1)] <- conj_matlm_beta(zero = z, y = y, X = X, V = V, U = NULL, Q0 = Q0, mu0 = NULL)
  out2[i+1, ] <- conj_matlm_beta(y = y, X = X, V = V, U = NULL, Q0 = Q1, mu0 = NULL)
  out3[i+1, as.vector(z == 1)] <- gs_matlm_beta(beta = out3[i, as.vector(z == 1)], zero = z, y = y, X = X, V = V, U = NULL, Q0 = Q0, mu0 = NULL)
  out4[i+1, ] <- gs_matlm_beta(beta = out4[i, ], y = y, X = X, V = V, U = NULL, Q0 = Q1, mu0 = NULL)
}

stopifnot(out[, c(2, 4)] == 0, out2[, c(2, 4)] == 0, out3[, c(2, 4)] == 0, out4[, c(2, 4)] == 0)

library(tidyverse)
abind::abind(out, out2, out3, out4, along = 0L) %>%
  reshape2::melt(c("method", "iter", "i")) %>%
  filter(iter > 9000, i != 2, i != 4) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ i, scales = "free") +
  geom_density()






