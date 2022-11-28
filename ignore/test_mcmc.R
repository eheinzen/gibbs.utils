mh_gamma <- function(curr, k, mean, tau, mult = 1) {
  lmean <- exp(mean + 0.5/tau)
  lvar <- (exp(1/tau) - 1)*lmean^2

  alpha <- (lmean^2 / lvar + k)*mult
  beta <- (lmean / lvar + 1)*mult

  proposal <- rgamma(1, shape = alpha, rate = beta)

  LL.prop <- dpois(k, proposal, log = TRUE) + dlnorm(proposal, meanlog = mean, sdlog = 1/sqrt(tau), log = TRUE)
  q.prop <- dgamma(proposal, shape = alpha, rate = beta, log = TRUE)

  LL.curr <- dpois(k, exp(curr), log = TRUE) + dlnorm(exp(curr), meanlog = mean, sdlog = 1/sqrt(tau), log = TRUE)
  q.curr <- dgamma(exp(curr), shape = alpha, rate = beta, log = TRUE)

  ratio <- (LL.prop - q.prop) - (LL.curr - q.curr)
  accept <- log(runif(1)) <= ratio
  if(accept) log(proposal) else curr
}

y <- matrix(1)
m <- matrix(0)
Q <- diag(1, 1)
out <- out2 <- out3 <- out4 <- out5 <- out6 <- out7 <- matrix(0, 10000, 1)
for(i in seq_len(10000-1)) {
  out[i+1, ] <- sample_pois_reg(out[i, , drop = FALSE], y, m, Q, method = "slice")
  out2[i+1, ] <- sample_pois_reg(out2[i, , drop = FALSE], y, m, Q, method = "mv ind quad")
  out3[i+1, ] <- mh_gamma(out3[i, ], y, m, as.vector(Q), mult = 1)
  out4[i+1, ] <- sample_pois_reg(out4[i, , drop = FALSE], y, m, Q, method = "mv quadratic")
  out5[i+1, ] <- sample_pois_reg(out5[i, , drop = FALSE], y, m, Q, method = "normal")
  out6[i+1, ] <- sample_pois_reg(out6[i, , drop = FALSE], y, m, as.vector(Q), method = "gamma")
  out7[i+1, ] <- sample_pois_reg(out7[i, , drop = FALSE], y, m, Q, method = "mv gamma")
}

library(tidyverse)
abind::abind(out, out2, out3, out4, out5, out6, out7, along = 0L) %>%
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
out <- out2 <- out4 <- out5 <- out7 <- matrix(0, 10000, 1)
for(i in seq_len(10000-1)) {
  out[i+1, ] <-  sample_binom_reg(out[i, , drop = FALSE],  k, n, m, Q, method = "slice")
  out2[i+1, ] <- sample_binom_reg(out2[i, , drop = FALSE], k, n, m, Q, method = "mv ind quad")
  out4[i+1, ] <- sample_binom_reg(out4[i, , drop = FALSE], k, n, m, Q, method = "mv quadratic")
  out5[i+1, ] <- sample_binom_reg(out5[i, , drop = FALSE], k, n, m, Q, method = "normal")
  out7[i+1, ] <- sample_binom_reg(out7[i, , drop = FALSE], k, n, m, Q, method = "mv beta")
}

abind::abind(out, out2, out4, out5, out7, along = 0L) %>%
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
out <- out2 <- array(0, c(10000, dim(p)))
for(i in seq_len(10000-1)) {
  out[i+1, , , ] <-  sample_multinom_reg(p = out[i, , , ], z = z, k = k, mean = mean, precision = Q)
  out2[i+1, , , ] <-  sample_multinom_reg2(p = out2[i, , , ], z = z, k = k, mean = mean, precision = Q)
}

stopifnot(out[, , , 1] == 0, out2[, , , 1] == 0)

library(tidyverse)
abind::abind(out, out2, along = 0L) %>%
  reshape2::melt(c("method", "iter", "r", "i", "j")) %>%
  filter(iter > 1000, j != 1) %>%
  ggplot(aes(x = value, color = factor(method))) +
  facet_wrap(~ r + i + j, scales = "free") +
  geom_density()
