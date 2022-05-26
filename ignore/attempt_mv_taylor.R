
qt_mv_binom_approx <- function(around, k, n, mean, Q) {
  ep <- exp(around)
  ep1 <- 1 + ep
  ep2 <- ep1*ep1

  # -H(x)
  newQ <- Q + diag(n*ep / ep2, length(around))
  newQ.inv <- chol_inv(newQ)
  # x - H^-1(x) g(x)
  newmean <- around + newQ.inv %*% (k - Q %*% (around - mean) - n*ep/ep1)

  list(newmean, newQ, newQ.inv)
}

binom_LL_mv <- function(p, k, n, mean, Q) {
  sum(k * p) - sum(n * log(1 + exp(p))) - 0.5*as.vector(t(p - mean) %*% Q %*% (p - mean))
}

qt_mv_binom <- function(p, k, n, mean, Q) {

  prop_params <- qt_mv_binom_approx(p, k, n, mean, Q);
  proposal <- chol_mvrnorm(1, prop_params[[1]], Sigma = prop_params[[3]]);

  orig_params <- qt_mv_binom_approx(proposal, k, n, mean, Q);

  ratio <- binom_LL_mv(proposal, k, n, mean, Q)
  ratio <- ratio - dmvnorm(proposal, prop_params[[1]], Q = prop_params[[2]], log = TRUE)

  ratio <- ratio - binom_LL_mv(p, k, n, mean, Q)
  ratio <- ratio + dmvnorm(p, orig_params[[1]], Q = orig_params[[2]], log = TRUE)
  accept <- accept_reject(ratio)
  # print(exp(ratio))
  # list(if(accept) proposal else p, accept, exp(ratio), proposal,  prop_params, orig_params);
  if(accept) proposal else p
}

doit <- function(p, k, n, mean, Q) {
  t(vapply(seq_len(nrow(p)), function(i) qt_mv_binom(p[i, ], k[i, ], n[i, ], mean[i, ], Q), numeric(ncol(p))))
}

set.seed(88)
means <- matrix(0, nrow = 24*8*30, ncol = 181)
n <- matrix(rpois(181*24*8*30, 50), ncol = 181)
k <- matrix(rbinom(length(n), size = n, prob = 0.5), ncol = 181)
Q <- diag(1, 181)
p <- matrix(0.5, nrow = 24*8*30, ncol = 181)
qt_mv_binom(p[1, ], k = k[1, ], n = n[1, ], mean = means[1, ], Q = Q) %>% str

t1 <- Sys.time(); doit(p, k, n, means, Q); print(t2 <- Sys.time() - t1)
t1 <- Sys.time(); ss_binom_reg(p, k, n, means, Q); print(t2 <- Sys.time() - t1)

library(tidyverse)
tibble(x = seq(-3, 3, by = 0.01), y = map_dbl(x, ~ binom_LL_mv(.x, 100, 1000, 0, diag(10, 1)))) %>%
  mutate(y = y - y[x == p[1]]) %>%
  plot(type = 'l')
params <- qt_mv_binom(0.5, k = 100, n = 1000, mean = 0, Q = diag(10, 1))
tibble(x = seq(-3, 3, by = 0.01), y = -0.5*(x - params[[5]][[1]][1])^2 * params[[5]][[2]][1]) %>%
  mutate(y = y - y[x == p[1]]) %>%
  lines(col = 2)

profvis::profvis(doit(p[1:100, ], k = k[1:100, ], n = n[1:100, ], mean = means[1:100, ], Q = Q))





library(tidyverse)
pois_LL_mv <- function(L, k, mean, Q) {
  sum((k * L - exp(L))[!is.na(k)]) - 0.5*as.vector(t(L - mean) %*% Q %*% (L - mean))
}
tibble(x = seq(0, 3, by = 0.01), y = map_dbl(x, ~ pois_LL_mv(.x, 1, 2, diag(100, 1)))) %>%
  mutate(y = y - y[x == 0]) %>%
  plot(type = 'l')
params <- mvqt_pois_approx(0, k = 1, mean = 2, Q = diag(100, 1))
tibble(x = seq(0, 3, by = 0.01), y = -0.5*(x - params[[1]][1])^2 * params[[2]][1]) %>%
  mutate(y = y - y[x == 0]) %>%
  lines(col = 2)


k <- matrix(c(1, NA), nrow = 1, ncol = 2)
prec <- diag(100, 2)
out3 <- out2 <- out1 <- matrix(0, 1000, 2)
for(i in 1:999) {
  out1[i+1, ] <- mh_pois_reg(out1[i, , drop = FALSE], k, mean = 2, precision = prec, proposal = "quad")
  out3[i+1, ] <- mh_pois_reg(out3[i, , drop = FALSE], k, mean = 2, precision = prec, proposal = "mv q")
  out2[i+1, ] <- ss_pois_reg(out2[i, , drop = FALSE], k, mean = 2, precision = prec)
}
hist(out3[-(1:2), ])
qt_pois_approx(0, 1, 2, 100)










