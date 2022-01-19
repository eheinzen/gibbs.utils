
k <- 1
n <- 100
tau <- 1
mu <- -5
p <- -6
p + (k + tau*mu - tau*p - n*exp(p)/(1 + exp(p)))/(tau + n*exp(p)/(1 + exp(p))^2)
(k + tau*mu + n*exp(p)*(p - 1 - exp(p))/(1 + exp(p))^2)/(tau + n*exp(p)/(1 + exp(p))^2)

qt_binom_approx(p, k, n, mu, tau)


LL <- function(x) k*x - n*log(1 + exp(x)) - 0.5*tau*(x - mu)^2

library(tidyverse)
f <- function() {
  approx <- qt_binom_approx(p, k = k, n = n, mean = mu, tau = tau)
  tibble(
    x = seq(-7, -3, by = 0.1),
    L = map_dbl(x, LL),
    L2 = map_dbl(x, dnorm, mean = approx[1], sd = approx[2], log = TRUE)
  ) %>%
    mutate(
      L = L - L[x == p],
      L2 = L2 - L2[x == p]
    )
}
blah <- f()
plot(y = blah$L, x = blah$x, type = 'l')
lines(y = blah$L2, x = blah$x, col = 2)


devtools::load_all()
out3 <- out2 <- out <- replace(numeric(1e4), TRUE, -6)
accept2 <- accept <- numeric(1e4)

for(i in 2:1e4) {
  tmp <- mh_binom_reg(p = out[i-1], k = 1, n = 100, mean = -5, precision = matrix(1), proposal = "quad")
  out[i] <- tmp
  accept[i] <- attr(tmp, "accept")
  tmp2 <- mh_binom_reg(p = out2[i-1], k = 1, n = 100, mean = -5, precision = 1, proposal = "normal")
  out2[i] <- tmp2
  accept2[i] <- attr(tmp2, "accept")
  out3[i] <- ss_binom_reg(p = out3[i-1], k = 1, n = 100, mean = -5, precision = 1)
}

tibble(
  quad = out,
  norm = out2,
  ss = out3
) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(x = value, color = name)) +
  geom_density()

plot(out, type = 'l')
table(accept)
table(accept2)



