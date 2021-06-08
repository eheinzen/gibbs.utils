
check_one_or_all <- function(x, len) {
  nm <- deparse(substitute(x))
  if(length(x) == 1) rep_len(x, len) else if(length(x) == len) x else stop(nm, " must be of length 1 or ", len)
}

my_inv <- function(x) {
  chol2inv(chol(x))
}


logit <- function(p) log(p) - log(1 - p)
expit <- function(x) 1 / (1 + exp(-x))

propose_01 <- function(p, scale, adjust) {
  logitp <- logit(p)
  sd <- sqrt(abs(logitp)+1)*scale*(0.8^adjust)
  sd <- ifelse(sd < 0.01, 0.01, sd)
  expit(stats::rnorm(length(logitp), logitp, sd))
}

propose_01_logdensity <- function(p, proposal, scale, adjust) {
  logitp <- logit(p)
  sd <- sqrt(abs(logitp)+1)*scale*(0.8^adjust)
  sd <- ifelse(sd < 0.01, 0.01, sd)
  dlogitnorm(proposal, mu = logitp, sd = sd, log = TRUE)
}
