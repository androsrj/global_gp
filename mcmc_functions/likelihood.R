library(mvtnorm)

### LOG LIKELIHOOD ###
logLik <- function(Sigma, mu) {
  dmvnorm(as.vector(Y), c(mu * J), Sigma, log = TRUE)
}

