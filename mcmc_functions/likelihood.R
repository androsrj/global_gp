library(mvtnorm)

### LOG LIKELIHOOD ###
logLik <- function(Sigma, beta) {
  dmvnorm(as.vector(Y), A * beta, Sigma, log = TRUE)
}

