library(mvtnorm)

### LOG LIKELIHOOD ###
logLik <- function(Sigma, beta) {
  dmvnorm(as.vector(newY), A %*% beta, Sigma, log = TRUE)
}

