library(mvtnorm)

### LOG LIKELIHOOD ###
logLik <- function(Sigma, beta) {
  beta.expanded <- beta[rep(1:nrow(beta), each = S), ]
  dmvnorm(as.vector(Y), rowSums(A * beta.expanded), Sigma, log = TRUE)
}

