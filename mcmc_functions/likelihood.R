library(mvtnorm)

### LOG LIKELIHOOD ###
logLik <- function(Sigma, beta) {
  beta.expanded <- beta[rep(1:nrow(beta), times = S), ]
  #dmvnorm(as.vector(Y), rep(1, S) %x% rowSums(X0 * beta), Sigma, log = TRUE)
  dmvnorm(as.vector(Y), rowSums(A * beta.expanded), Sigma, log = TRUE)
}

