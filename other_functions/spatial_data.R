library(fields)
library(splines2)
library(Matrix)

# Function to simulate spatial data
spatialData <- function(n, X, Z, K,
                        sigb2, thb, sigma2, theta, tau2, beta,
                        U = NULL, range = c(0, 10), dims = 2, eps=1e-6,
                        covariance = "exponential", intercept = TRUE) {
  
  # Sample the locations and put them into an n-by-dims matrix
  # Unless the coordinates are pre-supplied
  if (is.null(U)) {
    locations <- runif(n * dims, range[1], range[2])
    U <- matrix(locations, nrow = n, ncol = dims)
  }
  
  # Order U by sum of coordinates
  U <- U[order(rowSums(U)), ]
  
  # Compute the covariance matrix (symmetric)
  if (covariance == "exponential") {
    D <- rdist(U)
  } else if (covariance == "exp_squared") {
    D <- rdist(U)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  
  C <- lapply(1:K, function(k) {
    sigma2[k] * exp(-theta[k] * D)
  })
  
  # Covariance - h
  eta <- sapply(1:K, function(k) {
    t(rmvnorm(1, sigma = C[[k]]))
  })
  BF <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
  S <- nrow(Z)
  h <- c(sapply(1:S, \(i) rowSums(sapply(1:K, \(k) BF[i, k] * eta[ , k]))))
  #basis <- lapply(1:K, function(k) {
  #  Reduce("rbind", lapply(1:S, \(s) BF[s, k] * diag(n)))
  #})
  #B <- lapply(1:K, \(k) tcrossprod(basis[[k]] %*% C[[k]], basis[[k]]))
  B.eta <- lapply(1:S, \(s) Reduce("+", lapply(1:K, \(k) BF[s, k]^2 * C[[k]])))
  C.eta <- bdiag(B.eta)
  
  # Covariance - beta
  n <- nrow(X)
  if (intercept == TRUE) {
    X0 <- cbind(rep(1, n), X)
  } else {
    X0 <-  X
  }
  q <- ncol(X0)
  CB <- lapply(1:q, \(j) sigb2[j] * exp(-thb[j] * D))
  CXB <- Reduce("+", lapply(1:q, \(j) matrix(X0[ , j], nrow = n, ncol = n) * CB[[j]] *
                              matrix(X0[ , j], nrow = n, ncol = n, byrow = T)))
  B <- Reduce("cbind", lapply(1:q, \(j) t(rmvnorm(1, sigma = CB[[j]])))) + matrix(beta, nrow = n, ncol = q, byrow = TRUE)
  XB <- rep(1, S) %x% rowSums(X0 * B)
  
  # Final covariance matrix for Y
  Sigma <- diag(S) %x% CXB + C.eta + diag(rnorm(n * S, 0, sqrt(tau2)))
  #cat(min(eigen(Sigma)$values))
  
  # Generate Y
  Y <- t(rmvnorm(1, sigma = Sigma))
  
  # Return data
  return(list(X = X, Z = Z, Y = Y, B = B, h = as.vector(h), D = D, U = U, BF = BF))
}
