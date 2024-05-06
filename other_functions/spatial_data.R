library(fields)
library(splines2)

# Function to simulate spatial data
spatialData <- function(n, X, K,
                        sigma2, tau2, theta, beta,
                        S = NULL, range = c(0, 10), dims = 2, 
                        covariance = "exponential") {
  
  # Sample the locations and put them into an n-by-dims matrix
  # Unless the coordinates are pre-supplied
  if (is.null(S)) {
    locations <- runif(n * dims, range[1], range[2])
    S <- matrix(locations, nrow = n, ncol = dims)
  }
  
  # Order S by sum of coordinates
  S <- S[order(rowSums(S)), ]
  
  # Compute the covariance matrix (symmetric)
  if (covariance == "exponential") {
    D <- rdist(S)
  } else if (covariance == "exp_squared") {
    D <- rdist(S)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  #C <- sigma2 * exp(- theta * D)
  C <- lapply(1:K, function(k) {
    sigma2[k] * exp(- theta[k] * D)
  })
  
  # Sample h
  eta <- sapply(1:K, function(k) {
    t(rmvnorm(1, sigma = C[[k]]))
  })
  basis <- Bsplines_2D(X, df = c(sqrt(K), sqrt(K)))
  h <- rowSums(basis * eta)
  
  # Generate Y
  n <- nrow(X)
  #Y <- lapply(1:nSubj, \(i) X %*% beta + rep(Z[i], n) * gamma + W + eps)
  #X <- lapply(1:nSubj, \(i) cbind(rep(Z[i], n), X))
  Y <- X %*% beta + h + rnorm(n, 0, sqrt(tau2))
  
  # Return data
  return(list(X = X, Y = Y, h = as.vector(h), D = D, S = S))
}
