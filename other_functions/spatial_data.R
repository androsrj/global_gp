library(fields)
library(splines2)

# Function to simulate spatial data
spatialData <- function(n, X, Z, K,
                        sigf2, thf, sigma2, theta, tau2, beta,
                        U = NULL, range = c(0, 10), dims = 2, eps=1e-6,
                        covariance = "exponential") {
  
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
    DX <- rdist(X)
  } else if (covariance == "exp_squared") {
    D <- rdist(U)^2
    DX <- rdist(X)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  #C <- sigma2 * exp(- theta * D)
  C <- lapply(1:K, function(k) {
    sigma2[k] * exp(-theta[k] * D)
  })
  CX <- sigf2 * exp(-thf * DX)
  
  # Sample h
  eta <- sapply(1:K, function(k) {
    t(rmvnorm(1, sigma = C[[k]]))
  })
  basis <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
  S <- nrow(Z)
  h <- c(sapply(1:S, \(i) rowSums(sapply(1:K, \(k) basis[i, k] * eta[ , k]))))
  
  # Sample f
  f <- t(rmvnorm(1, sigma = matrix(1, S, S) %x% CX + diag(eps, S*n)))
  
  # Columns of ones combined with X matrix
  A <- rep(1, S) %x% cbind(matrix(1, nrow = n), X)
  
  # Generate Y
  Y <- A %*% beta + f + h + rnorm(n * S, 0, sqrt(tau2))
  
  # Return data
  return(list(X = X, Z = Z, Y = Y, h = as.vector(h), D = D, U = U))
}
