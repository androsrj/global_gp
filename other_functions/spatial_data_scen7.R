library(fields)

# Function to simulate spatial data
spatialData <- function(n, X, Z, 
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
  
  # Temporary matrix that combines U and Z (for sampling h later)
  S <- nrow(Z)
  UZ <- cbind(U[rep(1:n, times = S), ],
              Z[rep(1:S, each = n), ])
  
  # Compute the covariance matrix (symmetric)
  if (covariance == "exponential") {
    DZ <- rdist(UZ)
    DX <- rdist(X)
  } else if (covariance == "exp_squared") {
    DZ <- rdist(UZ)^2
    DX <- rdist(X)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  #C <- sigma2 * exp(- theta * D)
  CZ <- sigma2 * exp(-theta * DZ)
  CX <- sigf2 * exp(-thf * DX)
  
  # Sample h, which is a GP between both the locations U and global covariates Z
  h <- t(rmvnorm(1, sigma = CZ))
  
  # Sample f
  f <- t(rmvnorm(1, sigma = matrix(1, S, S) %x% CX + diag(eps, S*n)))
  
  # Columns of ones combined with X matrix
  A <- rep(1, S) %x% cbind(matrix(1, nrow = n), X)
  
  # Generate Y
  Y <- A %*% beta + f + h + rnorm(n * S, 0, sqrt(tau2))
  
  # Return data
  return(list(X = X, Z = Z, Y = Y, h = as.vector(h), D = D, U = U))
}
