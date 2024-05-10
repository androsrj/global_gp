library(fields)
library(splines2)

# Function to simulate spatial data
spatialData <- function(n, X, K,
                        sigma2, tau2, theta, mu,
                        U = NULL, range = c(0, 10), dims = 2, 
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
  } else if (covariance == "exp_squared") {
    D <- rdist(U)^2
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
  S <- nrow(X)
  h <- c(sapply(1:S, \(i) rowSums(sapply(1:K, \(k) basis[i, k] * eta[ , k]))))
  
  # Generate Y
  Y <- mu * matrix(1, nrow = n * S) + h + rnorm(n * S, 0, sqrt(tau2))
  
  # Return data
  return(list(X = X, Y = Y, h = as.vector(h), D = D, U = U))
}
