library(fields)

# Function to simulate spatial data
spatialData <- function(n, X, Z, 
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
  
  # Temporary matrix that combines U and Z (for sampling h later)
  S <- nrow(Z)
  UZ <- cbind(U[rep(1:n, times = S), ],
              Z[rep(1:S, each = n), ])
  
  # Compute the covariance matrix (symmetric)
  if (covariance == "exponential") {
    D <- rdist(U)
    DZ <- rdist(UZ)
  } else if (covariance == "exp_squared") {
    D <- rdist(U)^2
    DZ <- rdist(UZ)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  #C <- sigma2 * exp(- theta * D)
  CZ <- sigma2 * exp(-theta * DZ)
  
  # Sample h, which is a GP between both the locations U and global covariates Z
  h <- t(rmvnorm(1, sigma = CZ))
  
  # Sample beta
  n <- nrow(X)
  if (intercept == TRUE) {
    X0 <- cbind(rep(1, n), X)
  } else {
    X0 <-  X
  }
  q <- ncol(X0)
  DB <- lapply(1:q, \(j) matrix(X0[ , j], nrow = n, ncol = n) * 
                 (sigb2[j] * exp(-thb[j] * D)) *
                 matrix(X0[ , j], nrow = n, ncol = n, byrow = T))
  CB <- Reduce("+", DB)
  B <- t(rmvnorm(q, sigma = CB)) + matrix(beta, nrow = n, ncol = q, byrow = TRUE)
  XB <- rep(1, S) %x% (X0 %*% beta)
  
  # Generate Y
  Y <- XB + h + rnorm(n * S, 0, sqrt(tau2))
  
  # Return data
  return(list(X = X, Z = Z, Y = Y, B = B, h = as.vector(h), D = D, U = U))
}
