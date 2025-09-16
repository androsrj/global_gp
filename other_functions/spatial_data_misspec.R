library(fields)
library(Matrix)

# Function to simulate spatial data
spatialData <- function(n, X, Z, 
                        sigb2, thb, sigma2, theta, tau2, beta,
                        U = NULL, range = c(0, 10), dims = 2,
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
  CZ <- sigma2 * exp(-theta * DZ)
  
  # Sample h, which is a GP between both the locations U and global covariates Z
  #h <- t(rmvnorm(1, sigma = CZ))
  
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
  lon <- U[ , 1]
  lat <- U[ , 2]
  beta.surf <- cbind(
    lon - lat,
    lon + lat - 100,
    2 * lon - lat - 50
  )
  B <- Reduce("cbind", lapply(1:q, \(j) t(rmvnorm(1, sigma = CB[[j]])))) + beta.surf
  XB <- rep(1, S) %x% rowSums(X0 * B)
  
  # Final covariance matrix for Y
  Sigma <- CZ + tau2 * diag(n * S)
  
  # Generate Y
  Y <- t(rmvnorm(1, mean = XB, sigma = Sigma))
  
  # Return data
  return(list(X = X, Z = Z, Y = Y, B = B, D = D, U = U))
}
