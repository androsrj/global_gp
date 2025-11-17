library(fields)
library(splines2)
library(Matrix)

# Function to simulate spatial data
spatialData <- function(n, X, Z, K,
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
  
  # Compute the covariance matrix (symmetric)
  if (covariance == "exponential") {
    D <- rdist(U)
  } else if (covariance == "exp_squared") {
    D <- rdist(U)^2
  } else {
    stop("Covariance function must be either exponential or exp_squared.")
  }
  
  C <- lapply(1:K, function(k) {
    sigma2[k] * (1 + sqrt(3) * theta[k] * D) * exp(-sqrt(3) * theta[k] * D)
  })
  
  # Covariance - eta
  BF <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
  S <- nrow(Z)
  C.array <- simplify2array(C)  # array of dim n x n x K
  C.array <- aperm(C.array, c(3, 1, 2))  # now K x n x n
  W.list <- lapply(1:K, function(k) tcrossprod(BF[, k]))  # each S x S
  W.array <- simplify2array(W.list)  # S x S x K
  C.eta <- Reduce('+', lapply(1:K, function(k) kronecker(W.array[, , k], C.array[k, , ])))
  
  # Covariance - beta
  n <- nrow(X)
  if (intercept == TRUE) {
    X0 <- cbind(rep(1, n), X)
  } else {
    X0 <-  X
  }
  q <- ncol(X0)
  CB <- lapply(1:q, \(j) sigb2[j] * (1 + sqrt(3) * thb[j] * D) * exp(-sqrt(3) * thb[j] * D))
  CXB <- Reduce("+", lapply(1:q, \(j) matrix(X0[ , j], nrow = n, ncol = n) * CB[[j]] *
                              matrix(X0[ , j], nrow = n, ncol = n, byrow = T)))
  lon <- U[ , 1]
  lat <- U[ , 2]
  #beta.list <- vector("list", length = q)
  #for (j in 1:q) {
  #  beta.list[[j]] <- t(mvtnorm::rmvnorm(1, sigma = CB[[j]]))
  #}
  #beta.surf <- Reduce("cbind", beta.list)
  B <- Reduce("cbind", lapply(1:q, \(j) t(rmvnorm(1, sigma = CB[[j]]))))
  XB <- rep(1, S) %x% rowSums(X0 * B)
  
  # Final covariance matrix for Y
  #Sigma <- matrix(1, S, S) %x% CXB + C.eta + tau2 * diag(n * S)
  Sigma <- C.eta + tau2 * diag(n * S)
  
  # Generate Y
  Y <- t(rmvnorm(1, mean = XB, sigma = Sigma))
  
  # Return data
  return(list(X = X, Z = Z, Y = Y, B = B, D = D, U = U, BF = BF))
}
