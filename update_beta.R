library(Matrix)
library(fields)
library(mvtnorm)
source("other_functions/bsplines_2_3D.R")

obj <- readRDS("objects/small_scen1.RDS")[[1]]
scen <- 1
dir <- paste0("data/small/scen", scen, "/")
train <- readRDS(paste0(dir, "train.RDS"))
test <- readRDS(paste0(dir, "test.RDS"))
n.iter <- ncol(obj$paramSamples$sigma2)
#n.iter <- 25
K <- 9

# Response variable, Y
Y <- train$Y
n <- 100

# Predictors (local and global)
X <- cbind(rep(1, n), train$X)
p <- ncol(X)
Z <- train$Z
S <- nrow(Z)

# Distance matrix
D <- fields::rdist(train$U)

# Permutation matrix
mat <- matrix(c(rep(1:p, each = n), rep(1:n, times = p)), ncol = 2)
new.order <- order(mat[, 2], mat[, 1])
P <- diag(n*p)[new.order, ]

beta <- matrix(0, nrow = n*p, ncol = n.iter)
for (i in 1:n.iter) {
  
  # Parameter values
  sigma2 <- obj$paramSamples$sigma2[ , i]
  theta <- obj$paramSamples$theta[ , i]
  sigb2 <- obj$paramSamples$sigb2[ , i]
  thb <- obj$paramSamples$thb[ , i]
  tau2 <- obj$paramSamples$tau2[i]
  
  # Sample eta and calculate h
  BF <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
  C <- lapply(1:K, function(k) {
    sigma2[k] * exp(-theta[k] * D)
  })
  C.array <- simplify2array(C)  # array of dim n x n x K
  C.array <- aperm(C.array, c(3, 1, 2))  # now K x n x n
  W.list <- lapply(1:K, function(k) tcrossprod(BF[, k]))  # each S x S
  W.array <- simplify2array(W.list)  # S x S x K
  C.eta <- Reduce('+', lapply(1:K, function(k) kronecker(W.array[, , k], C.array[k, , ])))
  h <- c(mvtnorm::rmvnorm(1, sigma = C.eta))
  
  # Diagonal version of X
  X2 <- matrix(0, n, n * p)
  for (j in seq_len(n)) {
    cols <- ((j - 1) * p + 1):(j * p)  # slot for this row
    X2[j, cols] <- X[j, ]
  }
  
  # Big B
  Sigma.p.inv <- lapply(1:p, \(j) solve(sigb2[j] * exp(-thb[j] * D)))
  Sigma.inv <- as.matrix(bdiag(Sigma.p.inv))
  big.B <- solve(S * crossprod(X2 %*% P) / tau2 + Sigma.inv)
  
  # Little b
  little.b <- Reduce("+", lapply(1:S, function(s) {
    ind <- ((s-1)*n+1):(s*n)
    Y2 <- (Y - h)[ind, ]
    crossprod(X2 %*% P, Y2) / tau2
  }))
  
  # Sample beta
  beta[ , i] <- mvtnorm::rmvnorm(1, mean = big.B %*% little.b, sigma = big.B)
  
}

beta.est <- rowMeans(beta)







