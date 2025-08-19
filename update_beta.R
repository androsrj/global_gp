library(Matrix)
library(fields)
library(mvtnorm)
library(MBA)
source("other_functions/bsplines_2_3D.R")

obj <- readRDS("objects/small_scen1.RDS")[[1]]
scen <- 1
dir <- paste0("data/small/scen", scen, "/")
train <- readRDS(paste0(dir, "train.RDS"))
test <- readRDS(paste0(dir, "test.RDS"))
#n.iter <- ncol(obj$paramSamples$sigma2)
n.iter <- 100
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

# Mean of beta
beta.mean <- c(5, 2, -4)
beta.mean <- rep(0, 3)

beta <- matrix(0, nrow = n*p, ncol = n.iter)
y.preds <- matrix(0, nrow = n, ncol = n.iter)
for (i in 1:n.iter) {
  if (i %% 10 == 0) {
    cat(paste0("Beginning iteration ", i, ".\n"))
  }
  
  # Parameter values
  sigma2 <- obj$paramSamples$sigma2[ , i]
  theta <- obj$paramSamples$theta[ , i]
  sigb2 <- obj$paramSamples$sigb2[ , i]
  thb <- obj$paramSamples$thb[ , i]
  tau2 <- obj$paramSamples$tau2[i]
  #sigb2 <- seq(0.5, 1, length = p)
  #thb <- seq(0.05, 0.1, length = p)
  #sigma2 <- seq(5, 10, length = K)
  #theta <- seq(0.1, 0.5, length = K)
  #tau2 <- 0.2
  
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
    cols <- ((j - 1) * p + 1):(j * p)
    X2[j, cols] <- X[j, ]
  }
  
  # Big B
  Sigma.p.inv <- lapply(1:p, \(j) solve(sigb2[j] * exp(-thb[j] * D)))
  Sigma.inv <- as.matrix(bdiag(Sigma.p.inv))
  big.B <- solve(S * crossprod(X2 %*% P) / tau2 + Sigma.inv)
  
  # Little b
  little.b <- Reduce("+", lapply(1:S, function(s) {
    ind <- ((s-1)*n+1):(s*n)
    #Y2 <- (Y - h)[ind, ]
    Y2 <- Y[ind, ]
    crossprod(X2 %*% P, Y2) / tau2
  }))
  
  mu <- big.B %*% little.b
  
  # Sample beta
  beta[ , i] <- mvtnorm::rmvnorm(1, mean = mu, sigma = big.B)
  current.beta <- matrix(beta[ , i], ncol = 3)
  
  # Sample from posterior predictive.
  y.mean <- rowSums(X * current.beta)
  y.var <- tau2 * diag(n)
  y.preds[ , i] <- c(mvtnorm::rmvnorm(1, mean = y.mean, sigma = y.var))
}

beta.est <- rowMeans(beta)
beta.mat <- matrix(beta.est, ncol = 3)
colMeans(beta.mat)


for (j in 1:3) {
  x <- train$U[ , 1]
  y <- train$U[ , 2]
  
  # Plot 1: Estimated beta surface
  z <- beta.mat[ , j]
  data <- cbind(x, y, z)
  surf <- mba.surf(data, n = 100, no.X = 500, no.Y = 500)
  grid <- surf$xyz.est
  #fields::image.plot(grid$x, grid$y, grid$z,
  #                   main = paste0("Estimated Beta", j))
  #contour(grid$x, grid$y, grid$z, add = TRUE, col = "black")
  df <- data.frame(data)
  p1 <- ggplot(data = df, aes(x, y, colour = z)) + geom_point(size = 5)
  print(p1)
  
  # Plot 2: Actual beta surface
  z <- train$B[ , j]
  data <- cbind(x, y, z)
  surf <- mba.surf(data, n = 100, no.X = 500, no.Y = 500)
  grid <- surf$xyz.est
  #fields::image.plot(grid$x, grid$y, grid$z, 
  #                   main = paste0("Actual Beta", j))
  #contour(grid$x, grid$y, grid$z, add = TRUE, col = "black")
  df <- data.frame(data)
  p2 <- ggplot(data = df, aes(x, y, colour = z)) + geom_point(size = 5)
  print(p2)
}

preds <- rowMeans(y.preds)

# RMSE
sqrt(mean((Y[1:100,] - preds)^2))

# Compare to SD
sd(Y[1:100,])

# Width
lower <- apply(y.preds, 1, quantile, .025)
upper <- apply(y.preds, 1, quantile, .975)
mean(upper - lower)

# Coverage
mean(Y[1:100,] < upper & Y[1:100,] > lower)

