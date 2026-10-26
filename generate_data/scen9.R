source("../other_functions/bsplines_2_3D.R")
library(splines2)
library(Matrix)

rmvn <- function(n, mu=0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V),p))))
    stop("Dimension problem!")
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

# Parameter values
K <- 9
p <- 2
q <- 2
B <- as.matrix(c(0, 0, 0))
sigma.sq <- seq(0.5, 1.0, length = p + 1)
tau.sq <- 0.2
theta <- 0.1

# Generate data
set.seed(1)
nTrain <- 100
nTest <- 25
n <- nTrain + nTest
U <- cbind(runif(n, 0, 1), runif(n, 0, 1))
X <- cbind(rep(1, n), matrix(runif(n*p, 0, 10), nrow = n, ncol = p))
A <- t(bdiag(as.list(as.data.frame(t(X)))))
D <- as.matrix(dist(U))
C <- exp(-theta*D) %x% diag(sigma.sq)
w <- rmvn(1, rep(0,(p+1)*n), C)
mu <- as.vector(X %*% B + A %*% w)
Y <- matrix(rnorm(n, mu, sqrt(tau.sq)), ncol = 1)


# Training data
S <- 10
Z <- matrix(runif(q * S, 0, 100), ncol = q)
BF <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
YTrain <- matrix(Y[1:nTrain,], ncol=1)
for (s in 2:S) {
  y.s <- matrix(rnorm(nTrain, mu[1:nTrain], sqrt(tau.sq)), ncol = 1)
  YTrain <- rbind(YTrain, y.s)
}
train <- list(X = X[1:nTrain,], Y = YTrain, D = D[1:nTrain, 1:nTrain], 
              U = U[1:nTrain,], Z = Z, BF = BF)

# Testing data
test.ind <- (nTrain+1):n
STest <- 10
ZTest <- matrix(runif(q * STest, 0, 100), ncol = q)
BFTest <- Bsplines_2D(ZTest, df = c(sqrt(K), sqrt(K)))
YTest <- matrix(Y[test.ind,], ncol=1)
for (s in 2:STest) {
  y.s.test <- matrix(rnorm(nTest, mu[test.ind], sqrt(tau.sq)), ncol = 1)
  YTest <- rbind(YTest, y.s.test)
}
test <- list(X = X[test.ind,], Y = YTest, D = D[test.ind, test.ind], 
             U = U[test.ind,], Z = ZTest, BF = BFTest)

# Save
saveRDS(train, file = "../data/scen9/train.RDS")
saveRDS(test, file = "../data/scen9/test.RDS")
