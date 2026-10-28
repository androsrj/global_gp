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
B <- as.matrix(c(1, 3, 5))
sigma.sq <- seq(0.5, 1.0, length = p + 1)
tau.sq <- 0.2
theta <- 0.1

# Training data
set.seed(1)
n <- 100
S <- 10
Z <- matrix(runif(q * S, 0, 100), ncol = q)
BF <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
U <- cbind(runif(n, 0, 100), runif(n, 0, 100))
X <- cbind(rep(1, n), matrix(runif(n*p, 0, 10), nrow = n, ncol = p))
X.aug <- cbind(rep(1, S) %x% X, Z %x% rep(1, n))
A <- t(bdiag(as.list(as.data.frame(t(X)))))
D <- as.matrix(dist(U))
C <- exp(-theta*D) %x% diag(sigma.sq)
w <- rmvn(1, rep(0,(p+1)*n), C)
mu <- as.vector(X %*% B + A %*% w)
Y <- matrix(rnorm(n, mu, sqrt(tau.sq)), ncol = 1)
for (s in 2:S) {
  y.s <- matrix(rnorm(n, mu, sqrt(tau.sq)), ncol = 1)
  Y <- rbind(Y, y.s)
}
train <- list(X = X, Y = Y, D = D, U = U, Z = Z, BF = BF)

# Testing data
set.seed(1)
nTest <- 25
STest <- 10
ZTest <- matrix(runif(q * STest, 0, 100), ncol = q)
BFTest <- Bsplines_2D(ZTest, df = c(sqrt(K), sqrt(K)))
UTest <- cbind(runif(nTest, 0, 100), runif(nTest, 0, 100))
XTest <- cbind(rep(1, nTest), matrix(runif(nTest*p, 0, 10), nrow = nTest, ncol = p))
ATest <- t(bdiag(as.list(as.data.frame(t(XTest)))))
DTest <- as.matrix(dist(UTest))
CTest <- exp(-theta*DTest) %x% diag(sigma.sq)
wTest <- rmvn(1, rep(0,(p+1)*nTest), CTest)
muTest <- as.vector(XTest %*% B + ATest %*% wTest)
YTest <- matrix(rnorm(nTest, muTest, sqrt(tau.sq)), ncol = 1)
for (s in 2:STest) {
  y.s.test <- matrix(rnorm(nTest, muTest, sqrt(tau.sq)), ncol = 1)
  YTest <- rbind(YTest, y.s.test)
}
test <- list(X = XTest, Y = YTest, D = DTest, U = UTest, Z = ZTest, BF = BFTest)

# Save
saveRDS(train, file = "../data/scen9/train.RDS")
saveRDS(test, file = "../data/scen9/test.RDS")
