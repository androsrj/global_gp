source("other_functions/spatial_data.R")
source("other_functions/bsplines_2_3D.R")
mySeed <- 321

# Sample sizes
n <- 500
nTest <- 100
S <- 10
STest <- 5
K <- 9

# True parameter values
trueSigma2 <- seq(1, 5, length = K)
trueTau2 <- 0.2
trueTheta <- runif(K, 1, 10)
trueBeta <- c(8, -2)

# Generate training data
#set.seed(mySeed)
X <- cbind(matrix(1, nrow = n), 
           runif(n, -5, 5))
Z <- matrix(sort(runif(2 * S, -20, 20)), ncol = 2)
train <- spatialData(n = n, 
                     X = X,
                     Z = Z,
                     K = K,
                     sigma2 = trueSigma2, 
                     tau2 = trueTau2, 
                     theta = trueTheta, 
                     beta = trueBeta)
save(train, file = "data/train.RData")

# Generate testing data
#set.seed(mySeed)
XTest <- cbind(matrix(1, nrow = nTest), 
               runif(nTest, -5, 5))
ZTest <- matrix(sort(runif(2 * STest, 20, 20)), ncol = 2)
test <- spatialData(n = nTest, 
                    X = XTest, 
                    Z = ZTest,
                    K = K,
                    sigma2 = trueSigma2, 
                    tau2 = trueTau2, 
                    theta = trueTheta,
                    beta = trueBeta)
save(test, file = "data/test.RData")
save(trueTheta, file = "data/theta.RData")

trueSigma2
trueTheta
