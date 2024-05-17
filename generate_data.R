source("other_functions/spatial_data.R")
source("other_functions/bsplines_2_3D.R")
mySeed <- 321

# Sample sizes
n <- 100
nTest <- 10
S <- 5
STest <- 5
K <- 9

# True parameter values
trueSigma2 <- seq(1, 3, length = K)
trueTau2 <- 0.2
trueTheta <- runif(K, 1, 10)
trueMu <- 8

# Generate training data
set.seed(mySeed)
X <- matrix(sort(runif(2 * S)), ncol = 2)
train <- spatialData(n = n, 
                     X = X,
                     K = K,
                     sigma2 = trueSigma2, 
                     tau2 = trueTau2, 
                     theta = trueTheta, 
                     mu = trueMu)
save(train, file = "data/train.RData")

# Generate testing data
set.seed(mySeed)
X <- matrix(sort(runif(2 * STest)), ncol = 2)
test <- spatialData(n = nTest, 
                    X = X, 
                    K = K,
                    sigma2 = trueSigma2, 
                    tau2 = trueTau2, 
                    theta = trueTheta,
                    mu = trueMu)
save(test, file = "data/test.RData")
save(trueTheta, file = "data/theta.RData")

trueSigma2
trueTheta
