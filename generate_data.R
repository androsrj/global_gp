source("other_functions/spatial_data.R")
source("other_functions/bsplines_2_3D.R")
mySeed <- 200

# Sample sizes
n <- 500
nTest <- 10
nSubj <- 10
nTestSubj <- 5
K <- 9

# True parameter values
trueSigma2 <- 1:K
trueTau2 <- 0.2
trueTheta <- runif(K, 1, 10)
trueBeta <- c(2, -1)
p <- length(trueBeta)

# Generate training data
set.seed(mySeed)
X <- matrix(rnorm(n * length(trueBeta)), nrow = n, ncol = length(trueBeta))
X <- X[order(X[ , 1]), ]
train <- spatialData(n = n, 
                     X = X, 
                     K = K,
                     sigma2 = trueSigma2, 
                     tau2 = trueTau2, 
                     theta = trueTheta, 
                     beta = trueBeta)
save(train, file = "data/train.RData")

# Generate testing data
set.seed(mySeed)
X <- matrix(rnorm(nTest * length(trueBeta)), nrow = nTest, ncol = length(trueBeta))
X <- X[order(X[ , 1]), ]
test <- spatialData(n = nTest, 
                    X = X, 
                    K = K,
                    sigma2 = trueSigma2, 
                    tau2 = trueTau2, 
                    theta = trueTheta, 
                    beta = trueBeta)
save(test, file = "data/test.RData")
