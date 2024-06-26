source("other_functions/spatial_data.R")
source("other_functions/bsplines_2_3D.R")
mySeed <- 32158

# Sample sizes
n <- 1000
nTest <- 100
S <- 10
STest <- 10
K <- 9

# True parameter values
trueSigma2 <- seq(50, 100, length = K)
trueTau2 <- 0.2
trueTheta <- runif(K, 0.01, 0.1)
trueBeta <- c(3, -2)

# Generate training data
#set.seed(mySeed)
X <- cbind(matrix(1, nrow = n), 
           runif(n, -5, 5))
Z <- matrix(runif(2 * S, 0, 100), ncol = 2)
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
ZTest <- matrix(runif(2 * STest, 0, 100), ncol = 2)
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
