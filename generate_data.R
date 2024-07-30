source("other_functions/spatial_data.R")
source("other_functions/bsplines_2_3D.R")
mySeed <- 3215877

# Sample sizes
n <- 200
nTest <- 25
S <- 10
STest <- 5
K <- 9

# True parameter values
trueSigf2 <- 10
trueThf <- 1
trueSigma2 <- seq(50, 100, length = K)
trueTheta <- runif(K, 0.01, 0.1)
trueTau2 <- 0.2
trueBeta <- 3

# Generate training data
set.seed(mySeed)
X <- matrix(1, ncol=1, nrow = n)
Z <- matrix(runif(2 * S, 0, 100), ncol = 2)
train <- spatialData(n = n, 
                     X = X,
                     Z = Z,
                     K = K,
                     sigf2 = trueSigf2,
                     thf = trueThf,
                     sigma2 = trueSigma2, 
                     theta = trueTheta, 
                     tau2 = trueTau2, 
                     beta = trueBeta,
                     range = c(0, 100))
save(train, file = "data/train.RData")

set.seed(mySeed)
indexTest <- sample(n, nTest)
U <- train$U[indexTest, ]

# Generate testing data
#set.seed(mySeed)
XTest <- matrix(1, ncol=1, nrow = nTest)
ZTest <- matrix(runif(2 * STest, 0, 100), ncol = 2)
test <- spatialData(n = nTest, 
                    X = XTest, 
                    Z = ZTest,
                    K = K,
                    U = U,
                    sigf2 = trueSigf2,
                    thf = trueThf,
                    sigma2 = trueSigma2, 
                    theta = trueTheta, 
                    tau2 = trueTau2, 
                    beta = trueBeta,
                    range = c(0, 100))
test$index <- indexTest
save(test, file = "data/test.RData")
save(trueTheta, file = "data/theta.RData")

trueSigma2
trueTheta
