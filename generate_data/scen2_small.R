source("../other_functions/spatial_data.R")
source("../other_functions/bsplines_2_3D.R")
mySeed <- 4528
library(mvtnorm)

# Sample sizes
# Can have a "small" dataset with n = 100 and nTest = 25
# Then a "large" dataset with n = 500 and nTest = 100
n <- 100
nTest <- 25

# Number of subjects - can probably leave these alone
S <- 10
STest <- 10

# Number of BFE's and predictors - leave these alone
K <- 9
p <- 2

### True parameter values ###
# Need to play around with these

# Covariance parameters for beta
trueSigb2 <- seq(1, 3, length = p + 1)
trueThb <- seq(0.1, 0.2, length = p + 1)

# Covariance parameters for global covariates (each length K)
trueSigma2 <- seq(5, 10, length = K)
trueTheta <- seq(0.1, 0.5, length = K)

# Error variance
trueTau2 <- 0.2

##########################
# Generate training data #
set.seed(mySeed)
X <- matrix(runif(n*p, 0, 10), nrow = n, ncol = p)
Z <- matrix(runif(2 * S, 0, 100), ncol = 2)
train <- spatialData(n = n, 
                     X = X,
                     Z = Z,
                     K = K,
                     sigb2 = trueSigb2,
                     thb = trueThb,
                     sigma2 = trueSigma2, 
                     theta = trueTheta, 
                     tau2 = trueTau2, 
                     range = c(0, 100))
saveRDS(train, file = "../data/small/scen2/train.RDS")

set.seed(mySeed)
indexTest <- sample(n, nTest)
U <- train$U[indexTest, ]

# Generate testing data
set.seed(mySeed)
XTest <- matrix(runif(nTest*p, 0, 10), nrow = nTest, ncol = p)
ZTest <- matrix(runif(2 * STest, 0, 100), ncol = 2)
test <- spatialData(n = nTest, 
                    X = XTest, 
                    Z = ZTest,
                    K = K,
                    U = U,
                    sigb2 = trueSigb2,
                    thb = trueThb,
                    sigma2 = trueSigma2, 
                    theta = trueTheta, 
                    tau2 = trueTau2, 
                    range = c(0, 100))
sd(test$Y)
test$index <- indexTest
saveRDS(test, file = "../data/small/scen2/test.RDS")

