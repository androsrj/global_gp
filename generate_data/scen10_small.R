source("../other_functions/spatial_data_misspec.R")
mySeed <- 45213

# Sample sizes
# Can have a "small" dataset with n = 100 and nTest = 25
# Then a "large" dataset with n = 500 and nTest = 100
n <- 100
nTest <- 25

# Number of subjects - can probably leave these alone
S <- 10
STest <- 10

# Number of predictors
p <- 2

### True parameter values ###
# Need to play around with these

# Covariance parameters for local covariates
trueSigf2 <- 5
trueThf <- 1

# Covariance parameters for global covariates (each length K)
trueSigma2 <- 5
trueTheta <- 1

# Error variance
trueTau2 <- 0.2

# Regression coefficients
trueBeta <- c(1, 0.5, -1)



##########################
# Generate training data #
set.seed(mySeed)
X <- matrix(runif(n*p, 0, 10), nrow = n, ncol = p)
Z <- matrix(runif(2 * S, 0, 100), ncol = 2)
train <- spatialData(n = n, 
                     X = X,
                     Z = Z,
                     sigf2 = trueSigf2,
                     thf = trueThf,
                     sigma2 = trueSigma2, 
                     theta = trueTheta, 
                     tau2 = trueTau2, 
                     beta = trueBeta,
                     range = c(0, 100))
save(train, file = "../data/small/scen10/train.RData")

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
                    U = U,
                    sigf2 = trueSigf2,
                    thf = trueThf,
                    sigma2 = trueSigma2, 
                    theta = trueTheta, 
                    tau2 = trueTau2, 
                    beta = trueBeta,
                    range = c(0, 100))
test$index <- indexTest
save(test, file = "../data/small/scen10/test.RData")

