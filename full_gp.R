# SOURCES
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")

load("data/train.RData")
load("data/test.RData")
n <- nrow(train$X)
nTest <- nrow(test$X)
X <- train$X
Y <- train$Y
U <- train$U
D <- train$D
XTest <- test$X
YTest <- test$Y
UTest <- test$U
DTest <- test$D
K <- 9
propSD <- list(sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.05)

results <- mcmc(X = X, Y = Y, D = D,
                K = K,
                theta = runif(9, 0.5, 3),
                propSD = propSD,
                nIter = 2000, nBurn = 100, nThin=2,
                model = "full_gp")
results$posteriorMeans
results$acceptance
nSamples <- length(results$paramSamples[[3]])
plot(1:nSamples, results$paramSamples[[3]], type="l")

