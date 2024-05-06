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
S <- train$S
D <- train$D
XTest <- test$X
YTest <- test$Y
STest <- test$S
DTest <- test$D
propSD <- list(sigma2 = seq(0.05, 0.12, length = K),
               tau2 = 0.07)
K <- 9

results <- mcmc(X = X, Y = Y, D = D, S = S,
                K = K,
                theta = runif(9, 0.5, 3),
                propSD = propSD,
                nIter = 300, nBurn = 10, nThin=1,
                model = "full_gp",
                transform = FALSE)
results$posteriorMeans
results$acceptance
plot(1:290, results$paramSamples[[1]][1,])

