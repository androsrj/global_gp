# SOURCES
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")
library(fields)

nReps <- 10
size <- "small"
scen <- "scen10"
dir <- paste0("data/", size, "/", scen, "/")
load(paste0(dir, "train.RData"))
load(paste0(dir, "test.RData"))
n <- nrow(train$X)
nTest <- nrow(test$X)
X <- train$X
Z <- train$Z
Y <- train$Y
U <- train$U
D <- train$D
XTest <- test$X
ZTest <- test$Z
YTest <- test$Y
UTest <- test$U
DTest <- test$D
K <- 9
propSD <- list(sigf2 = 0.6,
               thf = 0.7,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(5, 10, length = K),
                 theta = rep(5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))

cat("Setup complete \n")
results <- vector("list", length = nReps)
for (i in 1:nReps) {
  results[[i]] <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                       starting = starting,
                       propSD = propSD,
                       nIter = 2000, nBurn = 1000, nThin=2,
                       model = "full_gp")
}

path <- paste0("objects/", size, "_", scen, ".RDS") 
saveRDS(results, file = path)

cat("Mean of Y: \n")
mean(train$Y)
cat("SD of Y: \n")
sd(train$Y)
cat("Sigma2: \n")
apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigma2), 1, mean)
cat("Theta: \n")
apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$theta), 1, mean)
cat("Sigma2_f: \n")
mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigf2))
cat("Theta_f: \n")
mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$thf))
cat("Beta: \n")
apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$beta), 1, mean)
cat("Tau2: \n")
mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$tau2))


nSamples <- length(results$paramSamples[[5]])
#plot(1:nSamples, results$paramSamples[[5]], type="l")
#saveRDS(results, file = "objects/global.RDS")

library(MBA)
library(fields)


results <- results[[1]]
STest <- nrow(test$Z)
rmse <- cvg <- width <- scores <- crps <- numeric(STest)
a <- .05
for (i in 1:STest) {
  truth <- YTest[(nTest*(i-1)+1):(i*nTest)]
  pred <- results$preds[2, (nTest*(i-1)+1):(i*nTest)]
  rmse[i] <- sqrt(mean((truth - pred)^2))
  lower <- results$preds[1, (nTest*(i-1)+1):(i*nTest)]
  upper <- results$preds[3, (nTest*(i-1)+1):(i*nTest)]
  cvg[i] <- mean(lower < truth & upper > truth)
  width[i] <- mean(upper - lower)
  scores[i] <- mean((upper - lower) + 
                      2/a * (lower - truth) * (truth < lower) + 
                      2/a * (truth - upper) * (truth > upper))
  predSamples <- t(results$predSamples[(nTest*(i-1)+1):(i*nTest), ])
  crps[i] <- mean(energy_score(truth, predSamples))
}

rmse
cat(paste0("Root MS error: ", round(mean(rmse), 3), "\n"))

cvg
cat(paste0("Mean coverage: ", round(mean(cvg), 3), "\n"))

width
cat(paste0("Mean width: ", round(mean(width), 3), "\n"))

scores
cat(paste0("Mean interval score: ", round(mean(scores), 3), "\n"))

crps
cat(paste0("Mean CRPS: ", round(mean(crps), 3), "\n"))
