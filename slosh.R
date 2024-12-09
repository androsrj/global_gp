# SOURCES
source("mcmc_functions/slosh_mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")
library(fields)
library(ggplot2)

load("data/flood_data.RData")

mySeed <- 1234
which.Z <- c(1:5)
n <- 100
nTest <- 25
S <- 10
STest <- 10

set.seed(mySeed)
which.storms <- sample(4000, S + STest)
train.storms <- which.storms[1:S]
test.storms <- which.storms[(S+1):(S+STest)]

coords.subset <- coords[coords$x < -74.81 & coords$x > -74.83 & coords$y < 39.08 & coords$y > 39.06, ]
ggplot(data = coords, aes(x=x, y=y, col=elev_meters)) + 
  geom_point() + 
  geom_vline(aes(xintercept = -74.81)) + 
  geom_vline(aes(xintercept = -74.83)) + 
  geom_hline(aes(yintercept = 39.06)) + 
  geom_hline(aes(yintercept = 39.08)) 
ggplot(data = coords.subset, aes(x=x, y=y, fill=elev_meters)) + 
  geom_tile()

source("coastlines.R")

set.seed(mySeed)
which.points <- sample(nrow(coords.subset), n + nTest)
train.index <- which.points[1:n]
test.index <- which.points[(n+1):(n+nTest)]

X <- cbind(coords.subset$elev_meters[train.index], 
           coords.subset$dist.east[train.index])
Z <- inputs[train.storms, which.Z]
Y <- matrix(c(t(as.matrix(out[train.storms, train.index]))), ncol = 1)
U <- coords[train.index, 1:2]
D <- fields::rdist(U)
XTest <- cbind(coords.subset$elev_meters[test.index], 
               coords.subset$dist.east[test.index])
ZTest <- inputs[test.storms, which.Z]
YTest <- matrix(c(t(as.matrix(out[test.storms, test.index]))), ncol = 1)
UTest <- coords[test.index, 1:2]
DTest <- fields::rdist(UTest)


K <- 9
propSD <- list(sigf2 = 0.6,
               thf = 20,
               sigma2 = seq(0.1, 0.25, length = K),
               tau2 = 0.4,
               theta = seq(0.1, 0.3, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 6,
                 thf = 5, 
                 tau2 = 0.1,
                 beta = rep(0,8))

cat("Setup complete \n")
results<- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
               starting = starting,
               propSD = propSD,
               nIter = 3000, nBurn = 2000, nThin=2,
               model = "full_gp")
saveRDS(results, file = "objects/slosh.RDS")

#theta
sd(YTest)
results$posteriorMeans
results$acceptance
nSamples <- length(results$paramSamples[[5]])
plot(1:nSamples, results$paramSamples[[5]], type="l")
saveRDS(results, file = "objects/global.RDS")

library(MBA)
library(fields)

STest <- nrow(ZTest)
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
