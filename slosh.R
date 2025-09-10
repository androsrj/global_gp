# SOURCES
#source("mcmc_functions/slosh_mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")
library(fields)
library(ggplot2)
library(Matrix)

load("data/slosh/flood_data.RData")

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

subsample <- coords$x < -74.82 & coords$x > -74.84 & coords$y < 39.07 & coords$y > 39.05
coords.subset <- coords[subsample, ]
out.subset <- out[ , subsample]
ggplot(data = coords, aes(x=x, y=y, col=elev_meters)) + 
  geom_point() + 
  geom_vline(aes(xintercept = -74.82)) + 
  geom_vline(aes(xintercept = -74.84)) + 
  geom_hline(aes(yintercept = 39.05)) + 
  geom_hline(aes(yintercept = 39.07)) 
ggplot(data = coords.subset, aes(x=x, y=y, fill=elev_meters)) + 
  geom_tile()

#source("coastlines.R")

set.seed(mySeed)
which.points <- sample(nrow(coords.subset), n + nTest)
train.index <- which.points[1:n]
test.index <- which.points[(n+1):(n+nTest)]

X <- cbind(coords.subset$elev_meters[train.index], 
           coords.subset$dist.east[train.index])
Z <- inputs[train.storms, which.Z]
Y <- matrix(c(t(as.matrix(out.subset[train.storms, train.index]))), ncol = 1)
U <- coords.subset[train.index, 1:2]
D <- fields::rdist(U)
XTest <- cbind(coords.subset$elev_meters[test.index], 
               coords.subset$dist.east[test.index])
ZTest <- inputs[test.storms, which.Z]
YTest <- matrix(c(t(as.matrix(out.subset[test.storms, test.index]))), ncol = 1)
UTest <- coords.subset[test.index, 1:2]
DTest <- fields::rdist(UTest)

flood.train <- list(X=X, Z=Z, Y=Y, U=U, D=D)
flood.test <- list(X=XTest, Z=ZTest, Y=YTest, U=UTest, D=DTest)
save(flood.train, flood.test, file = "data/slosh/flood_subset.RData")

K <- 9
#q <- ncol(X) + ncol(Z) + 1
q <- ncol(X) + 1
propSD <- list(sigb2 = seq(0.4, 0.6, length = q),
               thb = seq(0.3, 0.5, length = q),
               sigma2 = seq(0.4, 0.6, length = K),
               tau2 = 0.2,
               theta = seq(0.5, 0.8, length = K))
starting <- list(sigma2 = seq(0.01, 0.1, length = K),
                 theta = rep(0.5, K),
                 sigb2 = rep(1, q),
                 thb = rep(0.2, q), 
                 tau2 = 0.1,
                 beta = rep(0,7))

cat("Setup complete \n")
results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                starting = starting,
                propSD = propSD,
                nIter = 50, nBurn = 50, nThin = 2, nReport = 10,
                model = "full_gp")
saveRDS(results, file = "objects/slosh.RDS")

#theta
sd(YTest)
results$posteriorMeans
results$acceptance
nSamples <- length(results$paramSamples[[5]])

library(MBA)
library(fields)

STest <- nrow(ZTest)
rmse <- cvg <- width <- scores <- crps <- numeric(STest)
a <- .05
for (i in 1:STest) {
  truth <- YTest[(nTest*(i-1)+1):(i*nTest)]
  pred <- results$preds[2, (nTest*(i-1)+1):(i*nTest)]
  rmse[i] <- sqrt(mean((truth - pred)^2))
  lower <- pmax(0, results$preds[1, (nTest*(i-1)+1):(i*nTest)])
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