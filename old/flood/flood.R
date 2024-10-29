# SOURCES
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")

# Libraries
library(MBA)
library(fields)

# Data
load("data/flood_data.RData")

# Global covariates (Z)
STrain <- 10
STest <- 5
set.seed(123)
samp <- sample(1:nrow(inputs), STrain + STest)
stormsTrain <- samp[1:STrain]
stormsTest <- samp[(STrain+1):(STrain+STest)]
Z <- inputs[stormsTrain,3:4]
ZTest <- inputs[stormsTest,3:4]

# Local covariates (X)
nTest <- 20
#samp2 <- sample(1:nrow(coords))
#train <- samp2[1:n]
train <- which(coords$x > -74.86 & coords$x < -74.83 & coords$y > 39.15 & coords$y < 39.175)
#test <- samp2[(n+1):(n+nTest)]
n <- length(train)
test <- sample(train, nTest)
X <- matrix(c(rep(1, n), coords$elev_meters[train]), ncol=2)
XTest <- matrix(c(rep(1, nTest), coords$elev_meters[test]), ncol=2)

# Distance matrices (D)
U <- coords[train, 1:2]
UTest <- coords[test, 1:2]
D <- rdist(U)
DTest <- rdist(UTest)

# Response (Y)
Y <- matrix(as.matrix(out[stormsTrain, train]), ncol = 1)
YTest <- matrix(as.matrix(out[stormsTest, test]), ncol = 1)

# Other
K <- 9
propSD <- list(sigma2 = seq(0.05, 0.15, length = K),
               tau2 = 0.15)
#theta <- runif(9, 0.5, 3)
theta <- seq(10, 100, length = K)

results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                theta = theta,
                propSD = propSD,
                nIter = 500, nBurn = 200, nThin=2,
                model = "full_gp")

theta
results$posteriorMeans
results$acceptance
nSamples <- length(results$paramSamples[[3]])
plot(1:nSamples, results$paramSamples[[3]], type="l")
saveRDS(results, file = "objects/flood.RDS")

lwr <- min(c(YTest[1:(2*nTest)], results$preds[2,]))
upr <- max(c(YTest[1:(2*nTest)], results$preds[2,]))
lims <- c(lwr, upr)

pdf("figures/flood/subj1_true.pdf")
pred.surf <-  mba.surf(cbind(UTest, YTest[1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="True Surface, Storm 1", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/flood/subj1_global.pdf")
pred.surf <-  mba.surf(cbind(UTest, results$preds[2,1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="Global GP, Storm 1", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/flood/subj2_true.pdf")
pred.surf <-  mba.surf(cbind(UTest, YTest[(nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="True Surface, Storm 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/flood/subj2_global.pdf")
pred.surf <-  mba.surf(cbind(UTest, results$preds[2, (nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="Global GP, Storm 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

nTestSubj <- nrow(ZTest)
abs_error <- cvg <- width <- scores <- crps <- numeric(nTestSubj)
a <- .05
for (i in 1:nTestSubj) {
  truth <- YTest[(nTest*(i-1)+1):(i*nTest)]
  pred <- results$preds[2, (nTest*(i-1)+1):(i*nTest)]
  abs_error[i] <- mean(abs(truth - pred))
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

abs_error
cat(paste0("Mean absolute error: ", round(mean(abs_error), 3), "\n"))

cvg
cat(paste0("Mean coverage: ", round(mean(cvg), 3), "\n"))

width
cat(paste0("Mean width: ", round(mean(width), 3), "\n"))

scores
cat(paste0("Mean interval score: ", round(mean(scores), 3), "\n"))

crps
cat(paste0("Mean CRPS: ", round(mean(crps), 3), "\n"))

