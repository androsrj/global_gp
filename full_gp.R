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
load("data/theta.RData")
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
propSD <- list(sigma2 = seq(0.05, 0.15, length = K),
               tau2 = 0.15,
               theta = seq(0.1, 0.25, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 10,
                 thf = 1, 
                 tau2 = 0.1,
                 beta = 2)

results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                starting = starting,
                propSD = propSD,
                nIter = 100, nBurn = 100, nThin=2,
                model = "full_gp")

#theta
results$posteriorMeans
results$acceptance
nSamples <- length(results$paramSamples[[5]])
plot(1:nSamples, results$paramSamples[[5]], type="l")
saveRDS(results, file = "objects/global.RDS")

library(MBA)
library(fields)

lims <- c(-15, 15)

pdf("figures/subj1_true.pdf")
pred.surf <-  mba.surf(cbind(UTest, YTest[1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="True Surface, Subject 1", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/subj1_global.pdf")
pred.surf <-  mba.surf(cbind(UTest, results$preds[2,1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="Global GP, Subject 1", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/subj2_true.pdf")
pred.surf <-  mba.surf(cbind(UTest, YTest[(nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="True Surface, Subject 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/subj2_global.pdf")
pred.surf <-  mba.surf(cbind(UTest, results$preds[2, (nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="Global GP, Subject 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

STest <- nrow(test$Z)
abs_error <- cvg <- width <- scores <- crps <- numeric(STest)
a <- .05
for (i in 1:STest) {
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
