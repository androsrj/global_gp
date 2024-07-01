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
#Z <- train$Z
Z <- matrix(1, ncol = 2, nrow = nrow(train$Z))
Y <- train$Y
U <- train$U
D <- train$D
XTest <- test$X
#ZTest <- test$Z
ZTest <- matrix(1, ncol = 2, nrow = nrow(test$Z))
YTest <- test$Y
UTest <- test$U
DTest <- test$D
K <- 9
propSD <- list(sigma2 = seq(0.1, 1, length = K),
               tau2 = 0.15)
#theta <- runif(9, 0.5, 3)
theta <- trueTheta

results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                theta = theta,
                propSD = propSD,
                nIter = 400, nBurn = 100, nThin=2,
                model = "full_gp")

theta
results$posteriorMeans
results$acceptance
nSamples <- length(results$paramSamples[[3]])
plot(1:nSamples, results$paramSamples[[3]], type="l")

library(MBA)
library(fields)

pdf("figures/subj1.pdf")
pred.surf <-  mba.surf(cbind(UTest, YTest[1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="True Surface, Subject 1", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)

pred.surf <-  mba.surf(cbind(UTest, results$preds[2,1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="Predicted Surface, Subject 1", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

pdf("figures/subj2.pdf")
pred.surf <-  mba.surf(cbind(UTest, YTest[(nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="True Surface, Subject 2", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)

pred.surf <-  mba.surf(cbind(UTest, results$preds[2,(nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="Predicted Surface, Subject 2", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

nTestSubj <- nrow(test$Z)
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
