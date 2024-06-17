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
               tau2 = 0.15)
#theta <- runif(9, 0.5, 3)
theta <- trueTheta

results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                theta = theta,
                propSD = propSD,
                nIter = 150, nBurn = 100, nThin=2,
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

pred.surf <-  mba.surf(cbind(UTest, results$preds[2, (nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="Predicted Surface, Subject 2", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
dev.off()

