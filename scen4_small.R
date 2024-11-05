# SOURCES
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")

library(fields)
size <- "small"
scen <- "scen4"
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
               thf = 2,
               sigma2 = seq(0.1, 0.25, length = K),
               tau2 = 0.35,
               theta = seq(0.9, 1.5, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 1.5,
                 beta = c(2, 0, 0))

results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                starting = starting,
                propSD = propSD,
                nIter = 2000, nBurn = 2000, nThin=2,
                model = "full_gp")

path <- paste0("objects/", size, "_", scen, ".RDS") 
saveRDS(results, file = path)

#theta
mean(train$Y)
sd(train$Y)
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
