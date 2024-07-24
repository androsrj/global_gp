library(spBayes)

load("data/train.RData")
load("data/test.RData")
load("data/theta.RData")
source("other_functions/helper_functions.R")

d.max <- max(iDist(train$U))

r <- 2
n <- nrow(train$X)
nTest <- nrow(test$X)
priors <- list("phi.Unif"=list(rep(3/(0.75*d.max), r), rep(3/(0.001*d.max), r)),
               "sigma.sq.IG"=list(rep(2, r), rep(1, r)),
               "tau.sq.IG"=c(2, 1))

starting <- list("phi"=rep(3/(0.1*d.max), r), "sigma.sq"=rep(1, r), "tau.sq"=1)
tuning <- list("phi"=rep(0.1, r), "sigma.sq"=rep(0.05, r), "tau.sq"=0.1)
n.samples <- 500

m.3 <- spSVC(train$Y[1:n,] ~ train$X - 1, coords=train$U,
             starting=starting, svc.cols=c(1,2),
             tuning=tuning, priors=priors, cov.model="exponential",
             n.samples=n.samples, n.report=5000, n.omp.threads=4)

m.3 <- spRecover(m.3, start=floor(0.5*n.samples), thin=2,
                 n.omp.threads=4, verbose=FALSE)

STest <- nrow(test$Z)
abs_error <- cvg <- width <- scores <- crps <- numeric(STest)
a <- .05
for (i in 1:STest) {
  truth <- test$Y[(nTest*(i-1)+1):(nTest*i), ]
  m.3.pred <- spPredict(m.3, pred.covars=test$X,
                        pred.coords=test$U + rnorm(50, 0, 0.0001), thin=10,
                        joint=TRUE, n.omp.threads=4, verbose=FALSE)
  pred <- apply(m.3.pred$p.y.predictive.samples, 1, mean)
  abs_error[i] <- mean(abs(truth - pred))
  lower <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .025)
  upper <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .975)
  cvg[i] <- mean(lower < truth & upper > truth)
  width[i] <- mean(upper - lower)
  scores[i] <- mean((upper - lower) + 
                      2/a * (lower - truth) * (truth < lower) + 
                      2/a * (truth - upper) * (truth > upper))
  predSamples <- t(m.3.pred$p.y.predictive.samples)
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




