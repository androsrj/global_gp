# SOURCES
source("other_functions/helper_functions.R") 
library(fields)
library(ggplot2)
library(spBayes)

load("data//slosh/flood_subset.RData")

mySeed <- 1234
which.Z <- c(1:5)
n <- 100
nTest <- 25
S <- 10
STest <- 10

X <- flood.train$X
Z <- flood.train$Z
Y <- flood.train$Y
U <- flood.train$U
D <- flood.train$D
XTest <- flood.test$X
ZTest <- flood.test$Z
YTest <- flood.test$Y
UTest <- flood.test$U
DTest <- flood.test$D

d.max <- max(iDist(U))
r <- 2
n <- nrow(X)
nTest <- nrow(XTest)
priors <- list("phi.Unif"=list(rep(3/(0.75*d.max), r), rep(3/(0.001*d.max), r)),
               "sigma.sq.IG"=list(rep(2, r), rep(1, r)),
               "tau.sq.IG"=c(2, 1))

starting <- list("phi"=rep(3/(0.1*d.max), r), "sigma.sq"=rep(1, r), "tau.sq"=1)
tuning <- list("phi"=rep(0.1, r), "sigma.sq"=rep(0.05, r), "tau.sq"=0.1)
n.samples <- 5000

m.3 <- spSVC(Y[1:n,] ~ X, coords=U,
             starting=starting, svc.cols=c(1,2),
             tuning=tuning, priors=priors, cov.model="exponential",
             n.samples=n.samples, n.report=n.samples/10, n.omp.threads=4)

m.3 <- spRecover(m.3, start=floor(0.5*n.samples), thin=2,
                 n.omp.threads=4, verbose=FALSE)

rmse <- cvg <- len <- numeric(STest)
svc.preds <- matrix(0, nrow = STest, ncol = nTest)
for (k in 1:STest) {
  truth <- YTest[(nTest*(k-1)+1):(nTest*k), ]
  m.3.pred <- spPredict(m.3, pred.covars = cbind(rep(1, nTest), XTest),
                        pred.coords=UTest + rnorm(50, 0, 0.0001), thin=10,
                        joint=TRUE, n.omp.threads=4, verbose=FALSE)
  preds <- apply(m.3.pred$p.y.predictive.samples, 1, mean)
  rmse[k] <- sqrt(mean((truth - preds)^2))
  lower <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .025)
  upper <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .975)
  cvg[k] <- mean(lower < truth & upper > truth)
  len[k] <- mean(upper - lower)
  svc.preds[k, ] <- preds
}

cat(paste0("YTest Standard Dev: ", round(sd(YTest), 3), "\n"))

rmse
cat(paste0("Root MS error: ", round(mean(rmse), 3), "\n"))

cvg
cat(paste0("Mean coverage: ", round(mean(cvg), 3), "\n"))

len
cat(paste0("Mean width: ", round(mean(len), 3), "\n"))

# Save predictions to use in plots later
saveRDS(svc.preds, file = "objects/slosh_svc_preds.RDS")

