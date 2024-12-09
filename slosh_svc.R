# SOURCES
source("other_functions/helper_functions.R") 
library(fields)
library(ggplot2)
library(spBayes)

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
}

cat(paste0("YTest Standard Dev: ", round(sd(YTest), 3), "\n"))

rmse
cat(paste0("Root MS error: ", round(mean(rmse), 3), "\n"))

cvg
cat(paste0("Mean coverage: ", round(mean(cvg), 3), "\n"))

width
cat(paste0("Mean width: ", round(mean(width), 3), "\n"))