library(spBayes)
library(fields)

size <- "large"
scen <- "scen1"
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

fullX <- matrix(1, ncol=1, nrow=10) %x% train$X
fullCoords <- matrix(1, ncol=1, nrow=10) %x% train$U + rnorm(10000, 0, 1/10)
fullXTest <- matrix(1, ncol=1, nrow=10) %x% test$X

p <- 3
starting <- list("phi"=3/0.5, "sigma.sq"=10, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                 "phi.Unif"=c(0.5, 30), "sigma.sq.IG"=c(2, 2),
                 "tau.sq.IG"=c(2, 0.1))
cov.model <- "exponential"
n.report <- 10
n.samples <- 10
verbose <- TRUE
m.1 <- spLM(Y ~ fullX, coords=fullCoords, starting=starting,
            tuning=tuning, priors=priors.1, cov.model=cov.model,
            n.samples=n.samples, verbose=verbose, n.report=n.report)
m.1 <- spRecover(m.1, start=0.5*n.samples, thin=2)

testCoords <- matrix(1, ncol=1, nrow=nrow(test$Z)) %x% test$U + rnorm(500, 0, 1/10)
m.1.pred <- spPredict(m.1, pred.covars=cbind(matrix(1, ncol=1, nrow=nTest*nrow(test$Z)), fullXTest), 
                      pred.coords=testCoords,
                      start=0.5*n.samples)
y.hat <- apply(m.1.pred$p.y.predictive.samples, 1, mean)
quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}
y.quant <- apply(m.1.pred$p.y.predictive.samples, 1, quant)
save(m.1, m.1.pred, y.hat, y.quant, file = "objects/spLM.RData")

