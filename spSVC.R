library(spBayes)
library(coda)
obj <- readRDS("objects/small_scen1.RDS")
n <- 100

p.theta.samples <- as.mcmc(cbind(t(obj[[1]]$paramSamples$sigb2),
                                 obj[[1]]$paramSamples$tau2,
                                 t(obj[[1]]$paramSamples$thb)))
colnames(p.theta.samples) <- c("sigma.sq.(Intercept)", 
                               "sigma.sq.train$X1",
                               "sigma.sq.train$X2",
                               "tau.sq",
                               "phi.(Intercept)",
                               "phi.train$X1",
                               "phi.train$X2")
head(p.theta.samples)

acceptance <- 53.08
Y <- train$Y[1:n,]
X <- cbind(rep(1, n), train$X)
colnames(X) <- c("(Intercept)", "train$X1", "train$X2")
which.Z <- 1:3
Z <- X[ , which.Z]
colnames(Z) <- colnames(X)[which.Z]
center.scale <- FALSE
X.sc = scale(X)
coords <- train$U
cov.model <- "exponential"
nugget <- 1
beta.prior <- "flat"
beta.Norm <- 0
x.names <- colnames(X)
run.time <- m.3$run.time
K.diag <- 1
svc.cols <- 1:3
dropped.obs <- rep(FALSE, n)
n.samples <- nrow(p.theta.samples)

svc.obj <- list(p.theta.samples = p.theta.samples,
                acceptance = acceptance,
                Y = Y,
                X = X,
                Z = Z,
                center.scale = center.scale,
                X.sc = X.sc,
                coords = coords,
                cov.model = cov.model,
                nugget = nugget,
                beta.prior = beta.prior,
                beta.Norm = beta.Norm,
                x.names = x.names,
                run.time = run.time,
                K.diag = K.diag,
                svc.cols = svc.cols,
                dropped.obs = dropped.obs)
class(svc.obj) <- "spSVC"
coefs <- spRecover(svc.obj, start=floor(0.5*n.samples), thin=2,
                   n.omp.threads=4, verbose=FALSE)
coefs$p.beta.recover.samples
colMeans(coefs$p.beta.recover.samples)
dim(coefs$p.w.recover.samples)

beta.est <- colMeans(coefs$p.beta.recover.samples)
beta.est[1] + rowMeans(coefs$p.w.recover.samples)[1:100]
