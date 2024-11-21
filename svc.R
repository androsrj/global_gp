library(spBayes)

scen <- 1
nReps <- 10

path <- paste0("data/small/scen", scen) 
load(paste0(path, "/train.RData")) 
load(paste0(path, "/test.RData"))
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
n.samples <- 5000

results <- vector("list", nReps)
for (i in 1:nReps) {
  m.3 <- spSVC(train$Y[1:n,] ~ train$X, coords=train$U,
               starting=starting, svc.cols=c(1,2),
               tuning=tuning, priors=priors, cov.model="exponential",
               n.samples=n.samples, n.report=5000, n.omp.threads=4)
  
  m.3 <- spRecover(m.3, start=floor(0.5*n.samples), thin=2,
                   n.omp.threads=4, verbose=FALSE)
  
  results[[i]] <- m.3
  
}

saveRDS(results, paste0("objects/svc_scen", scen, ".RDS"))




