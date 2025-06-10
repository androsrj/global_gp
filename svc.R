library(spBayes)
source("other_functions/helper_functions.R")

run_svc <- function(scen, nReps) {
  d.max <- max(iDist(train$U))
  r <- 2
  n <<- nrow(train$X)
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
  return(paste0("Completed Scenario ", scen))
}

nReps <- 10

# Scenario 1
load("data/small/scen1/train.RData")
load("data/small/scen1/test.RData")
run_svc(1, nReps)


# Scenario 2
load("data/small/scen2/train.RData")
load("data/small/scen2/test.RData")
run_svc(2, nReps)


# Scenario 3
load("data/small/scen3/train.RData")
load("data/small/scen3/test.RData")
run_svc(3, nReps)


# Scenario 4
load("data/small/scen4/train.RData")
load("data/small/scen4/test.RData")
run_svc(4, nReps)


# Scenario 5
load("data/small/scen5/train.RData")
load("data/small/scen5/test.RData")
run_svc(5, nReps)


# Scenario 6
load("data/small/scen6/train.RData")
load("data/small/scen6/test.RData")
run_svc(6, nReps)


# Scenario 7
load("data/small/scen7/train.RData")
load("data/small/scen7/test.RData")
run_svc(7, nReps)


# Scenario 8
load("data/small/scen8/train.RData")
load("data/small/scen8/test.RData")
run_svc(8, nReps)


# Scenario 9
load("data/small/scen9/train.RData")
load("data/small/scen9/test.RData")
run_svc(9, nReps)


# Scenario 10
load("data/small/scen10/train.RData")
load("data/small/scen10/test.RData")
run_svc(10, nReps)


# Scenario 11
load("data/small/scen11/train.RData")
load("data/small/scen11/test.RData")
run_svc(11, nReps)

