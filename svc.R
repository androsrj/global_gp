library(spBayes)
source("other_functions/helper_functions.R")

run.svc <- function(scen, nReps) {
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
    m.3 <- spSVC(train$Y[1:nrow(train$X),] ~ train$X, coords=train$U,
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
train <- readRDS("data/small/scen1/train.RDS")
test <- readRDS("data/small/scen1/test.RDS")
run.svc(1, nReps)


# Scenario 2
train <- readRDS("data/small/scen2/train.RDS")
test <- readRDS("data/small/scen2/test.RDS")
run.svc(2, nReps)


# Scenario 3
train <- readRDS("data/small/scen3/train.RDS")
test <- readRDS("data/small/scen3/test.RDS")
run.svc(3, nReps)


# Scenario 4
train <- readRDS("data/small/scen4/train.RDS")
test <- readRDS("data/small/scen4/test.RDS")
run.svc(4, nReps)


# Scenario 5
train <- readRDS("data/small/scen5/train.RDS")
test <- readRDS("data/small/scen5/test.RDS")
run.svc(5, nReps)


# Scenario 6
train <- readRDS("data/small/scen6/train.RDS")
test <- readRDS("data/small/scen6/test.RDS")
run.svc(6, nReps)


# Scenario 7
train <- readRDS("data/small/scen7/train.RDS")
test <- readRDS("data/small/scen7/test.RDS")
run.svc(7, nReps)


# Scenario 8
train <- readRDS("data/small/scen8/train.RDS")
test <- readRDS("data/small/scen8/test.RDS")
run.svc(8, nReps)


# Scenario 9
train <- readRDS("data/small/scen9/train.RDS")
test <- readRDS("data/small/scen9/test.RDS")
run.svc(9, nReps)


# Scenario 10
train <- readRDS("data/small/scen10/train.RDS")
test <- readRDS("data/small/scen10/test.RDS")
run.svc(10, nReps)


# Scenario 11
train <- readRDS("data/small/scen11/train.RDS")
test <- readRDS("data/small/scen11/test.RDS")
run.svc(11, nReps)

