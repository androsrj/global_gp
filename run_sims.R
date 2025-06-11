# SOURCES
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")
library(MBA)
library(splines)
library(fields)
library(parallel) 
library(doParallel)
library(foreach) 

run.sims <- function(scen) {
  dir <- paste0("data/small/scen", scen, "/")
  train <- readRDS(paste0(dir, "train.RDS"))
  test <- readRDS(paste0(dir, "test.RDS"))
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
  
  propSD <- readRDS("objects/propSD.RDS")[[scen]]
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigf2 = 7,
                   thf = 2, 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  if (scen == 2 | scen == 8) {
    starting$thf <- 8
  }
  if (scen == 3 | scen == 9) {
    starting$sigf2 <- 16
  }
  if (scen == 4) {
    starting$tau2 <- 1.6
  }
  if (scen == 5 | scen == 10) {
    starting$sigma2 <- runif(K, 5, 10)
  } 
  if (scen == 6 | scen == 11) {
    starting$theta <- rep(2.5, K)
  }
  
  results <- vector("list", length = nReps)
  for (i in 1:nReps) {
    results[[i]] <- mcmc(X = X, Z = Z, Y = Y, D = D, 
			 XTest = XTest, ZTest = ZTest, 
			 YTest = YTest, DTest = DTest, K = K,
                         starting = starting,
                         propSD = propSD,
                         nIter = 200, nBurn = 100, nThin=2,
                         model = "full_gp")
  }
  
  path <- paste0("objects/small_scen", scen, ".RDS") 
  saveRDS(results, file = path)
  
  acc <- apply(sapply(1:nReps, \(i) unlist(results[[i]]$acceptance)), 1, mean)
  sigma2 <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigma2), 1, mean)
  theta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$theta), 1, mean)
  sigf2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigf2))
  thetaf <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$thf))
  beta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$beta), 1, mean)
  tau2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$tau2))
  
  #cat(paste0("Scenario ", scen, ":\n"))
  #cat("Avg. Acceptance: \n")
  #acc
  #cat("Sigma2: \n")
  #sigma2
  #cat("Theta: \n")
  #theta
  #cat("Sigma2_f: \n")
  #sigf2
  #cat("Theta_f: \n")
  #thetaf
  #cat("Beta: \n")
  #beta
  #cat("Tau2: \n")
  #tau2
}

nReps <- 10
nCores <- 11
cl <- makeCluster(nCores)
registerDoParallel(cl)
obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.sims(i)
stopCluster(cl)
