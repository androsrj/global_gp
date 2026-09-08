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
library(Matrix)
nReps <- nCores <- 10
set.seed(999)
which.scens <- 1:8

run.mcmc <- function(rep) {
  results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                  starting = starting,
                  propSD = propSD,
                  nIter = 2000, nBurn = 2000, nThin=2,
                  model = "full_gp")
  return(results)
}

##### SCENARIO 1 #####
scen <- 1
if (scen %in% which.scens) {
  dir <- paste0("data/scen", scen, "/")
  train <- readRDS(paste0(dir, "train.RDS"))
  test <- readRDS(paste0(dir, "test.RDS"))
  n <- nrow(train$X)
  nTest <- nrow(test$X)
  X <- train$X; XTest <- test$X
  Z <- train$Z; ZTest <- test$Z
  Y <- train$Y; YTest <- test$Y
  U <- train$U; UTest <- test$U
  D <- train$D; DTest <- test$D
  K <- 9
  q <- ncol(X) + 1
  propSD <- list(sigma2 = seq(0.2, 0.3, length = K),
                 theta = seq(0.7, 1.0, length = K),
                 sigb2 = seq(0.3, 0.4, length = q),
                 thb = seq(0.01, 0.1, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = rep(50, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.2, q),
                   thb = rep(0.2, q),
                   tau2 = 2)
  start.time <- Sys.time()
  run.mcmc(1)
  Sys.time() - start.time
}



