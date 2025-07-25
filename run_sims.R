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
nReps <- nCores <- 10
set.seed(999)
which.scens <- 2

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
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.3, 0.5, length = K),
                 theta = seq(0.7, 0.9, length = K),
                 sigb2 = seq(0.2, 0.3, length = q),
                 thb = seq(0.3, 0.4, length = q),
                 tau2 = 1.0)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(2, q),
                   thb = rep(0.5, q), 
                   tau2 = 2,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
}

##### SCENARIO 2 #####
scen <- 2
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.3, 0.5, length = K),
                 theta = seq(0.7, 0.9, length = K),
                 sigb2 = seq(0.2, 0.4, length = q),
                 thb = seq(0.5, 0.7, length = q),
                 tau2 = 0.9)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(0.2, q),
                   thb = rep(0.2, q), 
                   tau2 = 0.1,
                   beta = c(5, 2, -4))
  #cl <- makeCluster(nCores)
  #registerDoParallel(cl)
  #obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  #stopCluster(cl)
  obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  #acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  #acc
  print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  #obj[[1]]$posteriorMeans
  print(obj$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj$preds[2,] - test$Y)^2))
  print(rmse)
}

##### SCENARIO 3 #####
scen <- 3
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.1, 0.2, length = K),
                 theta = seq(0.2, 0.4, length = K),
                 sigb2 = seq(0.7, 0.9, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(17, q),
                   thb = rep(2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
}

##### SCENARIO 4 #####
scen <- 4
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.1, 0.25, length = K),
                 theta = seq(0.3, 0.6, length = K),
                 sigb2 = seq(0.7, 0.9, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.35)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(2, q),
                   thb = rep(0.2, q),  
                   tau2 = 1.5,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 5 #####
scen <- 5
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.2, 0.4, length = K),
                 theta = seq(0.4, 0.7, length = K),
                 sigb2 = seq(0.7, 0.9, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.35)
  starting <- list(sigma2 = runif(K, 5, 10),
                   theta = rep(.25, K),
                   sigb2 = rep(2, q),
                   thb = rep(0.2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 6 #####
scen <- 6
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.1, 0.2, length = K),
                 theta = seq(0.3, 0.6, length = K),
                 sigb2 = seq(0.7, 0.9, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.4)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(0.2, K),
                   sigb2 = rep(2, q),
                   thb = rep(0.2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 7 #####
scen <- 7
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.4, 0.6, length = K),
                 theta = seq(3, 4, length = K),
                 sigb2 = seq(1.3, 1.7, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.25)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(7, q),
                   thb = rep(2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 8 #####
scen <- 8
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.3, 0.5, length = K),
                 theta = seq(4, 5, length = K),
                 sigb2 = seq(0.5, 0.7, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.25)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(7, q),
                   thb = rep(2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 9 #####
scen <- 9
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.4, 0.6, length = K),
                 theta = seq(1.6, 2, length = K),
                 sigb2 = seq(1.2, 1.6, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.2)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(.25, K),
                   sigb2 = rep(17, q),
                   thb = rep(2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 10 #####
scen <- 10
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.4, 0.7, length = K),
                 theta = seq(2.5, 3, length = K),
                 sigb2 = seq(0.5, 0.7, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.25)
  starting <- list(sigma2 = runif(K, 5, 10),
                   theta = rep(.25, K),
                   sigb2 = rep(7, q),
                   thb = rep(2, q), 
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

##### SCENARIO 11 #####
scen <- 11
if (scen %in% which.scens) {
  dir <- paste0("data/small/scen", scen, "/")
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
  propSD <- list(sigma2 = seq(0.4, 0.6, length = K),
                 theta = seq(3, 4, length = K),
                 sigb2 = seq(1.5, 2.5, length = q),
                 thb = seq(0.3, 0.5, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = runif(K, 50, 100),
                   theta = rep(2, K),
                   sigb2 = rep(7, q),
                   thb = rep(2, q),  
                   tau2 = 0.3,
                   beta = c(0, 0, 0))
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields")) %dopar% run.mcmc(i)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  acc
  cat(paste0("And posterior means of: "))
  obj[[1]]$posteriorMeans
}

#acc <- apply(sapply(1:nReps, \(i) unlist(results[[i]]$acceptance)), 1, mean)
#sigma2 <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigma2), 1, mean)
#theta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$theta), 1, mean)
#sigf2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigf2))
#thetaf <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$thf))
#beta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$beta), 1, mean)
#tau2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$tau2))




