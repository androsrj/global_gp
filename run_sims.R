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
set.seed(9994)
which.scens <- 13

run.mcmc <- function(rep) {
  results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                  starting = starting,
                  propSD = propSD,
                  nIter = 5000, nBurn = 2000, nThin=2,
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
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.2, 0.3, length = K),
                 theta = seq(0.7, 1.2, length = K),
                 sigb2 = seq(0.3, 0.4, length = q),
                 thb = seq(0.4, 0.5, length = q),
                 tau2 = 0.25)
  starting <- list(sigma2 = rep(10, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(0.2, q), 
                   tau2 = 1.5)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.2, 0.4, length = K),
                 theta = seq(0.5, 0.7, length = K),
                 sigb2 = seq(0.5, 0.7, length = q),
                 thb = seq(0.9, 1.1, length = q),
                 tau2 = 0.2)
  starting <- list(sigma2 = rep(10, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(2, q), 
                   tau2 = 2)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.3, 0.5, length = K),
                 theta = seq(1.4, 1.6, length = K),
                 sigb2 = seq(0.4, 0.6, length = q),
                 thb = seq(1.2, 1.4, length = q),
                 tau2 = 0.2)
  starting <- list(sigma2 = runif(K, 1, 10),
                   theta = rep(.25, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(0.2, q),  
                   tau2 = 5)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  #obj <- run.mcmc(1)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  #print(obj$acceptance)
  print(acc)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.15, 0.2, length = K),
                 theta = seq(0.3, 0.5, length = K),
                 sigb2 = seq(0.4, 0.6, length = q),
                 thb = seq(2.2, 2.4, length = q),
                 tau2 = 0.5)
  starting <- list(sigma2 = rep(20, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.2, q),
                   thb = rep(0.2, q), 
                   tau2 = 2)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  #obj <- run.mcmc(1)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  #print(obj$acceptance)
  print(acc)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.2, 0.3, length = K),
                 theta = seq(0.5, 0.8, length = K),
                 sigb2 = seq(0.6, 0.8, length = q),
                 thb = seq(2.4, 2.7, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = runif(K, 1, 10),
                   theta = rep(2, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(0.2, q), 
                   tau2 = 2)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  #obj <- run.mcmc(1)
  stopCluster(cl)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  #print(obj$acceptance)
  print(acc)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.6, 0.9, length = K),
                 theta = seq(1, 1.5, length = K),
                 sigb2 = seq(0.3, 0.5, length = q),
                 thb = seq(0.01, 0.1, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = rep(50, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.2, q),
                   thb = rep(0.2, q),
                   tau2 = 2)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.7, 1.0, length = K),
                 theta = seq(0.7, 1.2, length = K),
                 sigb2 = seq(0.4, 0.5, length = q),
                 thb = seq(0.4, 0.5, length = q),
                 tau2 = 0.25)
  starting <- list(sigma2 = rep(10, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(0.2, q), 
                   tau2 = 1.5)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
                 theta = seq(0.7, 0.9, length = K),
                 sigb2 = seq(0.5, 0.7, length = q),
                 thb = seq(0.9, 1.1, length = q),
                 tau2 = 0.2)
  starting <- list(sigma2 = rep(10, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(2, q), 
                   tau2 = 2)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  propSD <- list(sigma2 = seq(0.5, 0.7, length = K),
                 theta = seq(1.4, 1.6, length = K),
                 sigb2 = seq(0.4, 0.6, length = q),
                 thb = seq(1.2, 1.4, length = q),
                 tau2 = 0.2)
  starting <- list(sigma2 = runif(K, 1, 10),
                   theta = rep(.25, K),
                   sigb2 = rep(0.1, q),
                   thb = rep(0.2, q),  
                   tau2 = 5)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
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
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
}


##### SCENARIO 12 #####
scen <- 12
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
  propSD <- list(sigma2 = seq(0.6, 0.8, length = K),
                 theta = seq(1.5, 2.0, length = K),
                 sigb2 = seq(0.3, 0.4, length = q),
                 thb = seq(0.01, 0.1, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = rep(50, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.2, q),
                   thb = rep(0.2, q),
                   tau2 = 1)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
}

##### SCENARIO 13 #####
scen <- 13
if (scen %in% which.scens) {
  source("mcmc_functions/smooth/mcmc_smooth.R") # Metropolis-Gibbs Sampler
  source("mcmc_functions/priors.R")
  source("mcmc_functions/jacobians.R")
  source("mcmc_functions/likelihood.R")
  source("mcmc_functions/smooth/posterior.R")
  source("mcmc_functions/smooth/helper_functions.R") # Other misc functions (not part of MCMC)
  source("other_functions/bsplines_2_3D.R")
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
  propSD <- list(sigma2 = seq(0.6, 0.8, length = K),
                 theta = seq(1.5, 2.0, length = K),
                 sigb2 = seq(0.3, 0.4, length = q),
                 thb = seq(0.01, 0.1, length = q),
                 tau2 = 0.3)
  starting <- list(sigma2 = rep(50, K),
                   theta = rep(.25, K),
                   sigb2 = rep(0.2, q),
                   thb = rep(0.2, q),
                   tau2 = 1)
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  obj <- foreach(i = 1:nCores, .packages = c("mvtnorm", "splines", "fields", "Matrix")) %dopar% run.mcmc(i)
  stopCluster(cl)
  #obj <- run.mcmc(1)
  saveRDS(obj, file = paste0("objects/small_scen", scen, ".RDS"))
  acc <- apply(sapply(1:nReps, \(i) unlist(obj[[i]]$acceptance)), 1, mean)
  cat(paste0("Finished Scenario ", scen, " with average acceptance of: "))
  print(acc)
  #print(obj$acceptance)
  cat(paste0("And posterior means of: "))
  print(obj[[1]]$posteriorMeans)
  cat("And RMSE of: ")
  rmse <- sqrt(mean((obj[[1]]$preds[2,] - test$Y)^2))
  print(rmse)
  cat("Compared to the test data SD of:")
  sd.test <- round(sd(test$Y), 3)
  print(sd.test)
  cat("Average prediction interval width of: ")
  print(mean(obj[[1]]$preds[3,] - obj[[1]]$preds[1,]))
  cat("With average coverage of: ")
  print(mean(obj[[1]]$preds[3,] > test$Y & obj[[1]]$preds[1,] < test$Y))
}
