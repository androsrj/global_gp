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

run.mcmc <- function(rep) {
  results <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                  starting = starting,
                  propSD = propSD,
                  nIter = 500, nBurn = 500, nThin=2,
                  model = "full_gp")
  return(results)
}

##### SCENARIO 1 #####
scen <- 1
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
               thb = seq(3, 5, length = q),
               tau2 = 0.4)
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

##### SCENARIO 2 #####
scen <- 2
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
               theta = seq(0.2, 0.4, length = K),
               sigb2 = seq(0.7, 0.9, length = q),
               thb = seq(3, 5, length = q),
               tau2 = 0.4)
starting <- list(sigma2 = runif(K, 50, 100),
                 theta = rep(.25, K),
                 sigb2 = rep(7, q),
                 thb = rep(8, q), 
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

##### SCENARIO 3 #####
scen <- 3
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
               thb = seq(3, 5, length = q),
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
acc

##### SCENARIO 4 #####
scen <- 4
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
               thb = seq(3, 5, length = q),
               tau2 = 0.35)
starting <- list(sigma2 = runif(K, 50, 100),
                 theta = rep(.25, K),
                 sigb2 = rep(7, q),
                 thb = rep(2, q),  
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

##### SCENARIO 5 #####
scen <- 5
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
               thb = seq(3, 5, length = q),
               tau2 = 0.35)
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

##### SCENARIO 6 #####
scen <- 6
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
               thb = seq(3, 5, length = q),
               tau2 = 0.4)
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

##### SCENARIO 7 #####
scen <- 7
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
               thb = seq(3, 5, length = q),
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

##### SCENARIO 8 #####
scen <- 8
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
               thb = seq(3, 5, length = q),
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

##### SCENARIO 9 #####
scen <- 9
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
               thb = seq(3, 5, length = q),
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

##### SCENARIO 10 #####
scen <- 10
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
               thb = seq(3, 5, length = q),
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

##### SCENARIO 11 #####
scen <- 11
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
               thb = seq(3, 5, length = q),
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

#acc <- apply(sapply(1:nReps, \(i) unlist(results[[i]]$acceptance)), 1, mean)
#sigma2 <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigma2), 1, mean)
#theta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$theta), 1, mean)
#sigf2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigf2))
#thetaf <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$thf))
#beta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$beta), 1, mean)
#tau2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$tau2))




