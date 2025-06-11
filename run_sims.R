# SOURCES
source("mcmc_functions/mcmc.R") # Metropolis-Gibbs Sampler
source("mcmc_functions/priors.R")
source("mcmc_functions/jacobians.R")
source("mcmc_functions/likelihood.R")
source("mcmc_functions/posterior.R")
source("other_functions/helper_functions.R") # Other misc functions (not part of MCMC)
source("other_functions/bsplines_2_3D.R")
library(MBA)
library(fields)

run.sims <- function(scen, nReps, starting, propSD) {
  dir <- paste0("data/", size, "/scen", scen, "/")
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
  K <- 9
  
  results <- vector("list", length = nReps)
  for (i in 1:nReps) {
    results[[i]] <- mcmc(X = X, Z = Z, Y = Y, D = D, K = K,
                         starting = starting,
                         propSD = propSD,
                         nIter = 2000, nBurn = 1000, nThin=2,
                         model = "full_gp")
  }
  
  path <- paste0("objects/", size, "_scen", scen, ".RDS") 
  saveRDS(results, file = path)
  
  acc <- apply(sapply(1:nReps, \(i) unlist(results[[i]]$acceptance)), 1, mean)
  sigma2 <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigma2), 1, mean)
  theta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$theta), 1, mean)
  sigf2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$sigf2))
  thetaf <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$thf))
  beta <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$beta), 1, mean)
  tau2 <- mean(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$tau2))
  
  return(list(acceptance = acc,
              sigma2 = sigma2, 
              theta = theta,
              sigf2 = sigf2,
              thetaf = thetaf,
              beta = beta,
              tau2 = tau2))
}

##### SCENARIO 1 #####
propSD <- list(sigf2 = 0.4,
               thf = 2,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.4,
               theta = seq(0.4, 0.7, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 6,
                 thf = 1.3, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen1 <- run.sims(1, 10, starting, propSD)
cat("Scenario 1: \n")
cat("Avg. Acceptance: \n")
scen1$acceptance
cat("Sigma2: \n")
scen1$sigma2
cat("Theta: \n")
scen1$theta
cat("Sigma2_f: \n")
scen1$sigf2
cat("Theta_f: \n")
scen1$thetaf
cat("Beta: \n")
scen1$beta
cat("Tau2: \n")
scen1$tau2

##### SCENARIO 2 #####
propSD <- list(sigf2 = 0.6,
               thf = 5,
               sigma2 = seq(0.1, 0.25, length = K),
               tau2 = 0.4,
               theta = seq(0.5, 0.8, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 6,
                 thf = 8, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen2 <- run.sims(2, 10, starting, propSD)
cat("Scenario 2: \n")
cat("Avg. Acceptance: \n")
scen2$acceptance
cat("Sigma2: \n")
scen2$sigma2
cat("Theta: \n")
scen2$theta
cat("Sigma2_f: \n")
scen2$sigf2
cat("Theta_f: \n")
scen2$thetaf
cat("Beta: \n")
scen2$beta
cat("Tau2: \n")
scen2$tau2

##### SCENARIO 3 #####
propSD <- list(sigf2 = 0.6,
               thf = 2,
               sigma2 = seq(0.1, 0.25, length = K),
               tau2 = 0.3,
               theta = seq(0.5, 0.8, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 15,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen3 <- run.sims(3, 10, starting, propSD)
cat("Scenario 3: \n")
cat("Avg. Acceptance: \n")
scen3$acceptance
cat("Sigma2: \n")
scen3$sigma2
cat("Theta: \n")
scen3$theta
cat("Sigma2_f: \n")
scen3$sigf2
cat("Theta_f: \n")
scen3$thetaf
cat("Beta: \n")
scen3$beta
cat("Tau2: \n")
scen3$tau2

##### SCENARIO 4 #####
propSD <- list(sigf2 = 0.6,
               thf = 2,
               sigma2 = seq(0.1, 0.25, length = K),
               tau2 = 0.35,
               theta = seq(0.9, 1.5, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 1.5,
                 beta = c(0, 0, 0))
scen4 <- run.sims(4, 10, starting, propSD)
cat("Scenario 4: \n")
cat("Avg. Acceptance: \n")
scen4$acceptance
cat("Sigma2: \n")
scen4$sigma2
cat("Theta: \n")
scen4$theta
cat("Sigma2_f: \n")
scen4$sigf2
cat("Theta_f: \n")
scen4$thetaf
cat("Beta: \n")
scen4$beta
cat("Tau2: \n")
scen4$tau2

##### SCENARIO 5 #####
propSD <- list(sigf2 = 0.6,
               thf = 1.5,
               sigma2 = seq(0.2, 0.4, length = K),
               tau2 = 0.35,
               theta = seq(0.9, 1.5, length = K))
starting <- list(sigma2 = seq(5, 10, length = K),
                 theta = rep(0.5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen5 <- run.sims(5, 10, starting, propSD)
cat("Scenario 5: \n")
cat("Avg. Acceptance: \n")
scen5$acceptance
cat("Sigma2: \n")
scen5$sigma2
cat("Theta: \n")
scen5$theta
cat("Sigma2_f: \n")
scen5$sigf2
cat("Theta_f: \n")
scen5$thetaf
cat("Beta: \n")
scen5$beta
cat("Tau2: \n")
scen5$tau2

##### SCENARIO 6 #####
propSD <- list(sigf2 = 0.6,
               thf = 0.7,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen6 <- run.sims(1, 10, starting, propSD)
cat("Scenario 6: \n")
cat("Avg. Acceptance: \n")
scen6$acceptance
cat("Sigma2: \n")
scen6$sigma2
cat("Theta: \n")
scen6$theta
cat("Sigma2_f: \n")
scen6$sigf2
cat("Theta_f: \n")
scen6$thetaf
cat("Beta: \n")
scen6$beta
cat("Tau2: \n")
scen6$tau2

##### SCENARIO 7 #####
propSD <- list(sigf2 = 0.6,
               thf = 0.7,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen7 <- run.sims(7, 10, starting, propSD)
cat("Scenario 7: \n")
cat("Avg. Acceptance: \n")
scen7$acceptance
cat("Sigma2: \n")
scen7$sigma2
cat("Theta: \n")
scen7$theta
cat("Sigma2_f: \n")
scen7$sigf2
cat("Theta_f: \n")
scen7$thetaf
cat("Beta: \n")
scen7$beta
cat("Tau2: \n")
scen7$tau2

##### SCENARIO 8 #####
propSD <- list(sigf2 = 0.6,
               thf = 1,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(5, K),
                 sigf2 = 4,
                 thf = 5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen8 <- run.sims(8, 10, starting, propSD)
cat("Scenario 8: \n")
cat("Avg. Acceptance: \n")
scen8$acceptance
cat("Sigma2: \n")
scen8$sigma2
cat("Theta: \n")
scen8$theta
cat("Sigma2_f: \n")
scen8$sigf2
cat("Theta_f: \n")
scen8$thetaf
cat("Beta: \n")
scen8$beta
cat("Tau2: \n")
scen8$tau2

##### SCENARIO 9 #####
propSD <- list(sigf2 = 0.6,
               thf = 0.7,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(5, K),
                 sigf2 = 10,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen9 <- run.sims(9, 10, starting, propSD)
cat("Scenario 9: \n")
cat("Avg. Acceptance: \n")
scen9$acceptance
cat("Sigma2: \n")
scen9$sigma2
cat("Theta: \n")
scen9$theta
cat("Sigma2_f: \n")
scen9$sigf2
cat("Theta_f: \n")
scen9$thetaf
cat("Beta: \n")
scen9$beta
cat("Tau2: \n")
scen9$tau2

##### SCENARIO 10 #####
propSD <- list(sigf2 = 0.6,
               thf = 0.7,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(5, 10, length = K),
                 theta = rep(5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen10 <- run.sims(10, 10, starting, propSD)
cat("Scenario 10: \n")
cat("Avg. Acceptance: \n")
scen10$acceptance
cat("Sigma2: \n")
scen10$sigma2
cat("Theta: \n")
scen10$theta
cat("Sigma2_f: \n")
scen10$sigf2
cat("Theta_f: \n")
scen10$thetaf
cat("Beta: \n")
scen10$beta
cat("Tau2: \n")
scen10$tau2

##### SCENARIO 11 #####
propSD <- list(sigf2 = 0.6,
               thf = 0.7,
               sigma2 = seq(0.1, 0.2, length = K),
               tau2 = 0.35,
               theta = seq(0.5, 0.9, length = K))
starting <- list(sigma2 = seq(50, 100, length = K),
                 theta = rep(5, K),
                 sigf2 = 4,
                 thf = 1.5, 
                 tau2 = 0.1,
                 beta = c(0, 0, 0))
scen11 <- run.sims(11, 10, starting, propSD)
cat("Scenario 11: \n")
cat("Avg. Acceptance: \n")
scen11$acceptance
cat("Sigma2: \n")
scen11$sigma2
cat("Theta: \n")
scen11$theta
cat("Sigma2_f: \n")
scen11$sigf2
cat("Theta_f: \n")
scen11$thetaf
cat("Beta: \n")
scen11$beta
cat("Tau2: \n")
scen11$tau2


