library(fields)
library(spBayes)
library(refund)
source("other_functions/helper_functions.R")
nScen <- 5
scenarios <- 7:(6+nScen)

# Which models do you want diagnostics for?
GGP <- TRUE
SVC <- TRUE
FOSR <- TRUE

#######################
###### Global GP ######
#######################
if (GGP == TRUE) {
  for (s in scenarios) {
    cat(paste0("Global GP, Scenario ", s, ": \n"))
    
    # Read in model results and original data
    results <- readRDS(paste0("objects/small_scen", s, ".RDS"))
    nReps <- length(results)
    train <- readRDS(paste0("data/small/scen", s, "/train.RDS"))
    test <- readRDS(paste0("data/small/scen", s, "/test.RDS"))
    
    # Std deviation of original data
    cat("SD of Y: \n")
    cat(round(sd(train$Y), 2))
    
    # Beta estimates and credible intervals
    beta.means <- apply(sapply(1:nReps, \(i) results[[i]]$posteriorMeans$beta), 1, mean)
    beta.lower <- apply(sapply(1:nReps, \(i) results[[i]]$credLower$beta), 1, mean)
    beta.upper <- apply(sapply(1:nReps, \(i) results[[i]]$credUpper$beta), 1, mean)
    beta.mat <- cbind(beta.means, beta.lower, beta.upper)
    colnames(beta.mat) <- c("Mean", "Lower", "Upper")
    rownames(beta.mat) <- c("Beta0", "Beta1", "Beta2")
    print(beta.mat)
    
    # RMSE and Coverage (averaged from all reps)
    STest <- nrow(test$Z)
    nTest <- nrow(test$X)
    a <- 0.05
    diagnostics <- sapply(1:nReps, function(j) {
      rmse.vec <- cvg.vec <- width.vec <- scores.vec <- crps.vec <- numeric(STest)
      for (i in 1:STest) {
        truth <- test$Y[(nTest*(i-1)+1):(i*nTest)]
        pred <- results[[j]]$preds[2, (nTest*(i-1)+1):(i*nTest)]
        rmse.vec[i] <- sqrt(mean((truth - pred)^2))
        lower <- results[[j]]$preds[1, (nTest*(i-1)+1):(i*nTest)]
        upper <- results[[j]]$preds[3, (nTest*(i-1)+1):(i*nTest)]
        cvg.vec[i] <- mean(lower < truth & upper > truth)
        #width.vec[i] <- mean(upper - lower)
        #scores.vec[i] <- mean((upper - lower) + 
        #                    2/a * (lower - truth) * (truth < lower) + 
        #                    2/a * (truth - upper) * (truth > upper))
        #predSamples <- t(results[[j]]$predSamples[(nTest*(i-1)+1):(i*nTest), ])
        #crps.vec[i] <- mean(energy_score(truth, predSamples))
      }
      c(rmse = mean(rmse.vec), cvg = mean(cvg.vec))
    })
    cat(paste0("RMSE: ", round(mean(diagnostics["rmse", ]), 3), "\n"))
    cat(paste0("Coverage: ", round(mean(diagnostics["cvg", ]), 3), "\n"))
  }
}


#######################
######### SVC #########
#######################
if (SVC == TRUE) {
  for (s in scenarios) {
    cat(paste0("SVC, Scenario ", s, ": \n"))
    
    # Read in model results and original data
    results <- readRDS(paste0("objects/svc_scen", s, ".RDS"))
    nReps <- length(results)
    train <- readRDS(paste0("data/small/scen", s, "/train.RDS"))
    test <- readRDS(paste0("data/small/scen", s, "/test.RDS"))
    
    # Std deviation of original data
    cat("SD of Y: \n")
    cat(round(sd(train$Y), 2))
    
    # Betas # Beta estimates and credible intervals
    beta.means <- apply(sapply(1:nReps, \(i) apply(results[[i]]$p.beta.recover.samples, 2, mean)), 1, mean)
    beta.lower <- apply(sapply(1:nReps, \(i) apply(results[[i]]$p.beta.recover.samples, 2, quantile, .025)), 1, mean)
    beta.upper <- apply(sapply(1:nReps, \(i) apply(results[[i]]$p.beta.recover.samples, 2, quantile, .975)), 1, mean)
    beta.mat <- cbind(beta.means, beta.lower, beta.upper)
    colnames(beta.mat) <- c("Mean", "Lower", "Upper")
    rownames(beta.mat) <- c("Beta0", "Beta1", "Beta2")
    print(beta.mat)
    
    # RMSE and Coverage (averaged from all reps)
    STest <- nrow(test$Z)
    nTest <- nrow(test$X)
    a <- 0.05
    diagnostics <- sapply(1:nReps, function(j) {
      rmse.vec <- cvg.vec <- width.vec <- scores.vec <- crps.vec <- numeric(STest)
      for (i in 1:STest) {
        truth <- test$Y[(nTest*(i-1)+1):(nTest*i), ]
        m.3.pred <- spPredict(results[[j]], pred.covars = cbind(rep(1, nTest), test$X),
                              pred.coords=test$U + rnorm(50, 0, 0.0001), thin=10,
                              joint=TRUE, n.omp.threads=4, verbose=FALSE)
        preds <- apply(m.3.pred$p.y.predictive.samples, 1, mean)
        rmse.vec[i] <- sqrt(mean((truth - preds)^2))
        lower <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .025)
        upper <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .975)
        cvg.vec[i] <- mean(lower < truth & upper > truth)
        #width.vec[i] <- mean(upper - lower)
      }
      c(rmse = mean(rmse.vec), cvg = mean(cvg.vec))
    })
    cat(paste0("RMSE: ", round(mean(diagnostics["rmse", ]), 3), "\n"))
    cat(paste0("Coverage: ", round(mean(diagnostics["cvg", ]), 3), "\n"))
  }
}



########################
######### FOSR #########
########################
if (FOSR == TRUE) {
  for (s in scenarios) {
    cat(paste0("FOSR, Scenario ", s, ": \n"))
    
    # Read in original data
    train <- readRDS(paste0("data/small/scen", s, "/train.RDS"))
    test <- readRDS(paste0("data/small/scen", s, "/test.RDS"))
    
    # Fit FOSR model
    n <- nrow(train$X)
    nTest <- nrow(test$X)
    S <- nrow(train$Z)
    STest <- nrow(test$Z)
    Y <- rbind(matrix(train$Y, nrow = n, ncol = S),
               matrix(test$Y, nrow = nTest, ncol = STest))
    X <- rbind(train$X, test$X)
    colnames(X) <- c("X1", "X2")
    dfl <- as.data.frame(X)
    dfl$Y <- Y
    train.index <- 1:n
    test.index <- (n+1):(n+nTest)
    fit <- bayes_fosr(data = dfl[train.index, ], Y ~ X1 + X2, est.method = "VB")
    
    # Std deviation of original data
    cat("SD of Y: \n")
    cat(round(sd(train$Y), 2))
    
    # Betas # Beta estimates and credible intervals
    beta.means <- apply(fit$beta.hat, 1, mean)
    beta.lower <- apply(fit$beta.LB, 1, mean)
    beta.upper <- apply(fit$beta.UB, 1, mean)
    beta.mat <- cbind(beta.means, beta.lower, beta.upper)
    colnames(beta.mat) <- c("Mean", "Lower", "Upper")
    rownames(beta.mat) <- c("Beta0", "Beta1", "Beta2")
    print(beta.mat)
    
    # Get model's predictions for test data
    preds <- predict(object = fit, newdata = dfl[test.index,])
    rmse <- sqrt(mean((preds - Y[test.index, ])^2))
    cat(paste0("RMSE: ", round(rmse, 3), "\n"))
  }
}




