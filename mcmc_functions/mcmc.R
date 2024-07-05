### FINAL MCMC SKETCHING FUNCTION ###
# Priors, Jacobians, likelihood, and posterior must already be sourced

mcmc <- function(X, Z, Y, D, K,
                 theta,
                 model = c("full_gp", "mpp", "sparse_gpp")[1],
                 propSD = c(0.1, 0.1), 
                 nIter = 1000, nBurn = 100, nThin = 2) {
  
  # Dimensions
  S <<- nrow(Z)
  n <<- length(Y) / S
  STest <<- nrow(ZTest)
  nTest <<- length(YTest) / STest
  J <- matrix(1, nrow = S, ncol = 1)
  JTest <- matrix(1, nrow = STest, ncol = 1)
  A <<- J %x% X
  ATest <<- JTest %x% XTest
  K <<- K
  p <<- ncol(X)
  
  # Save model type and theta globally
  model <<- model
  theta <<- theta
  BF <<- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
  basis <<- lapply(1:K, function(k) {
    Reduce("rbind", lapply(1:S, \(s) BF[s, k] * diag(n)))
  })
  BFTest <- Bsplines_2D(ZTest, df = c(sqrt(K), sqrt(K)))
  basisTest <<- lapply(1:K, function(k) {
    Reduce("rbind", lapply(1:STest, \(s) BFTest[s, k] * diag(nTest)))
  })
  
  # MCMC chain properties
  nIter <- nBurn + nIter # 15 to 20 thousand ideally
  
  # Tuning parameters for variance of each proposal distribution
  # Can be user-supplied
  sdSigma2 <- propSD[[1]]
  sdTau2 <- propSD[[2]]
  
  trSigma2 <- matrix(0, nrow = K, ncol = nIter)
  trTau2 <- numeric(nIter) # Transformed parameters
  beta <- matrix(0, nrow = p, ncol = nIter)
  acceptTau2 <- 0 # Track acceptance rates
  #acceptSigma2 <- rep(0, K)
  acceptSigma2 <- 0

  # Initial values of transformed parameters (except for beta, not transformed)
  trSigma2[, 1] <- rep(log(50), K)
  trTau2[1] <- log(0.1)
  beta[ , 1] <- 0
  
  # Base of covariance matrix for updating sigma2 and tau2 (only need to compute once)
  B <<- baseVariance(theta, D = D)
  Sigma <<- Reduce("+", lapply(1:K, \(k) exp(trSigma2[k, 1]) * B[[k]])) + exp(trTau2[1]) * diag(n * S)
    
  # Base of covariance matrix for predictions
  BTest <- lapply(1:K, \(k) tcrossprod(basisTest[[k]] %*% exp(-theta[k] * DTest), basisTest[[k]]))
  SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
    exp(trSigma2[k, 1]) * BTest[[k]]
  })) + exp(trTau2[1]) * diag(STest * nTest)
  
  # Initial predictions for test subjects
  YPreds <- matrix(data = NA, nrow = nTest * STest, ncol = nIter)
  YPreds[ , 1] <- t(rmvnorm(1, mean = ATest %*% beta[ , 1], sigma = SigmaTest))
  
  # Run Gibbs/Metropolis for one chain
  for (i in 2:nIter) {

    cat(paste0("Beginning iteration ", i, ".\n"))
    
    ### Metropolis update (sigma2) ###
    
    propTrSigma2 <- trSigma2[ , i - 1]
    #for (j in 1:K) {
      #lastTrSigma2 <- propTrSigma2
      #propTrSigma2[j] <- rnorm(1, mean = trSigma2[j , i - 1], sd = sdSigma2)
      #MHratio <- logRatioSigma2(propTrSigma2, 
      #                           lastTrSigma2, 
      #                           trTau2[i - 1],
      #                           beta[ , i - 1])
      
      propTrSigma2 <- rnorm(K, mean = trSigma2[ , i - 1], sd = sdSigma2)
      MHratio <<- logRatioSigma2(propTrSigma2, 
                                 trSigma2[, i - 1], 
                                 trTau2[i - 1],
                                 beta[ , i - 1])
      
      if(runif(1) < exp(MHratio)) {
        #trSigma2[j, i] <- propTrSigma2[j]
        #Sigma <<- SigmaProp
        #acceptSigma2[j] <- acceptSigma2[j] + 1
        trSigma2[, i] <- propTrSigma2
        Sigma <<- SigmaProp
        acceptSigma2 <- acceptSigma2 + 1
      } else {
        #trSigma2[j, i] <- trSigma2[j, i - 1]
        trSigma2[, i] <- trSigma2[, i - 1]
      }
    #}
    
    #cat("Sigma2 updated \n")
    
    ### Metropolis update (tau2) ###
    
    propTrTau2 <- rnorm(1, mean = trTau2[i - 1], sd = sdTau2)
    MHratio <- logRatioTau2(trSigma2[, i - 1], 
                            propTrTau2,
                            trTau2[i - 1],
                            beta[ , i - 1])
    #if (is.na(MHratio)) {
    #  break
    #}
    if (runif(1) < exp(MHratio)) { 
      trTau2[i] <- propTrTau2
      Sigma <<-  SigmaProp
      acceptTau2 <- acceptTau2 + 1
    } else {
      trTau2[i] <- trTau2[i - 1]
    }
    #cat("Tau2 updated \n")
    
    ### Gibbs update (beta) ###
    
    SigmaInv <- solve(Sigma)
    SigmaBeta <- solve(crossprod(A, SigmaInv %*% A) + 1)
    meanBeta <- SigmaBeta %*% crossprod(A, SigmaInv %*% Y)
    beta[ , i] <- t(rmvnorm(1, meanBeta, SigmaBeta))
    
    #cat("Beta updated \n")
    
    ### Posterior predictive sampling for test subjects ###
    SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
      exp(trSigma2[k, i]) * BTest[[k]]
    })) + exp(trTau2[i]) * diag(STest * nTest)
    YPreds[ , i] <- t(rmvnorm(1, mean = ATest %*% beta[ , i], sigma = SigmaTest))
  }
  #cat(paste0("MH Ratio is ", exp(MHratio), "\n"))
  #cat(paste0("Last TrTau2 was ", trTau2[i-1]), "\n")
  #cat(paste0("Proposed TrTau2 is ", propTrTau2), "\n")
  #cat(paste0("Beta is ", beta[ , i-1]), "\n")
  #return(list(prevTrSigma2 = trSigma2[,i-1], trSigma2 = trSigma2[,i]))
  
  # Acceptance rates (for Metropolis-sampled parameters)
  acceptance <- list(sigma2 = acceptSigma2 / nIter, 
                     tau2 = acceptTau2 / nIter)
  
  # Remove burn-in and perform thinning
  index <- seq(nBurn + 1, nIter, by = nThin)
  trSigma2 <- trSigma2[ , index]
  trTau2 <- trTau2[index]
  beta <- beta[ , index]
  YPreds <- YPreds[ , index]
  nSamples <- length(index)
  
  # Back-transform
  sigma2 <- exp(trSigma2)
  tau2 <- exp(trTau2)
  
  # Trace plots
  #pdf(paste0("../paper/figures/trace_plots/trace_plots_", model, ".pdf"))
  #plot(1:nSamples, sigma2, type = 'l', ylab = "Sigma2", main = "")
  #plot(1:nSamples, tau2, type = 'l', ylab = "Tau2", main = "")
  ##plot(1:nSamples, theta, type = 'l', ylab = "Trace Plot for theta")
  #plot(1:nSamples, beta[1,], type = 'l', ylab = "beta_0", main = "")
  #dev.off()  
  
  if (class(beta) == "numeric") {
    beta <- matrix(beta, nrow = 1)
  }
  
  # Posterior mean estimates (can be somewhat skewed because of back-transformations)
  posteriorMeans <- list(sigma2 = apply(sigma2, 1, mean),
                         tau2 = mean(tau2),
                         beta = apply(beta, 1, mean))
  
  # Posterior median estimates (more accurate)
  posteriorMedians <- list(sigma2 = apply(sigma2, 1, median),
                           tau2 = median(tau2),
                           beta = apply(beta, 1, median))
  
  # 95% credible interval bounds
  credLower <- list(sigma2 = apply(sigma2, 1, quantile, 0.025), 
                    tau2 = quantile(tau2, 0.025),
                    beta = apply(beta, 1, quantile, 0.025))
  credUpper <- list(sigma2 = apply(sigma2, 1, quantile, 0.975), 
                    tau2 = quantile(tau2, 0.975),
                    beta = apply(beta, 1, quantile, 0.975))
  
  # Posterior predictive results for test data
  #preds <- lapply(1:nTestSubj, function(j) {
  #  apply(YPreds[[j]], 1, quantile, c(0.025, 0.5, 0.975))
  #})
  preds <- apply(YPreds, 1, quantile, c(0.025, 0.5, 0.975))
  
  # Return results
  return(list(acceptance = acceptance, 
              posteriorMeans = posteriorMeans,
              posteriorMedians = posteriorMedians,
              credLower = credLower,
              credUpper = credUpper,
              preds = preds,
              predSamples = YPreds,
	      paramSamples = list(sigma2, tau2, beta)))
}


