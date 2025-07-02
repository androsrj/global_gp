### FINAL MCMC SKETCHING FUNCTION ###
# Priors, Jacobians, likelihood, and posterior must already be sourced

mcmc <- function(X, Z, Y, D, K,
                 starting,
                 propSD,
                 model = c("full_gp", "mpp", "sparse_gpp")[1],
                 nIter = 1000, nBurn = 100, nThin = 2) {
    
  # Dimensions
  S <<- nrow(Z)
  n <<- length(Y) / S
  STest <<- nrow(ZTest)
  nTest <<- length(YTest) / STest
  K <<- K
  #J <<- matrix(1, nrow = S, ncol = 1)
  #JTest <<- matrix(1, nrow = STest, ncol = 1)
  A <<- rep(1, S) %x% cbind(matrix(1, nrow = n, ncol = 1), X)
  ATest <<- rep(1, STest) %x% cbind(matrix(1, nrow = nTest, ncol = 1), XTest)
  p <<- ncol(X)
  X0 <<- cbind(rep(1, n), X)
  X0Test <<- cbind(rep(1, nTest), XTest)
  #DXFull <<- matrix(1, S, S) %x% rdist(scale(X))
  #DXTestFull <<- matrix(1, STest, STest) %x% rdist(scale(XTest))
  DB <- lapply(1:(p+1), \(j) matrix(X0[ , j], nrow = n, ncol = n) * 
                 (starting$sigb2[j] * exp(-starting$thb[j] * D)) *
                 matrix(X0[ , j], nrow = n, ncol = n, byrow = T))
  CBFull <<- matrix(1, S, S) %x% Reduce("+", DB)
  DBTest <- lapply(1:(p+1), \(j) matrix(X0Test[ , j], nrow = nTest, ncol = nTest) * 
                     (starting$sigb2[j] * exp(-starting$thb[j] * DTest)) *
                     matrix(X0Test[ , j], nrow = nTest, ncol = nTest, byrow = T))
  CBTestFull <<- matrix(1, STest, STest) %x% Reduce("+", DBTest)
  
  # Save model type and theta globally
  model <<- model
  #theta <<- theta
  BF <- Bsplines_2D(Z, df = c(sqrt(K), sqrt(K)))
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
  sdSigb2 <- propSD$sigb2
  sdThb <- propSD$thb
  sdSigma2 <- propSD$sigma2
  sdTheta <- propSD$theta
  sdTau2 <- propSD$tau2
  
  # Initialize vectors for MCMC
  trSigma2 <- trTheta <- matrix(0, nrow = K, ncol = nIter)
  trThb <- trSigb2 <- matrix(0, nrow = p + 1, ncol = nIter)
  trTau2 <- numeric(nIter) # Transformed parameters
  beta <- matrix(0, nrow = p+1, ncol = nIter)
  acceptTau2 <- 0 # Track acceptance rates
  acceptSigma2 <- 0
  acceptSigb2 <- 0
  acceptTheta <- 0
  acceptThb <- 0

  # Initial values of transformed parameters (except for beta, not transformed)
  trSigma2[, 1] <- log(starting$sigma2)
  trTheta[, 1] <- g(starting$theta)
  trSigb2[, 1] <- log(starting$sigb2)
  trThb[, 1] <- g(starting$thb)
  trTau2[1] <- log(starting$tau2)
  beta[ , 1] <- starting$beta
  
  # Base of covariance matrix for updating sigma2 and tau2
  B <<- baseVariance(theta = starting$theta, D = D)
  Sigma <<- Reduce("+", lapply(1:K, \(k) starting$sigma2[k] * B[[k]])) + 
    CBFull +
    #starting$sigb2 * exp(-starting$thb * DXFull) + 
    starting$tau2 * diag(n * S)
  
  # Base of covariance matrix for predictions
  BTest <- lapply(1:K, \(k) tcrossprod(basisTest[[k]] %*% exp(-starting$theta[k] * DTest), basisTest[[k]]))
  SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
    starting$sigma2[k] * BTest[[k]]
  })) + 
    CBTestFull +
    #starting$sigb2 * exp(-starting$thb * DXTestFull) + 
    starting$tau2 * diag(STest * nTest)
  
  # Initial predictions for test subjects
  YPreds <- matrix(data = NA, nrow = nTest * STest, ncol = nIter)
  YPreds[ , 1] <- t(rmvnorm(1, mean = ATest %*% beta[ , 1], sigma = SigmaTest))
  
  # Run Gibbs/Metropolis for one chain
  for (i in 2:nIter) {
    
    if (i %% 10 == 0) {
      cat(paste0("Beginning iteration ", i, ".\n"))
    }

    ### Metropolis update (sigma2_b) ###
    propTrSigb2 <- rnorm(p + 1, mean = trSigb2[ , i - 1], sd = sdSigb2)
    MHratio <- logRatioSigb2(propTrSigb2, 
                             trSigb2[ , i - 1],
                             trThb[ , i - 1],
                             trSigma2[ , i - 1],
                             trTheta[ , i - 1],
                             trTau2[i - 1],
                             beta[ , i - 1])
    #cat(paste0("MH Ratio is ", round(exp(MHratio), 2), "\n"))
    if(runif(1) < exp(MHratio)) {
      trSigb2[ , i] <- propTrSigb2
      Sigma <<- SigmaProp
      DB <<- DBNew
      CBFull <<- CBNew
      acceptSigb2 <- acceptSigb2 + 1
    } else {
      trSigb2[ , i] <- trSigb2[ , i - 1]
    }
    
    #cat(paste0("finished sigb2: ", round(exp(trSigb2[i]), 2), "\n"))
    #cat(paste0("Log likelihood is ", round(logLik(Sigma, beta[ , i-1]), 3), "\n"))
    
    ### Metropolis update (theta_b) ###
    propTrThb <- rnorm(p + 1, mean = trThb[ , i - 1], sd = sdThb)
    MHratio <- logRatioThb(propTrThb, 
                           trThb[ , i - 1],
                           trSigma2[ , i - 1],
                           trTheta[ , i - 1],
                           trSigb2[ , i],
                           trTau2[i - 1],
                           beta[ , i - 1])
    
    if(runif(1) < exp(MHratio)) {
      trThb[ , i] <- propTrThb
      Sigma <<- SigmaProp
      DB <<- DBNew
      CBFull <<- CBNew
      acceptThb <- acceptThb + 1
    } else {
      trThb[ , i] <- trThb[ , i - 1]
    }
    
    #cat(paste0("finished thb: ", round(gInv(trThb[i]), 2), "\n"))
    #cat(paste0("Log likelihood is ", round(logLik(Sigma, beta[ , i-1]), 3), "\n"))
    
    ### Metropolis update (sigma2) ###
    propTrSigma2 <- rnorm(K, mean = trSigma2[ , i - 1], sd = sdSigma2)
    MHratio <- logRatioSigma2(propTrSigma2, 
                               trSigma2[ , i - 1], 
                               trSigb2[i],
                               trThb[i],
                               trTheta[ , i - 1],
                               trTau2[i - 1],
                               beta[ , i - 1])
      
    if(runif(1) < exp(MHratio)) {
      trSigma2[, i] <- propTrSigma2
      Sigma <<- SigmaProp
      acceptSigma2 <- acceptSigma2 + 1
    } else {
      trSigma2[, i] <- trSigma2[, i - 1]
    }
    
    #cat("finished sigma2: ")
    #cat(round(exp(trSigma2[,i]), 2))
    #cat("\n")
    #cat(paste0("Log likelihood is ", round(logLik(Sigma, beta[ , i - 1]), 3), "\n"))
    
    ### Metropolis update (theta) ###
    propTrTheta <- rnorm(K, mean = trTheta[ , i - 1], sd = sdTheta)
    MHratio <- logRatioTheta(propTrTheta,
                              trTheta[ , i - 1],
                              trSigma2[ , i], 
                              trSigb2[i],
                              trThb[i],
                              trTau2[i - 1],
                              beta[ , i - 1])
    
    if(runif(1) < exp(MHratio)) {
      trTheta[, i] <- propTrTheta
      Sigma <<- SigmaProp
      B <<- BProp
      acceptTheta <- acceptTheta + 1
    } else {
      trTheta[, i] <- trTheta[, i - 1]
    }
    
    #cat("finished theta: ")
    #cat(round(gInv(trTheta[,i]), 2))
    #cat("\n")
    #cat(paste0("Log likelihood is ", round(logLik(Sigma, beta[ , i - 1]), 3), "\n"))
    
    ### Metropolis update (tau2) ###
    
    propTrTau2 <- rnorm(1, mean = trTau2[i - 1], sd = sdTau2)
    #cat(paste0("proposed tau2 is ", exp(propTrTau2)))
    MHratio <- logRatioTau2(propTrTau2,
                            trTau2[i - 1],
                            trSigb2[i], 
                            trThb[i], 
                            trSigma2[ , i], 
                            trTheta[ , i], 
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
    #cat(paste0("finished tau: ", round(exp(trTau2[i]), 2), "\n"))
    #cat(paste0("Log likelihood is ", round(logLik(Sigma, beta[ , i - 1]), 3), "\n"))
    
    ### Gibbs update (beta) ###
    
    SigmaInv <- solve(Sigma)
    SigmaBeta <- solve(crossprod(A, SigmaInv %*% A) + diag(p+1))
    meanBeta <- SigmaBeta %*% crossprod(A, SigmaInv %*% Y)
    beta[ , i] <- t(rmvnorm(1, meanBeta, SigmaBeta))
    #beta[i] <-  mean(Y)
    
    #cat(paste0("finished beta: ", round(beta[ , i], 2), "\n"))
    #cat(paste0("Log likelihood is ", round(logLik(Sigma, beta[ , i - 1]), 3), "\n"))
    #cat("beta updated \n")
    
    ### Posterior predictive sampling for test subjects ###
    #SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
    #  exp(trSigma2[k, i]) * BTest[[k]]
    #})) + exp(trTau2[i]) * diag(STest * nTest) + 
    #  exp(-gInv(trThb[i]) * DXTestFull)
    
    DBTest <- lapply(1:(p+1), \(j) matrix(X0Test[ , j], nrow = nTest, ncol = nTest) * 
                       (exp(trSigb2[j]) * exp(-gInv(trThb[j]) * DTest)) *
                       matrix(X0Test[ , j], nrow = nTest, ncol = nTest, byrow = T))
    CBTestFull <<- matrix(1, STest, STest) %x% Reduce("+", DBTest)
    SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
      exp(trSigma2[k, i]) * BTest[[k]]
    })) + CBTestFull + exp(trTau2[i]) * diag(STest * nTest)
    YPreds[ , i] <- t(rmvnorm(1, mean = ATest %*% beta[ , i], sigma = SigmaTest))
  }
  #cat(paste0("MH Ratio is ", exp(MHratio), "\n"))
  #cat(paste0("Last TrTau2 was ", trTau2[i-1]), "\n")
  #cat(paste0("Proposed TrTau2 is ", propTrTau2), "\n")
  #cat(paste0("beta is ", beta[ , i-1]), "\n")
  #return(list(prevTrSigma2 = trSigma2[,i-1], trSigma2 = trSigma2[,i]))
  
  # Acceptance rates (for Metropolis-sampled parameters)
  acceptance <- list(sigb2 = acceptSigb2 / nIter,
                     thb = acceptThb / nIter,
                     sigma2 = acceptSigma2 / nIter, 
                     theta = acceptTheta / nIter,
                     tau2 = acceptTau2 / nIter)
  
  # Remove burn-in and perform thinning
  index <- seq(nBurn + 1, nIter, by = nThin)
  trSigb2 <- trSigb2[index]
  trThb <- trThb[index]
  trSigma2 <- trSigma2[ , index]
  trTheta <- trTheta[ , index]
  trTau2 <- trTau2[index]
  beta <- beta[ , index]
  YPreds <- YPreds[ , index]
  nSamples <- length(index)
  
  # Back-transform
  sigb2 <- exp(trSigb2)
  thb <- gInv(trThb)
  sigma2 <- exp(trSigma2)
  theta <- gInv(trTheta)
  tau2 <- exp(trTau2)
  
  # Trace plots
  #pdf(paste0("../paper/figures/trace_plots/trace_plots_", model, ".pdf"))
  #plot(1:nSamples, sigma2, type = 'l', ylab = "Sigma2", main = "")
  #plot(1:nSamples, tau2, type = 'l', ylab = "Tau2", main = "")
  ##plot(1:nSamples, theta, type = 'l', ylab = "Trace Plot for theta")
  #plot(1:nSamples, beta[1,], type = 'l', ylab = "beta_0", main = "")
  #dev.off()  
  
  #if (class(beta) == "numeric") {
  #  beta <- matrix(beta, nrow = 1)
  #}
  
  # Posterior mean estimates (can be somewhat skewed because of back-transformations)
  posteriorMeans <- list(sigb2 = apply(sigb2, 1, mean),
                         thb = apply(thb, 1, mean),
                         sigma2 = apply(sigma2, 1, mean),
                         theta = apply(theta, 1, mean),
                         tau2 = mean(tau2),
                         beta = apply(beta, 1, mean))
  
  # Posterior median estimates (more accurate)
  posteriorMedians <- list(sigb2 = apply(sigb2, 1, median),
                           thb = apply(thb, 1, median),
                           sigma2 = apply(sigma2, 1, median),
                           theta = apply(theta, 1, median),
                           tau2 = median(tau2),
                           beta = apply(beta, 1, median))
  
  # 95% credible interval bounds
  credLower <- list(sigb2 = apply(sigb2, 1, quantile, 0.025),
                    thb = apply(thb, 1, quantile, 0.025),
                    sigma2 = apply(sigma2, 1, quantile, 0.025),
                    theta = apply(theta, 1, quantile, 0.025),
                    tau2 = quantile(tau2, 0.025),
                    beta = apply(beta, 1, quantile, .025))
  credUpper <- list(sigb2 = apply(sigb2, 1, quantile, 0.975),
                    thb = apply(thb, 1, quantile, 0.975),
                    sigma2 = apply(sigma2, 1, quantile, 0.975),
                    theta = apply(theta, 1, quantile, 0.975),
                    tau2 = quantile(tau2, 0.975),
                    beta = apply(beta, 1, quantile, .975))
  
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
              paramSamples = list(sigb2 = sigb2, 
                                  thb = thb, 
                                  sigma2 = sigma2, 
                                  theta = theta, 
                                  tau2 = tau2, 
                                  beta = beta)))
}


