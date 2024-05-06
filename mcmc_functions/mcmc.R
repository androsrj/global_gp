### FINAL MCMC SKETCHING FUNCTION ###
# Priors, Jacobians, likelihood, and posterior must already be sourced

mcmc <- function(X, Y, D, S, K,
                 theta,
                 model = c("full_gp", "mpp", "sparse_gpp")[1],
                 propSD = c(0.1, 0.1), 
                 nIter = 1000, nBurn = 100, nThin = 2,
                 mProp = 0.1, transform = TRUE) {
  
  # Save model type and theta globally
  model <<- model
  theta <<- theta
  basis <- Bsplines_2D(X, df = c(sqrt(K), sqrt(K)))
  basis <<- lapply(1:K, \(k) diag(basis[ , k]))
  basisTest <<- Bsplines_2D(XTest, df = c(sqrt(K), sqrt(K)))
  #basisTest <<- lapply(1:K, \(k) diag(basisTest[ , k]))
  
  #nSubj <- ifelse(class(X) == "list", length(X), 1)
  #subjs <<- 1:nSubj
  #nTestSubj <- length(test_subjects)
  
  # Dimensions
  n <- nrow(X)
  p <- ncol(X)
  
  # Generate phi and compress data, if desired
  if (transform == TRUE) {
    m <<- round(mProp * n)
    phi <<- matrix(rnorm(m * n, 0, 1 / sqrt(n)), nrow = m, ncol = n)
    newY <<- lapply(Y, \(y) phi %*% y) 
    newX <<- lapply(X, \(x) phi %*% x) 
  } else {
    m <<- n
    phi <<- diag(m)
    newY <<- Y
    newX <<- X
  }
  
  # MCMC chain properties
  nIter <- nBurn + nIter # 15 to 20 thousand ideally
  
  # Tuning parameters for variance of each proposal distribution
  # Can be user-supplied
  sdSigma2 <- propSD[[1]]
  sdTau2 <- propSD[[2]]
  #sdTheta <- propSD[3]
  
  #trSigma2 <- trTau2 <- trTheta <- numeric(nIter) # Transformed parameters
  trSigma2 <- matrix(0, nrow = K, ncol = nIter)
  trTau2 <- numeric(nIter) # Transformed parameters
  beta <- matrix(0, nrow = p, ncol = nIter) # Beta
  acceptSigma2 <- acceptTau2 <- 0 # Track acceptance rates
  
  # Initial values of transformed parameters (except for beta, not transformed)
  trSigma2[, 1] <- log(1:K)
  trTau2[1] <- log(0.2)
  beta[ , 1] <- rep(0, p)
  
  # Base of covariance matrix for updating sigma2 (only need to compute once)
  B <<- baseVariance(theta, D = D)
  #Sigma <<- exp(trSigma2[1]) * B + exp(trTau2[1]) * diag(m)
  Sigma <<- Reduce("+", lapply(1:K, \(k) exp(trSigma2[k, 1] * B[[k]]))) + exp(trTau2[1]) * diag(m)
    
  # Initial predictions for storm 1 (and non-transformed covariance matrix)
  BTest <<- lapply(1:K, function(k) {
    diag(basisTest[ , k]) %*% exp(-theta[k] * DTest) %*% diag(basisTest[ , k])
  })
  SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
    exp(trSigma2[k, 1]) * BTest[[k]]
  })) + exp(trTau2[1]) * diag(nTest)
  YPreds <- matrix(data = NA, nrow = nTest, ncol = nIter)
  YPreds[ , 1] <- t(rmvnorm(1, mean = as.vector(XTest %*% beta[ , 1]), sigma = SigmaTest))
  
  # Run Gibbs/Metropolis for one chain
  for (i in 2:nIter) {

    cat(paste0("Beginning iteration ", i, ".\n"))
    
    ### Metropolis update (sigma2) ###
    
    propTrSigma2 <- rnorm(K, mean = trSigma2[ , i - 1], sd = sdSigma2)
    MHratio <<- logRatioSigma2(propTrSigma2, 
                              trSigma2[, i - 1], 
                              trTau2[i - 1],
                              beta[ , i - 1])
    
    if(runif(1) < exp(MHratio)) {
      trSigma2[, i] <- propTrSigma2
      Sigma <<- SigmaProp
      acceptSigma2 <- acceptSigma2 + 1
    } else {
      trSigma2[, i] <- trSigma2[, i - 1]
    }
    
    ### Metropolis update (tau2) ###
    
    propTrTau2 <- rnorm(1, mean = trTau2[i - 1], sd = sdTau2)
    MHratio <- logRatioTau2(trSigma2[i - 1], 
                            propTrTau2,
                            trTau2[i - 1],
                            beta[ , i - 1])
    if (runif(1) < exp(MHratio)) { 
      trTau2[i] <- propTrTau2
      Sigma <<-  SigmaProp
      acceptTau2 <- acceptTau2 + 1
    } else {
      trTau2[i] <- trTau2[i - 1]
    }
    
    ### Gibbs update (beta) ###
    
    SigmaInv <- solve(Sigma)
    SigmaBeta <- (n / m) * t(newX) %*% SigmaInv %*% newX + diag(p)
    #SigmaBetaList <- lapply(newX, \(x) t(x) %*% SigmaInv %*% x)
    #SigmaBeta <- (n / m) * (Reduce("+", SigmaBetaList) + diag(p))
    SigmaBetaInv <- solve(SigmaBeta)
    meanBeta <- (n / m) * SigmaBetaInv %*% t(newX) %*% SigmaInv %*% newY
    #meanBetaList <- lapply(subjs, \(i) t(newX[[i]]) %*% SigmaInv %*% newY[[i]])
    #meanBeta <- (n / m) * SigmaBetaInv %*% Reduce("+", meanBetaList)
    beta[ , i] <- t(rmvnorm(1, meanBeta, SigmaBetaInv))
    
    ### Posterior predictive sampling for test subjects ###
    #SigmaTest <- exp(trSigma2[i]) * BTest + exp(trTau2[i]) * diag(nTest)
    SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
      exp(trSigma2[k, i]) * BTest[[k]]
    })) + exp(trTau2[i]) * diag(nTest)
    YPreds[ , i] <- t(rmvnorm(1, mean = as.vector(XTest %*% beta[ , i]), sigma = SigmaTest))
  }
  
  # Acceptance rates (for Metropolis-sampled parameters)
  acceptance <- c(sigma2 = acceptSigma2, 
                  tau2 = acceptTau2) / nIter
  
  # Remove burn-in and perform thinning
  index <- seq(nBurn + 1, nIter, by = nThin)
  trSigma2 <- trSigma2[ , index]
  trTau2 <- trTau2[index]
  #trTheta <- trTheta[index]
  beta <- beta[ , index]
  #YPreds <- lapply(1:nTestSubj, \(j) YPreds[[j]][ , index])
  YPreds <- YPreds[ , index]
  nSamples <- length(index)
  
  # Back-transform
  sigma2 <- exp(trSigma2)
  tau2 <- exp(trTau2)
  #theta <- fInv(trTheta)
  
  # Trace plots
  #pdf(paste0("../paper/figures/trace_plots/trace_plots_", model, ".pdf"))
  #plot(1:nSamples, sigma2, type = 'l', ylab = "Sigma2", main = "")
  #plot(1:nSamples, tau2, type = 'l', ylab = "Tau2", main = "")
  ##plot(1:nSamples, theta, type = 'l', ylab = "Trace Plot for theta")
  #plot(1:nSamples, beta[1, ], type = 'l', ylab = "Beta_1", main = "")
  #plot(1:nSamples, beta[p, ], type = 'l', ylab = "Beta_p", main = "")
  #dev.off()  
  
  # Posterior mean estimates (can be somewhat skewed because of back-transformations)
  posteriorMeans <- list(sigma2 = apply(sigma2, 1, mean),
                         tau2 = mean(tau2),
                         #theta = mean(theta),
                         beta = apply(beta, 1, mean))
  
  # Posterior median estimates (more accurate)
  posteriorMedians <- list(sigma2 = apply(sigma2, 1, median),
                           tau2 = median(tau2),
                           #theta = median(theta),
                           beta = apply(beta, 1, median))
  
  # 95% credible interval bounds
  credLower <- list(sigma2 = apply(sigma2, 1, quantile, 0.025), 
                    tau2 = quantile(tau2, 0.025),
                    #theta = quantile(theta, 0.025), 
                    beta = apply(beta, 1, quantile, 0.025))
  credUpper <- list(sigma2 = apply(sigma2, 1, quantile, 0.975), 
                    tau2 = quantile(tau2, 0.975),
                    #theta = quantile(theta, 0.975), 
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


