### LOG POSTERIORS NOT WRITTEN OUT FULLY HERE ###
### Only the components needed to update each parameter

logRatioSigma2 <- function(propTrSigma2, prevTrSigma2, trSigb2, trThb, trTheta, trTau2, beta) {
  propSigma2 <- exp(propTrSigma2)
  prevSigma2 <- exp(prevTrSigma2)
  SigmaProp <<- Sigma + Reduce("+", lapply(1:K, \(k) (propSigma2[k] - prevSigma2[k]) * B[[k]]))
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorSigma2(propSigma2)) - sum(logPriorSigma2(prevSigma2)) + # Log Priors
    sum(propTrSigma2) - sum(prevTrSigma2) # Jacobian is exp()
  # Taking the log and exp of trSigma2 cancels out
}

logRatioSigb2 <- function(propTrSigb2, prevTrSigb2, trThb, trSigma2, trTheta, trTau2, beta) {
  propSigb2 <- exp(propTrSigb2)
  prevSigb2 <- exp(prevTrSigb2)
  DBNew <<- lapply(1:(p+1), \(j) matrix(X0[ , j], nrow = n, ncol = n) * 
                 (propSigb2[j] * exp(-gInv(trThb[j]) * D)) *
                 matrix(X0[ , j], nrow = n, ncol = n, byrow = T))
  CBNew <<- diag(S) %x% Reduce("+", DBNew)
  SigmaProp <<- Sigma - CBFull + CBNew
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorSigma2(propSigb2)) - sum(logPriorSigma2(prevSigb2)) + # Log Priors
    sum(propTrSigb2) - sum(prevTrSigb2) # Jacobian is exp()
  # Taking the log and exp of trSigb2 cancels out
}

logRatioTau2 <- function(propTrTau2, prevTrTau2, trSigb2, trThb, trSigma2, trTheta, beta) {
  propTau2 <- exp(propTrTau2)
  prevTau2 <- exp(prevTrTau2)
  SigmaProp <<- Sigma + (propTau2 - prevTau2) * diag(n * S)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    logPriorTau2(propTau2) - logPriorTau2(prevTau2) + # Log Priors
    propTrTau2 - prevTrTau2 # Jacobian is exp()
  # Taking the log and exp of trTau2 cancels out
}

logRatioThb <- function(propTrThb, prevTrThb, trSigma2, trTheta, trSigb2, trTau2, beta) {
  propThb <- gInv(propTrThb)
  prevThb <- gInv(prevTrThb)
  #SigmaProp <<- Sigma - exp(trSigb2) * exp(-prevThb * DXFull) + exp(trSigb2) * exp(-propThb * DXFull)
  DBNew <<- lapply(1:(p+1), \(j) matrix(X0[ , j], nrow = n, ncol = n) * 
                 (exp(trSigb2)[j] * exp(-propThb[j] * D)) *
                 matrix(X0[ , j], nrow = n, ncol = n, byrow = T))
  CBNew <<- diag(S) %x% Reduce("+", DBNew)
  SigmaProp <<- Sigma - CBFull + CBNew
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorTheta(propThb)) - sum(logPriorTheta(prevThb)) + # Log Priors
    sum(logJac(propThb)) - sum(logJac(prevThb)) # Jacobians
}

logRatioTheta <- function(propTrTheta, prevTrTheta, trSigma2, trSigb2, trThb, trTau2, beta) {
  propTheta <- gInv(propTrTheta)
  prevTheta <- gInv(prevTrTheta)
  BProp <<- baseVariance(propTheta, D)
  SigmaProp <<- CBFull+
    Reduce("+", lapply(1:K, \(k) exp(trSigma2[k]) * BProp[[k]])) + 
    exp(trTau2) * diag(n * S)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorTheta(propTheta)) - sum(logPriorTheta(prevTheta)) + # Log Priors
    sum(logJac(propTheta)) - sum(logJac(prevTheta)) # Jacobians
}
