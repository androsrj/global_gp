### LOG POSTERIOR ###

logPost <- function(trSigf2, trThf, trSigma2, trTheta, trTau2, beta) {
  logLik(exp(trSigma2), exp(trTau2), theta, beta) + # Log Likelihood
    logPriorSigma2(exp(trSigf2)) + # Log Priors
    logPriorTheta(gInv(trThf)) + 
    logPriorSigma2(exp(trSigma2)) +
    logPriorTheta(gInv(trTheta))
    logPriorTau2(exp(trTau2)) + 
    logPriorbeta(beta) + 
    trTau2 + trSigma2 + trSigf2 + logJac(trThf) + logJac(trTheta) # Jacobians
  # Taking the log and exp of trTau2, trSigma2, and trSigf2 cancels out
}

logRatioSigma2 <- function(propTrSigma2, prevTrSigma2, trSigf2, trThf, trTheta, trTau2, beta) {
  propSigma2 <- exp(propTrSigma2)
  prevSigma2 <- exp(prevTrSigma2)
  SigmaProp <<- sparse(Sigma + Reduce("+", lapply(1:K, \(k) (propSigma2[k] - prevSigma2[k]) * B[[k]])), m = m)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorSigma2(propSigma2)) - sum(logPriorSigma2(prevSigma2)) + # Log Priors
    sum(propTrSigma2) - sum(prevTrSigma2) # Jacobian is exp()
  # Taking the log and exp of trSigma2 cancels out
}

logRatioSigf2 <- function(propTrSigf2, prevTrSigf2, trThf, trSigma2, trTheta, trTau2, beta) {
  propSigf2 <- exp(propTrSigf2)
  prevSigf2 <- exp(prevTrSigf2)
  SigmaProp <<- sparse(Sigma + (propSigf2 - prevSigf2) * exp(-gInv(trThf) * DXFull), m = m)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    logPriorSigma2(propSigf2) - logPriorSigma2(prevSigf2) + # Log Priors
    propTrSigf2 - prevTrSigf2 # Jacobian is exp()
  # Taking the log and exp of trSigf2 cancels out
}

logRatioTau2 <- function(propTrTau2, prevTrTau2, trSigf2, trThf, trSigma2, trTheta, beta) {
  propTau2 <- exp(propTrTau2)
  prevTau2 <- exp(prevTrTau2)
  SigmaProp <<- sparse(Sigma + (propTau2 - prevTau2) * diag(n * S), m = m)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    logPriorTau2(propTau2) - logPriorTau2(prevTau2) + # Log Priors
    propTrTau2 - prevTrTau2 # Jacobian is exp()
  # Taking the log and exp of trTau2 cancels out
}

logRatioThf <- function(propTrThf, prevTrThf, trSigma2, trTheta, trSigf2, trTau2, beta) {
  propThf <- gInv(propTrThf)
  prevThf <- gInv(prevTrThf)
  #SigmaProp <<- SigmaK + exp(trSigf2) * exp(-propThf * D)
  SigmaProp <<- sparse(Sigma - exp(trSigf2) * exp(-prevThf * DXFull) + exp(trSigf2) * exp(-propThf * DXFull), m = m)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    logPriorTheta(propThf) - logPriorTheta(prevThf) + # Log Priors
    logJac(propThf) - logJac(prevThf) # Jacobians
}

logRatioTheta <- function(propTrTheta, prevTrTheta, trSigma2, trSigf2, trThf, trTau2, beta) {
  propTheta <- gInv(propTrTheta)
  prevTheta <- gInv(prevTrTheta)
  BProp <<- baseVariance(propTheta, D)
  SigmaProp <<- exp(trSigf2) * exp(-gInv(trThf) * DXFull) +
    Reduce("+", lapply(1:K, \(k) exp(trSigma2[k]) * BProp[[k]])) + 
    exp(trTau2) * diag(n * S)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorTheta(propTheta)) - sum(logPriorTheta(prevTheta)) + # Log Priors
    sum(logJac(propTheta)) - sum(logJac(prevTheta)) # Jacobians
}
