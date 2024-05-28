### LOG POSTERIOR ###

logPost <- function(trSigma2, trTau2, theta, beta) {
  logLik(exp(trSigma2), exp(trTau2), theta, beta) + # Log Likelihood
    logPriorSigma2(exp(trSigma2)) + # Log Priors
    logPriorTau2(exp(trTau2)) + 
    logPriorbeta(beta) + 
    trTau2 + trSigma2 # + logJac(trTheta) # Jacobians
  # Taking the log and exp of trTau2 and trSigma2 cancels out
}

logRatioSigma2 <- function(propTrSigma2, prevTrSigma2, trTau2, beta) {
  propSigma2 <- exp(propTrSigma2)
  prevSigma2 <- exp(prevTrSigma2)
  SigmaProp <<- Sigma + Reduce("+", lapply(1:K, \(k) (propSigma2[k] - prevSigma2[k]) * B[[k]]))
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    sum(logPriorSigma2(propSigma2)) - sum(logPriorSigma2(prevSigma2)) + # Log Priors
    sum(propTrSigma2) - sum(prevTrSigma2) # Jacobian is exp()
  # Taking the log and exp of trSigma2 cancels out
}

logRatioTau2 <- function(trSigma2, propTrTau2, prevTrTau2, beta) {
  propTau2 <- exp(propTrTau2)
  prevTau2 <- exp(prevTrTau2)
  SigmaProp <<- Sigma + (propTau2 - prevTau2) * diag(n * S)
  
  logLik(SigmaProp, beta) - logLik(Sigma, beta) + # Log Likelihoods
    logPriorTau2(propTau2) - logPriorTau2(prevTau2) + # Log Priors
    propTrTau2 - prevTrTau2 # Jacobian is exp()
  # Taking the log and exp of trTau2 cancels out
}

