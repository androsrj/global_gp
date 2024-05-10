### Wasserstein mean for multiple Markov chains (multiple distributions) 
wasserstein <- function(results, time) {

  wassersteinAcc <- rowMeans(sapply(results, \(x) unlist(x$acceptance)))
  wassersteinMeans <- rowMeans(sapply(results, \(x) unlist(x$posteriorMedians)))
  wassersteinLower <- rowMeans(sapply(results, \(x) unlist(x$credLower)))
  wassersteinUpper <- rowMeans(sapply(results, \(x) unlist(x$credUpper)))
  
  predictions <- vector("list", nTestSubj)
  for (i in 1:nTestSubj) {
    predsList <- lapply(results, \(x) x$preds[[i]])
    predictions[[i]] <- Reduce("+", predsList) / length(predsList)
  }
    
  wassersteinResults <- list(acc = wassersteinAcc, 
                             means = wassersteinMeans, 
                             lower = wassersteinLower,
                             upper = wassersteinUpper, 
                             predictions = predictions,
                             time = time)
}

# Energy score (CRPS) calculation for predictions vs truth
energy_score = function(y, z){
  d = dim(z)
  n_samp = d[1]
  n = d[2]
  G = colMeans(sapply(1:n,function(i) sapply(1:n_samp, function(k) norm(z[k,i]-y[i],type='2'))))
  UQ = sapply(1:n,function(i) (1/(2*n_samp^2))*sum(as.matrix(distances::distances(z[,i]))))
  return(G - UQ)
}


# Calculates the base of the covariance matrix for likelihood function
baseVariance <- function(theta, D) {
  B <- lapply(1:K, \(k) tcrossprod(basis[[k]] %*% exp(-theta[k] * D), basis[[k]]))
  return(B)
  
}


# Convert list of cluster labels to vector of indices (for data subsetting)
list2Vec <- function(ls) {
  temp <- rep(seq_along(ls), lengths(ls))
  temp[unlist(ls)]
}

