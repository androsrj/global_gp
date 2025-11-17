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
var.eta <- function(sigma2, theta, D, BF) {
  C.list <- lapply(1:K, function(k) {
    sigma2[k] * (1 + sqrt(3) * theta[k] * D) * exp(-sqrt(3) * theta[k] * D)
  })
  C.array <- simplify2array(C.list)  # array of dim n x n x K
  C.array <- aperm(C.array, c(3, 1, 2))  # now K x n x n
  W.list <- lapply(1:K, function(k) tcrossprod(BF[, k]))  # each S x S
  W.array <- simplify2array(W.list)  # S x S x K
  C.eta <- Reduce('+', lapply(1:K, function(k) kronecker(W.array[, , k], C.array[k, , ])))
  return(C.eta)
}


# Convert list of cluster labels to vector of indices (for data subsetting)
list2Vec <- function(ls) {
  temp <- rep(seq_along(ls), lengths(ls))
  temp[unlist(ls)]
}

