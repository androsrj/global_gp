library(fields)
library(pracma)

# Function to select neighbors - need to take part of sparse() function
# and put it here.

neighbors <- function(C) {
  ###
}


# New sparsity function

sparse <- function(C, m) {
  
  #C <- exp(- theta * Dmat)
  n <- nrow(C)
  A <- matrix(0, n, n)
  D <- matrix(0, n, n)
  D[1, 1] <- C[1, 1]
  
  # Initial computation of A and D
  for(i in 1:(n - 1)) {
    A[i + 1, 1:i] <- solve(C[1:i, 1:i], C[1:i, i + 1])
    D[i + 1, i + 1] <- C[i + 1,i + 1] - dot(C[i + 1, 1:i], A[i + 1, 1:i])
  }
  
  # Sparse computaton of A and D
  N <- vector(mode = "list", length = n)
  for(i in 1:(n-1)) {
    neighbors <- which(sapply(1:i, \(j) (j <= i) & A[i + 1, j] != 0))
    if (i <= m) {
      N[[i + 1]] <- neighbors
    } else {
      N[[i + 1]] <- which(order(C[neighbors,1]) <= (m+1))[-1]
    }
    nn <- length(N[[i + 1]])
    A[i + 1, N[[i + 1]] ] <- C[i + 1, N[[i + 1]]] %*% 
      solve(C[ N[[i + 1]], N[[i + 1]] ] + tau2 * nn)
    D[i + 1, i + 1] <- C[i + 1, i + 1] + 
      tau2 - 
      C[i + 1, N[[i + 1]]] %*% 
      solve(C[N[[i + 1]], N[[i + 1]]] + tau2 * diag(nn)) %*% 
      C[N[[i + 1]], i + 1]
  }
  
  Ctilde <- solve(diag(n) - A) %*% D %*% solve(t(diag(n) - A))
  
  #if (!isSymmetric(Ctilde)) {
  #  Ctilde <- (Ctilde + t(Ctilde)) / 2
  #}
  
  CtildeInv <- t(diag(n) - A) %*% solve(D) %*% (diag(n) - A)
  
  return(list(Ctilde = Ctilde, CtildeInv = CtildeInv))
  #return(Ctilde)
}

S <- cbind(rnorm(100, 0, 10), rnorm(100, 0, 10))
tau2 <- 0.2
Dis <- fields::rdist(S)
C <- 2 * exp(-0.5 * Dis)
#Ctilde <- sparse(C, m = 10)$Ctilde
Ctilde2 <- sparse2(C, m = 10)$Ctilde
#mean(abs(Ctilde + tau2*diag(100) - Ctilde2))
#mean(abs(C - Ctilde2))


