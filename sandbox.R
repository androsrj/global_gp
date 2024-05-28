library(fields)
library(MBA)
library(mvtnorm)

source("generate_data.R")
theta <- trueTheta
sigma2 <- trueSigma2
tau2 <- 0.2
K <- 9

pred.surf <-  mba.surf(cbind(test$U, test$Y[1:10]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="True Surface, Subject 1", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)

BFTest <- Bsplines_2D(test$X, df = c(sqrt(K), sqrt(K)))
basisTest <<- lapply(1:K, function(k) {
  Reduce("rbind", lapply(1:STest, \(s) BFTest[s, k] * diag(nTest)))
})
DTest <- test$D
BTest <- lapply(1:K, \(k) tcrossprod(basisTest[[k]] %*% exp(-theta[k] * DTest), basisTest[[k]]))
SigmaTest <<- Reduce("+", lapply(1:K, function(k) {
  sigma2[k] * BTest[[k]]
})) + tau2 * diag(STest * nTest)
YPreds <- rmvnorm(1, mean = rep(8, nTest * STest), sigma = SigmaTest)
#YPreds <- apply(rmvnorm(1000, mean = rep(8, nTest * STest), sigma = SigmaTest), 2, mean)
pred.surf <-  mba.surf(cbind(test$U, YPreds[1:10]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="Predicted Surface, Subject 1", col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)

mean(abs(t(test$Y) - YPreds))
sd(test$Y)
