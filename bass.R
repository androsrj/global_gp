library(BASS)
source("other_functions/helper_functions.R")
load("data/train.RData")
load("data/test.RData")

S <- nrow(train$Z) 
STest <- nrow(test$Z)
n <- nrow(train$X)
nTest <- nrow(test$X)

Ybass <- matrix(c(train$Y), nrow = nrow(train$Z)) 
model <- bassPCA(train$Z, Ybass, n.pc = 3, n.cores = 1)
#model <- bassPCA(inputs[-stormsTest, ], out[-stormsTest, ], n.pc = 3, n.cores = 1)
predictions <- predict(model, test$Z)[ , , test$index]

bassPreds <- c(apply(predictions, 2:3, mean))
bassLower <- apply(predictions, 2:3, quantile, .025)
bassUpper <- apply(predictions, 2:3, quantile, .975)
bassCRPS <- mean(sapply(1:STest, function(i) {
  truth <- test$Y[((i-1)*nTest+1):(nTest*i), ]
  preds <- predictions[i, , ]
  mean(energy_score(truth, preds))
}))

#sqrt(mean((bassPreds - test$Y)^2))

nTestSubj <- nrow(test$Z)
abs_error <- cvg <- width <- scores <- numeric(nTestSubj)
a <- .05
for (i in 1:nTestSubj) {
  truth <- test$Y[(nTest*(i-1)+1):(i*nTest)]
  pred <- bassPreds[(nTest*(i-1)+1):(i*nTest)]
  abs_error[i] <- mean(abs(truth - pred))
  lower <- bassLower[(nTest*(i-1)+1):(i*nTest)]
  upper <- bassUpper[(nTest*(i-1)+1):(i*nTest)]
  cvg[i] <- mean(lower < truth & upper > truth)
  width[i] <- mean(upper - lower)
  scores[i] <- mean((upper - lower) +
                     2/a * (lower - truth) * (truth < lower) +
                     2/a * (truth - upper) * (truth > upper))
  #predSamples <- t(results$predSamples[(nTest*(i-1)+1):(i*nTest), ])
  #crps[i] <- mean(energy_score(truth, predSamples))
}

abs_error
cat(paste0("Mean absolute error: ", round(mean(abs_error), 3), "\n"))

cvg
cat(paste0("Mean coverage: ", round(mean(cvg), 3), "\n"))

width
cat(paste0("Mean width: ", round(mean(width), 3), "\n"))

scores
cat(paste0("Mean interval score: ", round(mean(scores), 3), "\n"))

bassCRPS
cat(paste0("Mean CRPS: ", round(mean(bassCRPS), 3), "\n"))
