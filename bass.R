library(BASS)
load("data/train.RData")
load("data/test.RData")

Ybass <- matrix(c(train$Y), nrow = nrow(train$Z))
model <- bassPCA(train$Z, Ybass, n.pc = 3, n.cores = 1)
#model <- bassPCA(inputs[-stormsTest, ], out[-stormsTest, ], n.pc = 3, n.cores = 1)
predictions <- predict(model, test$Z)[ , , test$index]

bassPreds <- c(apply(predictions, 2:3, mean))
bassLower <- apply(predictions, 2:3, quantile, .025)
bassUpper <- apply(predictions, 2:3, quantile, .975)
bassCRPS <- mean(sapply(1:nTestSubj, function(i) {
  truth <- out[stormsTest[i], indexTest]
  preds <- predictions[i, , ]
  mean(energy_score(truth, preds))
}))