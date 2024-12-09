# SOURCES
source("other_functions/helper_functions.R") 
library(fields)
library(ggplot2)
library(spBayes)
library(BASS)

load("data/flood_data.RData")

mySeed <- 1234
n <- 100
nTest <- 25
S <- 10
STest <- 10

set.seed(mySeed)
which.storms <- sample(4000, S + STest)
train.storms <- which.storms[1:S]
test.storms <- which.storms[(S+1):(S+STest)]

coords.subset <- coords[coords$x < -74.81 & coords$x > -74.83 & coords$y < 39.08 & coords$y > 39.06, ]
source("coastlines.R")

set.seed(mySeed)
which.points <- sample(nrow(coords.subset), n + nTest)
train.index <- which.points[1:n]
test.index <- which.points[(n+1):(n+nTest)]

inputs <- inputs[c(train.storms, test.storms), ]
out <- out[c(train.storms, test.storms), ]
gc()

strt <- Sys.time()
set.seed(mySeed)
model <- bassPCA(inputs[1:10, ], out[1:10, ], n.pc = 3, n.cores = 1)
predictions <- predict(model, inputs[11:20, ])[ , , test.index]

bassPreds <- apply(predictions, 2:3, mean)
bassLower <- apply(predictions, 2:3, quantile, .025)
bassUpper <- apply(predictions, 2:3, quantile, .975)
bassCRPS <- mean(sapply(1:STest, function(i) {
  truth <- out[(11:20)[i], test.index]
  preds <- predictions[i, , ]
  mean(energy_score(truth, preds))
}))
bassTime <- Sys.time() - strt

bass <- list(preds = bassPreds, 
             lower = bassLower, 
             upper = bassUpper,
             time = bassTime,
             crps = bassCRPS)
truth <- as.matrix(out[11:20, test.index])

rmse <- sqrt(mean((bassPreds - truth)^2))
cvg <- mean(bassLower < truth & bassUpper > truth)
len <- mean(bassUpper - bassLower)

cat(paste0("YTest Standard Dev: ", round(sd(truth), 3), "\n"))
cat(paste0("Root MS error: ", round(mean(rmse), 3), "\n"))
cat(paste0("Mean coverage: ", round(mean(cvg), 3), "\n"))
cat(paste0("Mean width: ", round(mean(len), 3), "\n"))
