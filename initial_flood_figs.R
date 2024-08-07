library(MBA)
library(fields)
library(geoR)
source('other_functions/helper_functions.R')
load("data/flood_data.RData")

# Global covariates (Z)
covars <- c(3, 4)
STrain <- 10
STest <- 5
#set.seed(123)
#samp <- sample(1:nrow(inputs), STrain + STest)
samp <- 51:65
stormsTrain <- samp[1:STrain]
stormsTest <- samp[(STrain+1):(STrain+STest)]
Z <- inputs[stormsTrain, covars]
ZTest <- inputs[stormsTest, covars]

# Local covariates (X)
n <- 500
nTest <- 20
samp2 <- sample(1:nrow(coords))
#train <- samp2[1:n]
train <- which(coords$x > -74.86 & coords$x < -74.83 & coords$y > 39.15 & coords$y < 39.175)
n <- length(train)
length(train)
#test <- samp2[(n+1):(n+nTest)]
test <- sample(train, nTest)
X <- matrix(c(rep(1, n), coords$elev_meters[train]), ncol=2)
XTest <- matrix(c(rep(1, nTest), coords$elev_meters[test]), ncol=2)

# Distance matrices (D)
U <- coords[train, 1:2]
UTest <- coords[test, 1:2]
D <- rdist(U)
DTest <- rdist(UTest)

# Response (Y)
Y <- matrix(as.matrix(out[stormsTrain, train]), ncol = 1)
YTest <- matrix(as.matrix(out[stormsTest, test]), ncol = 1)

par(mfrow=c(1,2))
# Variogram (locations)
#full_U <- as.matrix(U) %x% matrix(1, nrow = STrain, ncol = 1)
#vg_loc <- variog(coords = full_U, data = Y)
#vg_model <- lm(vg_loc$v ~ vg_loc$u)
#plot(vg_loc, main = "Locations")
#abline(a=coef(vg_model)[1], b=coef(vg_model)[2], 
#       col = "brown1", lwd=2, lty=2)

# Variogram (global covariates)
full_Z <- as.matrix(Z) %x% matrix(1, nrow = n, ncol = 1)
vg_global <- variog(coords = full_Z, data = Y)
vg_model2 <- lm(vg_global$v ~ vg_global$u)
plot(vg_global, main = "Training Set 2, Global Covariates")
abline(a=coef(vg_model2)[1], b=coef(vg_model2)[2], 
       col = "brown1", lwd=2, lty=2)


### BASS
library(BASS)
set.seed(92486)
model <- bassPCA(inputs[stormsTrain, ], out[stormsTrain, ], n.pc = 3, n.cores = 1)
predictions <- predict(model, inputs[stormsTest, ])[ , , test]
bassPreds <- apply(predictions, 2:3, mean)
# RMSE
sqrt(mean((c(bassPreds) - YTest)^2))
# CRPS
bassLower <- apply(predictions, 2:3, quantile, .025)
bassUpper <- apply(predictions, 2:3, quantile, .975)
bassCRPS <- mean(sapply(1:STest, function(i) {
  truth <- out[stormsTest[i], test]
  preds <- predictions[i, , ]
  mean(energy_score(truth, preds))
}))
bassCRPS

# Plot
pred.surf <-  mba.surf(cbind(UTest, bassPreds[(nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="BASS, Subject 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)

pred.surf <-  mba.surf(cbind(UTest, YTest[(nTest+1):(2*nTest), ]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="True Subject 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
