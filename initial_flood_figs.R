library(MBA)
library(fields)
library(geoR)
load("data/flood_data.RData")

# Global covariates (Z)
covars <- c(3, 4)
STrain <- 10
STest <- 5
set.seed(123)
samp <- sample(1:nrow(inputs), STrain + STest)
stormsTrain <- samp[1:STrain]
stormsTest <- samp[(STrain+1):(STrain+STest)]
Z <- inputs[stormsTrain, covars]
ZTest <- inputs[stormsTest, covars]

# Local covariates (X)
n <- 500
nTest <- 20
samp2 <- sample(1:nrow(coords))
#train <- samp2[1:n]
#train <- which(coords$x > -74.85 & coords$x < -74.83 & coords$y > 39.1 & coords$y < 39.125)
train <- which(coords$x > -74.85 & coords$x < -74.83 & coords$y > 39.075 & coords$y < 39.1)
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
full_U <- as.matrix(U) %x% matrix(1, nrow = STrain, ncol = 1)
vg_loc <- variog(coords = full_U, data = Y)
vg_model <- lm(vg_loc$v ~ vg_loc$u)
plot(vg_loc, main = "Locations")
abline(a=coef(vg_model)[1], b=coef(vg_model)[2], 
       col = "brown1", lwd=2, lty=2)

# Variogram (global covariates)
full_Z <- as.matrix(Z) %x% matrix(1, nrow = n, ncol = 1)
vg_global <- variog(coords = full_Z, data = Y)
vg_model2 <- lm(vg_global$v ~ vg_global$u)
plot(vg_global, main = "Global Covariates")
abline(a=coef(vg_model2)[1], b=coef(vg_model2)[2], 
       col = "brown1", lwd=2, lty=2)


