library(MBA)
library(fields)
load("data/flood_data.RData")

# Global covariates (Z)
STrain <- 10
STest <- 5
samp <- sample(1:nrow(inputs), STrain + STest)
stormsTrain <- samp[1:STrain]
stormsTest <- samp[(STrain+1):(STrain+STest)]
Z <- inputs[stormsTrain,3:4]
ZTest <- inputs[stormsTest,3:4]

# Local covariates (X)
n <- 500
nTest <- 20
samp2 <- sample(1:nrow(coords))
train <- samp2[1:n]
test <- samp2[(n+1):(n+nTest)]
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

# Storm 1
pred.surf <-  mba.surf(cbind(UTest, YTest[1:nTest]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="Storm 1", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)

# Storm 2
pred.surf <-  mba.surf(cbind(UTest, YTest[(nTest+1):(2*nTest)]), no.X=100, no.Y=100, extend=T)$xyz.est
image.plot(pred.surf, xaxs ="r", yaxs = "r", main="Storm 2", 
           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
contour(pred.surf, add=T)
