load("data/train.RData")
load("data/test.RData")
n <- nrow(train$X)
nTest <- nrow(test$X)
STrain <- nrow(train$Z)
STest <- nrow(test$Z)
library(geoR)

par(mfrow=c(1,2))
# Variogram (locations)
full_U <- train$U %x% matrix(1, nrow = STrain, ncol = 1)
vg_loc <- variog(coords = full_U, data = train$Y)
vg_model <- lm(vg_loc$v ~ vg_loc$u)
plot(vg_loc, main = "Locations")
abline(a=coef(vg_model)[1], b=coef(vg_model)[2], 
       col = "brown1", lwd=2, lty=2)

# Variogram (global covariates)
full_Z <- train$Z %x% matrix(1, nrow = n, ncol = 1)
vg_global <- variog(coords = full_Z, data = train$Y)
vg_model2 <- lm(vg_global$v ~ vg_global$u)
plot(vg_global, main = "Global Covariates")
abline(a=coef(vg_model2)[1], b=coef(vg_model2)[2], 
       col = "brown1", lwd=2, lty=2)

dists <- sv <- matrix(0, nrow = STrain, ncol = 10)
for (i in 1:STrain) {
  index <- seq(0, n*(STrain-1), by = n) + i
  vg <- variog(coords = train$Z, data = train$Y[index])
  dists[i, ] <- vg$u
  sv[i, ] <- vg$v
}
dists_agg <- apply(dists, 2, mean)
sv_agg <- apply(sv, 2, mean)
plot(dists_agg, sv_agg)
