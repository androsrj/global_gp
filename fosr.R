library(fields)
library(refund)
library(MBA)

scen <- "scen6"
size <- "small"

dir <- paste0("data/", size, "/", scen, "/")
load(paste0(dir, "train.RData"))
load(paste0(dir, "test.RData"))
n <- nrow(train$X)
nTest <- nrow(test$X)
S <- nrow(train$Z)
STest <- nrow(test$Z)

Y <- rbind(matrix(train$Y, nrow = n, ncol = S),
           matrix(test$Y, nrow = nTest, ncol = STest))
X <- rbind(train$X, test$X)
colnames(X) <- c("X1", "X2")
dfl <- as.data.frame(X)
dfl$Y <- Y

train.index <- 1:n
test.index <- (n+1):(n+nTest)

#fit <- fosr.vs(data = dfl[1:200,], formula = Y ~ ., method = "grMCP")
fit <- bayes_fosr(data = dfl[train.index, ], Y ~ X1 + X2, est.method = "VB")
#fit$beta.hat
apply(fit$beta.hat, 1, mean)
apply(fit$beta.LB, 1, mean)
apply(fit$beta.UB, 1, mean)

pred <- predict(object = fit, newdata = dfl[test.index,])
sqrt(mean((pred - Y[test.index, ])^2))
sd(test$Y)

#lims <- c(-15, 15)
#pdf("figures/subj1_fosr.pdf")
#pred.surf <-  mba.surf(cbind(test$U, pred[ , 1]), no.X=100, no.Y=100, extend=T)$xyz.est
#image.plot(pred.surf, xaxs ="r", yaxs = "r", zlim = lims, main="FOSR, Subject 1", 
#           cex.main = 1.5, col = hcl.colors(12, "YlOrRd", rev=TRUE))
#contour(pred.surf, add=T)
#dev.off()
