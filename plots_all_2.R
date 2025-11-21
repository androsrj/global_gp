library(ggplot2)
library(MBA)
library(fields)
library(latex2exp)

line.type <- 2
line.width <- 4
#which.scens <- c(1, 4, 7, 10)
which.scens <- 11:13
nScen <- length(which.scens)

# Beta0
pdf("figures/new/beta0_pt2.pdf", width = 5, height = 12)
layout.matrix <- matrix(c(1, 1, 
                          1, 1,
                          2, 5,
                          3, 6,
                          4, 7), nrow = 5, ncol = 2, byrow = TRUE)
layout(layout.matrix)

# Scens 1-10
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
nTest <- nrow(test$B)
beta0.true <- test$B[ , 1]

# True beta0
par(mar = c(5, 6.5, 7, 5))
mba.data.true <- data.frame(test$U, beta0.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_0$)"),
      cex.main = 2.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))

# GP
par(mar = c(2, 4, 3, 1))
for (i in which.scens) {
  path.gp <- paste0("objects/small_scen", i, ".RDS") 
  results.gp <- readRDS(path.gp)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta0.means <- results.gp$posteriorMeans$beta.test[1:nTest]
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = "fGP", 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
  mtext(paste0("Scenario ", i), side = 2, line = 2)
}

# SVC
par(mar = c(2, 2, 3, 3))
for (i in which.scens) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta.mu <- mean(results[[1]]$preds$p.beta.recover.samples[ , 1])
  w.mu <- apply(results[[1]]$preds$p.w.predictive.samples[1:nTest, ], 1, mean)
  beta0.means <- beta.mu + w.mu
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = "SVC", 
        col = tim.colors(64), cex.main = 1.5, 
        axes = FALSE)
}

dev.off()



# Beta1
pdf("figures/new/beta1_pt2.pdf", width = 5, height = 12)
layout.matrix <- matrix(c(1, 1, 
                          1, 1,
                          2, 5,
                          3, 6,
                          4, 7), nrow = 5, ncol = 2, byrow = TRUE)
layout(layout.matrix)

# Scens 1-10
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
nTest <- nrow(test$B)
beta1.true <- test$B[ , 2]

# True beta0
par(mar = c(5, 6.5, 7, 5))
mba.data.true <- data.frame(test$U, beta1.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_1$)"),
      cex.main = 2.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))

# GP
par(mar = c(2, 4, 3, 1))
for (i in which.scens) {
  path.gp <- paste0("objects/small_scen", i, ".RDS") 
  results.gp <- readRDS(path.gp)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta1.means <- results.gp$posteriorMeans$beta.test[(nTest+1):(2*nTest)]
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = "fGP", 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
  mtext(paste0("Scenario ", i), side = 2, line = 2)
}

# SVC
par(mar = c(2, 2, 3, 3))
for (i in which.scens) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta.mu <- mean(results[[1]]$preds$p.beta.recover.samples[ , 2])
  w.mu <- apply(results[[1]]$preds$p.w.predictive.samples[(nTest+1):(2*nTest), ], 1, mean)
  beta1.means <- beta.mu + w.mu
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = "SVC", 
        col = tim.colors(64), cex.main = 1.5, 
        axes = FALSE)
}

dev.off()




# Beta2
pdf("figures/new/beta2_pt2.pdf", width = 5, height = 12)
layout.matrix <- matrix(c(1, 1, 
                          1, 1,
                          2, 5,
                          3, 6,
                          4, 7), nrow = 5, ncol = 2, byrow = TRUE)
layout(layout.matrix)

# Scens 1-10
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
nTest <- nrow(test$B)
beta2.true <- test$B[ , 3]

# True beta0
par(mar = c(5, 6.5, 7, 5))
mba.data.true <- data.frame(test$U, beta2.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_2$)"),
      cex.main = 2.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))

# GP
par(mar = c(2, 4, 3, 1))
for (i in which.scens) {
  path.gp <- paste0("objects/small_scen", i, ".RDS") 
  results.gp <- readRDS(path.gp)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta2.means <- results.gp$posteriorMeans$beta.test[(2*nTest+1):(3*nTest)]
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = "fGP", 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
  mtext(paste0("Scenario ", i), side = 2, line = 2)
}

# SVC
par(mar = c(2, 2, 3, 3))
for (i in which.scens) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta.mu <- mean(results[[1]]$preds$p.beta.recover.samples[ , 3])
  w.mu <- apply(results[[1]]$preds$p.w.predictive.samples[(2*nTest+1):(3*nTest), ], 1, mean)
  beta2.means <- beta.mu + w.mu
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = "SVC", 
        col = tim.colors(64), cex.main = 1.5, 
        axes = FALSE)
}

dev.off()
