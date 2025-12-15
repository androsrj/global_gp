library(ggplot2)
library(MBA)
library(fields)
library(latex2exp)

line.type <- 2
line.width <- 4
which.scens <- c(1:3, 7:9, 11:13)
nScen <- length(which.scens)


#################################
############# Beta0 #############
#################################

pdf("figures/new/beta0.pdf", width = 5, height = 16)
par(mar = c(2, 3, 3, 1), oma = c(0, 0, 3, 0), mfrow = c(nScen, 2))

for (i in which.scens) {
  
  # Data for given scenario
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta0.true <- test$B[ , 1]
  
  # True beta0 for given scenario
  mba.data.true <- data.frame(test$U, beta0.true)
  mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp.true$xyz.est, main = paste0("Scenario ", i, " "),
        cex.main = 1.25, col = tim.colors(64))
  
  # Model results for given scenario
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta0.means <- results$posteriorMeans$beta.test[1:nTest]
  
  # Estimated beta0 for given scenario
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.25,
        axes = FALSE)
  
  # Label each column of plots
  mtext(TeX("True surface of $\\beta_0$"), 
        side = 3, line = 1, outer = TRUE, 
        at = 0.27, cex = 1.1)
  mtext(TeX("Estimated surface of $\\beta_0$"), 
        side = 3, line = 1, outer = TRUE, 
        at = 0.77, cex = 1.1)
}
dev.off()


#################################
############# Beta1 #############
#################################

pdf("figures/new/beta1.pdf", width = 5, height = 16)
par(mar = c(2, 3, 3, 1), oma = c(0, 0, 3, 0), mfrow = c(nScen, 2))

for (i in which.scens) {
  
  # Data for given scenario
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta1.true <- test$B[ , 2]
  
  # True beta1 for given scenario
  mba.data.true <- data.frame(test$U, beta1.true)
  mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp.true$xyz.est, main = paste0("Scenario ", i, " "),
        cex.main = 1.25, col = tim.colors(64))
  
  # Model results for given scenario
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta1.means <- results$posteriorMeans$beta.test[(nTest+1):(2*nTest)]
  
  # Estimated beta1 for given scenario
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.25,
        axes = FALSE)
  
  # Label each column of plots
  mtext(TeX("True surface of $\\beta_1$"), 
        side = 3, line = 1, outer = TRUE, 
        at = 0.27, cex = 1.1)
  mtext(TeX("Estimated surface of $\\beta_1$"), 
        side = 3, line = 1, outer = TRUE, 
        at = 0.77, cex = 1.1)
}
dev.off()


#################################
############# Beta2 #############
#################################

pdf("figures/new/beta2.pdf", width = 5, height = 16)
par(mar = c(2, 3, 3, 1), oma = c(0, 0, 3, 0), mfrow = c(nScen, 2))

for (i in which.scens) {
  
  # Data for given scenario
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta2.true <- test$B[ , 3]
  
  # True beta2 for given scenario
  mba.data.true <- data.frame(test$U, beta2.true)
  mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp.true$xyz.est, main = paste0("Scenario ", i, " "),
        cex.main = 1.25, col = tim.colors(64))
  
  # Model results for given scenario
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta2.means <- results$posteriorMeans$beta.test[(2*nTest+1):(3*nTest)]
  
  # Estimated beta2 for given scenario
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.25,
        axes = FALSE)
  
  # Label each column of plots
  mtext(TeX("True surface of $\\beta_2$"), 
        side = 3, line = 1, outer = TRUE, 
        at = 0.27, cex = 1.1)
  mtext(TeX("Estimated surface of $\\beta_2$"), 
        side = 3, line = 1, outer = TRUE, 
        at = 0.77, cex = 1.1)
}
dev.off()


#### Posterior distributions of tau2 for each scenario
# Density plots for tau2
pdf("figures/new/tau2.pdf", width = 10, height = 8)
par(mfrow = c(3, 3), mar = c(4, 3, 3, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in which.scens) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  tau2.samples <- results$paramSamples$tau2
  if (i == 4 | i == 10 | i == 12) {
    true_tau2 <- 2
  } else {
    true_tau2 = 0.2
  }
  hist(tau2.samples, 
       xlab = paste0("Scenario ", i),
       main = "", ylab = "", cex.lab = 1.75)
  abline(v = true_tau2, lty = line.type, lwd = line.width, col = "blue")
  mtext(TeX("Posterior Samples of $\\tau^2$"), side = 3, line = -1, outer = TRUE, cex = 1.5)
}
dev.off()


