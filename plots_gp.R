library(ggplot2)
library(MBA)
library(fields)
library(latex2exp)

nScen <- 10
nReps <- 10
line.type <- 2
line.width <- 4
nTest <- 25

# Surface plots for true beta surfaces
test <- readRDS(paste0("data/small/scen", 1, "/test.RDS"))
beta0.true <- test$U[ , 1] - test$U[ , 2]
beta1.true <- test$U[ , 1] + test$U[ , 2] - 100
beta2.true <- 2 * test$U[ , 1] - test$U[ , 2] - 50

pdf("figures/gp/beta_true.pdf", width = 10, height = 3)
par(mfrow = c(1,3), mar = c(5, 5, 4, 8) + 0.2)

# Beta0
mba.data.true <- data.frame(test$U, beta0.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_0$)"),
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           axis.args = list(cex.axis = 1.5))

# Beta1
mba.data.true <- data.frame(test$U, beta1.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_1$)"),
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           axis.args = list(cex.axis = 1.5))

# Beta2
mba.data.true <- data.frame(test$U, beta2.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_2$)"),
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           axis.args = list(cex.axis = 1.5))

dev.off()

# Surface plots for beta0
pdf("figures/gp/beta0_gp.pdf", width = 10, height = 4)
par(mfrow = c(2,5), mar = c(3, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  beta0.means <- results$posteriorMeans$beta.test[1:nTest]
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scen. ", i), 
        col = tim.colors(64), cex.main = 1.5)
  mtext(TeX("$\\beta_0"), side = 3, line = 1, outer = TRUE, cex = 1.5)
}
dev.off()

# Surface plots for beta1
pdf("figures/gp/beta1_gp.pdf", width = 10, height = 4)
par(mfrow = c(2,5), mar = c(3, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  beta1.means <- results$posteriorMeans$beta.test[(nTest+1):(2*nTest)]
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scen. ", i), 
        col = tim.colors(64), cex.main = 1.5)
  mtext(TeX("$\\beta_1"), side = 3, line = 1, outer = TRUE, cex = 1.5)
}
dev.off()

# Surface plots for beta2
pdf("figures/gp/beta2_gp.pdf", width = 10, height = 4)
par(mfrow = c(2,5), mar = c(3, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  beta2.means <- results$posteriorMeans$beta.test[(2*nTest+1):(3*nTest)]
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scen. ", i), 
        col = tim.colors(64), cex.main = 1.5)
  mtext(TeX("$\\beta_2"), side = 3, line = 1, outer = TRUE, cex = 1.5)
}
dev.off()

# Density plots for tau2
pdf("figures/gp/tau2_gp.pdf", width = 10, height = 4)
par(mfrow = c(2,5))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  tau2_samples <- results$paramSamples$tau2
  if (i == 4) {
    true_tau2 <- 2
  } else {
    true_tau2 = 0.2
  }
  hist(tau2_samples, 
       xlab = paste0("Scenario ", i),
       main = "", ylab = "", cex.lab = 1.75)
  abline(v = true_tau2, lty = line.type, lwd = line.width, col = "blue")
  mtext(TeX("$\\tau^2$"), side = 3, line = -2.5, outer = TRUE, cex = 1.5)
}
dev.off()



### BOXPLOTS FOR PREDICTIVE DIAGNOSTICS ###

# Get diagnostics (rmse, coverage, length)
rmse <- cvg <- len <- c()
std.dev <- rep(0, nScen)
for (i in 1:nScen) { 
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  preds <- lapply(1:nReps, \(j) results[[j]]$preds[2,])
  lower <- lapply(1:nReps, \(j) results[[j]]$preds[1,])
  upper <- lapply(1:nReps, \(j) results[[j]]$preds[3,])
  rmse <- c(rmse, sapply(1:nReps, \(j) sqrt(mean((preds[[j]] - test$Y)^2))))
  cvg <- c(cvg, sapply(1:nReps, \(j) mean(lower[[j]] < test$Y & upper[[j]] > test$Y)))
  len <- c(len, sapply(1:nReps, \(j) mean(upper[[j]] - lower[[j]])))
  std.dev[i] <- sd(test$Y)
}

# Organize data
df <- data.frame(RMSE = rmse, 
                 Coverage = cvg,
                 Length = len,
                 Scenario = factor(rep(1:nScen, each = nReps)),
                 Rep = rep(1:nReps, nScen))

# Interval Coverage
ggplot(df, aes(x = Scenario, y = Coverage, group = Scenario, fill = Scenario)) + 
  geom_boxplot() +
  theme_bw() +
  labs(title = "Interval Coverage (95%) for Global GP",
       x = "", y = "")
ggsave("figures/gp/coverage_gp.pdf")

# Interval Length
ggplot(df, aes(x = Scenario, y = Length, group = Scenario, fill = Scenario)) + 
  geom_boxplot() +
  theme_bw() + 
  labs(title = "Interval Length (95%) for Global GP",
       x = "", y = "")
ggsave("figures/gp/length_gp.pdf")

# Table for RMSE, coverage, and length
avg.rmse <- aggregate(data = df, RMSE ~ Scenario, mean)
avg.cvg <- aggregate(data = df, Coverage ~ Scenario, mean)[, 2]
avg.length <- aggregate(data = df, Length ~ Scenario, mean)[, 2]
cbind(avg.rmse, std.dev, avg.cvg, avg.length)

# Scenario 11

# Read in data and model results
path <- paste0("objects/small_scen", 11, ".RDS") 
results <- readRDS(path)[[1]]
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
beta0.true <- test$B[ , 1]
beta1.true <- test$B[ , 2]
beta2.true <- test$B[ , 3]
beta0.means <- results$posteriorMeans$beta.test[1:nTest]
beta1.means <- results$posteriorMeans$beta.test[(nTest+1):(2*nTest)]
beta2.means <- results$posteriorMeans$beta.test[(2*nTest+1):(3*nTest)]

# Make legend limits for plots
expand.by <- 1.2
beta0.lims <- c(min(c(beta0.true, beta0.means)), max(c(beta0.true, beta0.means))) * expand.by
beta1.lims <- c(min(c(beta1.true, beta1.means)), max(c(beta1.true, beta1.means))) * expand.by
beta2.lims <- c(min(c(beta2.true, beta2.means)), max(c(beta2.true, beta2.means))) * expand.by

# Surface plots for true beta surfaces
pdf("figures/gp/beta_true_11.pdf", width = 10, height = 3)
par(mfrow = c(1,3), mar = c(5, 5, 4, 8) + 0.2)

# Beta0
mba.data.true <- data.frame(test$U, beta0.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_0$)"),
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           zlim = beta0.lims, axis.args = list(cex.axis = 1.5))

# Beta1
mba.data.true <- data.frame(test$U, beta1.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_1$)"),
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           zlim = beta1.lims, axis.args = list(cex.axis = 1.5))

# Beta2
mba.data.true <- data.frame(test$U, beta2.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp.true$xyz.est, main = TeX("True surface ($\\beta_2$)"),
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           zlim = beta2.lims, axis.args = list(cex.axis = 1.5))

dev.off()

# Plot beta estimated surfaces for scen 11
pdf("figures/gp/beta_scen11.pdf", width = 10, height = 3)
par(mfrow = c(1,3), mar = c(5, 5, 4, 8) + 0.2)

# Beta0
mba.data <- data.frame(test$U, beta0.means)
mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp$xyz.est, main = TeX("Estimated Surface $\\beta_0"), 
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           zlim = beta0.lims, col = tim.colors(64), cex.main = 1.5)

# Beta1
mba.data <- data.frame(test$U, beta1.means)
mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp$xyz.est, main = TeX("Estimated Surface $\\beta_1"), 
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           zlim = beta1.lims, col = tim.colors(64), cex.main = 1.5)

# Beta2
mba.data <- data.frame(test$U, beta2.means)
mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
image.plot(mba.interp$xyz.est, main = TeX("Estimated Surface $\\beta_2"), 
           cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
           zlim = beta2.lims, col = tim.colors(64), cex.main = 1.5)

dev.off()


