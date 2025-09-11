library(ggplot2)
library(spBayes)
library(MBA)
library(fields)
library(latex2exp)
nScen <- 6
nReps <- 10
line.type <- 2
line.width <- 4
nTest <- 25

# Surface plots for beta0
pdf("figures/svc/beta0_svc.pdf", width = 8, height = 6)
par(mfrow = c(2,3), mar = c(3, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  beta.mu <- mean(results[[1]]$preds$p.beta.recover.samples[ , 1])
  w.mu <- apply(results[[1]]$preds$p.w.predictive.samples[1:nTest, ], 1, mean)
  beta0.means <- beta.mu + w.mu
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scen. ", i), 
        col = tim.colors(64), cex.main = 1.5)
  mtext(TeX("$\\beta_0"), side = 3, line = 1, outer = TRUE, cex = 1.5)
}
dev.off()

# Surface plots for beta1
pdf("figures/svc/beta1_svc.pdf", width = 8, height = 6)
par(mfrow = c(2,3), mar = c(3, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  beta.mu <- mean(results[[1]]$preds$p.beta.recover.samples[ , 2])
  w.mu <- apply(results[[1]]$preds$p.w.predictive.samples[(nTest+1):(2*nTest), ], 1, mean)
  beta1.means <- beta.mu + w.mu
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scen. ", i), 
        col = tim.colors(64), cex.main = 1.5)
  mtext(TeX("$\\beta_1"), side = 3, line = 1, outer = TRUE, cex = 1.5)
}
dev.off()

# Surface plots for beta2
pdf("figures/svc/beta2_svc.pdf", width = 8, height = 6)
par(mfrow = c(2,3), mar = c(3, 4, 2, 2) + 0.1, oma = c(0, 0, 4, 0))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  beta.mu <- mean(results[[1]]$preds$p.beta.recover.samples[ , 3])
  w.mu <- apply(results[[1]]$preds$p.w.predictive.samples[(2*nTest+1):(3*nTest), ], 1, mean)
  beta2.means <- beta.mu + w.mu
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scen. ", i), 
        col = tim.colors(64), cex.main = 1.5)
  mtext(TeX("$\\beta_2"), side = 3, line = 1, outer = TRUE, cex = 1.5)
}
dev.off()

# Density plots for tau2
pdf("figures/svc/tau2_svc.pdf", width = 8, height = 6)
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  tau2_samples <- results[[1]]$model$p.theta.samples[ , 4]
  if (i == 4) {
    true_tau2 <- 2
  } else {
    true_tau2 <- 0.2
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
for (i in 1:nScen) { 
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  STest <- nrow(test$Z)
  nTest <- nrow(test$X)
  a <- .05
  for (j in 1:nReps) {
    rmse_vec <- cvg_vec <- len_vec <- numeric(STest)
    for (k in 1:STest) {
      truth <- test$Y[(nTest*(k-1)+1):(nTest*k), ]
      pred.samples <- results[[j]]$preds$p.y.predictive.samples
      preds <- apply(pred.samples, 1, mean)
      rmse_vec[k] <- sqrt(mean((truth - preds)^2))
      lower <- apply(pred.samples, 1, quantile, .025)
      upper <- apply(pred.samples, 1, quantile, .975)
      cvg_vec[k] <- mean(lower < truth & upper > truth)
      len_vec[k] <- mean(upper - lower)
    }
    rmse_vec <- mean(rmse_vec)
    cvg_vec <- mean(cvg_vec)
    len_vec <- mean(len_vec)
    
    rmse <- c(rmse, rmse_vec)
    cvg <- c(cvg, cvg_vec)
    len <- c(len, len_vec)
  }
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
  scale_y_continuous(limits = c(0.75, 0.95)) + 
  labs(title = "Interval Coverage (95%) for SVC",
       x = "", y = "")
ggsave("figures/svc/coverage_svc.pdf")

# Interval Length
ggplot(df, aes(x = Scenario, y = Length, group = Scenario, fill = Scenario)) + 
  geom_boxplot() +
  theme_bw() + 
  scale_y_continuous(limits = c(150, 225)) + 
  labs(title = "Interval Length (95%) for SVC",
       x = "", y = "")
ggsave("figures/svc/length_svc.pdf")
