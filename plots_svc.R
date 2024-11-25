library(ggplot2)
library(spBayes)
nScen <- 6
nReps <- 10
line.type <- 2
line.width <- 4

# Density plots for beta0
pdf("figures/svc/beta0_svc.pdf")
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  beta0_samples <- results[[1]]$p.beta.recover.samples[, 1]
  hist(beta0_samples, 
       xlab = paste0("Scenario ", i),
       main = "",
       xlim = c(0, 6),
       breaks = 10)
  abline(v = 1, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Beta_0 Samples", side = 3, line = - 2, outer = TRUE)
}
dev.off()

# Density plots for beta1
pdf("figures/svc/beta1_svc.pdf")
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  beta1_samples <- results[[1]]$p.beta.recover.samples[, 2]
  hist(beta1_samples, 
       xlab = paste0("Scenario ", i),
       main = "",
       xlim = c(-0.5, 1))
  abline(v = 0.5, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Beta_1 Samples", side = 3, line = - 2, outer = TRUE)
}
dev.off()

# Density plots for beta2
pdf("figures/svc/beta2_svc.pdf")
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  beta2_samples <- results[[1]]$p.beta.recover.samples[, 3]
  hist(beta2_samples, 
       xlab = paste0("Scenario ", i),
       main = "",
       xlim = c(-1.5, -0.5),
       breaks = 10)
  abline(v = -1, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Beta_2 Samples", side = 3, line = - 2, outer = TRUE)
}
dev.off()

# Density plots for tau2
pdf("figures/svc/tau2_svc.pdf")
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  tau2_samples <- results[[1]]$p.theta.samples[, 3]
  if (i == 4) {
    true_tau2 <- 2
  } else {
    true_tau2 = 0.2
  }
  hist(tau2_samples, 
       xlab = paste0("Scenario ", i),
       main = "")
  abline(v = true_tau2, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Tau2 Samples", side = 3, line = - 2, outer = TRUE)
}
dev.off()

### BOXPLOTS FOR PREDICTIVE DIAGNOSTICS ###

# Get diagnostics (rmse, coverage, length)
rmse <- cvg <- len <- c()
for (i in 1:nScen) { 
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  load(paste0("data/small/scen", i, "/test.RData"))
  STest <- nrow(test$Z)
  nTest <- nrow(test$X)
  a <- .05
  for (j in 1:nReps) {
    rmse_vec <- cvg_vec <- len_vec <- numeric(STest)
    for (k in 1:STest) {
      truth <- test$Y[(nTest*(k-1)+1):(nTest*k), ]
      m.3.pred <- spPredict(results[[j]], pred.covars = cbind(rep(1, nTest), test$X),
                            pred.coords=test$U + rnorm(50, 0, 0.0001), thin=10,
                            joint=TRUE, n.omp.threads=4, verbose=FALSE)
      preds <- apply(m.3.pred$p.y.predictive.samples, 1, mean)
      rmse_vec[k] <- sqrt(mean((truth - preds)^2))
      lower <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .025)
      upper <- apply(m.3.pred$p.y.predictive.samples, 1, quantile, .975)
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
  labs(title = "Interval Coverage (95%)",
       x = "", y = "")
ggsave("figures/svc/coverage_svc.pdf")

# Interval Length
ggplot(df, aes(x = Scenario, y = Length, group = Scenario, fill = Scenario)) + 
  geom_boxplot() +
  theme_bw() + 
  labs(title = "Interval Length (95%)",
       x = "", y = "")
ggsave("figures/svc/length_svc.pdf")
