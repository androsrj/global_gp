library(ggplot2)
nScen <- 5
nReps <- 10
line.type <- 2
line.width <- 4

# Density plots for beta0
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  beta0_samples <- results[[1]]$paramSamples$beta[1,]
  hist(beta0_samples, 
       xlab = paste0("Scenario ", i),
       main = "",
       xlim = c(-2.5, 4),
       breaks = 10)
  abline(v = 1, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Beta_0 Samples", side = 3, line = - 2, outer = TRUE)
}

# Density plots for beta1
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  beta1_samples <- results[[1]]$paramSamples$beta[2,]
  hist(beta1_samples, 
       xlab = paste0("Scenario ", i),
       main = "",
       xlim = c(-1.5, 2))
  abline(v = 0.5, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Beta_1 Samples", side = 3, line = - 2, outer = TRUE)
}

# Density plots for beta2
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  beta2_samples <- results[[1]]$paramSamples$beta[3,]
  hist(beta2_samples, 
       xlab = paste0("Scenario ", i),
       main = "",
       xlim = c(-2, 1.5),
       breaks = 10)
  abline(v = -1, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Beta_2 Samples", side = 3, line = - 2, outer = TRUE)
}

# Density plots for theta_f
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  thf_samples <- results[[3]]$paramSamples$thf
  if (i == 2) {
    true_thf <- 10
  } else {
    true_thf = 1
  }
  hist(thf_samples, 
       xlab = paste0("Scenario ", i),
       main = "")
  abline(v = true_thf, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Theta_f Samples", side = 3, line = - 2, outer = TRUE)
}


# Density plots for sigf2
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  sigf2_samples <- results[[2]]$paramSamples$sigf2
  if (i == 3) {
    true_sigf2 <- 20
  } else {
    true_sigf2 <- 5
  }
  hist(sigf2_samples, 
       xlab = paste0("Scenario ", i),
       main = "")
  abline(v = true_sigf2, lty = line.type, lwd = line.width, col = "skyblue4")
  mtext("Sigma2_f Samples", side = 3, line = - 2, outer = TRUE)
}

# Density plots for tau2
par(mfrow = c(2,3))
for (i in 1:nScen) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  tau2_samples <- results[[1]]$paramSamples$tau2
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


### BOXPLOTS FOR PREDICTIVE DIAGNOSTICS ###

# Get diagnostics (rmse, coverage, length)
rmse <- cvg <- len <- c()
for (i in 1:nScen) { 
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  load(paste0("data/small/scen", i, "/test.RData"))
  preds <- lapply(1:nReps, \(j) results[[j]]$preds[2,])
  lower <- lapply(1:nReps, \(j) results[[j]]$preds[1,])
  upper <- lapply(1:nReps, \(j) results[[j]]$preds[3,])
  rmse <- c(rmse, sapply(1:nReps, \(j) sqrt(mean((preds[[j]] - test$Y)^2))))
  cvg <- c(cvg, sapply(1:nReps, \(j) mean(lower[[j]] < test$Y & upper[[j]] > test$Y)))
  len <- c(len, sapply(1:nReps, \(j) mean(upper[[j]] - lower[[j]])))
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

# Interval Length
ggplot(df, aes(x = Scenario, y = Length, group = Scenario, fill = Scenario)) + 
  geom_boxplot() +
  theme_bw() + 
  labs(title = "Interval Length (95%)",
       x = "", y = "")

