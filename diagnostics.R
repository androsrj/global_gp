library(ggplot2)
library(scoringRules)
library(scoringutils)
nScen <- 9
nReps <- 10

##########################################
##### Predictive diagnostics for fGP #####
##########################################
rmse <- cvg <- len <- score <- crps <- c()
std.dev <- rep(0, nScen)
for (i in 1:nScen) { 
  path <- paste0("objects/fgp_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  preds <- lapply(1:nReps, \(j) results[[j]]$preds[2,])
  lower <- lapply(1:nReps, \(j) results[[j]]$preds[1,])
  upper <- lapply(1:nReps, \(j) results[[j]]$preds[3,])
  rmse <- c(rmse, sapply(1:nReps, \(j) sqrt(mean((preds[[j]] - test$Y)^2))))
  cvg <- c(cvg, sapply(1:nReps, \(j) mean(lower[[j]] < test$Y & upper[[j]] > test$Y)))
  len <- c(len, sapply(1:nReps, \(j) mean(upper[[j]] - lower[[j]])))
  score <- c(score, sapply(1:nReps, \(j) mean(scoringutils:::interval_score(c(test$Y), lower[[j]], upper[[j]], interval_range = rep(95, nrow(test$Y))))))
  pred.samples <- do.call(cbind, lapply(1:nReps, \(j) results[[j]]$predSamples))
  crps <- c(crps, mean(crps_sample(c(test$Y), pred.samples)))
  std.dev[i] <- sd(test$Y)
}

# Organize data
df <- data.frame(RMSE = rmse, 
                 Coverage = cvg,
                 Length = len,
		 Score = score,
		 CRPS = crps,
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
avg.rmse <- aggregate(data = df, RMSE ~ Scenario, mean)[ , 2]
avg.cvg <- aggregate(data = df, Coverage ~ Scenario, mean)[ , 2]
avg.length <- aggregate(data = df, Length ~ Scenario, mean)[ , 2]
avg.score <- aggregate(data = df, Score ~ Scenario, mean)[ , 2]
avg.crps <- aggregate(data = df, CRPS ~ Scenario, mean)[ , 2]
cat("Predictive diagnostics for fGP: \n")
data.frame(Scen = 1:nScen, std.dev, avg.rmse, avg.cvg, avg.length, avg.score, avg.crps)


##########################################
##### Predictive diagnostics for SVC #####
##########################################
rmse <- cvg <- len <- score <- crps <- c()
std.dev <- rep(0, nScen)
for (i in 1:nScen) { 
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)
  test <- readRDS(paste0("data/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  STest <- nrow(test$Z)
  nTest <- nrow(test$X)
  a <- .05
  for (j in 1:nReps) {
    rmse_vec <- cvg_vec <- len_vec <- score_vec <- crps_vec <- numeric(STest)
    for (k in 1:STest) {
      truth <- test$Y[(nTest*(k-1)+1):(nTest*k), ]
      pred.samples <- results[[j]]$preds$p.y.predictive.samples
      preds <- apply(pred.samples, 1, mean)
      rmse_vec[k] <- sqrt(mean((truth - preds)^2))
      lower <- apply(pred.samples, 1, quantile, .025)
      upper <- apply(pred.samples, 1, quantile, .975)
      cvg_vec[k] <- mean(lower < truth & upper > truth)
      len_vec[k] <- mean(upper - lower)
      score_vec[k] <- mean(scoringutils:::interval_score(truth, lower, upper, interval_range = rep(95, length(truth))))
      crps_vec[k] <- mean(crps_sample(truth, pred.samples))
    }
    rmse_vec <- mean(rmse_vec)
    cvg_vec <- mean(cvg_vec)
    len_vec <- mean(len_vec)
    score_vec <- mean(score_vec)
    crps_vec <- mean(crps_vec)
    
    rmse <- c(rmse, rmse_vec)
    cvg <- c(cvg, cvg_vec)
    len <- c(len, len_vec)
    score <- c(score, score_vec)
    crps <- c(crps, crps_vec)
  }
  std.dev[i] <- sd(test$Y)
}

# Organize data
df <- data.frame(RMSE = rmse, 
                 Coverage = cvg,
                 Length = len,
		 Score = score,
		 CRPS = crps,
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

# Table for RMSE, coverage, and length
avg.rmse <- aggregate(data = df, RMSE ~ Scenario, mean)[ , 2]
avg.cvg <- aggregate(data = df, Coverage ~ Scenario, mean)[ , 2]
avg.length <- aggregate(data = df, Length ~ Scenario, mean)[ , 2]
avg.score <- aggregate(data = df, Score ~ Scenario, mean)[ , 2]
avg.crps <- aggregate(data = df, CRPS ~ Scenario, mean)[ , 2]
cat("Diagnostics for SVC: \n")
data.frame(Scen = 1:nScen, std.dev, avg.rmse, avg.cvg, avg.length, avg.score, avg.crps)

