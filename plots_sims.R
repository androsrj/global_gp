library(ggplot2)
library(MBA)
library(fields)
library(latex2exp)

line.type <- 2
line.width <- 4
which.scens <- c(1:8)
nScen <- length(which.scens)
#subj <- 1

####################################
######## Surface plots of Y ########
####################################

for (i in which.scens) {
  
  # True Y for given scenario
  pdf(paste0("figures/sims/y_true_scen", i, ".pdf"), width = 5, height = 3)
  par(mar = c(2, 2, 1, 1))
  test <- readRDS(paste0("data/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  mba.data.true <- data.frame(test$U, test$Y[1:nTest, ])
  mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
  image.plot(mba.interp.true$xyz.est, main = "",
        cex.main = 1.25, col = tim.colors(64))
  
  dev.off()
  
  # Estimated Y for given scenario
  pdf(paste0("figures/sims/y_pred_scen", i, ".pdf"), width = 5, height = 3)
  par(mar = c(2, 2, 1, 1))
  path <- paste0("objects/fgp_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  y.hat <- results$preds[2, 1:nTest]
  mba.data <- data.frame(test$U, y.hat)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image.plot(mba.interp$xyz.est, main = "", 
        col = tim.colors(64), cex.main = 1.25)
  dev.off()
}


#### Posterior distributions of tau2 for each scenario
# Density plots for tau2
for (i in which.scens) {
  pdf(paste0("figures/sims/tau2_scen", i, ".pdf"), width = 5, height = 4)
  par(mar = c(2, 2, 1, 1))
  path <- paste0("objects/fgp_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  tau2.samples <- results$paramSamples$tau2
  if (i == 2) {
    true_tau2 <- 2
  } else {
    true_tau2 = 0.2
  }
  hist(tau2.samples, 
       #xlab = paste0("Scenario ", i),
       xlab = "",
       main = "", ylab = "", cex.lab = 1.75)
  abline(v = true_tau2, lty = line.type, lwd = line.width, col = "blue")
  #mtext(TeX("Posterior Samples of $\\tau^2$"), side = 3, line = -1, outer = TRUE, cex = 1.5)
  dev.off()
}


