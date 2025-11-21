library(ggplot2)
library(MBA)
library(fields)
library(latex2exp)

load("data/slosh/flood_subset.RData")
STest <- nrow(flood.test$Z)
nTest <-  nrow(flood.test$U)
for (storm in 1:STest) {
  
  # Test index for storm s
  index <- (nTest*(storm-1)+1):(nTest*storm)
  true.y <- flood.test$Y[index, ]
  
  # Read in model predictions
  slosh.gp <- readRDS("objects/slosh.RDS")
  preds.gp <- slosh.gp$preds[2, index]
  slosh.svc <- readRDS("objects/slosh_svc_preds.RDS")
  preds.svc <- slosh.svc[storm, ]
  slosh.fosr <- readRDS("objects/slosh_fosr_preds.RDS")
  preds.fosr <- unname(slosh.fosr[ , storm])
  all.data <- c(true.y, preds.gp, preds.svc, preds.fosr)
  
  color.lims <- 1.03 * c(min(all.data), max(all.data))
  
  # Start PDF
  pdf(paste0("figures/slosh/storm", storm, ".pdf"), width = 10, height = 7)
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 4), oma = c(0, 2, 0, 2))
  
  # Actual water levels
  mba.data.true <- data.frame(flood.test$U, true.y)
  mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
  image.plot(mba.interp.true$xyz.est, main = paste0("True Water Levels - Storm ", storm),
             cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
             axis.args = list(cex.axis = 1.5),
             zlim = color.lims)
  
  # GP estimated water levels
  mba.data.gp <- data.frame(flood.test$U, preds.gp)
  mba.interp.gp <- mba.surf(mba.data.gp, no.X=100, no.Y=100, extend=TRUE)
  image.plot(mba.interp.gp$xyz.est, main = paste0("fGP Estimated Water Levels - Storm ", storm),
             cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
             axis.args = list(cex.axis = 1.5),
             zlim = color.lims)
  
  # SVC estimated water levels
  mba.data.svc <- data.frame(flood.test$U, preds.svc)
  mba.interp.svc <- mba.surf(mba.data.svc, no.X=100, no.Y=100, extend=TRUE)
  image.plot(mba.interp.svc$xyz.est, main = paste0("SVC Estimated Water Levels - Storm ", storm),
             cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
             axis.args = list(cex.axis = 1.5),
             zlim = color.lims)
  
  # FOSR estimated water levels
  mba.data.fosr <- data.frame(flood.test$U, preds.fosr)
  mba.interp.fosr <- mba.surf(mba.data.fosr, no.X=100, no.Y=100, extend=TRUE)
  image.plot(mba.interp.fosr$xyz.est, main = paste0("FOSR Estimated Water Levels - Storm ", storm),
             cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
             axis.args = list(cex.axis = 1.5),
             zlim = color.lims)
  
  # Close PDF
  dev.off()
}

