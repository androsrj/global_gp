nScen <- 6
rmse <- cvg <- len <- c()
for (i in 1:nScen) { 
  path <- paste0("objects/svc_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  load(paste0("data/small/scen", i, "/test.RData"))
  STest <- nrow(test$Z)
  nTest <- nrow(test$X)
  
  # Betas
  cat(paste0("Point estimates: ", apply(results$p.beta.recover.samples, 2, mean), "\n"))
  cat(paste0("Lower: ", apply(results$p.beta.recover.samples, 2, quantile, .025), "\n"))
  cat(paste0("Upper: ", apply(results$p.beta.recover.samples, 2, quantile, .975), "\n"))
  
  rmse_vec <- cvg_vec <- len_vec <- numeric(STest)
  for (k in 1:STest) {
    truth <- test$Y[(nTest*(k-1)+1):(nTest*k), ]
    m.3.pred <- spPredict(results, pred.covars = cbind(rep(1, nTest), test$X),
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
rmse
cvg
