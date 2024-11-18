
nReps <- 10
for (i in 1:5) { 
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)
  load(paste0("data/small/scen", i, "/test.RData"))
  preds <- lapply(1:nReps, \(j) results[[j]]$preds[2,])
  lower <- lapply(1:nReps, \(j) results[[j]]$preds[1,])
  upper <- lapply(1:nReps, \(j) results[[j]]$preds[3,])
  paramMeans <- lapply(1:nReps, \(j) results[[j]]$posteriorMeans)
  rmse <- sapply(1:nReps, \(j) sqrt(mean((preds[[j]] - test$Y)^2)))
  cvg <- sapply(1:nReps, \(j) mean(lower[[j]] < test$Y & upper[[j]] > test$Y))
}
