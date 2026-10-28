library(GWmodel)
library(scoringRules)
library(scoringutils)
scens <- 1:9
nScen <- length(scens)

results <- data.frame(
  scen = scens,
  sd = numeric(nScen),
  rmse = numeric(nScen),
  coverage = numeric(nScen),
  length = numeric(nScen),
  score = numeric(nScen)
)
for (i in scens) {
  cat(paste0("Beginning scenario ", i, "\n"))
  train <- readRDS(paste0("../data/scen", i, "/train.RDS"))
  test <- readRDS(paste0("../data/scen", i,"/test.RDS"))
  nSubj <- nrow(train$Z)
  n <- nrow(train$X)
  nTest <- nrow(test$X)

  if (i == 9) {
    train$X <- train$X[ , -1]
    test$X <- test$X[ , -1]
  }
  
  rmse <- cvg <- len <- score <- numeric(nSubj)
  for (j in 1:nSubj) {
    train.index <- ((j-1)*n+1):(j*n)
    test.index <- ((j-1)*nTest+1):(j*nTest)
    
    # Compile training data
    gwr.train <- data.frame(train$X, train$U, train$Y[train.index, ])
    colnames(gwr.train) <- c("X1", "X2", "U1", "U2", "Y")
    coordinates(gwr.train) <- ~ U1 + U2
    
    # Compile testing data
    gwr.test <- data.frame(test$X, test$U, test$Y[test.index, ])
    colnames(gwr.test) <- c("X1", "X2", "U1", "U2", "Y")
    coordinates(gwr.test) <- ~ U1 + U2
    
    # Fit model
    out <- capture.output(
      bw <- bw.gwr(Y ~ ., data = gwr.train, approach = "AICc")
    )
    gwr.model <- gwr.basic(Y ~ ., data = gwr.train, bw = bw)
    
    # Check training error
    #preds.train <- gwr.model$SDF@data$yhat
    #sqrt(mean((train$Y[1:n, ] - preds.train)^2))
    #sd(train$Y)
    
    # Predict on test data
    preds.gwr <- gwr.predict(
      formula = Y ~ .,
      data = gwr.train,
      predictdata = gwr.test,
      bw = gwr.model$GW.arguments$bw,
      kernel = gwr.model$GW.arguments$kernel,
      adaptive = gwr.model$GW.arguments$adaptive
    )
    preds.test <- preds.gwr$SDF@data$prediction
    upper <- preds.test + qnorm(.975) * sqrt(preds.gwr$SDF@data$prediction_var)
    lower <- preds.test - qnorm(.975) * sqrt(preds.gwr$SDF@data$prediction_var)
    
    # Check testing error and other diagnostics
    y.true <- test$Y[test.index, ]
    rmse[j] <- sqrt(mean((y.true - preds.test)^2))
    cvg[j] <- mean(y.true < upper & y.true > lower)
    len[j] <- mean(upper - lower)
    score[j] <- mean(scoringutils:::interval_score(y.true, lower, upper, interval_range = rep(95, length(y.true))))
    
    # Check RMSE for estimating beta surfaces on test data
    #beta0.est <- preds.gwr$SDF@data$Intercept_coef
    #beta1.est <- preds.gwr$SDF@data$X1_coef
    #beta2.est <- preds.gwr$SDF@data$X2_coef
    
    # Beta0
    #sqrt(mean((test$B[ , 1] - beta0.est)^2))
    #sd(test$B[ , 1])
    
    # Beta1
    #sqrt(mean((test$B[ , 2] - beta1.est)^2))
    #sd(test$B[ , 2])
    
    # Beta2
    #sqrt(mean((test$B[ , 3] - beta2.est)^2))
    #sd(test$B[ , 3])
  }
  
  if (i != 9) {
    results[i, ] <- c(i, round(sd(test$Y), 2), round(mean(rmse), 2), 
                      round(mean(cvg), 2), round(mean(len), 2), round(mean(score), 2))
  } else {
    results[i, ] <- c(i, round(sd(test$Y[1:nTest, ]), 2), round(rmse[1], 2), 
                      round(cvg[1], 2), round(len[1], 2), round(score[1], 2))
  }

}
results




