library(spBayes)

train <- readRDS("../data/scen4/train.RDS")
test <- readRDS("../data/scen4/test.RDS")

coords <- train$U
y <- train$Y[1:100, ]
X <- data.frame(train$X)
mm <- model.matrix(data = X, ~ X1 + X2)
p <- ncol(mm)
x1 <- X[,1]
x2 <- X[,2]

# Priors
priors <- list(
  beta.Norm = list(rep(0, p), diag(1000, p)),
  sigma.sq.IG = c(2, 1),
  tau.sq.IG = c(2, 1),
  phi.Unif = c(3 / max(dist(coords)),
               3 / min(dist(coords)))
)

# Starting values
starting <- list(
  beta = rep(0, p),
  sigma.sq = var(y),
  tau.sq = var(y) / 10,
  phi = mean(priors$phi.Unif)
)

# Tuning parameters
tuning <- list(
  sigma.sq = 0.1,
  tau.sq = 0.1,
  phi = 0.1
)

# Fit model
gp.fit <- spLM(
  y ~ x1 + x2,
  coords = coords,
  starting = starting,
  tuning = tuning,
  priors = priors,
  cov.model = "exponential",
  n.samples = 10000,
  verbose = TRUE
)

# Recover posterior samples
gp.fit <- spRecover(
  gp.fit,
  start = 5000,
  thin = 5,
  verbose = FALSE
)

# Prediction
coords_test <- test$U
XTest <- data.frame(test$X)
mmTest <- model.matrix(data = XTest, ~ X1 + X2)
pred <- spPredict(
  gp.fit,
  pred.coords = coords_test,
  pred.covars = mmTest,
  start = 5000,
  thin = 5,
  verbose = FALSE
)
point.preds <- apply(pred$p.y.predictive.samples, 1, mean)
truth <- test$Y[1:25,]
sqrt(mean((point.preds - truth)^2))
sd(truth)

