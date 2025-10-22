library(randomForest)
train <- readRDS("../data/small/scen1/train.RDS")
test <- readRDS("../data/small/scen1/test.RDS")
nSubj <- nrow(train$Z)
n <- nrow(train$X)
nTest <- nrow(test$X)

# Compile training data
x.train <- matrix(rep(train$X, nSubj), byrow = TRUE, ncol = 2)
u.train <- matrix(rep(train$U, nSubj), byrow = TRUE, ncol = 2)
z.train <- matrix(rep(train$Z, n), byrow = TRUE, ncol = 2)
rf.train <- data.frame(x.train, u.train, z.train, train$Y)
colnames(rf.train) <- c("X1", "X2", "U1", "U2", "Z1", "Z2", "Y")

# Compile testing data
x.test <- matrix(rep(test$X, nSubj), byrow = TRUE, ncol = 2)
u.test <- matrix(rep(test$U, nSubj), byrow = TRUE, ncol = 2)
z.test <- matrix(rep(test$Z, nTest), byrow = TRUE, ncol = 2)
rf.test <- data.frame(x.test, u.test, z.test, test$Y)
colnames(rf.test) <- c("X1", "X2", "U1", "U2", "Z1", "Z2", "Y")

# Fit model
rf.model <- randomForest(Y ~ ., data = rf.train, ntree = 1000, mtry = 2)
print(rf.model)

# Get out-of-sample predictions
pred <- predict(rf.model, rf.test)

# Test error
sqrt(mean((test$Y - pred)^2))
sd(test$Y)

# Train error
sqrt(mean((train$Y - rf.model$predicted)^2))
sd(train$Y)








