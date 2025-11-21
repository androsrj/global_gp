library(ggplot2)
library(MBA)
library(fields)
library(latex2exp)

line.type <- 2
line.width <- 4
which.scens1 <- c(1, 4, 7)
which.scens2 <- 11:13
nScen <- length(which.scens1) + length(which.scens2)


#################################
############# Beta0 #############
#################################

pdf("figures/new/beta0.pdf", width = 5.5, height = 12)
layout.matrix <- matrix(c(1, 5, 
                          2, 6,
                          3, 7,
                          4, 8), nrow = 4, ncol = 2, byrow = TRUE)
layout(layout.matrix)

# Data for scens 1-10
test <- readRDS(paste0("data/small/scen", 1, "/test.RDS"))
nTest <- nrow(test$B)
beta0.true <- test$U[ , 1] - test$U[ , 2]

# True beta0, scens 1-10
par(mar = c(3, 3, 5.5, 2) + 0.1)
mba.data.true <- data.frame(test$U, beta0.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = "(Scenarios 1 - 10)",
      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))
mtext(TeX("True surface of $\\beta_0$"), side = 3, line = 3.5)

# Est. beta, scens 1, 4, 7
par(mar = c(3, 3, 3, 2) + 0.1)
for (i in which.scens1) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta0.means <- results$posteriorMeans$beta.test[1:nTest]
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
}

# Data for scens 11-13
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
nTest <- nrow(test$B)
beta0.true <- test$B[ , 1]

# True beta0, scens 11-13
par(mar = c(3, 3, 5.5, 2) + 0.1)
mba.data.true <- data.frame(test$U, beta0.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = "(Scenarios 11 - 13)",
      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))
mtext(TeX("True surface of $\\beta_0$"), side = 3, line = 3.5)

# Est. beta, scens 11-13
par(mar = c(3, 3, 3, 2) + 0.1)
for (i in which.scens2) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta0.means <- results$posteriorMeans$beta.test[1:nTest]
  mba.data <- data.frame(test$U, beta0.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
}

dev.off()




#################################
############# Beta1 #############
#################################

pdf("figures/new/beta1.pdf", width = 5.5, height = 12)
layout.matrix <- matrix(c(1, 5, 
                          2, 6,
                          3, 7,
                          4, 8), nrow = 4, ncol = 2, byrow = TRUE)
layout(layout.matrix)

# Data for scens 1-10
test <- readRDS(paste0("data/small/scen", 1, "/test.RDS"))
nTest <- nrow(test$B)
beta1.true <- test$U[ , 1] + test$U[ , 2] - 100

# True beta1, scens 1-10
par(mar = c(3, 3, 5.5, 2) + 0.1)
mba.data.true <- data.frame(test$U, beta1.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = "(Scenarios 1 - 10)",
      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))
mtext(TeX("True surface of $\\beta_1$"), side = 3, line = 3.5)

# Est. beta, scens 1, 4, 7
par(mar = c(3, 3, 3, 2) + 0.1)
for (i in which.scens1) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta1.means <- results$posteriorMeans$beta.test[(nTest+1):(2*nTest)]
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
}

# Data for scens 11-13
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
nTest <- nrow(test$B)
beta1.true <- test$B[ , 2]

# True beta1, scens 11-13
par(mar = c(3, 3, 5.5, 2) + 0.1)
mba.data.true <- data.frame(test$U, beta1.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = "(Scenarios 11 - 13)",
      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))
mtext(TeX("True surface of $\\beta_1$"), side = 3, line = 3.5)

# Est. beta, scens 11-13
par(mar = c(3, 3, 3, 2) + 0.1)
for (i in which.scens2) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta1.means <- results$posteriorMeans$beta.test[(nTest+1):(2*nTest)]
  mba.data <- data.frame(test$U, beta1.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
}

dev.off()





#################################
############# Beta2 #############
#################################

pdf("figures/new/beta2.pdf", width = 5.5, height = 12)
layout.matrix <- matrix(c(1, 5, 
                          2, 6,
                          3, 7,
                          4, 8), nrow = 4, ncol = 2, byrow = TRUE)
layout(layout.matrix)

# Data for scens 1-10
test <- readRDS(paste0("data/small/scen", 1, "/test.RDS"))
nTest <- nrow(test$B)
beta2.true <- 2 * test$U[ , 1] - test$U[ , 2] - 50

# True beta2, scens 1-10
par(mar = c(3, 3, 5.5, 2) + 0.1)
mba.data.true <- data.frame(test$U, beta2.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = "(Scenarios 1 - 10)",
      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))
mtext(TeX("True surface of $\\beta_2$"), side = 3, line = 3.5)

# Est. beta, scens 1, 4, 7
par(mar = c(3, 3, 3, 2) + 0.1)
for (i in which.scens1) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta2.means <- results$posteriorMeans$beta.test[(2*nTest+1):(3*nTest)]
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
}

# Data for scens 11-13
test <- readRDS(paste0("data/small/scen", 11, "/test.RDS"))
nTest <- nrow(test$B)
beta2.true <- test$B[ , 3]

# True beta0, scens 11-13
par(mar = c(3, 3, 5.5, 2) + 0.1)
mba.data.true <- data.frame(test$U, beta2.true)
mba.interp.true <- mba.surf(mba.data.true, no.X=100, no.Y=100, extend=TRUE)
image(mba.interp.true$xyz.est, main = "(Scenarios 11 - 13)",
      cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, 
      col = tim.colors(64))
mtext(TeX("True surface of $\\beta_2$"), side = 3, line = 3.5)

# Est. beta, scens 11-13
par(mar = c(3, 3, 3, 2) + 0.1)
for (i in which.scens2) {
  path <- paste0("objects/small_scen", i, ".RDS") 
  results <- readRDS(path)[[1]]
  test <- readRDS(paste0("data/small/scen", i, "/test.RDS"))
  nTest <- nrow(test$B)
  beta2.means <- results$posteriorMeans$beta.test[(2*nTest+1):(3*nTest)]
  mba.data <- data.frame(test$U, beta2.means)
  mba.interp <- mba.surf(mba.data, no.X=100, no.Y=100, extend=TRUE)
  image(mba.interp$xyz.est, main = paste0("Scenario ", i), 
        col = tim.colors(64), cex.main = 1.5,
        axes = FALSE)
}

dev.off()
