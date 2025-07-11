### JACOBIAN AND TRANSFORMATIONS (for theta) ###

# Transformation for theta, where trTheta = log((theta - a) / (b - theta))
g <- function(x, a = 0.01, b = 20) {
  log((a - x) / (x - b))
}

# Inverse transformation for theta
gInv <- function(trTheta, a = 0.01, b = 20) {
  (b * exp(trTheta) + a) / (1 + exp(trTheta))
}

# Log-Jacobian for theta, (log-derivative of gInv function above)
logJac <- function(trTheta, a = 0.01, b = 20) {
  # log( (b - a) * exp(trTheta) / (1 + exp(trTheta))^2 ) # or simplify, as below
  log(b - a) + trTheta - 2 * log(1 + exp(trTheta))
}
