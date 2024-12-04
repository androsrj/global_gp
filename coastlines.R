dist2d <- function(a, b, c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  d
} 

west <- function(x) {
  1.7*x + 166.45
}
east <- function(x) {
  1.25*x + 132.5
}

x1.west <- -75
y1.west <- west(x1.west)

x2.west <- -74.8
y2.west <- west(x2.west)

x1.east <- -74.9
y1.east <- east(x1.east)

x2.east <- -74.6
y2.east <- east(x2.east)

#coords.subset$dist.west <- 0
coords.subset$dist.east <- 0
for (i in 1:nrow(coords.subset)) {
  loc <- c(coords.subset$x[i], coords.subset$y[i])
  
  #start <- c(x1.west, y1.west)
  #end <- c(x2.west, y2.west)
  #coords.subset$dist.west[i] <- dist2d(loc, start, end)
  
  start <- c(x1.east, y1.east)
  end <- c(x2.east, y2.east)
  coords.subset$dist.east[i] <- dist2d(loc, start, end)
}
