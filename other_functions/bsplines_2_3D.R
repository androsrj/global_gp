library(splines)

Bsplines_2D = function(X, df = NULL, knots = NULL, degree = 3)
{
  # df = a 2 vector giving number of degrees of freedom in each direction, if NULL, chosen as length(knots) + degree
  # knots = list of internal brakpoints in each direction. if NULL chosen as quantiles
  # eg knots = list(c(1,2,3),c(1,2,3))
  # eg df = c(5,5)
  # atleast one of df or knots is needed.
  
  basis_x = bs(X[,1], df[1], knots[[1]], degree, Boundary.knots = c(min(X) - 1, max(X) + 1))
  basis_y = bs(X[,2], df[2], knots[[2]], degree, Boundary.knots = c(min(X) - 1, max(X) + 1))
  
  prodbasis = function(bx,by) as.numeric(outer(bx,by))
  basis_xy = matrix(0, nrow(basis_x), ncol(basis_x) * ncol(basis_y))
  for(i in 1:nrow(basis_x))
  {
    basis_xy[i,] = prodbasis(basis_x[i,],basis_y[i,])
  }
  return(basis_xy)
}

Bsplines_3D = function(X, df = NULL, knots = NULL, degree = 3)
{
  # df = a 3 vector giving number of degrees of freedom in each direction, if NULL, chosen as length(knots) + degree
  # knots = list of internal brakpoints in each direction. if NULL chosen as quantiles
  # eg knots = list(c(1,2,3),c(1,2,3),c(1,2,3))
  # eg df = c(5,5,5)
  # atleast one of df or knots is needed.
  
  basis_x = bs(X[,1], df[1], knots[[1]], degree)
  basis_y = bs(X[,2], df[2], knots[[2]], degree)
  basis_z = bs(X[,3], df[3], knots[[3]], degree)
  
  prodbasis = function(bx,by) as.numeric(outer(bx,by))
  basis_xyz = matrix(0, nrow(basis_x), ncol(basis_x) * ncol(basis_y)* ncol(basis_z))
  for(i in 1:nrow(basis_x))
  {
    basis_xyz[i,] = prodbasis(prodbasis(basis_x[i,],basis_y[i,]),basis_z[i,])
  }
  return(basis_xyz)
}


# 2D example
#X = matrix(runif(2000),1000,2)
#basis = Bsplines_2D(X, df = c(7,7))
#z = apply(X,1,function(x) sin(x[1]) + tan(x[2]^2))
#lmod = lm(z ~ basis)
#summary(lmod)

# 3D example
#X = matrix(runif(15000),5000,3)
#basis = Bsplines_3D(X, df = c(7,7,7))
#z = apply(X,1,function(x) sin(x[1]) + tan(x[2]^2) + cos(prod(x)))
#lmod = lm(z ~ basis)
#summary(lmod)
