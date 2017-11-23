GPDsample <- function(n,scale,shape,shapeV,scaleV,shapescaleV){
  library(lhs)
  
  X <- randomLHS(n,2) 
  params <- matrix(0, nrow=n, ncol=2)

  rho = shapescaleV/(sqrt(shapeV)*sqrt(scaleV)) # parameter correlation
  
  # Assume bivariate normally distributed parameter uncertainty 
  params[,1] <- qnorm( X[,1], mean=0, sd=1 ) * sqrt( shapeV ) + shape
  params[,2] <- qnorm( X[,2], mean=0, sd=1 ) * sqrt( scaleV ) * sqrt( 1 - rho^2 ) +
    (scale + rho * sqrt( scaleV ) / sqrt( shapeV ) * (params[,1] - shape) )
  
  params[,2] <- mapply( function(x) max( .0001, x ), params[,2] ) # remove negative

 return( params )
  
}