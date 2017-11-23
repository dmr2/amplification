GPDNExceed <- function(z,lambda,MHHW,scale,shape){
  
  # Function "GPDNExceed"
  #
  # Calculate log of the number of exceedances of z from a
  # Poisson-Generalized Pareto Distribution with Poisson mean
  # lambda and the specified shape and scale factor. Assumes
  # the exceedance threshold has already been removed from z.
  #
  # This function ensures that values returned stay within the
  # support of the GPD.
  #
  # For values of z below zero, the function will treat as though
  # z = 0 unless MHHW is specified. If MHHW is specified
  # value, exceedances below zero will be assumed to fall on a
  # Gumbel distribution between lambda exceedances at zero and a
  # specified value at z = MHHW(1) < 0. If MHHW(2) exists, it is
  # the value of exceedances at z = MHHW(1); otherwise, will default
  # to 365.25/2.
  
  # Reference:
  # M.K. Buchanan, R.E. Kopp, M. Oppenheimer, and C. Tebaldi. (2016).
  # Allowances for evolving coastal flood risk under uncertain local
  # sea-level rise. Climatic Change.
  
  # History
  # 1/25/2017 (DMR): Function created
  #

  exponent <- function(a, pow) (abs(a)^pow)*sign(a)
  
  z0 <- z

  z <- pmax(z,0) # > 0 only
  
  if ( shape == 0 ){ 
    
    result <- mapply( function(x) lambda*exp(-1*x/scale), z)
    
  }else if( shape < 0 ){ 
    
  #  z <- mapply( function(x) min(x, .99999* -scale/shape), z) # values stay within the support of the GPD
    z <- pmin(z, .99999* -scale/shape) # values stay within the support of the GPD
    result <- mapply( function(x) lambda* exponent( 1 + (shape*x /scale), -1/shape), z)
    
  }else{
    
    result <- mapply( function(x) lambda* exponent( 1 + (shape*x /scale), -1/shape), z)
    
  }

  # Put flood heights below 99.9th threshold and greater than MHHW on a Gumbel distribution
  if( length(MHHW) >= 1 ){
    
   # z0 <- mapply( function(x) max(-MHHW[1], x), z0)
   z0 <- pmax(MHHW[1], z0)
    indx <- which( z0 < 0 )
    
    if( length( MHHW ) >= 2 ){
      
      exceedMHHW <- MHHW[2]
      
    }else{
      
      exceedMHHW = 365.25/2
    }
    
    #result[indx] <- lambda*exponent( exceedMHHW/lambda, (z0[indx]/ -MHHW[1]) )
    #print(lambda*exponent( exceedMHHW/lambda, (-0.4165/ MHHW[1]) ))
    result[indx] <- lambda*exponent( exceedMHHW/lambda, (z0[indx]/ MHHW[1]) )
  }

#  if( sum(is.nan(result) ) > 0 ){
    
 #   stop( paste("result:",result,"poisson:",lambda,"mu:",MHHW,"scale:",scale,"shape:",shape,sep=" ") )
    
 # }else{
    
    return(result)
    
 # }
  
}