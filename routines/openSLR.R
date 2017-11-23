getRSLsamps <- function(fil){
  
  # Open relative sea-level rise Monte Carlo samples (K14 format)
  
  # Output:
  #   - samples[ MC RSL samples (in meters), Years ]
  #   - Years
  
  
  print(" ")
  print( paste( "Opening SL MC samples:", fil) )
  print(" ")
  
  x = read.table( paste(fil, sep=""), skip=1, sep="\t", header=FALSE)
  
  years = x[,1] # Years corresponding to SLR
  
  iii <- as.data.frame( t(x) ) # Transpose to [SAMPLES,YEARS]
  iii <- iii[2:10001,]/1000 # millimeters to meters
  
  # remove samples that may be physically implausible (i.e., the 99.9th percentile)
  q99 <- apply( iii, 2, quantile, probs=c(0.999), na.rm=T)
  samples <- matrix( NaN, nrow=10000, ncol=length(years) )
  
  for(t in 1:length( years )){
    samples[,t] <- pmin(iii[,t],q99[t])
  }
  
  return ( list(samples=samples, years=years) )
  
}