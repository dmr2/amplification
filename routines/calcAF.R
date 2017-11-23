get_AF <- function(rp,z,Ne_hist,histCurveSamps,slrsamps) {
  
  q <- 1/rp
  
  # The height of the historical flood with frequency q
  hist_height <- z[which.min(abs(Ne_hist-q))] 
  
  mi <- NULL
  ElogNprojshiftuncht <- NULL
  
  randi <- sample(1:1000, 10000, replace=T) # random GPD sample
  
  for(iii in 1:length(slrsamps)){
    mi <- which.min(abs(z-(hist_height-slrsamps[iii]))) # the 
    ElogNprojshiftuncht[iii] = log(histCurveSamps[randi[iii],mi])
  }
  
  af = exp(ElogNprojshiftuncht-log(q))
  
  return ( list(AF=af,HistHeight=hist_height)  )
}