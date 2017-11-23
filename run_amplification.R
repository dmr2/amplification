#!/opt/local/bin/R
# run_amplification.R
# 
# last updated: Fri Sep  1 03:01:48 EDT 2017
# 
#  Calculate flood amplification factors for
#  given return periods and selected years for
#  multiple global mean surface temperature
#  scenarios.
#
#  Inputs:
#   - a file with GPD parameters for each site
#   - Monte Carlo sea-level rise projections 
#
#  Outputs:
#   - writes amplification factors to a table

rm(list=ls(all=TRUE))
setwd("/Users/dmr/Dropbox/IPCC Sea Level/amplification")

source("GPDsample.R")
source("GPDNExceed.R")
source("openSLR.R")
source("calcAF.R")


# Location of local sea level rise Monte Carlo samples
dir = "/Users/dmr/Dropbox/IPCC\ Sea\ Level/slr_samples"

targ_years = seq(2010,2150,10) 

scenarios = c("1.5C","2.0C","2.5C")
slab = c("1p5degree","2p0degree","2p5degree")
slab_out = c("1p5degree","2p0degree","2p5degree")

rp_list = c(10,20,50,100,500)

# Open GPD parameters and historical flood data
dat <- read.csv("GPDfits_uhawaii_projectLSLfmt.tsv",header=T,sep="\t")
sites <- dat$Site # tide gauge sites

AFexp <- array(NaN,dim=c(length(sites),length(targ_years),
                      length(scenarios),length(rp_list)))
AFq17 <- array(NaN,dim=c(length(sites),length(targ_years),
                         length(scenarios),length(rp_list)))
AFq83 <- array(NaN,dim=c(length(sites),length(targ_years),
                         length(scenarios),length(rp_list)))

histFloodHeight <- array(NaN,dim=c(length(sites),length(targ_years),
                                  length(scenarios),length(rp_list)))
FutureFreq <- array(NaN,dim=c(length(sites),length(targ_years),
                                   length(scenarios),length(rp_list)))
slr <- array(NaN,dim=c(length(sites),length(scenarios),length(targ_years),3))


for(j in 1:length(sites)){
#for(j in which(sites=="The Battery, NY"):which(sites=="The Battery, NY")){ # Only process "The Battery, NYC"

  site <- dat$Site[j] # tide gauge site name
  scale <- dat$Scale50[j] # median scale parameter
  shape <- dat$Shape50[j] # median shape parameter
  UHid <- dat$UHAWAII_ID[j]
  threshold <- dat$Q99[j] # GPD threshold
  lambda <- dat$Lambda[j] # mean Poisson arrival rate of threshold
  shapescaleV <- dat$Vscaleshape[j] # covariance of GPD scale and shape parameters
  shapeV <- dat$Vshape[j] # variance of shape parameter
  scaleV <- dat$Vscale[j] # variance of scale parameter
  gauge <- dat$PSMSL_ID[j] # Tide gauge ID
  basin <- dat$Basin[j]
  
 # Account for GPD parameter uncertainty by making draws from a
 # bivariate normal distribution using Latin hypercube sampling
  GPD <- GPDsample(1000, scale, shape, shapeV, scaleV, shapescaleV)

  z <- seq(0,10,.01) # some flood heights (meters above tide gauge MHHW)
  
  # Expected historical flood height curve (No SLR) (GPD uncertainty)
  histCurveSamps <- matrix(NaN, nrow=length(GPD[,2]), ncol=length(z))
  for(iii in 1:length(GPD[,2]) ){
    histCurveSamps[iii,] <- GPDNExceed(z-threshold,lambda,-threshold,GPD[iii,2],GPD[iii,1])
  }
  
  Ne_hist <- apply(histCurveSamps,2,mean,na.rm=T) 
  
  for( s in 1:length(scenarios) ){
  
    # Get sea level rise Monte Carlo Samples
    fil <- paste( dir,"/LSLproj_MC_",gauge,"_",slab[s],".tsv",sep="")
    if (!file.exists(fil)){
      fil <- paste( dir,"/LSLprojMC_",gauge,"_",slab[s],".tsv",sep="")
    }
    SLR <- getRSLsamps( fil )
    
    SLRMC <- SLR$samples
    years <- SLR$years 
  
    for( t in 1:length(targ_years) ){
      
      print(paste("[",dat$Site[j],"] Calculating AFs for year: ",targ_years[t]," Scenario: ", scenarios[s]))
      indx = which(years==targ_years[t])
      
      slr[j,s,t,] <- quantile(SLRMC[,indx],probs=c(.17,.5,.83))
      
      for (r in 1:length(rp_list)){
        
       print(paste("   ",rp_list[r]))
       af <- get_AF(rp_list[r],z,Ne_hist,histCurveSamps,SLRMC[,indx])

       AFexp[j,t,s,r] <- mean(af$AF)
       FutureFreq[j,t,s,r] <- AFexp[j,t,s,r] * 1/rp_list[r]
       
       # Mask AF if dominated by tidal events (i.e., every other day)
       limit <- (365.25/2)/(1/rp_list[r])
       if( AFexp[j,t,s,r] >= limit ){
         AFexp[j,t,s,r] = -999.999
         FutureFreq[j,t,s,r] <- -999.999
       }
       
       #AFq17[j,t,s,r] <- quantile(af$AF,probs = c(.05,0.16666,.5,0.8333,.95))
       #AFq83[j,t,s,r] <- quantile(af$AF,probs = 0.83333)
       histFloodHeight[j,t,s,r] <- af$HistHeight
      }
      
    } # each target year
  } # each scenario
}# each tide gauge

# Write all AFs for all sites to a table
for( s in 1:length(scenarios) ){
  for( t in 1:length(targ_years) ){
    
    # put the data we want into data frames for writing to disk
    
    df <- data.frame(Site=dat$Site,Country=dat$Country,Region=dat$Region,
                     ID=dat$PSMSL_ID,UHAWAII_ID=dat$UHAWAII_ID,
                     Lat=dat$Lat,Lon=dat$Lon,
                     Scenario=rep(slab_out[s],length(dat$Site)),
                     Year=rep(targ_years[t],length(dat$Site)),
                     SLR_17=slr[,s,t,1],SLR_50=slr[,s,t,2],SLR_83=slr[,s,t,3],
                     GPD_threshold=dat$Q99,
                     Hist_Freq_10=rep(.1,length(dat$Site)),Hist_Height_10=histFloodHeight[,t,s,1],
                     Fut_Freq_10=FutureFreq[,t,s,1],AF_10=AFexp[,t,s,1],
                     Hist_Freq_20=rep(.05,length(dat$Site)),Hist_Height_20=histFloodHeight[,t,s,2],
                     Fut_Freq_20=FutureFreq[,t,s,2],AF_20=AFexp[,t,s,2],
                     Hist_Freq_50=rep(.02,length(dat$Site)),Hist_Height_50=histFloodHeight[,t,s,3],
                     Fut_Freq_50=FutureFreq[,t,s,3],AF_50=AFexp[,t,s,3],
                     Hist_Freq_100=rep(.01,length(dat$Site)),Hist_Height_100=histFloodHeight[,t,s,4],
                     Fut_Freq_100=FutureFreq[,t,s,4],AF_100=AFexp[,t,s,4],
                     Hist_Freq_500=rep(.002,length(dat$Site)),Hist_Height_500=histFloodHeight[,t,s,5],
                     Fut_Freq_500=FutureFreq[,t,s,5],AF_500=AFexp[,t,s,5])
                    
    outf <- paste("output/af_stab_",slab_out[s],"_sealevel_",targ_years[t],".csv",sep="")
    write.table(df,outf,sep=",",row.names = FALSE)
  }
}
