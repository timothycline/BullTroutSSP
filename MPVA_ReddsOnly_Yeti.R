#Multipopulation viability analysis using only bull trout redd data.
library(here)
library(dplyr)
library(rjags)
library(R2jags)
library(DescTools)
library(MCMCvis)



## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
## Revert to 'sequential' setup of PSOCK cluster in RStudio Console on macOS and R 4.0.0
# if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
#     Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
#   parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
# }

#Read in Redd data
Redds<-read.csv(here('Data','Raw','USGS_Harmonized_BLT_Redd_Data_BLTPatchID.csv'),header=T,stringsAsFactors = F) %>% 
  filter(Year != 999) %>% 
  mutate(SiteID = paste(gsub(" ", "", CoreArea, fixed = TRUE),gsub(" ", "", LocalPopulation, fixed = TRUE),gsub(" ", "", Waterbody, fixed = TRUE),sep='.')) %>%
  filter(!(CoreArea %in% c('Bitterroot River','West Fork Bitterroot River')))
  #filter(Kovach_ID != '' & CoreArea2015 != 'CANADA' & !is.na(LocalPopulation2015))
Redds$Redds[Redds$Redds==0] <- -999 #note true zeros

#Just the Site level metadata organization
Redds_Meta <- read.csv(here('Data','ReddSiteMetadata.csv'),header=T,stringsAsFactors = F)

#Build out Redd time series by MT_ID
Redd_X <- xtabs(Redds ~ SiteID + Year,data=Redds) #build time series table automatically adds zeros where a count was not performed in a year
Redd_X[Redd_X == 0] <- NA #rewrite those zeros as NA values
Redd_X[Redd_X == -999] <- 0 #add back in the true zeros

Redd_X <- Redd_X[rowSums(!is.na(Redd_X))>=5,]
#nrow(Redd_X)

#unique(Redd_X)

# pdf('LocalPop_ReddTimeSeries.pdf')
# for(i in 1:nrow(Redd_X)){
#   MT1<-row.names(Redd_X)[i]
#   CA1<-Redds$CoreArea[match(MT1,Redds$SiteID)]
#   LP1 <- Redds$LocalPopulation[match(MT1,Redds$SiteID)]
#   plot(as.numeric(colnames(Redd_X)),Redd_X[i,],type='o',pch=16,main=paste(CA1,MT1,sep='.'),ylab='Redd count')
# }
# dev.off()

# library(MARSS)
# m1 <- MARSS(Redd_X,model=list(m=1,R='diagonal and equal'),form='dfa')
# plot(as.numeric(colnames(Redd_X)),m1$states[1,],type='l',xlim=c(1990,2020),ylim=c(-4,4))
# points(as.numeric(colnames(Redd_X)),m1$states[1,]+m1$states.se[1,],type='l',lty=3)
# points(as.numeric(colnames(Redd_X)),m1$states[1,]-m1$states.se[1,],type='l',lty=3)
# abline(h=0)
# 
# Browns<-read.csv('~/Downloads/BrownTroutTrend_Ripped.csv')
# points(Browns[,1]-2,Browns[,2],type='l',col='darkred',lwd=3)


PatchCovars <- read.csv(here('Data','CovariateData','BLT_Patches_US_Metrics_Oct28_21.csv'),header=T,stringsAsFactors=F) %>% 
  filter(Subpatch_Flag == 998) #FullPatches Only
#997 a smaller patch the corresponds to where the fish actually spawn
#998 the full patch
#head(PatchCovars)


#JAGS models structure for MPVA
modelScript.name <- "MPVA_OBSERR_POIS_BLT.txt"
jagsscript <- cat("
model {
    
    for(p in 1:nPop){
    
      
      ## N ~ Initialize N first year
      N1[p] ~ dgamma(0.001, 0.001)
      N[p,startYear[p]-1] ~ dpois(N1[p]);
      
      #N1[p] ~ dunif(1,maxN1[p]);
      #lN[p,startYear[p]-1] <- log(N1[p]);
      #N[p,1] ~ dpois(N1[p]);
      
      
      #PROCESS MODEL
      for(t in startYear[p]:(endYear[p]+nYearPred)){
        #lN[p,t] ~ dnorm( (lN[p,t-1] + r[p,t] + phi[p,t]*exp(lN[p,t-1])) , pow(sigmaR[p],-2) );
        N[p,t] ~ dpois( min(1e6, max(0, N[p,t-1]) ) * exp(R[p,t]) );
        
        R[p,t] ~ dnorm(r[p,t] + phi[p,t] * max(0, N[p,t-1]) / extent[p], pow(sigmaR[LP[p]],-2));
      
        #R[p,t] ~ dnorm(r[p,t] + phi[p,t] * max(0, N[p,t-1]), pow(sigmaR[p],-2));
        #R[i,t] ~ dnorm(r[p,t] + phi[p,t] * max(0, N[p,t-1] + reintro[p,t-1]) / extent[p], pow(sigmaR[p],-2))
        
        r[p,t] <- b0r[p];
        phi[p,t] <- b0phi[p];
        
      }
      
      ## sigmaR[i] ~ environmental stochasticity in population i (standard deviation in realized growth rate)
      
      
      ##POP PARAMETERS FROM LOCAL POPS
        ## REGRESSION PARAMETERS FOR 
        #b0r[p] ~ dnorm( LP_b0r[LP[p]] , pow(sigma_LP_b0r[LP[p]],-2) );
        b0r[p] <- LP_b0r[LP[p]];
        
      ## REGRESSION PARAMETERS FOR CARRYING CAPACITY
        #b0phi[p] ~ dt(0, pow(1,-2), 1) T(,0);
        #b0phi[p] ~ dnorm( LP_b0phi[LP[p]] , pow(sigma_LP_b0phi[LP[p]],-2) ) T(,0);
        b0phi[p] <- LP_b0phi[LP[p]]
        
      #OBSERVATION MODEL  
      for(t in startYear[p]:endYear[p]){
        A[p,t] ~ dbinom(N[p,t],p);
        F[p,t] ~ dpois(FCR * d[p,t]);
        ReddCounts[p,t] ~ dpois(A[p,t] + F[p,t]);
      }
      
    }
    
    #PRIORS FOR OBS MODEL
    FCR ~ dgamma(0.001,0.001); #FALSE COUNT RATE
    p ~ dbeta(8,2); #PROB OF OBSERVING A REDD
    
    #PRIORS FOR PROCESS MODEL
    #musig ~ dunif(0,10);
    #sigsig ~ dunif(0,10);
    
    #LOCAL POP PARAMETERS FROM CORE AREA
    for(lp in 1:nLocalPops){
      LP_b0r[lp] ~ dnorm( CA_b0r[CA[lp]], pow(sigma_CA_b0r[CA[lp]],-2));
      LP_b0phi[lp] ~ dnorm( CA_b0phi[CA[lp]] , pow(sigma_CA_b0phi[CA[lp]],-2)) T(,0);
      
      sigmaR[lp] ~ dt(0, pow(1,-2), 1) T(0.01,);
      #sigma_LP_b0r[lp] ~ dt(0, pow(1,-2), 1) T(0,);
      #sigma_LP_b0phi[lp] ~ dt(0, pow(1,-2), 1) T(0,);
    }
    
    ##PRIORS FOR CORE AREAS
    for(ca in 1:nCoreAreas){
      CA_b0r[ca] ~ dunif(-2,2);
      CA_b0phi[ca] ~ dunif(-1,0);
      
      sigma_CA_b0r[ca] ~ dt(0, pow(1,-2), 1) T(0,);
      sigma_CA_b0phi[ca] ~ dt(0, pow(1,-2), 1) T(0,);
    }
    
}  
", 
file = here('Analysis','MPVA',modelScript.name))


NPRED<-0 #REFERS TO FORWARD PREDICTIONS
n_in_full <- cbind(matrix(NA,nrow=nrow(Redd_X),ncol=1),as.matrix(Redd_X),matrix(NA,nrow=nrow(Redd_X),ncol=NPRED))
#
n_in <- n_in_full[sample.int(nrow(n_in_full),10,replace=F),]

StartYears <- apply(n_in,1,FUN=function(x){return(which(!is.na(x))[1])})
EndYears <- rep(ncol(Redd_X)+1,nrow(Redd_X))

LocalPops <- Redds_Meta$LocalPopulation[match(row.names(n_in),Redds_Meta$SiteID)]
CoreAreas <- Redds_Meta$CoreArea[match(row.names(n_in),Redds_Meta$SiteID)]

LocalPops_Num <- as.numeric(factor(LocalPops))
CoreArea_Num <- as.numeric(factor(CoreAreas))

extent_in <- rep(1,nrow(n_in))

#Distance Surveyed in a p x t matrix
d_in <- matrix(1000,nrow=nrow(n_in),ncol=ncol(n_in))
jags.data <- list(ReddCounts = n_in,
                  d = d_in,
                  extent = extent_in,
                  startYear=StartYears,
                  endYear = EndYears,
                  nYearPred = NPRED,
                  nPop=nrow(n_in),
                  nLocalPops=length(LocalPops),
                  nCoreAreas=length(CoreAreas),
                  LP = LocalPops_Num,
                  CA = CoreArea_Num)
jags.params <- c("N",'p','FCR',
                 'LP_b0r','LP_b0phi',
                 'CA_b0r','CA_b0phi',
                 'sigma_CA_b0r','sigma_CA_b0phi',
                 'sigmaR')
# nIter <- 5000
# mod_lm <- R2jags::jags(jags.data, parameters.to.save = jags.params,
#           model.file = here('Analysis','MPVA',modelScript.name), n.chains = 3, n.burnin = nIter/2, n.thin = 10,
#           n.iter = nIter)
mod_lm <- jags.parallel(jags.data, parameters.to.save = jags.params,
                       model.file = here('Analysis','MPVA',modelScript.name), n.chains = 10, n.burnin = 50000, n.thin = 10,
                       n.iter = 100000)

saveRDS(mod_lm,file=here('mod_lm.RDS'))
#jm <- rjags::jags.model(here('Analysis','MPVA',modelScript.name),n.chains=3,data=jags.data)
#update(jm,2000)
#s<-coda.samples(jm, jags.params,n.iter=2000)

# Rhats<-mod_lm$BUGSoutput$summary[,'Rhat']
# Rhats[which(Rhats>1.1)]
# plot(mod_lm)
# 
# apply(mod_lm$BUGSoutput$sims.list$N,c(2,3),median)
# dim(mod_lm$BUGSoutput$sims.list$N)
# Nsummary<-MCMCsummary(s,params='N')
# 
# b0r_summary <- MCMCsummary(s,params=c('CA_b0r'))
# 
# Nmeans <- matrix(Nsummary$mean,nrow=nrow(Redd_X),ncol=ncol(Redd_X)+5)
# Nl95 <- matrix(Nsummary$`2.5%`,nrow=nrow(Redd_X),ncol=ncol(Redd_X)+5)
# Nh95 <- matrix(Nsummary$`97.5%`,nrow=nrow(Redd_X),ncol=ncol(Redd_X)+5)
# 
# Nsummary
# 
# #Npredsim<-mod_lm$BUGSoutput$sims.list$N
# #Npred<-apply(Npredsim,c(2,3),FUN=mean)
# 
# #names(mod_lm$BUGSoutput$sims.list)
# 
# plot(as.numeric(names(Redd_X[5,])),Redd_X[5,],pch=16,xlim=c(1979,2025))
# points(1979:2025,Nmeans[5,],type='l',col='darkorange')
# 
# plot(Redd_X[5,])
# points(Nmeans[5,],type='l')
# 
# s$trace$b0phi
