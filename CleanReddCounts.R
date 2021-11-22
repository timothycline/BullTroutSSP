library(here)
library(dplyr)
library(tidyr)

#Read in harmonized redd file
HarmonyRedds <- read.csv(here('Data','Raw','USGS_Harmonized_BLT_Redd_Data_BLTPatchID.csv'),header=T,stringsAsFactors = F) %>% 
  filter(Year != 999) %>% 
  mutate(SiteID = paste(gsub(" ", "", CoreArea, fixed = TRUE),gsub(" ", "", LocalPopulation, fixed = TRUE),gsub(" ", "", Waterbody, fixed = TRUE),sep='.')) %>%
  filter(!(CoreArea %in% c('Bitterroot River','West Fork Bitterroot River')))
HarmonyRedds[1:20,] #check

unique(HarmonyRedds$CoreArea) %>% sort()
#Pull out just the unique sites, etc
#Keep record of if any Kovach_ID's
ReddSiteMetadata <- HarmonyRedds %>% group_by(RecoveryUnit,CoreArea,LocalPopulation,BLTPatchID,SiteID) %>% summarize(Kovach_ID = paste(unique(na.omit(Kovach_ID)),collapse='.'))
ReddSiteMetadata
write.table(ReddSiteMetadata,file=here('Data','ReddSiteMetadata.csv'),row.names=F,sep=',')

#Add a column for data included in kovach
HarmonyRedds <- HarmonyRedds %>% mutate(Kovach_YN = ifelse(is.na(Kovach_ID),'N','Y'))
HarmonyRedds[1:20,] #check

#Extract Max Year included in kovach
maxKovachYear <- max(HarmonyRedds %>% filter(Kovach_YN == 'Y') %>% pull(Year))
maxKovachYear

#Pull the unique KovachID identifiers
KovachIDs <- HarmonyRedds %>% pull(Kovach_ID) %>% unique() %>% na.omit()
KovachIDs #check

#list of kovach IDs with all sub MT_ID within them
GroupedIDs<-lapply(KovachIDs, FUN=function(x){HarmonyRedds %>% filter(Kovach_ID==x) %>% pull(SiteID) %>% unique()})
names(GroupedIDs) <- KovachIDs

Nsubpops<-lapply(GroupedIDs,FUN=length) # list of how many subpops in a kovach ID
length(Nsubpops) #of kovach IDs
length(which(Nsubpops>1)) # number of kovach IDs with more than 1 subpop
sum(unlist(Nsubpops[which(Nsubpops>1)])) # number of subpops collapsed into the 24 kovach IDs

### If you want aggregated counts by KOVACH ID
if(F){
#for each kovach ID, subset the data, and aggregate if necessary
HarmonyRedds_Aggregated <- lapply(GroupedIDs,FUN=function(x){
    ind<-which(names(GroupedIDs)=='Kovach_58')  
    x <- GroupedIDs[[ind]]  
    subsub <- HarmonyRedds %>% filter(BLTPatchID %in% x) %>% mutate(Kovach_ID = unique(na.omit(.$Kovach_ID))) #fill all rows with the proper kovach ID
    subsub$Redds_Numeric[subsub$Redds_Numeric==0] <- -999 #write -999 for TRUE zero values to avoid omitting them in the next step
    xtab<-subsub %>% xtabs(Redds_Numeric ~ MT_ID + Year,.) #xtabs builds a nice table, but puts 0's for NAs, but its ok because true zeros are labeled above as -999
    xtab[xtab==0] <- NA #fill xtabs 0's with NA
    xtab[xtab==-999] <- 0 #replace the temporary -999 with true 0's
    
    #Return years that will be omitted because of missing counts in some sites
    OmittedNew <- apply(xtab,2,FUN=function(x){return(sum(is.na(x)))})
    OmittedNew <- names(OmittedNew)[which(OmittedNew>0)] %>% as.numeric()
    
    totalRedds <- colSums(xtab) #aggregate counts by year, if a site is missing it automatically replaces with NA
    
    #some counts were omitted from kovach either due to bad counts or missing counts within (not easy to distinguish)
    OmittedOld <- subsub %>% filter(Year < 2014 & Kovach_YN =='N') %>% pull(Year)
    
    OmittedTotal <- c(OmittedOld,OmittedNew) %>% unique() %>% sort()
    
    totalRedds <- totalRedds[which(!(names(totalRedds) %in% OmittedTotal))]
    
    otp <- data.frame(RecoveryUnit=rep(subsub$RecoveryUnit[1],length(totalRedds)),
                      CoreArea=rep(subsub$CoreArea[1],length(totalRedds)),
                      LocalPopulation = rep(subsub$LocalPopulation[1],length(totalRedds)),
                      KovachID = rep(subsub$Kovach_ID[1],length(totalRedds)),
                      MT_ID = rep(paste(x,collapse='.'),length(totalRedds)),
                      Year = as.numeric(names(totalRedds)),
                      Redds_Numeric = totalRedds)
    otp <- otp %>% na.omit()
    return(otp)
})

#returns the years omitted from each Kovach ID
HarmonyRedds_Omitted <- lapply(GroupedIDs,FUN=function(x){
  subsub <- HarmonyRedds %>% filter(MT_ID %in% x) %>% mutate(Kovach_ID = unique(na.omit(.$Kovach_ID))) #fill all rows with the proper kovach ID
    subsub$Redds_Numeric[subsub$Redds_Numeric==0] <- -999 #write -999 for TRUE zero values to avoid omitting them in the next step
    xtab<-subsub %>% xtabs(Redds_Numeric ~ MT_ID + Year,.) #xtabs builds a nice table, but puts 0's for NAs, but its ok because true zeros are labeled above as -999
    xtab[xtab==0] <- NA #fill xtabs 0's with NA
    xtab[xtab==-999] <- 0 #replace the temporary -999 with true 0's
    
    #Return years that will be omitted because of missing counts in some sites
    OmittedNew <- apply(xtab,2,FUN=function(x){return(sum(is.na(x)))})
    OmittedNew <- names(OmittedNew)[which(OmittedNew>0)] %>% as.numeric()
    
    totalRedds <- colSums(xtab) #aggregate counts by year, if a site is missing it automatically replaces with NA
    
    #some counts were omitted from kovach either due to bad counts or missing counts within (not easy to distinguish)
    OmittedOld <- subsub %>% filter(Year < 2014 & Kovach_YN =='N') %>% pull(Year)
    
    OmittedTotal <- c(OmittedOld,OmittedNew) %>% unique() %>% sort()
    
    return(OmittedTotal)
    
})

HarmonyRebuild <- HarmonyRedds_Aggregated %>% bind_rows(.id=NULL)

write.table(HarmonyRebuild,file=here('Data','Modified','KovachAggregated_LocalPopulation_ReddData_Through2020.csv'),row.names=F,sep=',')

write.table(HarmonyRebuild %>% pivot_wider(names_from = Year,values_from=Redds_Numeric,names_sort=T),file=here('Data','Modified','KovachAggregated_LocalPopulation_ReddData_Through2020_Wide.csv'),row.names=F,sep=',')

}

##Use all data other than what is filtered at the top of this script
if(T){
  ReddTable <- xtabs(HarmonyRedds)
}
