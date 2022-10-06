#Reading in files output by analysis loop, manipulating data into summary data frames, file pathways made relative to '.Data/' need to be checked#
library(dplyr)
library(plyr)
library(sp)
library(rgdal)
require(terra)
require(raster)
require(ggmap)
require(maptools)
require(rgeos)
require(ecospat)
require(scales)
library(ade4)
library(reshape2)
library(tidyr)



##COLATING AND PLOTTING ECOSPAT ANALYSIS OUTPUTS ##

SimShiftNA<-list.files(path = ".Data//Europe_ebba_species/Niche_Analyses_summaries",
                       pattern = "Niche_similarity_shift_NA_test.Rds",
                       recursive = TRUE,
                       full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(SimShiftNA))
file_names <- substr(file_names,1,nchar(file_names)-35)
#Extract data
SimShiftNA <- lapply(SimShiftNA, function(i){
  readRDS(i)})
#give names from file name
names(SimShiftNA)<-file_names

SimShiftInt<-list.files(path = ".Data//Europe_ebba_species/Niche_Analyses_summaries",
                       pattern = "Niche_similarity_shift_int_test.Rds",
                       recursive = TRUE,
                       full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(SimShiftInt))
file_names <- substr(file_names,1,nchar(file_names)-36)
#Extract data
SimShiftInt <- lapply(SimShiftInt, function(i){
  readRDS(i)})
#give names from file name
names(SimShiftInt)<-file_names


SimConsNA<-list.files(path = ".Data//Europe_ebba_species/Niche_Analyses_summaries",
                      pattern = "Niche_similarity_conservatism_NA_test.Rds",
                      recursive = TRUE,
                      full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(SimConsNA))
file_names <- substr(file_names,1,nchar(file_names)-42)
#Extract data
SimConsNA <- lapply(SimConsNA, function(i){
  readRDS(i)})
#give names from file name
names(SimConsNA)<-file_names

SimConsInt<-list.files(path = ".Data//Europe_ebba_species/Niche_Analyses_summaries",
                      pattern = "Niche_similarity_conservatism_int_test.Rds",
                      recursive = TRUE,
                      full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(SimConsInt))
file_names <- substr(file_names,1,nchar(file_names)-43)
#Extract data
SimConsInt <- lapply(SimConsInt, function(i){
  readRDS(i)})
#give names from file name
names(SimConsInt)<-file_names


TestMetrics<-list(SimShiftNA, SimShiftInt, SimConsNA, SimConsInt)
names(TestMetrics)<-c("SimShiftNA","SimShiftInt", "SimConsNA", "SimConsInt")



niche.dyn<-list.files(path = ".Data//Europe_ebba_species/",
                      pattern = "niche.dyn.Rds",
                      recursive = TRUE,
                      full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(niche.dyn))
file_names <- substr(file_names,1,nchar(file_names)-14)
#Extract data
niche.dyn <- lapply(niche.dyn, function(i){
  readRDS(i)})
#give names from file name
names(niche.dyn)<-file_names


HistRangeCent<-list.files(path = ".Data//Europe_ebba_species/",
                          pattern = "HistRangeCent.Rds",
                          recursive = TRUE,
                          full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(HistRangeCent))
file_names<-substr(file_names,1,nchar(file_names)-18)
#Extract data
HistRangeCent <- lapply(HistRangeCent, function(i){
  readRDS(i)})
#give names from file name
names(HistRangeCent)<-file_names


CurrentRangeCent<-list.files(path = ".Data//Europe_ebba_species/",
                             pattern = "CurrentRangeCent.Rds",
                             recursive = TRUE,
                             full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(CurrentRangeCent))
file_names<-substr(file_names,1,nchar(file_names)-21)
#Extract data
CurrentRangeCent <- lapply(CurrentRangeCent, function(i){
  readRDS(i)})
#give names from file name
names(CurrentRangeCent)<-file_names

RangeCentShift<-list.files(path = ".Data//Europe_ebba_species/",
                           pattern = "RangeCentShift.Rds",
                           recursive = TRUE,
                           full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(RangeCentShift))
file_names<-substr(file_names,1,nchar(file_names)-19)
#Extract data
RangeCentShift <- lapply(RangeCentShift, function(i){
  readRDS(i)})
#give names from file name
names(RangeCentShift)<-file_names


HistRangeSize<-list.files(path = ".Data//Europe_ebba_species/",
                          pattern = "HistRangeSize.Rds",
                          recursive = TRUE,
                          full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(HistRangeSize))
file_names<-substr(file_names,1,nchar(file_names)-18)
#Extract data
HistRangeSize <- lapply(HistRangeSize, function(i){
  readRDS(i)})
#give names from file name
names(HistRangeSize)<-file_names

CurrentRangeSize<-list.files(path = ".Data//Europe_ebba_species/",
                             pattern = "CurrentRangeSize.Rds",
                             recursive = TRUE,
                             full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(CurrentRangeSize))
file_names<-substr(file_names,1,nchar(file_names)-21)
#Extract data
CurrentRangeSize <- lapply(CurrentRangeSize, function(i){
  readRDS(i)})
#give names from file name
names(CurrentRangeSize)<-file_names

RangeExpansion<-list.files(path = ".Data//Europe_ebba_species/",
                           pattern = "RangeExpansion.Rds",
                           recursive = TRUE,
                           full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(RangeExpansion))
file_names<-substr(file_names,1,nchar(file_names)-19)
#Extract data
RangeExpansion <- lapply(RangeExpansion, function(i){
  readRDS(i)})
#give names from file name
names(RangeExpansion)<-file_names

RangeContraction<-list.files(path = ".Data//Europe_ebba_species/",
                             pattern = "RangeContraction.Rds",
                             recursive = TRUE,
                             full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(RangeContraction))
file_names<-substr(file_names,1,nchar(file_names)-21)
#Extract data
RangeContraction <- lapply(RangeContraction, function(i){
  readRDS(i)})
#give names from file name
names(RangeContraction)<-file_names

RangeStability<-list.files(path = ".Data//Europe_ebba_species/",
                           pattern = "RangeStability.Rds",
                           recursive = TRUE,
                           full.names = TRUE)
#retrieve file names for objects
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(RangeStability))
file_names<-substr(file_names,1,nchar(file_names)-19)
#Extract data
RangeStability <- lapply(RangeStability, function(i){
  readRDS(i)})
#give names from file name
names(RangeStability)<-file_names


####Grouping, averaging, comparing and plotting outputs####
##coerce loop lists into datframe to simplify data handling##
df<-melt(TestMetrics)#Coerce to datframe
df<-subset(df,L3!="sim")#remove simulation results data
df1<-dcast(subset(df,L3!="obs"),L1 + L2 ~ L3)#extract columns of p-values and rows of species names
df<-dcast(subset(df,L3=="obs"), L1 + L2 ~ L4)#extract rows of species and columns of output values D, expansion, stability, unfilling and I

df<-merge(df,df1,by=c("L1","L2"))

names(df)[names(df)=="L1"]<-"Test"
names(df)[names(df)=="L2"]<-"species"



write.csv(df,".Data//Europe_ebba_species/Niche_Analyses_summaries/All_test_data.csv")


###ANNA's proportioning method###

# empty df to store data in
rel_niche_dynamics <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(rel_niche_dynamics) <- c("species", "test", "metric", "proportion", "total")

df1<-subset(df,substr(df$Test, 1, 8) == "SimShift")

specs <- unique(df1$species)

for(spec in specs){
  
    
    inter <- subset(df1, species == spec & substr(df1$Test, nchar(df1$Test)-2, nchar(df1$Test)) == "Int")
    whole <- subset(df1, species == spec & substr(df1$Test, nchar(df1$Test)-1, nchar(df1$Test)) == "NA")
    
    e_a <- inter$expansion
    s_a <- inter$stability
    u_a <- inter$unfilling
    
    e <- whole$expansion
    s <- whole$stability
    u <- whole$unfilling
    
    
    # native niche (non-analogue)
    s_x <- 1 - u
    
    # propotion non-analogue invaded niche stability and native stability
    x <- s / s_x
    
    # total niche (non-analogue)
    t <- u * x + s + e 
    
    # proportion of unfilling, stability and expansion towards total non-analogue niche
    u_rel <- u * x / t
    s_rel <- s / t
    e_rel <- e / t
    
    # native niche (analogue)
    s_y <- 1 - u_a
    
    # propotion analogue invaded niche stability and native stability
    y <- s_a / s_y
    
    # total niche (analogue)
    t_a <- u_a * y + s_a + e_a 
    
    # proportion of unfilling, stability and expansion towards total analogue niche
    s_rel_a1 <- s_a / t_a
    e_rel_a1 <- e_a / t_a
    u_rel_a1 <- u_a * y / t_a
    
    
    i <- s_rel / s_rel_a1
    
    # define target values 
    s_rel_a <- s_rel_a1 * i # stability (%)
    u_rel_a <- u_rel_a1 * i # unfilling (%)
    e_rel_a <- e_rel_a1 * i # expansion (%)
    
    a <- u_rel - u_rel_a # abandonment (%)
    p <- e_rel - e_rel_a # pioneering (%)
    
    results <- c(s_rel_a, u_rel_a, e_rel_a, a, p)
    
    #  to check whether values add up to 1
    total <- a + u_rel_a + s_rel_a + e_rel_a + p
    
    
    # add info to df
    rel_niche_dynamics <- rbind(rel_niche_dynamics,
                                data.frame(species = rep(spec, 5), test=rep(substr(inter$Test, 1, nchar(inter$Test)-3), 5),
                                           metric = fin_metrics, percentage = results,
                                           total = rep(total, 5)))
    
    
 } # end loop over species

#Cwrite output to file
write.csv(rel_niche_dynamics,".Data//Europe_ebba_species/Niche_Analyses_summaries/relative_niche_indeces.csv")

#Test types determine the p values for D, I, stability, unfilling etc. but actual values remain unchanged, can examine a bit below, but no impact on following figures, so can drop 3/4 tests before moving on


##Add Niche centroids and shifts to analyses dataframe##
dyndf<-melt(niche.dyn)#Coerce to dataframe
dyndf<-subset(dyndf,!L2%in%c("dyn","dynamic.index.w"))#isolate niche centroid data
dyndf<-dcast(dyndf, L1 ~ L2)
names(dyndf)[names(dyndf)=="L1"]<-"species"
df<-subset(df,df$Test=="SimShiftNA")##Restroct to one test type to remove duplication
df<-merge(df,dyndf,by='species',all=T)

###Plot trends in buffer sizes across species for all key output metrics###
bufferdf<-df
bufferdf$Buffer<-substr(bufferdf$species,nchar(bufferdf$species)-4,nchar(bufferdf$species)-2)
bufferdf<-ddply(bufferdf, .(Buffer), summarize, mean.D=mean(D), D.sd=sd(D), mean.I=mean(I), I.sd=sd(I), mean.stability=mean(stability), stability.sd=sd(stability), mean.unfilling= mean(unfilling), unfilling.sd= sd(unfilling),mean.expansion =mean(expansion), expansion.sd= sd(expansion),mean.niche.centroid.shift=mean(niche.centroid.shift),niche.centroid.shift.sd=sd(niche.centroid.shift),mean.p.D=mean(p.D), p.D.sd=sd(p.D), mean.p.I=mean(p.I), p.I.sd=sd(p.I), mean.p.stability=mean(p.stability), p.stability.sd=sd(p.stability), mean.p.unfilling= mean(p.unfilling), p.unfilling.sd= sd(p.unfilling),mean.p.expansion =mean(p.expansion), p.expansion.sd= sd(p.expansion))



##Add historic range centroids to analyses dataframe##
HistCentdf<-melt(HistRangeCent)
HistCentdf<-dcast(HistCentdf, L1 ~ Var2)
names(HistCentdf)<-c("species","Historic.x","Historic.y")
#df<-merge(df,HistCentdf,by='species') -With buffeer loops, species names not matching

##Add current range centroids to analyses dataframe##
CurrentCentdf<-melt(CurrentRangeCent)
CurrentCentdf<-dcast(CurrentCentdf, L1 ~ Var2)
names(CurrentCentdf)<-c("species","Current.x",'Current.y')
df1<-merge(HistCentdf,CurrentCentdf,by='species')

##Add Range centroid shift to analyses dataframe##
RangeShiftdf<-melt(RangeCentShift)
RangeShiftdf<-dcast(RangeShiftdf, L1 ~ Var2)
names(RangeShiftdf)[names(RangeShiftdf)=="L1"]<-"species"
df1<-merge(df1,RangeShiftdf,by='species')

##Add range size, stability, expansion and contraction to dataframe##
HistRangedf<-melt(HistRangeSize)
names(HistRangedf)<-c('Historic_range/ncells',"species")
df1<-merge(df1,HistRangedf,by='species') 

CurrentRangedf<-melt(CurrentRangeSize)
names(CurrentRangedf)<-c('Current_range/ncells',"species")
df1<-merge(df1,CurrentRangedf,by='species')  

RangeExpansiondf<-melt(RangeExpansion)
names(RangeExpansiondf)<-c('Range_Expansion',"species")
df1<-merge(df1,RangeExpansiondf,by='species')  

RangeContractiondf<-melt(RangeContraction)
names(RangeContractiondf)<-c('Range_Contraction',"species")
df1<-merge(df1,RangeContractiondf,by='species')  

RangeStabilitydf<-melt(RangeStability)
names(RangeStabilitydf)<-c('Range_Stability',"species")
df1<-merge(df1,RangeStabilitydf,by='species')  

##Calculate gross range change to add
df1$RangeChange_ncells<-df1$`Current_range/ncells`-df1$`Historic_range/ncells`
summary(df1$RangeChange_ncells)

#Extract and merge 500km buffer only for example plotting#
df<-df[grep("500km$", df$species), ]
df$species<-substr(df$species,1,nchar(df$species)-13)
df<-merge(df,df1,by='species',all=T)

##Write summary data to csv for reference
write.csv(df,".Data//Europe_ebba_species/Niche_Analyses_summaries/All_Species_Analysis_data.csv")

df$RangeShiftMagnitude<-df$RangeShiftMagnitude/1000#Convert metres to km for range shift

####At this point need to define groups of species to analyse, whether taxonomic or functional, then these can be summarized to compare niche and range shifts across these groups (though analyses above may also need refining)####
df$species

df$Groups<-round(seq(from=1,to=5,length.out=nrow(df)),digits=0)##Add arbitrary group ID's for playing around with, as binomials are alphabetised, generated ID's where adjascent rows more likely to be in the same group, so groups are slighly less arbitrary and likley to contain species in the same genus

GrpSmry<-ddply(df, .(Groups), summarize, mean.D=mean(D), D.sd=sd(D), mean.I=mean(I), I.sd=sd(I), mean.stability=mean(stability), stability.sd=sd(stability), mean.unfilling= mean(unfilling), unfilling.sd= sd(unfilling),mean.expansion =mean(expansion), expansion.sd= sd(expansion),mean.niche.centroid.shift=mean(niche.centroid.shift),niche.centroid.shift.sd=sd(niche.centroid.shift),east.shift=mean(EastShift),east.shift.sd= sd(EastShift), north.shift= mean(NorthShift),north.shift.sd= sd(NorthShift), range.shift.mag=mean(RangeShiftMagnitude), range.shift.mag.sd=sd(RangeShiftMagnitude))