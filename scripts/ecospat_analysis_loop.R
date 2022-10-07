## PARELLEL PROCESSING LOOP FOR ECOSPAT ANALYSES OF BREEDING BIRDS ##
# File pathways made relative to "./Data" need to be checked #

##required packages##
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
library(foreach)
library(doParallel)
registerDoParallel(5) # Set number of cores for parallel processing

# Read in EBBA species data inputs #
EBBA1 <- readOGR("./Data/input_data/EBCC/", "Filtered_EBBA1")
EBBA2 <- readOGR("./Data/input_data/EBCC/", "Filtered_EBBA1")

#Read in  climate data#
Bioclim80s <- brick("./Data/climate_inputs/Bioclim80s.tif")
Bioclim10s <- brick("./Data/climate_inputs/Bioclim10s.tif")

#atribute bioclim names#
bioclim <- c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3","bio4", "bio5","bio6", "bio7","bio8", "bio9")
#,"c3ann","c3nfx","c3per","c4ann","c4per","pastr","primf","primn","range","secdf","secdn","secma","secmb","urban")

names(Bioclim80s) <- bioclim
names(Bioclim10s) <- bioclim


####LARGE ANALYSIS LOOP####
foreach (i= unique(EBBA1$species), .packages= c("sp", "rgdal", "terra", "raster", "ggmap", "maptools", "rgeos", "ecospat", "scales", "ade4", "tidyr", "plyr", "dplyr")) %dopar% {

#Isolate single species to simplify explorations#
Sp_EBBA2 <- subset(EBBA2, EBBA2$brdlf__==i)
Sp_EBBA1 <- subset(EBBA1, EBBA1$species==i)

###Calculate range centroid shifts (Eastward and Northward components and magnitude) for all species, in metres (CRS units)- Could elaborate with polygon intersects to look at area change?###
HistRangeCent <- coordinates(gCentroid(Sp_EBBA1))
CurrentRangeCent <- coordinates(gCentroid(Sp_EBBA2))

EastShift <- CurrentRangeCent[,1] - HistRangeCent[,1]
NorthShift <- CurrentRangeCent[,2] - HistRangeCent[,2]

RangeShiftMagnitude <- sqrt((EastShift^2) + (EastShift^2))
RangeCentShift <- cbind(EastShift, NorthShift, RangeShiftMagnitude)

saveRDS(HistRangeCent, file = paste("./Data/Europe_ebba_species/", i, "_HistRangeCent.Rds", sep=""))

saveRDS(CurrentRangeCent, file = paste("./Data/Europe_ebba_species/", i, "_CurrentRangeCent.Rds", sep=""))

saveRDS(RangeCentShift, file = paste("./Data/Europe_ebba_species/", i, "_RangeCentShift.Rds", sep=""))

##Calculate range overlap,contraction & expansion etc. based on overlapping occurance cells in EBBA1 and 2##
EBBA1_Cells <- as./Data.frame(terra::intersect(Sp_EBBA1,EBBA2_Grid))[c('species','cell50x50')]

EBBA1_Cells$species <- "pres"
names(EBBA1_Cells) <- c("EBBA1_Status","EBBA_Grid_Cell")
HistRangeSize <- nrow(EBBA1_Cells)

saveRDS(HistRangeSize, file = paste("./Data/Europe_ebba_species/", i, "_HistRangeSize.Rds", sep=""))

EBBA2_Cells <- as./Data.frame(terra::intersect(Sp_EBBA2,EBBA2_Grid))[c('brdlf__','cell50x50')]

EBBA2_Cells$brdlf__ <- "pres"
names(EBBA2_Cells) <- c("EBBA2_Status","EBBA_Grid_Cell")
CurrentRangeSize <- nrow(EBBA2_Cells)

saveRDS(CurrentRangeSize, file = paste("./Data/Europe_ebba_species/", i, "_CurrentRangeSize.Rds", sep=""))

EBBA_Cells <- merge(EBBA1_Cells,EBBA2_Cells,by="EBBA_Grid_Cell",all=T)

EBBA_Cells <- replace_na(EBBA_Cells,list(EBBA1_Status='abs',EBBA2_Status='abs'))

RangeExpansion <- nrow(EBBA_Cells[which(EBBA_Cells$EBBA1_Status == 'abs'
                                           & EBBA_Cells$EBBA2_Status =='pres'),]) /nrow(EBBA2_Cells) #presented as a proportion of the recent range, to compare with niche expansion from ecospat, which uses the recent data.


saveRDS(RangeExpansion, file = paste("./Data/Europe_ebba_species/", i, "_RangeExpansion.Rds", sep=""))

RangeContraction <- nrow(EBBA_Cells[which(EBBA_Cells$EBBA1_Status == 'pres'
                                             & EBBA_Cells$EBBA2_Status =='abs'),]) /nrow(EBBA1_Cells) #presented as a proportion of the historic range, to compare with niche unfilling from ecospat, which uses the historic data.

saveRDS(RangeContraction,file = paste("./Data/Europe_ebba_species/",i,"_RangeContraction.Rds",sep=""))

RangeStability <- nrow(EBBA_Cells[which(EBBA_Cells$EBBA1_Status == 'pres'
                                           & EBBA_Cells$EBBA2_Status =='pres'),]) /nrow(EBBA2_Cells) #presented as a proportion of the recent range, to compare with niche stability from ecospat, which uses the recent data.

saveRDS(RangeStability,file = paste("./Data/Europe_ebba_species/",i,"_RangeStability.Rds",sep=""))

write.csv(EBBA_Cells,paste("./Data/Europe_ebba_species/EBBA_pres_Abs_comparison_",i,".csv",sep=""))

#Trial different buffer extents and impacts,start 300-700km, in 100km steps (error running at 900km)

for(j in seq(from=300000,to=700000,by=100000)){

  presence_buffer80s <- raster::buffer(Sp_EBBA1,width=j,disolve=T,byid=T)
  bg_buf80s <- extent(presence_buffer80s)
  
  presence_buffer10s <- raster::buffer(Sp_EBBA2,width=j,disolve=T,byid=T)
  bg_buf10s <- extent(presence_buffer10s)
  
  spBioclim80s <- crop(x=Bioclim80s,y=bg_buf80s)
  spBioclim10s <- crop(x=Bioclim10s,y=bg_buf10s)
  
  spBioclim80s <- mask(x=spBioclim80s,mask=presence_buffer80s)
  spBioclim10s <- mask(x=spBioclim10s,mask=presence_buffer10s)

  #Extract Species-climate data and enviro background data as dataframesand remove NA's in extraction.
  EBBA1Enviros <- (terra::extract(spBioclim80s,Sp_EBBA1))
  EBBA1Enviros <- EBBA1Enviros[!rowSums(!is.finite(EBBA1Enviros)),]
  
  EBBA2Enviros <- (terra::extract(spBioclim10s,Sp_EBBA2))
  EBBA2Enviros <- EBBA2Enviros[!rowSums(!is.finite(EBBA2Enviros)),]
  
   
  Background80s <- as.matrix(spBioclim80s)#For background sample containing all cells in buffer.
  Background80s <- Background80s[!rowSums(!is.finite(Background80s)),]

  Background10s <- as.matrix(spBioclim10s)
  Background10s <- Background10s[!rowSums(!is.finite(Background10s)),]


  #Create Environmental PCA object
  pca.cal <- dudi.pca(rbind(Background10s,Background80s), scannf= FALSE, nf= 2)


  # Plot variable contribution with "ecospat.plot.contrib()"
  tiff(filename = paste("./Data/Europe_ebba_species/", i, " Buffer ", j/1000, "km environmental variable contribution.tiff", sep=""), width= 1000, height= 1000, res= 125)
  
ecospat.plot.contrib(contrib = pca.cal$co, eigen = pca.cal$eig) # The correlation circle indicates the contribution of original predictors to the PCA axes.

dev.off()

### Run Ecospat Analyses ###

# calculate the PCA scores
historic.scores <- suprow(pca.cal,EBBA1Enviros[,1:19])$li
recent.scores <- suprow(pca.cal, EBBA2Enviros[,1:19])$li
historic.env.scores <- suprow(pca.cal, Background80s[,1:19])$li
recent.env.scores <- suprow(pca.cal, Background10s[,1:19])$li
climate.allcell.scores <- pca.cal$li

# Calculate the occurence density grids
grid.clim.historic <- ecospat.grid.clim.dyn(climate.allcell.scores, historic.env.scores, historic.scores, R= 100)
grid.clim.recent <- ecospat.grid.clim.dyn(climate.allcell.scores, recent.env.scores, recent.scores, R= 100)

# Calculate niche overlap
# Compute Schnoer´s D, index of niche overlap
D.overlap <- ecospat.niche.overlap(grid.clim.historic, grid.clim.recent, cor = TRUE)$D

# Niche equivalency test according to Warren et al. 
# It is recommended to use at least 1000 replications for the similarity test.

# version for niche shift- Add option for Analogue (in addition to Non-Analogue) climates
shift.test.NA <- ecospat.niche.similarity.test(grid.clim.historic, grid.clim.recent, rep = 1500, overlap.alternative = "lower", expansion.alternative = "higher", stability.alternative = "lower", unfilling.alternative = "higher", intersection = NA, rand.type = 2)

saveRDS(shift.test.NA,file = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_similarity_shift_NA_test.Rds",sep=""))

shift.test.int <- ecospat.niche.similarity.test(grid.clim.historic, grid.clim.recent, rep = 1500, overlap.alternative = "lower", expansion.alternative = "higher", stability.alternative = "lower", unfilling.alternative = "higher", intersection = 0, rand.type = 2)

saveRDS(shift.test.int,file = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_similarity_shift_int_test.Rds",sep=""))


# version for niche conservatism
cons.test.NA <- ecospat.niche.similarity.test(grid.clim.historic, grid.clim.recent, rep = 1500, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower", intersection = NA, rand.type = 2)

saveRDS(cons.test.NA,file = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_similarity_conservatism_NA_test.Rds",sep=""))

cons.test.int <- ecospat.niche.similarity.test(grid.clim.historic, grid.clim.recent, rep = 1500, overlap.alternative = "higher", expansion.alternative = "lower", stability.alternative = "higher", unfilling.alternative = "lower", intersection = 0, rand.type = 2)

saveRDS(cons.test.int,file = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_similarity_conservatism_int_test.Rds",sep=""))

# plot similarity tests
tiff(filename = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_similarity_NA_test.tiff",sep=""))

ecospat.plot.overlap.test(shift.test.NA, "D", "Niche Similarity")

dev.off()

tiff(filename = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_similarity_int_test.tiff",sep=""))

ecospat.plot.overlap.test(shift.test.int, "D", "Niche Similarity")

dev.off()

tiff(filename = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_conservatism_NA_test.tiff",sep=""))

ecospat.plot.overlap.test(cons.test.NA, "D", "Niche Similarity")

dev.off()

tiff(filename = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_conservatism_int_test.tiff",sep=""))

ecospat.plot.overlap.test(cons.test.int, "D", "Niche Similarity")

dev.off()

# delimiting niche categories and quantifying niche dynamics in analogue climate 
niche.dyn <- ecospat.niche.dyn.index(grid.clim.historic, grid.clim.recent, intersection = NA)#run with 0=only overlapping climates, NA= all available climates

# visualizing niche categories, niche dynamics and climate analogy between ranges
tiff(filename = paste("./Data/Europe_ebba_species/Niche_Analyses_summaries/",i,"_Buffer_",j/1000,"km_Niche_and_enviro_space.tiff",sep=""),width=1000,height=1000,res=150)

ecospat.plot.niche.dyn(grid.clim.historic, grid.clim.recent, quant = 0.25, interest = 2, title = "Niche overlap", name.axis1 = "PC1", name.axis2 = "PC2",
                       col.exp="green", col.unf = "red", col.stab = "blue", colZ1 = "blue", colZ2="green")

ecospat.shift.centroids(historic.scores, recent.scores, climate.allcell.scores, climate.allcell.scores)

#legend(3, 7.5, legend = c("expansion", "stability", "unfilling"), col = c("green", "blue", "red"), pch = 15, cex = 1.5, pt.cex = 2)

dev.off()

##Extract niche-space centroids##
niche.dyn$historic.niche.centroid.x <- mean(historic.scores$Axis1)
niche.dyn$historic.niche.centroid.y <- mean(historic.scores$Axis2)

niche.dyn$recent.niche.centroid.x <- mean(recent.scores$Axis1)
niche.dyn$recent.niche.centroid.y <- mean(recent.scores$Axis2)

##Calculate magnitude of niche centroid shift with some Pythagoras##
niche.dyn$niche.centroid.shift <- sqrt(((niche.dyn$recent.niche.centroid.x - niche.dyn$historic.niche.centroid.x)^2) + ((niche.dyn$recent.niche.centroid.y - niche.dyn$historic.niche.centroid.y)^2))

niche.dyn %>% saveRDS(file = paste("./Data/Europe_ebba_species/", i, "_Buffer_", j/1000, "km_niche.dyn.Rds", sep=""))

}}

