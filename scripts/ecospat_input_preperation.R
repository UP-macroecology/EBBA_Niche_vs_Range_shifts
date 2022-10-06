#Prepares raster climate layers to be uses as inputs for ecospat_analysis_loop

##File pathways made relative to ./Data need to be checked.

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
library(dplyr)
library(sf)
library(magrittr)
library(reshape2)

##PREPARING EBBA DATA FOR ECOSPAT ANALYSES ##



#####readRDS bioclim data####
Bioclim80s <- stack(list.files(file.path("./Data/Europe_Bioclimate_1981-90_2009-2018/1981-1990"), pattern = ".tif$", full.names = TRUE))

Bioclim10s <- stack(list.files(file.path("./Data/Europe_Bioclimate_1981-90_2009-2018/2009-2018"), pattern = ".tif$", full.names = TRUE))

##Suplement bioclims with LUH2 data from middle years - Investigated but discared at time of departure. Long-term averages should be used instead of midpoint years if LUH2 use is picked up again.
#Maybe remove secma and secmb? as these are not in gridcell fraction units as other Luh2 layers.
LandCover85 <- stack(list.files(file.path("./Data/LUH2_DAT"), pattern = "1985_luh2.50km.tif", full.names = TRUE))

#Bioclim80s <- stack (Bioclim80s,LandCover85)

LandCover15 <- stack(list.files(file.path("./Data/LUH2_DAT"), pattern = "2015_luh2.50km.tif", full.names = TRUE))

#Bioclim10s <- stack (Bioclim10s,LandCover15)



#make layer names consistent
bioclim<-c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", "bio16", "bio17", "bio18", "bio19", "bio2", "bio3","bio4", "bio5","bio6", "bio7","bio8", "bio9")
#,"c3ann","c3nfx","c3per","c4ann","c4per","pastr","primf","primn","range","secdf","secdn","secma","secmb","urban")
names(Bioclim80s)<-bioclim
names(Bioclim10s)<-bioclim

#Read in mask for eastern European countries to drop from analyses
MaskCountries<-readOGR("./Data","Non-comparable_region_mask")


#and the inverse for cropping sp dat
KeepCountries<-readOGR("./Data","Comparable_region_mask")

#Use Masks to remove data from Eastern Europe and the Canaries
Bioclim80s<-terra::mask(Bioclim80s,MaskCountries,inverse=T)
Bioclim80s<-terra::mask(Bioclim80s,KeepCountries,inverse=F,touches=T)##Removes rough edges around coastlines left over from previous mask, also all other raster cells that had a centroid not over land ... maybe benfit, but also maybe problem?

Bioclim10s<-terra::mask(Bioclim10s,MaskCountries,inverse=T)
Bioclim10s<-terra::mask(Bioclim10s,KeepCountries,inverse=F,touches=T)##Removes rough edges around coastlines left over from previous mask, also all other raster cells that had a centroid not over land ... maybe benfit, but also maybe problem?

#Mask EBBA Grid to get total retained cells#
RetainedEBBA_Grid<-terra::intersect(EBBA2_Grid,KeepCountries)
EBBAcells<-nrow(RetainedEBBA_Grid)

#write raster bricks to file to be read in by analysis loop
writeRaster(Bioclim80s,"./Data/climate_inputs/Bioclim80s.tif",overwrite=T)
writeRaster(Bioclim10s,"./Data/climate_inputs/Bioclim10s.tif",overwrite=T)
