####Generating species list and preliminary species inputs for requesting aligned data from raw EBBA1 and EBBA2, file pathways made relative to .Data/ need to be checked####
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

##Read in EBBA data##- UNLESS EDITTING CONVERSIONS, SKIP!
EBBA1<-read.table("./Data/input_data/EBCC/EBBA1/EBBA1_EPSG3035-Coords.csv",sep=",",header=T)
EBBA1<-subset(EBBA1,EBBA1$issue!="COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID")#Remove absence points from EBBA1

EBBA2<-read.table("./Data/input_data/EBCC/EBBA2/ebba2_data_occurrence_50km.csv",sep=";", header=TRUE)

#Remove missing coordinate#
EBBA1<-EBBA1[!is.na(EBBA1$ETRS89x),]

#Identify coordinates in EBBA1 data#
coordinates(EBBA1)<-~ETRS89x+ETRS89y

#Adding spatial data to EBBA2 from shapefile layer##
EBBA2_Grid<-readOGR(dsn="./Data/input_data/EBCC/EBBA2/ebba2_grid50x50_v1/ebba2_grid50x50_v1.shp",layer='ebba2_grid50x50_v1')
EBBA_CRS<-proj4string(EBBA2_Grid)
EBBA2_Grid<-gCentroid(EBBA2_Grid, byid = T,id=EBBA2_Grid$cell50x50)
EBBA2_Grid<-as.data.frame(EBBA2_Grid@coords)
EBBA2_Grid$cell50x50<-row.names(EBBA2_Grid)
head(EBBA2_Grid)

EBBA2<-merge(EBBA2_Grid,EBBA2,by='cell50x50',All=T)
coordinates(EBBA2)<-~x+y

# projection (equal area projection for Europe):
albers_projection <- "ESRI:102013" # attention: "ESRI" not "EPSG"
CRS(albers_projection)#to copy proj4string, unsure of direct extraction

## Give CRS to EBBAs ##
proj4string(EBBA1)<-EBBA_CRS
proj4string(EBBA2)<-EBBA_CRS

##Align CRS across datasets to match the equal area CHELSA climate reprojection##
EBBA1<-spTransform(EBBA1,CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"))

EBBA2<-spTransform(EBBA2,CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"))

EBBA2_Grid<-spTransform(EBBA2_Grid,CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"))

EBBA1<-intersect(EBBA1,KeepCountries)

EBBA2<-intersect(EBBA2,KeepCountries)


#write clipped data to file to speed up analyses in future
writeOGR(EBBA1,"./Data/input_data/EBCC/","Spatially_Clipped_EBBA1",driver="ESRI Shapefile",overwrite=T)
writeOGR(EBBA2,"./Data/input_data/EBCC/","Spatially_Clipped_EBBA2",driver="ESRI Shapefile",overwrite=T)

#Read in cropped EBBA data for analyses#
EBBA1<-readOGR("./Data/input_data/EBCC/","Spatially_Clipped_EBBA1")

EBBA2<-readOGR("./Data/input_data/EBCC/","Spatially_Clipped_EBBA2")

##

#Unify taxonomy between EBBA's based on EBBA methods chapter-Removing all species from Table 6 unless necesary taxon changes clear!
#change 3 species in EBBA1
EBBA1[which(EBBA1$species=="Acanthis hornemanni"),'species']<- "Acanthis flammea"

EBBA1[which(EBBA1$species=="Oceanodroma castro"),'species']<- "Hydrobates castro"

EBBA1[which(EBBA1$species=="Parus lugubris"),'species']<- "Poecile lugubris"

EBBA1[which(EBBA1$species=="Regulus ignicapillus"),'species']<- "Regulus ignicapilla"

EBBA1[which(EBBA1$species=="Serinus citrinella"),'species']<- "Carduelis citrinella"

#Remove non-comparible species from EBBA2 (and some from EBBA1)
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Branta hutchinsii"),'brdlf__']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Calonectris diomedea"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Calonectris borealis"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Calonectris borealis"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Hydrobates monteiroi"),'brdlf__']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Larus cachinnas"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Larus michahellis"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Larus michahellis"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Picus viridis"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Picus sharpei"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Picus viridis"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Lanius excubitor"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Lanius meridionalis"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Lanius excubitor"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Poecile hyrcanus"),'brdlf__']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Oenanthe xanthoprymna"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Oenanthe chrysopygia"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Oenanthe xanthoprymna"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Oenanthe xanthoprymna"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Oenanthe chrysopygia"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Oenanthe xanthoprymna"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Phylloscopus collybita"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Phylloscopus tristis"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Phylloscopus ibericus"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Phylloscopus collybita"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Phylloscopus bonelli"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Phylloscopus orientalis"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Phylloscopus bonelli"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Phylloscopus nitidus"),'brdlf__']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Sylvia hortensis"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Sylvia crassirostris"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Sylvia hortensis"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Sylvia sarda"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Sylvia balearica"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Sylvia sarda"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Sylvia cantillans"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Sylvia subalpina"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Sylvia cantillans"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Iduna pallida"),'brdlf__']
EBBA2<-EBBA2[which(EBBA2$brdlf__!="Iduna opaca"),'brdlf__']

EBBA1<-EBBA1[which(EBBA1$species!="Hippolais pallida"),'species']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Regulus  madeirensis"),'brdlf__']

EBBA2<-EBBA2[which(EBBA2$brdlf__!="Carduelis corsicana"),'brdlf__']

###Aligned taxons as much as possible above, matching species from those remaining below, keeping all species that have more than 20 records in both EBBAs###

#Get unique species names and counts, then remove all species with less than 5 records
EBBA2_Sp<-as.data.frame(table(EBBA2$brdlf__))#539 species
EBBA2_Sp<-as.character(subset(EBBA2_Sp,EBBA2_Sp$Freq>19)$Var1)

#Identify species dropped from EBBA1
EBBA1_Dropped<-subset(EBBA1,!(species %in% EBBA2_Sp))
EBBA1_Dropped<-unique(EBBA1_Dropped$species)#70 species dropped, can investigate more later

#subset EBBA1 to inlcude ONLY species with records in EBBA2
EBBA1<-subset(EBBA1,species %in% EBBA2_Sp)
EBBA1_keptSp<-unique(EBBA1$species)#379 species remaining

#Identify species to keep in EBBA2
EBBA1_Sp<-as.data.frame(table(EBBA1$species))#484 species
EBBA1_Sp<-as.character(subset(EBBA1_Sp,EBBA1_Sp$Freq>19)$Var1)

#Identify species dropped from EBBA2
EBBA2_Dropped<-subset(EBBA2,!(brdlf__ %in% EBBA1_Sp))
EBBA2_Dropped<-unique(EBBA2_Dropped$brdlf__)#160 species dropped, can investigate more later

#subset EBBA1 to inlcude ONLY species with records in EBBA2
EBBA2<-subset(EBBA2,brdlf__ %in% EBBA1_Sp)
EBBA2_keptSp<-unique(EBBA2$brdlf__)#379 species remaining


#write clipped data to file to speed up analyses in future
#writeOGR(EBBA1,"./Data/input_data/EBCC/","Taxon_matched_EBBA1",driver="ESRI Shapefile",overwrite=T)
#writeOGR(EBBA2,"./Data/input_data/EBCC/","Taxon_matched_EBBA2",driver="ESRI Shapefile",overwrite=T)

#Read in cropped EBBA data for analyses#
EBBA1<-readOGR("./Data/input_data/EBCC/","Taxon_matched_EBBA1")

EBBA2<-readOGR("./Data/input_data/EBCC/","Taxon_matched_EBBA2")

##Exporting species distibutions for examination and identification of remaining issues##
for(i in unique(EBBA1$species)){
  
  #Isolate single species#
  Sp_EBBA2<-subset(EBBA2,EBBA2$brdlf__==i)
  Sp_EBBA1<-subset(EBBA1,EBBA1$species==i)
  
  #Plot species distributions from EBBA1 and 2, exporting so I can go through and identify species with unrealistic distributions and/or seabirds with coastal distributions that are not suited to these analyses#
  tiff(filename = paste("./Data/EBBA_distribution_maps/",i,"_EBBA_map.tiff",sep=""),width=1000,height=1000,res=125)
  plot(Sp_EBBA1,col=alpha('red',1),pch=1,cex=0.75,main=paste(i,"red=EBBA1", "blue=EBBA2",sep=" "),xlim=c(xmin(EBBA1), xmax(EBBA1)), ylim=c(ymin(EBBA1), ymax(EBBA1)))
  plot(Sp_EBBA2,col=alpha('blue',1),pch=16,cex=0.5,add=T)
  plot(KeepCountries,add=T)
  dev.off()
}

#Removing seabird/coastal breeding species#- By eyeballing maps, will need a more robust and citable method eventually
SeaBirds<-read.csv2("C:/Users/James Hunter/OneDrive - Universität Potsdam/Desktop/Hunter_NichevsRangeShift_2022/input_data/BirdFuncDat.txt",header=T,sep="\t")#BirdFuncDat comes from the EltonTraits database
SeaBirds<-subset(SeaBirds,PelagicSpecialist==1)
SeaBirds<-SeaBirds$Scientific

EBBA1<-EBBA1[which(!EBBA1$species %in% SeaBirds),'species']
EBBA2<-EBBA2[which(!EBBA2$brdlf__ %in% SeaBirds),'brdlf__']

##Filtering out the most common birds##
EBBA2_common<-as.data.frame(table(EBBA2$brdlf__))
EBBA2_common<-as.character(subset(EBBA2_common,EBBA2_common$Freq>(0.9*EBBAcells))$Var1)#common if occur in >90% EBBA cells, removes 20 species

EBBA2<- as.data.frame(EBBA2) %>% group_by(brdlf__) %>% filter(n()<(0.9*EBBAcells)) %>% ungroup()#Use EBBA2 to filter as less confident presences in EBBA1 are acurate

Sp_List<-unique(EBBA2$brdlf__) #Now 326 species
Sp_List <- sub(" ", "_", Sp_List)#replace space with underscore to match IUCN file names

##Filter Out species with <5% range in EBBA comparible cells. Range from IUCN Red List gloabl polgons##
#List IUCN species polygons# 
#!!Make sure polygons are unzipped in folder before attempting!!#
target_path <- file.path(".Data//IUCN_range_polys/All_shapefiles")

IUCN_species <- list.files(target_path, full.names = F, pattern = ".shp")
IUCN_species <- substr(IUCN_species,1,nchar(IUCN_species)-13)

#Match IUCN names to EBBA species and identify unmatched EBBA species to harmonize#
Matched_sp<-IUCN_species[which(IUCN_species %in% Sp_List)]

unmatched_sp<-Sp_List[which(!Sp_List %in% IUCN_species)]##Examine these names and find suitable replacement for harmonization

#Harmonisation notes (changed IUCN file names to match those used in EBBA)
# - Carduelis_chloris -> Chloris_chloris
# - Carduelis_flammea -> Acanthis_flammea
# - Carduelis_spinus -> Spinus_spinus
# - Carduelis_cannabina -> Linaria_cannabina
# - Milaria_calandra -> Emberiza_calandra
# - Parus_ater -> Periparus_ater
# - Hirundo_rupestris -> Ptyonoprogne_rupestris
# - Hirundo_daurica -> Cecropis_daurica
# - Parus_caeruleus -> Cyanistes_caeruleus
# - Parus_cristatus -> Lophophanes_cristatus
# - parus_palustris -> Poecile_palustris
# - Carduelis_flavirostris -> Linaria_flavirostris
# - Parus_montanus -> Poecile_montanus
# - Parus_lugubris -> Poecile_lugubris
# - Parus_cinctus -> Poecile_cinctus
# - Sturnus_roseus -> Pastor_roseus
# - Hippolais_caligata -> Iduna_caligata
#!! Passer_italiae not mapped by IUCN!! But is European, so should be retained in species list#

#List IUCN files with full names#
IUCN_polys <- list.files(target_path, full.names = T, pattern = ".shp")
#Restrict IUCN polygons to just species included in the species list
IUCN_polys<-IUCN_polys[grep(paste(Sp_List, collapse="|"),IUCN_polys)]

##Read in EBBA polygons and restrict to retained region
EBBA2_poly<-readOGR("./Data/input_data/EBCC/EBBA2/ebba2_grid50x50_v1/ebba2_grid50x50_v1.shp")%>%
  spTransform(CRS('+proj=laea +lat_0=10 +lon_0=-81 +ellps=WGS84 +units=m +no_defs'))%>%
  gBuffer(byid=TRUE, width=0)

KeepCountries<-readOGR(".Data/","Comparable_region_mask")%>%
  spTransform(CRS('+proj=laea +lat_0=10 +lon_0=-81 +ellps=WGS84 +units=m +no_defs'))%>%
  gBuffer(byid=TRUE, width=0)#Put into global CRS and buffer to avoid gIntersect error

EBBA2_poly<-gIntersection(EBBA2_poly,KeepCountries)%>%
  gBuffer(byid=TRUE, width=0)



EBBA_range_cover<-list()

for(i in IUCN_polys[1:306]){
  glob_range<-readOGR(i) %>% 
    spTransform(CRS('+proj=laea +lat_0=10 +lon_0=-81 +ellps=WGS84 +units=m +no_defs'))%>%
    gBuffer(byid=TRUE, width=0)
  EBBA_range<-gIntersection(glob_range,EBBA2_poly)
  if (is.null(EBBA_range)){
    EBBA_range_cover[i]<-0
  } else {
    EBBA_range_cover[i]<-round(area(EBBA_range)/sum(area(glob_range))*100,2)
  }}

#Coerce list into table with species names, then subset by %
RangeCoverdf<-melt(EBBA_range_cover)
RangeCoverdf<-RangeCoverdf[,c(2,1)]
names(RangeCoverdf)<-c("species","range_cover_percent")
RangeCoverdf$species<-substr(RangeCoverdf$species,77,nchar(RangeCoverdf$species)-13)

#Restrict to finalised species list
write.csv(RangeCoverdf,"./Data/IUCN_poly_covered.csv")

hist(RangeCoverdf$range_cover_percent, main="Range within EBBA Region", xlab="Percentage range covered")

Sp_list<-subset(RangeCoverdf,range_cover_percent>5)$species #Likely threshold not defensible, method needs reconsidering

write.csv(Sp_List,"./Data/Change_Data_request_species_list.csv")

writeOGR(EBBA1,"./Data/input_data/EBCC/","Filtered_EBBA1",driver="ESRI Shapefile",overwrite=T)
writeOGR(EBBA2,"./Data/input_data/EBCC/","Filtered_EBBA2",driver="ESRI Shapefile",overwrite=T)