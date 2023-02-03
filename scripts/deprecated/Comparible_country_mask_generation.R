# xx not necessary any more!

##Generate templates to mask Out Areas not suitable for comparison between EBBA1 and 2 according to EBBA2 methods chapter (https://ico.ams3.digitaloceanspaces.com/ebba2/final_website/docs/EBBA2%20book%20chapter%20Methods.pdf).

#Regions removed from EBBA grid: Cyprus, Asian Turkey, Canary Islands, Russia, Kahzakstan, Georgia, Armenia and Azerbaijan.

#File pathways made relative to ./Data need to be checked.

library(sp)

#read in global boundary shapefile
Worldcountries <- readOGR("./Data/Countries geojson/countries.geojson")

#transform CRS
Worldcountries <- spTransform(Worldcountries,CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0
+y_0=0 +ellps=intl +units=m +no_defs"))

#Areas to exclude
MaskCountries <- subset(Worldcountries,ADMIN %in% c("Northern Cyprus","Cyprus","Akrotiri Sovereign Base Area","Turkey","Russia", "Kazakhstan", "Georgia", "Armenia", "Azerbaijan"))

#Areas to keep (or not exclude)
KeepCountries <- subset(Worldcountries,!(ADMIN %in% c("Northern Cyprus","Cyprus","Akrotiri Sovereign Base Area","Turkey","Russia", "Kazakhstan", "Georgia", "Armenia", "Azerbaijan",'Iran')))

#Canary Islands classed as 'Spain', so take some manipulation to isolate and remove
Canaries <- subset(Worldcountries,ADMIN %in% "Spain")
Canaries <- terra::crop(Canaries,extent(c(-18.1672257149999, -14, 27.6422386740001,32)))
MaskCountries <- spRbind(MaskCountries,Canaries)

plot(MaskCountries)

#Convert to proj cooordinates
MaskCountries <- spTransform(MaskCountries, CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"))


#write country mask to file
writeOGR(MaskCountries, "./Data", "Non-comparable_region_mask", driver= "ESRI Shapefile", overwrite= T)

writeOGR(KeepCountries, "./Data", "Comparable_region_mask", driver= "ESRI Shapefile", overwrite= T)


# adapted James code: ----
library(sf)
library(terra)
transfer <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "EBBA_niche_vs_range_shift", "Data")

# EBBA2 methods chapter:
# "We decided to restrict the change map to a well defined, continuous geographical area based on coverage in EBBA1. 
# Cyprus, the Asian part of Turkey and the Canary Islands were not included at all in EBBA1. 
# Moreover, a very large area including most of the European parts of Russia and Kazakhstan, Georgia, Armenia and Azerbaijan was only partially surveyd at the time.
# Consequently, this vast area in eastern Europe was excluded from the geographical area used in the change maps.
# More specifically, 50-km squares that lie totally or mostly (>70% of their area) within the above mentioned areas
# were not included."

#read in global boundary shapefile
Worldcountries <- read_sf(file.path(transfer, "Countries geojson", "countries.geojson"))

# prepared shapefiles with regions to keep in QGIS (removed isolated islands of Spain and Portugal (Canaries, Azores), Overseas France...)
keep_countries <- read_sf(file.path(transfer, "Comparable_countries.shp"))
plot(st_geometry(keep_countries))
# keep_countries_names <- c("Albania", "Andorra", "Austria", "Belgium", "Bulgaria", "Bosnia and Herzegovina",
#                           "Belarus", "Switzerland", "Czech Republic", "Germany", "Denmark", "Spain", "Estonia", 
#                           "Finland", "France", "United Kingdom", "Greece", "Croatia", "Hungary", "Ireland",
#                           "Iceland", "Italy", "Kosovo", "Liechtenstein", "Lithuania", "Luxembourg", "Latvia",
#                           "Moldova", "Macedonia", "Montenegro", "Netherlands", "Norway", "Poland", "Portugal",
#                           "Romania", "San Marino", "Republic of Serbia", "Slovakia", "Slovenia", "Sweden", "Ukraine", "Vatican")
# 
# # Areas to keep:
# KeepCountries <- subset(Worldcountries, ADMIN %in% keep_countries_names)
