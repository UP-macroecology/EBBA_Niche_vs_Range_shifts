##Generate templates to mask Out Areas not suitable for comparison between EBBA1 and 2 according to EBBA2 methods chapter (https://ico.ams3.digitaloceanspaces.com/ebba2/final_website/docs/EBBA2%20book%20chapter%20Methods.pdf): Cyprus,Asian Turkey, Canary Islands, Russia, Kahzakstan, Georgia, Armenia and Azerbaijan.

#File pathways made relative to ./Data need to be checked.

library(sp)

#Limmit global Download to national boundary layers for mask
Worldcountries<-readOGR("./Data/Countries geojson/countries.geojson")
#transform World countries to use as mask and tidy up coastlines
Worldcountries<-spTransform(Worldcountries,CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0
+y_0=0 +ellps=intl +units=m +no_defs"))
#Areas to exclude
MaskCountries<-subset(Worldcountries,ADMIN %in% c("Northern Cyprus","Cyprus","Akrotiri Sovereign Base Area","Turkey","Russia", "Kazakhstan", "Georgia", "Armenia", "Azerbaijan"))
#Areas to keep (or not exclude)
KeepCountries<-subset(Worldcountries,!(ADMIN %in% c("Northern Cyprus","Cyprus","Akrotiri Sovereign Base Area","Turkey","Russia", "Kazakhstan", "Georgia", "Armenia", "Azerbaijan",'Iran')))

Canaries<-subset(Worldcountries,ADMIN %in% "Spain")

plot(Canaries)

ext(Canaries)

Canaries<-terra::crop(Canaries,extent(c(-18.1672257149999, -14, 27.6422386740001,32)))

plot(Canaries)

MaskCountries<-spRbind(MaskCountries,Canaries)

plot(MaskCountries)

#Convert to proj cooordinates
MaskCountries<-spTransform(MaskCountries,CRS("+proj=aea +lat_0=30 +lon_0=10 +lat_1=43 +lat_2=62 +x_0=0
+y_0=0 +ellps=intl +units=m +no_defs"))


#write country mask to file
writeOGR(MaskCountries,"./Data","Non-comparable_region_mask",driver="ESRI Shapefile",overwrite=T)

writeOGR(KeepCountries,"./Data","Comparable_region_mask",driver="ESRI Shapefile",overwrite=T)
