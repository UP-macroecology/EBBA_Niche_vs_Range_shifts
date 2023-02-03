# Filter species to use in analysis, preparation for step 5: 

# (Step 5: Exclude species whose climatic niche is not well covered in Europe.
# To assess the climatic niche a species uses throughout its range we use Chelsa data.)

# in this script, Chelsa data of EBBA1 period (1981 - 1990) are projected to a global equal area projection (Interrupted Goode Homolosine)

# notes:
# - this script was run from the cluster, file paths may need to be updated
# - projecting the rasters yields slightly different results depending on the GDAL version used:
# 1) cluster "default": R version 4.0.4 (2021-02-15) - GEOS 3.9.0, GDAL 3.2.2, PROJ 7.2.1: "GDAL Error 1: Too many points failed to transform, unable to compute output bounds."
# 2) cluster R version 4.1.0 (access on ecoc9 by running "module load R/4.1.0-foss-2021a" and install all necessary packages) - GEOS 3.9.1, GDAL 3.3.0, PROJ 8.0.1: "GDAL Error 1: PROJ: igh: Invalid latitude" (may be due to PROJ version)
# 3) R version 4.1.2 on my laptop and GEOS 3.10.2, GDAL 3.4.1, PROJ 7.2.1: neither warnings nor errors
# -> I used option 2) since 2) and 3) yield same results, I did not get any messages for 3) and 1) duplicates Iceland and parts of Greenland

library(gdalUtilities)

datashare <- file.path("/mnt","ibb_share","zurell","envidat","biophysical","CHELSA_V2","global") 
chelsa_tifs <- list.files(datashare, full.names = FALSE, 
                          pattern = paste0("(", paste(1981:1990, collapse = "|"), ")_V.2.1.tif")) # 480

# folder to store reprojected CHELSA data: 
chelsa_birdlife_path <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "EBBA_Niche_vs_Range_shifts", "Data", "Chelsa_projected")
# create folder if it doesn't exist yet:
if(!dir.exists(chelsa_birdlife_path)){
  dir.create(chelsa_birdlife_path, recursive = TRUE)
}

# create file paths for the reprojected data:
names <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_50km.tif")
names <- file.path(chelsa_birdlife_path, names)

# global equal area projection used by SoilGrids (Homolosine projection applied to the WGS84 datum)
# from: https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection
homolosine <- 'PROJCS["Homolosine", 
                     GEOGCS["WGS 84", 
                            DATUM["WGS_1984", 
                                  SPHEROID["WGS 84",6378137,298.257223563, 
                                           AUTHORITY["EPSG","7030"]], 
                                  AUTHORITY["EPSG","6326"]], 
                            PRIMEM["Greenwich",0, 
                                   AUTHORITY["EPSG","8901"]], 
                            UNIT["degree",0.0174532925199433, 
                                 AUTHORITY["EPSG","9122"]], 
                            AUTHORITY["EPSG","4326"]], 
                     PROJECTION["Interrupted_Goode_Homolosine"], 
                     UNIT["Meter",1]]'

# reproject downloaded CHELSA layers:
options(warn = 1) # default warn = 0; 1 = warnings are printed as they occur

for(i in 1:length(chelsa_tifs)){
  
  print(paste(i, "of", length(chelsa_tifs)))
  
  # reproject Chelsa data:
  gdalUtilities::gdalwarp(srcfile = file.path(datashare, 
                                              chelsa_tifs[i]),
                          dstfile = names[i],
                          overwrite = TRUE,
                          tr = c(50000, 50000), # target resolution
                          r = "bilinear", # resampling method
                          t_srs = homolosine)
}

# could be speeded up by parallelising the code