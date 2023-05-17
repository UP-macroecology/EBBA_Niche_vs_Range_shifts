# Bioclimatic variables used for niche analysis:

# notes:
# - run from cluster
# - before running the script run "module load R/4.1.0-foss-2021a" on ecoc9 to ensure that R version 4.1.0 (2021-05-18) and GDAL 3.3.0 are used

library(gdalUtilities)
library(terra)
library(doParallel)
library(dplyr)
library(raster)
library(stringr)

# ------------------------ #
#       Set-up:         ####
# ------------------------ #

# file paths:
datashare_EBCC <- file.path("/mnt","ibb_share","zurell", "biodat", "distribution", "EBCC")
datashare_Chelsa <- file.path("/mnt","ibb_share","zurell","envidat","biophysical","CHELSA_V2","global") 

# folder to store reprojected Chelsa data:
chelsa_EPSG3035 <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Chelsa_EPSG3035")
if(!dir.exists(chelsa_EPSG3035)){dir.create(chelsa_EPSG3035, recursive = TRUE)}

# register cores for parallel computation:
registerDoParallel(cores = 8)
getDoParWorkers() # check registered number of cores

# ----------------------------------------- #
#   Project Chelsa data of EBBA 2 area:  ####
# ----------------------------------------- #

# Lambert azimuthal equal-area projection (ETRS89-extended / LAEA Europe, EPSG:3035)
# historic (1981-1990) and recent time period (2009-2018)

chelsa_tifs <- list.files(datashare_Chelsa, full.names = FALSE, 
                          pattern = paste0("(", paste(c(1981:1990, 2009:2018), collapse = "|"), ")_V.2.1.tif")) # 960

# names for reprojected Chelsa files:
names <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_EPSG3035_50km.tif")

# CRS: Lambert azimuthal equal-area projection:
lambert_projection <- "EPSG:3035"


## create a mask: ----

# to create the mask, I use the EBBA 2 spatial grid and not only the cells that can be compared across EBBA 1 and EBBA 2
# to avoid holes resulting from non-comparable cells

# rasterize shapefile:
gdalUtilities::gdal_rasterize(src_datasource = file.path(datashare_EBCC, "EBBA2", "ebba2_grid50x50_v1", "ebba2_grid50x50_v1.shp"),
                              dst_filename = file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "EBBA2_area.tif"),
                              burn = 1,
                              at = TRUE, # ALL_TOUCHED rasterization: all pixels touched by polygons get value 1
                              tr = c(50000, 50000), # target resolution (same unit as src_datasource)
                              a_nodata = -99999)

EBBA_mask <- rast(file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "EBBA2_area.tif"))
mask_ext <- ext(EBBA_mask)


## project Chelsa layers: ----

foreach(s = 1:length(chelsa_tifs), 
        .packages = c("gdalUtilities") , 
        .verbose = TRUE) %dopar% {
  
  # reproject Chelsa data:
  gdalUtilities::gdalwarp(srcfile = file.path(datashare_Chelsa, chelsa_tifs[s]),
                          dstfile = file.path(chelsa_EPSG3035, names[s]),
                          overwrite = TRUE,
                          tr = c(50000, 50000), # target resolution
                          r = "bilinear", # resampling method
                          t_srs = lambert_projection,
                          te = c(mask_ext[1], mask_ext[3], mask_ext[2], mask_ext[4]))
  
  # mask reprojected Chelsa data:
  rast(file.path(chelsa_EPSG3035, names[s])) %>% 
    mask(mask = EBBA_mask) %>%
    writeRaster(file.path(chelsa_EPSG3035, names[s]), overwrite = TRUE)
} 


# ----------------------------------------- #
#   Calculate bioclimatic variables:     ####
# ----------------------------------------- #

# once for historic period (1981-1990), once for recent period (2009-2018)

vars <- c("pr", "tas", "tasmin", "tasmax")
months <- str_pad(1:12, width = 2, pad = "0")

# load mask using raster package, this is required by dismo::biovars())
EBBA_mask <- raster(file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "EBBA2_area.tif"))

# folder to store bioclim rasters for historic time period:
bioclim_1981_1990 <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Bioclim_1981_1990")
if(!dir.exists(bioclim_1981_1990)){dir.create(bioclim_1981_1990, recursive = TRUE)}

# folder to store bioclim rasters for recent time period:
bioclim_2009_2018 <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Bioclim_2009_2018")
if(!dir.exists(bioclim_2009_2018)){dir.create(bioclim_2009_2018, recursive = TRUE)} 


# loop over time periods:
for(t in 1:2){

  print(paste("time period", t))
  
  if(t == 1){
    years <- 1981:1990
    output_folder <- bioclim_1981_1990
  } else {
      years <- 2009:2018
      output_folder <- bioclim_2009_2018}

  # masked Chelsa files corresponding to respective time period:
  chelsa_masked_files <- list.files(chelsa_EPSG3035, 
                                    pattern = paste0("(", paste(years, collapse = "|"), ")",".*","tif$"), 
                                    full.names = TRUE)
  
  # create a spatial brick template using the mask and 
  # create bricks to store the month-wise means across the selected years for each bioclim variable
  out <- brick(EBBA_mask, values = FALSE)
  pr_mean <- tasmin_mean <- tasmax_mean <- out
  
  # loop over months:
  for(j in as.numeric(months)){
    
    print(paste("month", j))
    
    # calculate precipitation mean for the current month (across all years)
    pr_mean[[j]] <- chelsa_masked_files[which(grepl(paste0("pr_", months[j]),
                                                    chelsa_masked_files))] %>% 
      stack %>% 
      mean
    
    # calculate tasmin mean for the current month (across all years)
    tasmin_mean[[j]] <- chelsa_masked_files[which(grepl(paste0("tasmin_", months[j]), 
                                                        chelsa_masked_files))] %>% 
      stack %>% 
      mean
    
    # calculate tasmax mean for the current month (across all years)
    tasmax_mean[[j]] <- chelsa_masked_files[which(grepl(paste0("tasmax_", months[j]), chelsa_masked_files))] %>% 
      stack %>% 
      mean
  }
  
  # calculate bioclimatic variables:
  biovars <- dismo::biovars(prec = pr_mean,  
                            tmin = tasmin_mean, 
                            tmax = tasmax_mean)
  
  biovars_rast <- terra::rast(biovars) # convert to terra object
  
  # save tifs:
  terra::writeRaster(biovars_rast,
              filename = file.path(output_folder, 
                                   paste0("CHELSA_", names(biovars), "_", min(years), "_", max(years), "_", "50km.tif")), 
              overwrite = TRUE)
}