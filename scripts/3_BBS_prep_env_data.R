# Bioclimatic variables used for niche analysis:

# notes:
# - run from cluster
# - before running the first part of script run "module load R/4.1.0-foss-2021a" on ecoc9 to ensure that R version 4.1.0 (2021-05-18) and GDAL 3.3.0 are used
# - but: calculation of bioclimatic variables only works when not running "module load R/4.1.0-foss-2021a" on ecoc9 beforehand -> did it in two steps

library(sf)
library(gdalUtilities)
library(terra)
library(doParallel)
library(dplyr)
library(raster)
library(stringr)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")

datashare_Chelsa <- file.path("/mnt","ibb_share","zurell","envidat","biophysical","CHELSA_V2","global") 
#datashare_Chelsa <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "envidat", "biophysical", "CHELSA_V2", "global") 

# folder to store reprojected Chelsa data:
chelsa_Albers_proj <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Chelsa_Albers_proj")
if(!dir.exists(chelsa_Albers_proj)){dir.create(chelsa_Albers_proj, recursive = TRUE)}

# register cores for parallel computation:
registerDoParallel(cores = 5)
getDoParWorkers() # check registered number of cores


# --------------------------------------------- #
#   Project Chelsa data of conterminous US:  ####
# --------------------------------------------- #

# historic and recent time periods:

chelsa_tifs <- list.files(datashare_Chelsa, full.names = FALSE,
                          pattern = paste0("(", paste(c(1980:1983, 1987:1990, 2015:2018), collapse = "|"), ")_V.2.1.tif"))

# names for reprojected Chelsa files:
names <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_Albers_1km.tif")

# CRS: Albers Equal Area projection:
albers_projection <- "ESRI:102003" # coordinate reference system


## create a mask: ----

# conterminous US, output of  2_3_BBS_species_filtering_5_climatic_niche_analysis.R:
# project with Albers projection:
read_sf(file.path(data_dir, "conterminousUS.shp")) %>%
  st_transform(crs = albers_projection) %>%
  st_write(obj = ., file.path(data_dir, "conterminousUS_Alb_proj.shp"), append = FALSE)

# rasterize shapefile:
gdalUtilities::gdal_rasterize(src_datasource = file.path(data_dir, "conterminousUS_Alb_proj.shp"),
                              dst_filename = file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "conterminousUS_1km.tif"),
                              burn = 1,
                              at = TRUE, # ALL_TOUCHED rasterization: all pixels touched by polygons get value 1
                              tr = c(1000, 1000), # 1 km target resolution (same unit as src_datasource)
                              a_nodata = -99999)


## project Chelsa layers: ----

foreach(s = 1:length(chelsa_tifs),
        .packages = c("gdalUtilities", "terra") ,
        .verbose = TRUE,
        .inorder = FALSE,
        .errorhandling = "remove") %dopar% {

  contUS_mask <- rast(file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "conterminousUS_1km.tif")) # read inside foreach since SpatRasters cannot be imported
  mask_ext <- ext(contUS_mask)

  # reproject Chelsa data:
  gdalUtilities::gdalwarp(srcfile = file.path(datashare_Chelsa,
                                              chelsa_tifs[s]),
                          dstfile = file.path(chelsa_Albers_proj, names[s]),
                          overwrite = TRUE,
                          tr = c(1000, 1000), # target resolution
                          r = "bilinear", # resampling method
                          t_srs = albers_projection,
                          te = c(mask_ext[1], mask_ext[3], mask_ext[2], mask_ext[4]))

  # mask reprojected Chelsa data:
  terra::rast(file.path(chelsa_Albers_proj, names[s])) %>%
    terra::mask(mask = contUS_mask) %>%
    terra::writeRaster(file.path(chelsa_Albers_proj, names[s]), overwrite = TRUE)

  names[s] # use file paths to avoid getting an error when trying to return SpatRasters
        }


# ----------------------------------------- #
#   Calculate bioclimatic variables:     ####
# ----------------------------------------- #

# works when not running "module load R/4.1.0-foss-2021a" on ecoc9 beforehand

time_periods <- list(1980:1983, 1987:1990, 2015:2018)

vars <- c("pr", "tas", "tasmin", "tasmax")
months <- str_pad(1:12, width = 2, pad = "0")

# load mask using raster package, this is required by dismo::biovars())
contUS_mask <- raster(file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "conterminousUS_1km.tif"))
#contUS_mask <- raster(file.path(data_dir, "conterminousUS_1km.tif"))

# loop over time periods:
for(t in 1:length(time_periods)){

  print(paste("time period", t))

  years <- time_periods[[t]]

  output_folder <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", paste0("Bioclim_", min(years), "_", max(years)))
  #output_folder <- file.path(data_dir, paste0("Bioclim_", min(years), "_", max(years)))
  if(!dir.exists(output_folder)){dir.create(output_folder, recursive = TRUE)}

  # masked Chelsa files corresponding to respective time period:
  chelsa_masked_files <- list.files(chelsa_Albers_proj,
                                    pattern = paste0("(", paste(years, collapse = "|"), ")",".*","tif$"),
                                    full.names = TRUE)

  # create a spatial brick template using the mask and
  # create bricks to store the month-wise means across the selected years for each bioclim variable
  out <- brick(contUS_mask, values = FALSE)
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
                                   paste0("CHELSA_", names(biovars), "_", min(years), "_", max(years), "_", "1km.tif")),
              overwrite = TRUE)
}