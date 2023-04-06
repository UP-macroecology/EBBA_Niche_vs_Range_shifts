# Filter species to use in analysis, step 5: 

# Step 5: Exclude species whose climatic niche is not well covered in Europe:
# we compare the European distribution data (EBBA2) with the Birdlife range maps (of 2022),
# of the Birdlife range maps we use the area where a species is considered extant throughout the year (for resident species) 
# or extant during the breeding season (for migrants).
# To assess the climatic niche a species uses throughout its range we use Chelsa data.

# notes:
# - parts of the script (extracting shapefiles and masking Chelsa data) were run from the cluster, file paths may need to be updated
# - did not run everything on the cluster because there were errors that did not occur on my laptop -> package versions? need to check whether they still occur xx
# - ecospat dependencies ‘biomod2’, ‘randomForest’ depend on R >= 4.1.0
# -> before running on the cluster run "module load R/4.1.0-foss-2021a" on ecoc9 to ensure that R version 4.1.0 (2021-05-18) is used

library(dplyr)
library(sf)
library(doParallel)
library(gdalUtilities)
library(terra)
library(ecospat)
library(ade4)
library(raster)
library(dismo) # requires spatial data formats of raster package

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")
#data_dir <- file.path("Data")

# EBBA change data:
datashare_EBCC <- file.path("/mnt", "ibb_share", "zurell", "biodat", "distribution", "EBCC", "EBBA_change")
#datashare_EBCC <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "EBCC", "EBBA_change")

# Bird life range maps:
datashare_Birdlife <- file.path("/mnt", "ibb_share", "zurell", "biodat", "distribution", "Birdlife", "BOTW_2022")
#datashare_Birdlife <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "Birdlife", "BOTW_2022")

# Chelsa data:
datashare_Chelsa <- file.path("/mnt", "ibb_share", "zurell", "envidat", "biophysical", "CHELSA_V2", "global") 
#datashare_Chelsa <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "envidat", "biophysical", "CHELSA_V2", "global") 

# folder where projected CHELSA data are stored:
chelsa_birdlife_path <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA")

# species left after filtering step 4 (pelagic specialists, rare and very common species 
# and species occurring only in one atlas period excluded)
load(file = file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_prelim.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
species_filtered <- sort(sub(" ", "_", unique(EBBA1_prep$species)))

# output folder for niche dynamics plots Birdlife range vs. EBBA area:
EBBA_BL_plot_dir <- file.path(data_dir, "plots", "EBBA2_BL22_niche_dyn")
if(!dir.exists(EBBA_BL_plot_dir)){
  dir.create(EBBA_BL_plot_dir, recursive = TRUE)
}

# register cores for parallel computation:
registerDoParallel(cores = 2)
#getDoParWorkers() # check registered number of cores


# ------------------------------------------------------- #
#    1. Create mask to extract relevant Chelsa data:   ####
# ------------------------------------------------------- #

# create mask that covers Birdlife range polygons of all species left after filtering step 4 and EBBA grid
# (resolution = 50 km, as EBBAs)

## extract shapefiles from Birdlife 2022 geodatabase:: -------------------------

# only parts of range where species is considered extant (presence = 1) and which is used either throughout the whole year (seasonal = 1) or during the breeding season (seasonal = 2)

# account for taxonomic changes between EBBA (older) and BL range maps (newer):
# use taxonomy of EBBA:
spec_name_change_df <- data.frame("EBBA_name" = "Psittacula krameri", "BL_name" = "Alexandrinus krameri")
spec_name_change_df[2,] <- c("Sylvia communis", "Curruca communis")
spec_name_change_df[3,] <- c("Sylvia conspicillata", "Curruca conspicillata")
spec_name_change_df[4,] <- c("Sylvia curruca", "Curruca curruca")
spec_name_change_df[5,] <- c("Sylvia hortensis", "Curruca hortensis")
spec_name_change_df[6,] <- c("Sylvia melanocephala", "Curruca melanocephala")
spec_name_change_df[7,] <- c("Sylvia nisoria", "Curruca nisoria")
spec_name_change_df[8,] <- c("Sylvia ruppeli", "Curruca ruppeli")
spec_name_change_df[9,] <- c("Sylvia undata", "Curruca undata")

# output folder for Birdlife range shapefiles:
if(!dir.exists(file.path(data_dir, "Birdlife_ranges_EBBA", "Shapefiles_2022"))){
  dir.create(file.path(data_dir, "Birdlife_ranges_EBBA", "Shapefiles_2022"), recursive = TRUE)
}

# extract shapefiles:
foreach(s = 1:length(species_filtered),
        .packages = c("gdalUtilities"), 
        .verbose = TRUE,
        .errorhandling = "remove",
        .inorder = FALSE) %dopar% {
          
          spec <- sub("_", " ", species_filtered[s])
          
          # here adaptations when BL uses more modern taxonomy than EBBA:
          if(spec %in% spec_name_change_df$EBBA_name){
            spec <- spec_name_change_df$BL_name[which(spec_name_change_df$EBBA_name == spec)]
          }
          
          gdalUtilities::ogr2ogr(src_datasource_name = file.path(datashare_Birdlife, "BOTW.gdb"),
                                 layer = "All_Species",
                                 where = paste0("sci_name = '", spec, "' AND (SEASONAL = '1' OR SEASONAL = '2') AND PRESENCE = '1'"), # SEASONAL = 1: resident throughout the year, SEASONAL = 2: breeding season
                                 dst_datasource_name = file.path(data_dir, "Birdlife_ranges_EBBA", "Shapefiles_2022", paste0(species_filtered[s], ".shp")),
                                 overwrite = TRUE)
        }


## rasterize and project Birdlife shapefiles: ----------------------------------

# for each species: convert Birdlife range polygons to raster, 
# project raster and change resolution to 50 km:

# global equal area projection (Interrupted Goode Homolosine)
# used by SoilGrids (Homolosine projection applied to the WGS84 datum)
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


# list shapefiles:
BL_range_shps <- list.files(file.path(data_dir, "Birdlife_ranges_EBBA", "Shapefiles_2022"), full.names = TRUE, pattern = ".shp")

# output folder for Birdlife range rasters:
if(!dir.exists(file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022"))){
  dir.create(file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022"), recursive = TRUE)
}

foreach(i = 1:length(BL_range_shps),
        .packages = c("gdalUtilities"), 
        .verbose = TRUE,
        .errorhandling = "remove",
        .inorder = FALSE) %dopar% {
          

          # rasterize polygon:
          gdalUtilities::gdal_rasterize(src_datasource = file.path(BL_range_shps[i]),
                                        dst_filename = file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022", paste0(species_filtered[i], "_WGS84.tif")),
                                        burn = 1,
                                        tr = c(0.1, 0.1), # target resolution in degrees (same unit as src_datasource) (0.1° enough since final resolution will be 50 km)
                                        a_nodata = -99999) # value for cells with missing data
          
          # project and change resolution to 50 km (as EBBA):
          gdalUtilities::gdalwarp(srcfile = file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022", paste0(species_filtered[i], "_WGS84.tif")), 
                                  dstfile = file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022", paste0(species_filtered[i], "_50km.tif")),
                                  overwrite = TRUE,
                                  tr = c(50000, 50000), # target resolution in meters (same as unit of target srs)
                                  r = "max", # resampling method
                                  t_srs = homolosine # target spatial reference
          )
        }


## rasterize and project EBBA change grid: -------------------------------------

# rasterize shapefile:
gdalUtilities::gdal_rasterize(src_datasource = file.path(datashare_EBCC, "EBBA_change", "ebba2_grid50x50_change_v1.shp"),
                              dst_filename = file.path(data_dir, "EBBA_change_50km.tif"),
                              burn = 1,
                              tr = c(50000, 50000), # target resolution in degrees (same unit as src_datasource)
                              a_nodata = -99999) # value for cells with missing data

# project and change resolution to 50 km:
gdalUtilities::gdalwarp(srcfile = file.path(data_dir, "EBBA_change_50km.tif"), 
                        dstfile = file.path(data_dir, "EBBA_change_50km_homolosine.tif"), 
                        overwrite = TRUE,
                        tr = c(50000, 50000), # target resolution in meters (same as unit of target srs)
                        r = "max", # resampling method
                        t_srs = homolosine # target spatial reference
)


## mask: overlay Birdlife ranges of all species and EBBA grid: -----------------


# list rasters:
BL_range_tifs <- list.files(file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022"), full.names = TRUE, pattern = "_50km.tif$")

# Birdlife ranges of all species:
ranges <- lapply(BL_range_tifs, rast)

# add EBBA grid:
ranges <- append(ranges, rast(file.path(data_dir, "EBBA_change_50km_homolosine.tif")))

# overlay rasters, add 100 km buffer and resample to match projected Chelsa rasters:
# (buffer to make sure that even after resampling mask covers raster grid points of species ranges)
# (resampling aligns raster origins, necessary for masking)

ranges_combined <- do.call(terra::merge, ranges) %>% 
  buffer(width = 100000) %>% # Warning: "[merge] rasters did not align and were resampled" -> fine, aligns ranges of different species
  resample(y = rast(file.path(data_dir, "Chelsa", "CHELSA_tasmin_11_2014_V.2.1_50km.tif")),
           method = "near") # projected Chelsa raster as template

# save mask:
writeRaster(ranges_combined, 
            filename = file.path(data_dir, "Birdlife_ranges_mask_300323.tif"),
            overwrite = TRUE,
            NAflag = FALSE)


## exploration: compare EBBA1 range with range maps from Birdlife: -------------

library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)

world <- ne_countries(scale = "small", returnclass = "sf")

EBBA1_taxunif_sf <- read_sf(file.path(data_dir, "EBBA1_comparable_harmonized.shp")) # output of 1_prep_EBBA_data.R
EBBA_change_grid <- read_sf(file.path(datashare_EBCC, "EBBA_change", "ebba2_grid50x50_change_v1.shp")) 

# for which species should distributions be plotted:
spec = 3

# Birdlife range map:
BL_range_sf <- read_sf(BL_range_shps[spec])

# EBBA distribution:
EBBA1_range_sf <- EBBA1_taxunif_sf %>%
  st_drop_geometry %>% 
  filter(species == sub("_", " ", species_filtered[spec])) %>%
  right_join(EBBA_change_grid, by = c("cell50x50" = "cell50x50_")) %>%
  st_as_sf()

ggplot(world) +
  geom_sf(color = "gray50") +
  geom_sf(data = st_transform(BL_range_sf, crs = st_crs(world)), fill = "red", alpha = 0.5) +
  geom_sf(data = st_transform(EBBA1_range_sf, crs = st_crs(world)), fill = "blue", colour = NA, alpha = 0.5) +
  ggtitle(sub("_", " ", species_filtered[spec]))

# ------------------------------- #
#     2. Mask Chelsa data      ####
# ------------------------------- #

# names of original Chelsa files:
chelsa_tifs <- list.files(datashare_Chelsa, full.names = FALSE, pattern = paste0("(", paste(2009:2018, collapse = "|"), ")_V.2.1.tif")) # 480

# names of projected Chelsa files:
names_proj <- list.files(file.path(chelsa_birdlife_path, "Chelsa_projected"), full.names = TRUE)

# folder to store masked CHELSA data:
chelsa_masked_path <- file.path(data_dir, "Chelsa", "Chelsa_masked_global")
# create folder if it doesn't exist yet:
if(!dir.exists(chelsa_masked_path)){
  dir.create(chelsa_masked_path, recursive = TRUE)
}

# create file paths for the reprojected and masked data:
names_masked <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_50km_masked.tif")
names_masked <- file.path(chelsa_masked_path, names_masked)

# mask:
ranges_mask <- rast(file.path(data_dir, "Birdlife_ranges_mask_300323.tif"))

# mask reprojected rasters:
for(i in 1:length(chelsa_tifs)){
  
  print(i)
  
  names_proj[i] %>%
    rast %>%
    mask(mask = ranges_mask) %>%
    writeRaster(names_masked[i], overwrite = TRUE)
}


# ----------------------------------------------------------- #
#     3. Calculate bioclim variables from Chelsa data:     ####
# ----------------------------------------------------------- #

vars <- c("pr", "tas", "tasmin", "tasmax")
months <- stringr::str_pad(1:12, width = 2, pad = "0")
years <- 2009:2018

# files masked Chelsa data:
chelsa_masked_files <- list.files(chelsa_masked_path, pattern = "tif$", full.names = TRUE)

# load mask using raster package, this is required by dismo::biovars())
ranges_mask <- raster(file.path(data_dir, "Birdlife_ranges_mask_230323.tif"))

# create a spatial brick template using the mask and 
# create bricks to store the month-wise means across the selected years for each bioclim variable
out <- brick(ranges_mask, values = FALSE)
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

# calculate rasters of bioclim variables:
biovars <- dismo::biovars(prec = pr_mean,  
                          tmin = tasmin_mean, 
                          tmax = tasmax_mean)

biovars_rast <- rast(biovars) # convert to terra object

# save tifs:
# create folder if it doesn't exist yet:
bioclim_folder <- file.path(data_dir, "Bioclim_EBBA2")
if(!dir.exists(bioclim_folder)){
  dir.create(bioclim_folder, recursive = TRUE)
}
writeRaster(biovars_rast,
            filename = file.path(bioclim_folder, 
                                 paste0("CHELSA_", names(biovars_rast), "_", min(years), "_", max(years), "_", "50km.tif")), 
            overwrite = TRUE)


# --------------------------------------- #
#     4. Climatic niche analyses:      ####
# --------------------------------------- #

# calculate stability index: how much of the species global climatic niche is covered in Europe

# read comparable EBBA1 data (output of 1_prep_EBBA_data.R):
EBBA2_comp_sf <- st_read(file.path(data_dir, "EBBA2_comparable.shp")) # output of 1_prep_EBBA_data.R

# files of rasterized and projected Birdlife ranges:
BL_range_tifs

# loop over species:
stability_df <- foreach(i = 1:length(species_filtered),
                          .combine = rbind,
                          .packages = c("ecospat", "ade4", "dplyr", "sf", "terra"), 
                          .verbose = TRUE,
                          .errorhandling = "remove") %dopar% {
                            
                            # data frame to store results:
                            stability_spec_df <- data.frame("species" = sub("_", " ", species_filtered[i]),
                                                       "stability" = NA,
                                                       "niche_breadth_zcor" = NA,
                                                       "niche_breadth_Z" = NA)
                            
                            # read tifs with bioclim variables:
                            # need to be loaded within foreach since SpatRasters and SpatVectors are non-exportable objects
                            biovars_rast <- rast(file.path(bioclim_folder, 
                                                           paste0("CHELSA_", paste0("bio", 1:19), "_", min(years), "_", max(years), "_", "50km.tif")))
                            
                            # EBBA2 occurrence point of species i:
                            
                            EBBA2_spec_sf <- EBBA2_comp_sf %>%
                              st_transform(crs = homolosine) %>% # transform to Homolosine projection to match bioclim variables
                              filter(species == sub("_", " ", species_filtered[i])) %>% 
                              vect # convert to terra object
                            
                            # add bioclim variables to occurrences:
                            EBBA2_spec_vars_df <- cbind(values(EBBA2_spec_sf[, c("species", "cell50x50")]),
                                                        terra::extract(biovars_rast, EBBA2_spec_sf))
                            
                            # Birdlife range of species i:
                            
                            BL_file <- grep(pattern = species_filtered[i], BL_range_tifs, value = TRUE)
                            
                            BL_pts_spec <- BL_file %>% 
                              rast %>% # load raster
                              as.points # convert raster to points (cell centroids)
                            
                            # add bioclim variables to occurrences:
                            BL_spec_vars_df <- cbind(data.frame(species = sub("_", " ", species_filtered[i]),
                                                                cell50x50 = "BL"),
                                                     terra::extract(biovars_rast, BL_pts_spec)) %>% 
                              filter(complete.cases(.)) # for few occurrence points "on the edge of the globe" no Chelsa data are available; because raster grid points are used as occurrence points and not cell centroids, and resampling is done, few points end up outside of global Chelsa range)
                            
                            # assess climate niche by using the first 2 PCA axes:
                            # calibrating the PCA in the whole study area, including both ranges:
                            pca.env <- dudi.pca(rbind(EBBA2_spec_vars_df, BL_spec_vars_df)[, paste0("bio", 1:19)],
                                                scannf = FALSE,
                                                nf = 2) # number of axes
                            
                            # predict the scores on the PCA axes:
                            # EBBA2 used as z1 (corresponds to native distribution in tutorials)
                            # Birdlife range used as z2 (corresponds to invasive distribution in tutorials)
                            scores.globclim <- pca.env$li # PCA scores for EBBA + Birdlife distribution
                            scores.clim.EBBA <- suprow(pca.env, EBBA2_spec_vars_df[, paste0("bio", 1:19)])$li # PCA scores for EBBA distribution
                            scores.clim.BL <- suprow(pca.env, BL_spec_vars_df[, paste0("bio", 1:19)])$li # PCA scores for Birdlife distribution
                            
                            # calculate the Occurrence Densities Grid for EBBA and Birdlife distribution:
                            
                            # EBBA2 distribution area:
                            grid.clim.EBBA <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                                                    glob1 = scores.clim.EBBA, 
                                                                    sp = scores.clim.EBBA, # same for background and occurrences, should be fine for our purpose
                                                                    R = 100, # grid resolution
                                                                    th.sp = 0)
                            # Birdlife range area:
                            grid.clim.BL <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                                                  glob1 = scores.clim.BL, 
                                                                  sp = scores.clim.BL, # same for background and occurrences, should be fine for our purpose
                                                                  R = 100, 
                                                                  th.sp = 0)
                            
                            stability_spec_df$niche_breadth_zcor <- vegan::diversity(as.vector(grid.clim.BL$z.cor)) #xx extract for later trait analysis
                            stability_spec_df$niche_breadth_Z <- vegan::diversity(as.vector(grid.clim.BL$Z)) #xx extract for later trait analysis
                            
                            # assess niche dynamics between EBBA2 (as z1) and Birdlife range map (as z2) 
                            # => the stability index shows how much of the species global climatic niche 
                            # is covered in Europe:
                            
                            EBBA_BL_niche_dyn <- ecospat.niche.dyn.index(grid.clim.EBBA, 
                                                                         grid.clim.BL,
                                                                         intersection = NA) # analysis performed on EBBA + Birdlife extent
                            stability_spec_df$stability <- EBBA_BL_niche_dyn$dynamic.index.w['stability']
                            
                            # save plot of niche dynamics:
                            jpeg(filename = file.path(EBBA_BL_plot_dir, paste0(species_filtered[i], ".jpg")),
                                 quality = 100, height = 920, width = 920
                            )
                            
                            ecospat.plot.niche.dyn(grid.clim.EBBA, 
                                                   grid.clim.BL, 
                                                   quant = 0.1,
                                                   interest = 2, # 1 = EBBA density, 2 = BL density
                                                   title = paste(species_filtered[i], "stability:", round(stability_spec_df$stability,2)), 
                                                   name.axis1 = "PC1",
                                                   name.axis2 = "PC2")
                            dev.off()
                            # blue = stability
                            # red = expansion
                            # green = unfilling
                            
                            stability_spec_df
                          }

stability_df %>% 
  arrange(-stability) %>%  View

# save species stability values :
write.csv(stability_df, 
          file = file.path(data_dir, "species_stability_EBBA2_BL22_030423.csv"),
          row.names = FALSE)

# select for now 150 species with highest coverage of climatic niche in Europe:
stability_df_selection <- stability_df %>% 
  arrange(-stability) %>% 
  slice_head(n = 150)

# save species selection:
write.csv(stability_df_selection, 
          file = file.path(data_dir, "species_selection_EBBA2_BL22_030423.csv"),
          row.names = FALSE)

# plots regarding stability:
plot(sort(stability_df$stability, decreasing = TRUE),
     ylab = "stability", xlab = "number of species", las = 1)
abline(h = 0.7, col = "red")
abline(h = 0.5, col = "blue")