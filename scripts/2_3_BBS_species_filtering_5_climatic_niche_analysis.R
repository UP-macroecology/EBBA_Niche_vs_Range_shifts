# Filter BBS species to use in analysis, step 5: 
# Exclude species whose climatic niche is not well covered in the conterminous US:
# We compare the global Birdlife range and the part of the range covered by the conterminous US.
# To assess the climatic niche we use Chelsa data.

# notes:
# - parts of the script (extracting shapefiles and masking Chelsa data) were run from the cluster, file paths may need to be updated
# - did not run everything on the cluster because there were errors that did not occur on my laptop -> package versions?
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
#data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")
data_dir <- file.path("data", "BBS_analysis")

# Birdlife range maps:
#datashare_Birdlife <- file.path("/mnt", "ibb_share", "zurell", "biodat", "distribution", "Birdlife", "BOTW_2022")
datashare_Birdlife <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "Birdlife", "BOTW_2022")

# Chelsa data:
#datashare_Chelsa <- file.path("/mnt", "ibb_share", "zurell", "envidat", "biophysical", "CHELSA_V2", "global") 
datashare_Chelsa <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data","envidat","biophysical","CHELSA_V2","global") 

# folder with projected CHELSA data (output of 2_2_species_filtering_5_project_Chelsa.R):
chelsa_proj_path <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA")

# output folder for niche dynamics plots Birdlife range vs. conterminous US area:
#US_BL_plot_dir <- file.path(data_dir, "plots", "contUS_vs_BL_niche_dyn")
US_BL_plot_dir <- file.path("plots", "coverage_climatic_niche", "contUS_vs_BL_niche_dyn")
if(!dir.exists(US_BL_plot_dir)){dir.create(US_BL_plot_dir, recursive = TRUE)}

# load species left after filter steps 1-4, for both versions of historic time period:
load(file = file.path(data_dir, "BBS_prep_steps1-4_hist81-83.RData")) # output of 2_1_BBS_species_filtering_1-4.R
species_filtered_V1 <- sort(sub(" ", "_", unique(hist_prep_df$species))) # 212
load(file = file.path(data_dir, "BBS_prep_steps1-4_hist96-98.RData")) # output of 2_1_BBS_species_filtering_1-4.R
species_filtered_V2 <- sort(sub(" ", "_", unique(hist_prep_df$species))) # 312
species_filtered <- sort(unique(c(species_filtered_V1, species_filtered_V2))) # 313

# register cores for parallel computation:
registerDoParallel(cores = 2)

# ------------------------------------------------------- #
#    1. Create mask to extract relevant Chelsa data:   ####
# ------------------------------------------------------- #

# create mask that covers Birdlife range polygons of all species left after filtering step 4
# (for species filtering we use same resolution as for EBBA data: 50 km)

## extract shapefiles from Birdlife 2022 geodatabase: --------------------------

# only parts of range where species is considered extant (presence = 1) and which 
# is used either throughout the whole year (seasonal = 1) or during the breeding season (seasonal = 2)

# account for taxonomic changes between BBS and BL range maps:
# (suitable resources:
# https://www.iucnredlist.org/
# https://explorer.natureserve.org/Search)
spec_name_change_df <- data.frame("BBS_name" = "Antigone canadensis", "BL_name" = "Grus canadensis")
spec_name_change_df[2,] <- c("Aphelocoma woodhouseii", "Aphelocoma californica")
spec_name_change_df[3,] <- c("Centronyx bairdii", "Passerculus bairdii")
spec_name_change_df[4,] <- c("Centronyx henslowii", "Passerculus henslowii")
spec_name_change_df[5,] <- c("Chroicocephalus philadelphia", "Larus philadelphia")
spec_name_change_df[6,] <- c("Coccothraustes vespertinus", "Hesperiphona vespertina")
spec_name_change_df[7,] <- c("Colaptes auratus cafer", "Colaptes cafer")
spec_name_change_df[8,] <- c("Dryobates albolarvatus", "Leuconotopicus albolarvatus")
spec_name_change_df[9,] <- c("Himantopus mexicanus", "Himantopus himantopus")
spec_name_change_df[10,] <- c("Icterus bullockii", "Icterus bullockiorum")
spec_name_change_df[11,] <- c("Junco hyemalis caniceps", "Junco hyemalis")
spec_name_change_df[12,] <- c("Junco hyemalis hyemalis", "Junco hyemalis")
spec_name_change_df[13,] <- c("Junco hyemalis oreganus", "Junco hyemalis")
spec_name_change_df[14,] <- c("Leucophaeus atricilla", "Larus atricilla")
spec_name_change_df[15,] <- c("Leucophaeus pipixcan", "Larus pipixcan")
spec_name_change_df[16,] <- c("Phalacrocorax auritus", "Nannopterum auritus")
spec_name_change_df[17,] <- c("Phalaropus tricolor", "Steganopus tricolor")
spec_name_change_df[18,] <- c("Pica nuttalli", "Pica nutalli")
spec_name_change_df[19,] <- c("Regulus calendula", "Corthylio calendula")
spec_name_change_df[20,] <- c("Setophaga coronata audoboni", "Setophaga auduboni")
spec_name_change_df[21,] <- c("Setophaga coronata coronata", "Setophaga coronata")
spec_name_change_df[22,] <- c("Butorides virescens", "Butorides striata")
spec_name_change_df[23,] <- c("Dryobates villosus", "Leuconotopicus villosus")
spec_name_change_df[24,] <- c("Dryocopus pileatus", "Hylatomus pileatus")
spec_name_change_df[25,] <- c("Colaptes auratus auratus", "Colaptes auratus")

# output folder for Birdlife range shapefiles:
if(!dir.exists(file.path(data_dir, "Birdlife_ranges_BBS", "Shapefiles_2022"))){dir.create(file.path(data_dir, "Birdlife_ranges_BBS", "Shapefiles_2022"), recursive = TRUE)}

# shapefile extraction:
foreach(s = 1:length(species_filtered),
        .packages = c("gdalUtilities"), 
        .verbose = TRUE,
        .errorhandling = "remove",
        .inorder = FALSE) %dopar% {
          
          spec <- sub("_", " ", species_filtered[s])
          
          # adjust species name if BL uses other taxonomy than BBS:
          if(spec %in% spec_name_change_df$BBS_name){
            spec <- spec_name_change_df$BL_name[which(spec_name_change_df$BBS_name == spec)]
          }
          # extract shapefile:
          gdalUtilities::ogr2ogr(src_datasource_name = file.path(datashare_Birdlife, "BOTW.gdb"),
                                 layer = "All_Species",
                                 where = paste0("sci_name = '", spec, "' AND (SEASONAL = '1' OR SEASONAL = '2') AND PRESENCE = '1'"), # SEASONAL = 1: resident throughout the year, SEASONAL = 2: breeding season
                                 dst_datasource_name = file.path(data_dir, "Birdlife_ranges_BBS", "Shapefiles_2022", paste0(species_filtered[s], ".shp")),
                                 overwrite = TRUE)
        }


## rasterize and project Birdlife shapefiles: ----------------------------------

# for each species: convert Birdlife range polygons to raster, 
# project raster and change resolution to 50 km (same as EBBA)

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

# shapefiles:
BL_range_shps <- list.files(file.path(data_dir, "Birdlife_ranges_BBS", "Shapefiles_2022"), full.names = TRUE, pattern = ".shp")

# output folder for Birdlife range rasters:
if(!dir.exists(file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022"))){dir.create(file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022"), recursive = TRUE)}

# calculate rasters:
foreach(i = 1:length(species_filtered),
        .packages = c("gdalUtilities"), 
        .verbose = TRUE,
        .errorhandling = "remove",
        .inorder = FALSE) %dopar% {
          
          # rasterize polygon:
          gdalUtilities::gdal_rasterize(src_datasource = grep(pattern = species_filtered[i], BL_range_shps, value = TRUE),
                                        dst_filename = file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022", paste0(species_filtered[i], "_WGS84.tif")),
                                        burn = 1,
                                        tr = c(0.1, 0.1), # target resolution in degrees (same unit as src_datasource)
                                        a_nodata = -99999) # value for cells with missing data
          
          # project and change resolution to 50 km:
          gdalUtilities::gdalwarp(srcfile = file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022", paste0(species_filtered[i], "_WGS84.tif")), 
                                  dstfile = file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022", paste0(species_filtered[i], "_50km.tif")),
                                  overwrite = TRUE,
                                  tr = c(50000, 50000), # target resolution in meters (same as unit of target srs)
                                  r = "max", # resampling method
                                  t_srs = homolosine # target spatial reference
          )
          unlink(file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022", paste0(species_filtered[i], "_WGS84.tif"))) # delete unprojected raster
        }


## mask: overlay Birdlife ranges of all species: -------------------------------

# range raster files:
BL_range_tifs <- list.files(file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022"), full.names = TRUE, pattern = "_50km.tif$")

# load as list:
ranges <- lapply(BL_range_tifs, rast)

# overlay rasters, add 100 km buffer and resample to match projected Chelsa rasters:
# (buffer to make sure that even after resampling mask covers raster grid points of species ranges)
# (resampling aligns raster origins, necessary for masking)

# in case of fatal R errors while merging all rasters in one step, merge in multiple steps:
ranges_combined1 <- do.call(terra::merge, ranges[1:200]) %>% 
  buffer(width = 100000) %>% # Warning: "[merge] rasters did not align and were resampled" -> fine, aligns ranges of different species
  resample(y = rast(file.path(data_dir, "Chelsa", "CHELSA_tasmin_11_2014_V.2.1_50km.tif")), # projected Chelsa raster as template, output of 2_2_species_filtering_5_project_Chelsa.R
           method = "near") # 
ranges_combined2 <- do.call(terra::merge, ranges[201:312]) %>% 
  buffer(width = 100000) %>%
  resample(y = rast(file.path(data_dir, "Chelsa", "CHELSA_tasmin_11_2014_V.2.1_50km.tif")),
           method = "near")
ranges_combined <- do.call(terra::merge, list(ranges_combined1,ranges_combined2))

# save mask:
writeRaster(ranges_combined, 
            filename = file.path(data_dir, "BBS_Birdlife_ranges_mask.tif"),
            overwrite = TRUE, NAflag = FALSE)


## exploration: compare BBS records with range maps from Birdlife: -------------

world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")

# for which species should range be plotted:
spec = 1

# Birdlife range map:
BL_range_sf <- read_sf(BL_range_shps[grep(BL_range_shps, pattern = species_filtered[spec])])

# BBS records historic period:
BBS_spec_occ <- st_read(file.path(data_dir, "BBS_historic_centr_proj_hist81-83.shp")) %>%  # output of 1_BBS_prep_data.R
  filter(species == sub("_", " ", species_filtered[spec]) & pres == 1)

ggplot2::ggplot(world) +
  ggplot2::geom_sf(color = "gray50") +
  ggplot2::geom_sf(data = st_transform(BL_range_sf, crs = st_crs(world)), fill = "red", alpha = 0.5) +
  ggplot2::geom_sf(data = st_transform(BBS_spec_occ, crs = st_crs(world)), colour = "blue") +
  ggplot2::ggtitle(sub("_", " ", species_filtered[spec]))


# ------------------------------- #
#     2. Mask Chelsa data      ####
# ------------------------------- #

# Chelsa files on datashare:
chelsa_tifs <- list.files(datashare_Chelsa, full.names = FALSE, 
                          pattern = paste0("(", paste(2015:2018, collapse = "|"), ")_V.2.1.tif"))

# projected Chelsa files:
names_proj <- list.files(file.path(chelsa_proj_path, "Chelsa_projected"), full.names = TRUE)

# folder to store masked Chelsa data:
chelsa_masked_path <- file.path(data_dir, "Chelsa", "Chelsa_masked_global")
# create folder if it doesn't exist yet:
if(!dir.exists(chelsa_masked_path)){dir.create(chelsa_masked_path, recursive = TRUE)}

# file paths for the masked data:
names_masked <- file.path(chelsa_masked_path, 
                          paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_50km_masked.tif"))

# load mask:
ranges_mask <- rast(file.path(data_dir, "BBS_Birdlife_ranges_mask.tif"))

# mask projected rasters:
for(i in 1:length(chelsa_tifs)){
  
  print(i)
  
  names_proj[i] %>%
    rast %>%
    terra::mask(mask = ranges_mask) %>%
    terra::writeRaster(names_masked[i], overwrite = TRUE)
}


# --------------------------------------------------------------- #
#     3. Calculate bioclimatic variables from Chelsa data:     ####
# --------------------------------------------------------------- #

vars <- c("pr", "tas", "tasmin", "tasmax")
months <- stringr::str_pad(1:12, width = 2, pad = "0")
years <- 2015:2018

# masked Chelsa data:
chelsa_masked_files <- list.files(chelsa_masked_path, pattern = "tif$", full.names = TRUE)

# load mask using raster package, this is required by dismo::biovars())
ranges_mask <- raster(file.path(data_dir, "BBS_Birdlife_ranges_mask.tif"))

# spatial brick template, bricks to store the month-wise means across the selected years for each bioclimatic variable:
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
bioclim_folder <- file.path(data_dir, paste0("Bioclim_global_", min(years), "_", max(years)))
if(!dir.exists(bioclim_folder)){dir.create(bioclim_folder, recursive = TRUE)}

writeRaster(biovars_rast,
            filename = file.path(bioclim_folder, 
                                 paste0("CHELSA_", names(biovars), "_", min(years), "_", max(years), "_", "50km.tif")), 
            overwrite = TRUE)


# --------------------------------------- #
#     4. Climatic niche analyses:      ####
# --------------------------------------- #

# quantify how much of the global climatic niche of a species is covered by the conterminous US, 
# compare global range to part of the global range located in the conterminous US
# (don't use BBS occurrence records because BBS does not cover conterminous US as well as EBBA data cover Europe)
# use stability index as a measure to quantify niche coverage

# outline of conterminous US:

# load conterminous US and save as shapefile:
st_write(obj = spData::us_states, file.path(data_dir, "conterminousUS.shp"))

# rasterize and transform in the same way as BL shapefiles -> same accuracy for comparison of global range and conterminous US
gdalUtilities::gdal_rasterize(src_datasource = file.path(data_dir, "conterminousUS.shp"),
                              dst_filename = file.path(data_dir, "contUS.tif"),
                              burn = 1,
                              tr = c(0.1, 0.1),
                              a_nodata = -99999) 
gdalUtilities::gdalwarp(srcfile = file.path(data_dir, "contUS.tif"), 
                        dstfile = file.path(data_dir, "contUS_50km.tif"),
                        overwrite = TRUE,
                        tr = c(50000, 50000),
                        r = "max",
                        t_srs = homolosine)
unlink(file.path(data_dir, "contUS.tif"))

# convert raster back to polygon to apply spatial filter in a later step:
contUS_rast_poly <- rast(file.path(data_dir, "contUS_50km.tif")) %>% 
  as.polygons %>% 
  st_as_sf

# files of rasterized and projected Birdlife ranges:
BL_range_tifs

# data frame for PCA eigenvalues
BBS_global_niche_PCA <- data.frame(species=species_filtered,PCA_percent=NA)

# loop over species:
stability_df <- foreach(i = 1:length(species_filtered),
                        .combine = rbind,
                        .packages = c("ecospat", "ade4", "dplyr", "sf", "terra"), 
                        .verbose = TRUE,
                        .errorhandling = "remove") %dopar% {
                          
                          # data frame to store results:
                          stability_spec_df <- data.frame("species" = sub("_", " ", species_filtered[i]),
                                                          "stability" = NA)
                          
                          # read tifs with bioclim variables:
                          # need to be loaded within foreach since SpatRasters and SpatVectors are non-exportable objects
                          biovars_rast <- rast(file.path(bioclim_folder, 
                                                         paste0("CHELSA_", paste0("bio", 1:19), "_", min(years), "_", max(years), "_", "50km.tif")))
                          
                          # Birdlife range of species i; raster cell centroids as occurrence points:
                          BL_pts <- grep(pattern = species_filtered[i], BL_range_tifs, value = TRUE) %>% 
                            rast %>%
                            as.points %>% # convert raster to points (cell centroids)
                            st_as_sf
                          
                          # add bioclim variables:
                          BL_vars_df <- cbind(data.frame(species = sub("_", " ", species_filtered[i])),
                                                   terra::extract(biovars_rast, BL_pts)) %>% 
                            filter(complete.cases(.)) # in (few) cases where occurrence points may be "on the edge of the globe" so no Chelsa data are available
                          
                          # Birdlife cell centroids within the conterminous US:
                          contUS_pts <- BL_pts %>% 
                            st_join(y = contUS_rast_poly, join  = st_intersects, left = FALSE)
                          
                          # add bioclim variables:
                          contUS_vars_df <- cbind(data.frame(species = sub("_", " ", species_filtered[i])),
                                                   terra::extract(biovars_rast, contUS_pts)) %>% 
                            filter(complete.cases(.))

                          # assess climate niche by using the first 2 PCA axes:
                          # calibrating the PCA in the whole study area, including both ranges:
                          pca.env <- dudi.pca(rbind(contUS_vars_df, BL_vars_df)[, paste0("bio", 1:19)],
                                              scannf = FALSE,
                                              nf = 2) # number of axes
                          
                          # How much climate variation explained by first two axes:
                          BBS_global_niche_PCA[i,2] = sum(pca.env$eig[1:2]/sum( pca.env$eig))
                          
                          # predict the scores on the PCA axes:
                          # occurrences within conterminous US used as z1 (corresponds to native distribution in tutorials)
                          # Birdlife range used as z2 (corresponds to invasive distribution in tutorials)
                          scores.globclim <- pca.env$li # PCA scores for cont. US + Birdlife distribution
                          scores.clim.contUS <- suprow(pca.env, contUS_vars_df[, paste0("bio", 1:19)])$li # PCA scores for cont. US distribution
                          scores.clim.BL <- suprow(pca.env, BL_vars_df[, paste0("bio", 1:19)])$li # PCA scores for Birdlife distribution
                          
                          # calculate the Occurrence Densities Grid for cont. US and Birdlife distribution:
                          
                          # cont. US distribution area:
                          grid.clim.contUS <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                                                  glob1 = scores.clim.contUS, 
                                                                  sp = scores.clim.contUS,
                                                                  R = 100, # grid resolution
                                                                  th.sp = 0)
                          # Birdlife range area:
                          grid.clim.BL <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                                                glob1 = scores.clim.BL, 
                                                                sp = scores.clim.BL,
                                                                R = 100, 
                                                                th.sp = 0)

                          # assess niche dynamics between conterminous US range (as z1) and Birdlife range map (as z2) 
                          # => the stability index shows how much of the species global climatic niche 
                          # is covered by conterminous US:
                          
                          niche_dyn <- ecospat.niche.dyn.index(grid.clim.contUS,
                                                               grid.clim.BL,
                                                               intersection = NA) # analysis performed on cont. US + Birdlife extent
                          stability_spec_df$stability <- niche_dyn$dynamic.index.w['stability']
                          
                          # save plot of niche dynamics:
                          jpeg(filename = file.path(US_BL_plot_dir, paste0(species_filtered[i], ".jpg")),
                               quality = 100, height = 920, width = 920
                          )
                          
                          ecospat.plot.niche.dyn(grid.clim.contUS, 
                                                 grid.clim.BL, 
                                                 quant = 0.1,
                                                 interest = 2, # 1 = cont. US density, 2 = BL density
                                                 title = paste(species_filtered[i], "stability:", round(stability_spec_df$stability, 2)), 
                                                 name.axis1 = "PC1",
                                                 name.axis2 = "PC2")
                          dev.off()
                          # blue = stability
                          # red = expansion
                          # green = unfilling
                          
                          stability_spec_df
                        }

stability_df %>% 
  arrange(-stability) # 309 species
# 4 missing species:
# Anas_crecca: BL range does not cover conterminous US 
# Cistothorus_platensis: BL range does not cover conterminous US
# Columba_livia: BL range does not cover conterminous US
# Passerella_iliaca: BL breeding range does not cover conterminous US

# save species stability values :
write.csv(stability_df, 
          file = file.path(data_dir, "species_stability_contUS_BL22.csv"),
          row.names = FALSE)

write.csv(BBS_global_niche_PCA, 
          file = file = file.path(data_dir, "BBS_global_niche_PCA.csv"), 
          row.names=F)

# plots regarding stability:
plot(sort(stability_df$stability, decreasing = TRUE),
     ylab = "stability", xlab = "number of species", las = 1)
abline(h = 0.7, col = "red")
abline(h = 0.5, col = "blue")
