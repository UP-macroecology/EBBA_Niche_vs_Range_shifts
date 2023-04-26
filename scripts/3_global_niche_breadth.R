# determine niche breadth of global Birdlife breeding range using the Shannon index
# for EBBA / BBS species left after 4 filtering steps (2_1_EBBA/BBS_species_filtering_1-4.R)

library(dplyr)
library(terra)
library(doParallel)
library(ecospat)
library(ade4)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

data_dir <- file.path("Data")

dataset <- "EBBA"
#dataset <- "BBS"

# register cores for parallel computation:
registerDoParallel(cores = 2)
#getDoParWorkers() # check registered number of cores


# ------------------------------ #
#          Load data:         ####
# ------------------------------ #

if(dataset == "EBBA"){
  
  years <- 2009:2018
  
  bioclim_folder <- file.path(data_dir, "Bioclim_EBBA2")
  
  # output folder for niche plots:
  plots_dir <- file.path("plots", "EBBA_BL_niche_breadth")
  if(!dir.exists(plots_dir)){
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # species left after filtering step 4:
  load(file = file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_prelim.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
  species_filtered <- sort(sub(" ", "_", unique(EBBA1_prep$species)))
  
  # files of rasterized and projected Birdlife ranges:
  BL_range_tifs <- list.files(file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022"), full.names = TRUE, pattern = "_50km.tif$")
  
  
}

if(dataset == "BBS"){
  
  years <- 2016:2018
  
  bioclim_folder <- file.path(data_dir, paste0("Bioclim_raster_BBS_", min(years), "_", max(years)))
  
  # output folder for niche plots:
  plots_dir <- file.path("plots", "BBS_BL_niche_breadth")
  if(!dir.exists(plots_dir)){
    dir.create(plots_dir, recursive = TRUE)
  }
  
  # load species left after filter steps 1-4, for both versions of historic time period:
  load(file = file.path(data_dir, "BBS_prep_steps1-4_V1.RData")) # output of 2_1_BBS_species_filtering_1-4.R
  species_filtered_V1 <- sort(sub(" ", "_", unique(hist_prep_df$species)))
  load(file = file.path(data_dir, "BBS_prep_steps1-4_V2.RData")) # output of 2_1_BBS_species_filtering_1-4.R
  species_filtered_V2 <- sort(sub(" ", "_", unique(hist_prep_df$species)))
  species_filtered <- sort(unique(c(species_filtered_V1, species_filtered_V2)))
  
  # files of rasterized and projected Birdlife ranges:
  BL_range_tifs <- list.files(file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022"), full.names = TRUE, pattern = "_50km.tif$")
  
}


# ------------------------------ #
#       Niche breadth :       ####
# ------------------------------ #

# niche breadth is quantified using the Shannon index of the occurrence density grid corrected for environmental prevalence
# to be comparable across species, the environmental background must be the same:

# bioclim variables:
biovars_rast <- rast(file.path(bioclim_folder, 
                               paste0("CHELSA_", paste0("bio", 1:19), "_", min(years), "_", max(years), "_", "50km.tif")))

# global raster points (= environmental background, absences and presences):
BL_global_points <- biovars_rast[[1]]%>% 
  as.points

# add bioclim variables to raster points:
BL_vars_df <- terra::extract(biovars_rast, BL_global_points)

# assess climate niche by using the first 2 PCA axes:
# calibrating the PCA in the whole study area:
pca.env <- dudi.pca(BL_vars_df[, paste0("bio", 1:19)], scannf = FALSE,
                    nf = 2) # number of axes

# predict the scores on the PCA axes:
scores.globclim <- pca.env$li # PCA scores for global raster points


# loop over species:
niche_breadth_df <- foreach(i = 1:length(species_filtered),
                        .combine = rbind,
                        .packages = c("ecospat", "ade4", "dplyr","terra"), 
                        .verbose = TRUE,
                        .errorhandling = "remove") %dopar% {
                          
                          
                          # data frame to store results:
                          nb_spec_df <- data.frame("species" = sub("_", " ", species_filtered[i]),
                                                   "niche_breadth_zcor" = NA)
                          
                          # Birdlife range of species i:
                          BL_pts_spec <- grep(pattern = species_filtered[i], BL_range_tifs, value = TRUE) %>% 
                            rast %>% # load raster
                            as.points # convert raster to points (cell centroids)
                          
                          #plot(BL_global_points)
                          #plot(BL_pts_spec, add = TRUE, col = "red")
                          
                          # read tifs with bioclim variables:
                          # need to be loaded within foreach since SpatRasters and SpatVectors are non-exportable objects
                          biovars_rast <- rast(file.path(bioclim_folder, 
                                                         paste0("CHELSA_", paste0("bio", 1:19), "_", min(years), "_", max(years), "_", "50km.tif")))
                          
                          # add bioclim variables to occurrences:
                          BL_spec_vars_df <- cbind(data.frame(species = sub("_", " ", species_filtered[i])),
                                                   terra::extract(biovars_rast, BL_pts_spec)) %>% 
                            filter(complete.cases(.)) # should be complete, just to be sure
                      
                          # PCA scores for species occurrences:
                          scores.sp <- suprow(pca.env, BL_spec_vars_df[, paste0("bio", 1:19)])$li # PCA scores for BL species distribution
                          
                          # calculate the Occurrence Densities Grid for Birdlife distribution:
                          
                          # Birdlife range area:
                          grid.clim.BL <- ecospat.grid.clim.dyn(glob = scores.globclim, # background pixels of the whole study area
                                                                glob1 = scores.globclim, # same for background pixels of the species
                                                                sp = scores.sp, # environmental values for the occurrences of the species
                                                                R = 100, 
                                                                th.sp = 0)
                          
                          nb_spec_df$niche_breadth_zcor <- vegan::diversity(as.vector(grid.clim.BL$z.cor), index = "shannon")
                          
                          # save plot of niche dynamics:
                          jpeg(filename = file.path(plots_dir, paste0(species_filtered[i], ".jpg")),
                               quality = 100, height = 920, width = 920)
                          
                          ecospat.plot.niche(grid.clim.BL, title = paste(species_filtered[i], round(nb_spec_df$niche_breadth_zcor,2)), 
                                                   name.axis1 = "PC1", name.axis2 = "PC2", cor = TRUE)
                          dev.off()

                          nb_spec_df
                        }


niche_breadth_df
# save resulting data frame:
write.csv(niche_breadth_df, file = file.path(data_dir, paste0(dataset, "_niche_breadth.csv")), row.names = FALSE)


# explorations: ---

EBBA_nb <- read.csv(file.path(data_dir, "EBBA_niche_breadth.csv"))
BBS_nb <- read.csv(file.path(data_dir, "BBS_niche_breadth.csv"))

hist(EBBA_nb$niche_breadth_zcor)
hist(BBS_nb$niche_breadth_zcor)

# some species are included in EBBA as well as in BBS analysis,
# niche breadth differs because of different years on which bioclim variables are based on
EBBA_BBS_specs <- EBBA_nb %>% 
  inner_join(BBS_nb, by = c("species"))
# 22 species
cor.test(EBBA_BBS_specs$niche_breadth_zcor.x, EBBA_BBS_specs$niche_breadth_zcor.y) # highly correlated
