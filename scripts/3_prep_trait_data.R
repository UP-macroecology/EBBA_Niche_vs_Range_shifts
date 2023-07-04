# Assemble trait information for European and North American bird species:
# 1) traits from AVONET data
# 2) calculate global climatic niche breadth of Birdlife range polygons based on bioclimatic variables

library(dplyr)
library(terra)
library(doParallel)
library(ecospat)
library(ade4)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

#dataset <- "EBBA"
dataset <- "BBS"

# project data:
if(dataset == "EBBA"){
  data_dir <- file.path("data", "EBBA_analysis")
} else {
  data_dir <- file.path("data", "BBS_analysis") 
}

# register cores for parallel computation:
registerDoParallel(cores = 2)
#getDoParWorkers() # check registered number of cores


# ------------------------------------- #
#           AVONET trait data:       ####
# ------------------------------------- #

# read AVONET data:

# Tobias et al. 2022: AVONET: morphological, ecological and geographical data for all birds.
# I downloaded the AVONET dataset from https://figshare.com/s/b990722d72a26b5bfead -> AVONET Supplementary dataset 1
# (link in publication)
# (I didn't use the package "traitdata" to access AVONET because the package seems to use 
# the data provided to reproduce the analyses of Tobias at al. 2022 instead of the full AVONET data set
# (e.g. species names only available in one unspecified format while main AVONET data set provides 3 formats))

# sheet with taxonomy according to Birdlife International:
avonet_BL <- readxl::read_xlsx(path = file.path("data", "AVONET Supplementary dataset 1.xlsx"),
                               sheet = "AVONET1_BirdLife")

# join AVONET data to selected species:
if(dataset == "EBBA"){
  
  # final species selection:
  load(file = file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_final.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
  EBBA_species_final <- sort(unique(EBBA1_prep$species))

  # AVONET data for these species:
  EBBA_avonet <- subset(avonet_BL, avonet_BL$Species1 %in% EBBA_species_final)
  # all selected EBBA species are found in AVONET
  write.csv(EBBA_avonet, file = file.path(data_dir, "AVONET_EBBA_species.csv"),
            row.names = FALSE)
  
} else {
  
  # selected species:
  BBS_species_final <- read.csv(file = file.path(data_dir, "species_stability_contUS_BL22.csv")) %>% 
    filter(stability >= 0.5) %>% 
    pull(species)
  
  # for some species BBS uses different names than BL:
  # (used also in 2_3_BBS_species_filtering_5_climatic_niche_analysis.R)
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
  
  # match AVONET data either on BBS species name or on BL species name:
  
  BBS_sel_species_df <- bind_rows( # join on either BBS name or BL name and bind rows
    
    # join on Birdlife species name:
    data.frame("BBS_species" = BBS_species_final) %>% 
      left_join(spec_name_change_df, by = c("BBS_species" = "BBS_name")) %>% 
      inner_join(avonet_BL, by = c("BL_name" = "Species1")),
    
    # join on BBS species name:
    data.frame("BBS_species" = BBS_species_final) %>% 
      left_join(spec_name_change_df, by = c("BBS_species" = "BBS_name")) %>% 
      inner_join(avonet_BL, by = c("BBS_species" = "Species1"))
  ) %>% 
    mutate(BL_name = ifelse(is.na(BL_name), BBS_species, BL_name)) %>% 
    arrange(BBS_species)
  
  write.csv(BBS_sel_species_df, file = file.path(data_dir, "AVONET_BBS_species.csv"),
            row.names = FALSE)
  
}


# --------------------------------- #
#        global niche breadth    ####
# --------------------------------- #

# quantify niche breadth of global Birdlife breeding range using the Shannon index


## load data:  -----------------------------------------------------------------

if(dataset == "EBBA"){
  
  years <- 2012:2017
  
  bioclim_folder <- file.path(data_dir, paste0("Bioclim_global_", min(years), "_", max(years)))
  
  # output folder for niche plots:
  plots_dir <- file.path("plots", "niche_breadth_species", "EBBA_BL_niche_breadth")
  if(!dir.exists(plots_dir)){dir.create(plots_dir, recursive = TRUE)}
  
  # species left after filtering step 5:
  load(file = file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_final.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
  species_filtered <- sort(sub(" ", "_", unique(EBBA1_prep$species)))
  
  # files of rasterized and projected Birdlife ranges:
  BL_range_tifs <- list.files(file.path(data_dir, "Birdlife_ranges_EBBA", "Raster_2022"), full.names = TRUE, pattern = "_50km.tif$")
  
}

if(dataset == "BBS"){
  
  years <- 2015:2018
  # years <- 2016:2018
  
  bioclim_folder <- file.path(data_dir, paste0("Bioclim_global_", min(years), "_", max(years)))
  
  # output folder for niche plots:
  plots_dir <- file.path("plots", "niche_breadth_species", "BBS_BL_niche_breadth")
  if(!dir.exists(plots_dir)){dir.create(plots_dir, recursive = TRUE)}
  
  # selected species (both versions of historic time period):
  species_filtered <- read.csv(file = file.path(data_dir, "species_stability_contUS_BL22.csv")) %>% 
    filter(stability >= 0.5) %>% 
    pull(species) %>%
    sub(" ", "_", .)
  
  # files of rasterized and projected Birdlife ranges:
  BL_range_tifs <- list.files(file.path(data_dir, "Birdlife_ranges_BBS", "Raster_2022"), full.names = TRUE, pattern = "_50km.tif$")
  
}


## calculate niche breadth: ------------------------------------------------

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

# save resulting data frame:
write.csv(niche_breadth_df, file = file.path(data_dir, paste0(dataset, "_niche_breadth.csv")), row.names = FALSE)


## explorations: ---------------------------------------------------------------

EBBA_nb <- read.csv(file.path("data", "EBBA_analysis", "EBBA_niche_breadth.csv"))
BBS_nb <- read.csv(file.path("data", "BBS_analysis", "BBS_niche_breadth.csv"))

hist(EBBA_nb$niche_breadth_zcor)
hist(BBS_nb$niche_breadth_zcor)

# some species are included in EBBA as well as in BBS analysis,
# niche breadth differs because of different years on which bioclim variables are based on
EBBA_BBS_specs <- EBBA_nb %>% 
  inner_join(BBS_nb, by = c("species"))

cor.test(EBBA_BBS_specs$niche_breadth_zcor.x, EBBA_BBS_specs$niche_breadth_zcor.y) # highly correlated