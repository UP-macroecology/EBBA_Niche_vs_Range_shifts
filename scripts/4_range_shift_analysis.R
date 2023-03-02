# Range overlap between time periods: 

# similar to niche analysis, but instead of considering two dimensional PCA space derived from climatic conditions,
# analysis is conducted in geographic space (axes = X and Y coordinates instead of PC1 and PC2)

# note: 
# - ecospat dependencies ‘biomod2’, ‘randomForest’ depend on R >= 4.1.0
# -> before running the script on the cluster run "module load R/4.1.0-foss-2021a" on ecoc9 to ensure that R version 4.1.0 (2021-05-18) is used

# outputs:
# for each species:
# - central range coordinate of historic and recent time period (centroid)
# - Euclidian distance between range centroids of historic and recent time period
# - north-south shift and east-west shift
# - size of historic and recent range (number of cells)
# - range overlap (Schoener's D)
# - results of range shift tests
# - results of range conservatism tests
# - range dynamic metrics: range stability, expansion and unfilling between historic and recent time period
# - standardised range dynamic metrics: stability, expansion, unfilling
# - plots of range dynamics

library(sf)
library(dplyr)
library(terra)
library(doParallel)
library(ecospat)
library(ade4)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# paths to data:

#EBBA_data_dir <- file.path("Data")
EBBA_data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")

# should environmental background be species-specific (500 km buffer around presences) or the same for all species (whole EBBA area):
env_background_species_specific <- TRUE # set to FALSE to use whole EBBA area as environmental background

# folder for range dynamics plots:
if(env_background_species_specific){
  plots_dir <- file.path(EBBA_data_dir, "plots", "EBBA_1_2_range_dyn_env_species_spec")
} else {
  plots_dir <- file.path(EBBA_data_dir, "plots", "EBBA_1_2_range_dyn_env_EBBA")
}
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# output file:
if(env_background_species_specific){
  results_file <- file.path(EBBA_data_dir, "range_shift_results_env_species_spec.csv")
} else {
  results_file <- file.path(EBBA_data_dir, "range_shift_results_env_EBBA.csv")
}


# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# comparable EBBA cells:
EBBA1_cells <- read_sf(file.path(EBBA_data_dir, "EBBA1_comparable.shp")) %>% # output of 1_prep_EBBA_data.R
  select(-species) %>% 
  distinct(cell50x50, .keep_all = TRUE) # keep geometry; would be the same when using EBBA2

# species selection: 
sel_species <- read.csv(file.path(EBBA_data_dir, "EBBA_niche_range_shifts_species_selection.csv"))$species %>% # output of 2_3_species_filtering_5_climatic_niche_analysis.R
  sort # alphabetically sorted (to avoid confusion with indices later)

# EBBA 1, only selected species:
EBBA1 <- read_sf(file.path(EBBA_data_dir, "EBBA1_comparable_harmonized.shp")) %>% # output of 1_prep_EBBA_data.R
  filter(species %in% sel_species)

# EBBA 2, only selected species:
EBBA2 <- read_sf(file.path(EBBA_data_dir, "EBBA2_comparable_harmonized.shp")) %>% # output of 1_prep_EBBA_data.R
  filter(species %in% sel_species)


# ------------------------------------------------ #
#       Range shift between time periods:     ####
# ------------------------------------------------ #

# data frame to store results:
range_shift_df <- data.frame("species" = sel_species,
                             "centroid_hist_X" = NA,
                             "centroid_hist_Y" = NA,
                             "centroid_rec_X" = NA,
                             "centroid_rec_Y" = NA,
                             "centroid_dist" = NA,
                             "centroid_NS_shift" = NA, # north-south shift
                             "centroid_EW_shift" = NA, # east-west shift
                             "n_cells_hist" = NA,
                             "n_cells_rec" = NA
                             )

# central range coordinate historic time of each species:
centroids_hist_sf <- EBBA1 %>%
  group_by(species) %>% 
  summarize(geometry = st_union(geometry), .groups = "keep") %>% 
  st_centroid

range_shift_df[, c("centroid_hist_X", "centroid_hist_Y")] <- st_coordinates(centroids_hist_sf)

# check:
#plot(st_geometry(EBBA1[which(EBBA1$"species" == sel_species[96]),]))
#plot(st_geometry(centroids_hist_sf[96,]), col = "red", pch = 19, add = TRUE)

# central range coordinate recent time of each species:
centroids_rec_sf <- EBBA2 %>%
  group_by(species) %>% 
  summarize(geometry = st_union(geometry), .groups = "keep") %>% 
  st_centroid

range_shift_df[, c("centroid_rec_X", "centroid_rec_Y")] <- st_coordinates(centroids_rec_sf)

# euclidian distance between range centroids of historic and recent time period:
range_shift_df$centroid_dist <- st_distance(centroids_hist_sf, centroids_rec_sf, by_element = TRUE) 

# north-south shift:
range_shift_df$centroid_NS_shift <- units::set_units(range_shift_df$centroid_rec_Y - range_shift_df$centroid_hist_Y, "m") # positive values = north shift, negative values = south shift

# east-west shift:
range_shift_df$centroid_EW_shift <- units::set_units(range_shift_df$centroid_rec_X - range_shift_df$centroid_hist_X, "m") # positive values = east shift, negative values = west shift

# number of cells of historic range:
range_shift_df$n_cells_hist <- EBBA1 %>% 
  st_drop_geometry %>% 
  group_by(species) %>% 
  summarize(n_cells = n(), .groups = "keep") %>% 
  pull(n_cells)

# number of cells of recent range:
range_shift_df$n_cells_rec <- EBBA2 %>% 
  st_drop_geometry %>% 
  group_by(species) %>% 
  summarize(n_cells = n(), .groups = "keep") %>% 
  pull(n_cells)


# ------------------------------------------------ #
#       Range overlap between time periods:     ####
# ------------------------------------------------ #

## mask to restrict occurrence density grids to land area: ----

# by first buffering comparable cell centroids and then using ALL_TOUCHED rasterization (at = TRUE), 
# we get a raster that approx. matches the comparable EBBA area and that has no holes due to single non-comparable cells

# (buffered cell centroids rasterized using the same resolution and extent as that of the occurrence density grid
# calculated with ecospat.grid.clim.dyn(..., R = 110) when the whole EBBA area is considered as the environmental background; 
# to get these values, I first calculated the grid without masking using the function ecospat.grid.clim.dyn() as below 
# but with geomask = NULL, and extracted the resolution and extent of the resulting grid)

EBBA1_cells %>% 
  st_buffer(dist = 25000) %>% # buffer to get polygons instead of points for rasterization
  st_write(file.path(EBBA_data_dir, "EBBA_comparable_buffer25km.shp"), append = FALSE)

gdalUtilities::gdal_rasterize(src_datasource = file.path(EBBA_data_dir, "EBBA_comparable_buffer25km.shp"),
                              dst_filename = file.path(EBBA_data_dir, "EBBA_mask_ecospat.tif"),
                              burn = 1,
                              at = TRUE, # ALL_TOUCHED rasterization: all pixels touched by polygons get value 1
                              tr = c(53558.6, 48762.57), # target resolution, values from occurrence density grid
                              a_nodata = -99999,
                              te = c(750722.8, 1240212, 6642168, 6604095)) # extent of occurrence density grid

# when using species specific backgrounds, the occurrence density grids of the individual species differ in terms of
# extent and resolution because ecospat only allows to specify the number of grid cells (argument R)
# by further converting the raster to a SpatialPolygons object we can use it to mask the occurrence density grids of all species
# independently of their resolution and extent:

EBBA_mask_sp <- rast(file.path(EBBA_data_dir, "EBBA_mask_ecospat.tif")) %>% 
  as.polygons %>% 
  st_as_sf %>% 
  as_Spatial # same as EBBA_mask, but as SpatialPolygons object


## environmental background (if environmental background = whole EBBA area for all species): ----

# coordinates of all comparable EBBA cells:
EBBA_cell_coords <- as.data.frame(st_coordinates(EBBA1_cells))


# loop over species:

# register cores for parallel computation:
registerDoParallel(cores = 10)
#getDoParWorkers() # check registered number of cores

range_ecospat_df <- foreach(s = 1:length(sel_species),
                            .combine = rbind,
                            .packages = c("ecospat", "ade4", "dplyr", "sf"),
                            .verbose = TRUE) %dopar% {
                              
                              # data frame to store results for single species:
                              range_ecospat_s <- data.frame("species" = sel_species[s],
                                                            # range overlap:
                                                            "D" = NA, # Schoener's D
                                                            # results of similarity test when testing for range shift:
                                                            "shift_p_D_NA" = NA, # NA = non-analogue conditions
                                                            "shift_p_E_NA" = NA,
                                                            "shift_p_S_NA" = NA,
                                                            "shift_p_U_NA" = NA,
                                                            "shift_p_D_A" = NA, # A = analogue conditions
                                                            "shift_p_E_A" = NA,
                                                            "shift_p_S_A" = NA,
                                                            "shift_p_U_A" = NA,
                                                            # results of similarity test when testing for range conservatism:
                                                            "cons_p_D_NA" = NA,
                                                            "cons_p_E_NA" = NA,
                                                            "cons_p_S_NA" = NA,
                                                            "cons_p_U_NA" = NA,
                                                            "cons_p_D_A" = NA,
                                                            "cons_p_E_A" = NA,
                                                            "cons_p_S_A" = NA,
                                                            "cons_p_U_A" = NA,
                                                            # range dynamic metrics:
                                                            "range_stability_NA" = NA,
                                                            "range_expansion_NA" = NA,
                                                            "range_unfilling_NA" = NA,
                                                            "range_stability_A" = NA,
                                                            "range_expansion_A" = NA,
                                                            "range_unfilling_A" = NA)
                              
                              ## coordinates of presence / absence points: -----
                              
                              # EBBA 1 occurrence density grid:
                              EBBA1_presence <- EBBA1 %>% 
                                filter(species == sel_species[s]) %>% 
                                st_coordinates %>% 
                                as.data.frame
                               
                              if(env_background_species_specific){
                                
                                # EBBA 1 environmental background:
                                EBBA1_spec_buffer <- EBBA1 %>% 
                                  filter(species == sel_species[s]) %>% 
                                  st_buffer(dist = 500000) # buffer to determine cells used as absences
                                
                                EBBA1_env <- EBBA1_cells %>% 
                                  st_filter(EBBA1_spec_buffer, .pred = st_within) %>% 
                                  st_coordinates %>% 
                                  as.data.frame
                                
                              }
                              
                              # EBBA 2 occurrence density grid:
                              EBBA2_presence <- EBBA2 %>% 
                                filter(species == sel_species[s]) %>% 
                                st_coordinates %>% 
                                as.data.frame
                              
                              if(env_background_species_specific){
                                
                                # EBBA 2 environmental background:
                                EBBA2_spec_buffer <- EBBA2 %>%
                                  filter(species == sel_species[s]) %>% 
                                  st_buffer(dist = 500000) # buffer to determine cells used as absences
                                
                                EBBA2_env <- EBBA1_cells %>% 
                                  st_filter(EBBA2_spec_buffer, .pred = st_within) %>% 
                                  st_coordinates %>% 
                                  as.data.frame
                                }
    
                              
                              ## occurrence densities grids for EBBA1 and EBBA2 distribution: ------
                              
                              if(env_background_species_specific){
                                
                                glob <- rbind(EBBA1_env, EBBA2_env) # coordinates instead of PCA scores
                                glob1_EBBA1 <- EBBA1_env
                                glob1_EBBA2 <- EBBA2_env
                                
                                } else {
                                  
                                  glob <- EBBA_cell_coords
                                  glob1_EBBA1 <- EBBA_cell_coords
                                  glob1_EBBA2 <- EBBA_cell_coords
                                }
                              
                              grid.range.EBBA1 <- ecospat.grid.clim.dyn(glob = glob, # whole study area, both time periods
                                                                        glob1 = glob1_EBBA1, # env. background
                                                                        sp = EBBA1_presence, # species occurrence
                                                                        R = 110, # grid resolution, results together with extended extent in spatial resolution of 53 x 49 km when using whole EBBA area as background
                                                                        geomask = EBBA_mask_sp,
                                                                        extend.extent = c(-2e+05, 2e+05, -2e+05, 2e+05)) # relative to original extent
                              
                              # check:
                              #ecospat.plot.niche(z = grid.range.EBBA1, name.axis1 = "X", name.axis2 = "Y", title = sel_species[s])
                              #points(EBBA_cell_coords, col = "blue", cex = 0.3)
                              #points(EBBA1_presence, col = "red", cex = 0.3)
                  
                              # gridding EBBA2 range:
                              grid.range.EBBA2 <- ecospat.grid.clim.dyn(glob = glob, # whole study area, both time periods
                                                                        glob1 = glob1_EBBA2, # env. background
                                                                        sp = EBBA2_presence, # species occurrence
                                                                        R = 110, # grid resolution, results together with extended extent in spatial resolution of 53 x 49 km when using whole EBBA area as background
                                                                        geomask = EBBA_mask_sp,
                                                                        extend.extent = c(-2e+05, 2e+05, -2e+05, 2e+05)) # relative to original extent
                              
                              # check:
                              #ecospat.plot.niche(z = grid.range.EBBA2, name.axis1 = "X", name.axis2 = "Y", title = sel_species[s])
                              #points(EBBA_cell_coords, col = "blue", cex = 0.3)
                              #points(EBBA2_presence, col = "red", cex = 0.3)
                              
                              
                              # range overlap: Schoener's D:
                              range_ecospat_s$D <- ecospat.niche.overlap(grid.range.EBBA1, grid.range.EBBA2, cor = TRUE)$D # correct occurrence densities by prevalence of the environments
                              
                              ## range similarity tests: -----------------------
                              
                              ### tests for range shifting: ----
                              
                              # non-analogue conditions:
                              sim_test_shift_NA <- ecospat.niche.similarity.test(grid.range.EBBA1,
                                                                                 grid.range.EBBA2,
                                                                                 rep = 1000,
                                                                                 overlap.alternative = "lower", # test for niche shifting / divergence: niche overlap is lower than random
                                                                                 expansion.alternative = "higher", # test for niche shifting: expansion higher than random
                                                                                 stability.alternative = "lower", # test for niche shifting: stability lower than random
                                                                                 unfilling.alternative = "higher", # test for niche shifting: unfilling higher than random
                                                                                 intersection = NA, 
                                                                                 rand.type = 2)
                              
                              range_ecospat_s$shift_p_D_NA <- sim_test_shift_NA$p.D
                              range_ecospat_s$shift_p_E_NA <- sim_test_shift_NA$p.expansion
                              range_ecospat_s$shift_p_S_NA <- sim_test_shift_NA$p.stability
                              range_ecospat_s$shift_p_U_NA <- sim_test_shift_NA$p.unfilling
                              
                              # plot results of similarity test:
                              #ecospat.plot.overlap.test(sim_test_shift_NA, "D", "Similarity")
                              
                              # analogue conditions:
                              sim_test_shift_A <- ecospat.niche.similarity.test(grid.range.EBBA1,
                                                                                grid.range.EBBA2,
                                                                                rep = 1000,
                                                                                overlap.alternative = "lower", # test for niche shifting / divergence: niche overlap is lower than random
                                                                                expansion.alternative = "higher", # test for niche shifting: expansion higher than random
                                                                                stability.alternative = "lower", # test for niche shifting: stability lower than random
                                                                                unfilling.alternative = "higher", # test for niche shifting: unfilling higher than random
                                                                                intersection = 0, 
                                                                                rand.type = 2)
                              
                              range_ecospat_s$shift_p_D_A <- sim_test_shift_A$p.D
                              range_ecospat_s$shift_p_E_A <- sim_test_shift_A$p.expansion
                              range_ecospat_s$shift_p_S_A <- sim_test_shift_A$p.stability
                              range_ecospat_s$shift_p_U_A <- sim_test_shift_A$p.unfilling
                              
                              
                              ### tests for range conservatism: ----
                              
                              # non-analogue conditions:
                              sim_test_cons_NA <- ecospat.niche.similarity.test(grid.range.EBBA1, 
                                                                                grid.range.EBBA2, 
                                                                                rep = 1000, 
                                                                                overlap.alternative = "higher", # test for niche conservatism: niche overlap is higher than random
                                                                                expansion.alternative = "lower", # test for niche conservatism: expansion lower than random
                                                                                stability.alternative = "higher", # test for niche conservatism: stability higher than random
                                                                                unfilling.alternative = "lower", # test for niche conservatism: unfilling lower than random
                                                                                intersection = NA, 
                                                                                rand.type = 2)
                              
                              range_ecospat_s$cons_p_D_NA <- sim_test_cons_NA$p.D
                              range_ecospat_s$cons_p_E_NA <- sim_test_cons_NA$p.expansion
                              range_ecospat_s$cons_p_S_NA <- sim_test_cons_NA$p.stability
                              range_ecospat_s$cons_p_U_NA <- sim_test_cons_NA$p.unfilling
                              
                              # plot results of similarity test:
                              #ecospat.plot.overlap.test(sim_test_cons_NA, "D", "Similarity")
                              
                              # analogue conditions:
                              sim_test_cons_A <- ecospat.niche.similarity.test(grid.range.EBBA1, 
                                                                               grid.range.EBBA2, 
                                                                               rep = 1000, 
                                                                               overlap.alternative = "higher", # test for niche conservatism: niche overlap is higher than random
                                                                               expansion.alternative = "lower", # test for niche conservatism: expansion lower than random
                                                                               stability.alternative = "higher", # test for niche conservatism: stability higher than random
                                                                               unfilling.alternative = "lower", # test for niche conservatism: unfilling lower than random
                                                                               intersection = 0, 
                                                                               rand.type = 2)
                              
                              range_ecospat_s$cons_p_D_A <- sim_test_cons_A$p.D
                              range_ecospat_s$cons_p_E_A <- sim_test_cons_A$p.expansion
                              range_ecospat_s$cons_p_S_A <- sim_test_cons_A$p.stability
                              range_ecospat_s$cons_p_U_A <- sim_test_cons_A$p.unfilling
                              
                              
                              ## indices for range dynamics: -------------------
                   
                              # between EBBA1 (as z1) and EBBA2 (as z2)
                              
                              # non-analogue conditions:
                              EBBA12_range_dyn_NA <- ecospat.niche.dyn.index(grid.range.EBBA1, 
                                                                             grid.range.EBBA2,
                                                                             intersection = NA) # analysis performed on both ranges of EBBA1 and EBBA2
                              
                              range_ecospat_s$range_stability_NA <- EBBA12_range_dyn_NA$dynamic.index.w['stability']
                              range_ecospat_s$range_expansion_NA <- EBBA12_range_dyn_NA$dynamic.index.w['expansion']
                              range_ecospat_s$range_unfilling_NA <- EBBA12_range_dyn_NA$dynamic.index.w['unfilling']
                              
                              # analogue conditions:
                              EBBA12_range_dyn_A <- ecospat.niche.dyn.index(grid.range.EBBA1, 
                                                                            grid.range.EBBA2,
                                                                            intersection = 0) # analysis performed on intersection of ranges of EBBA1 and EBBA2
                              
                              range_ecospat_s$range_stability_A <- EBBA12_range_dyn_A$dynamic.index.w['stability']
                              range_ecospat_s$range_expansion_A <- EBBA12_range_dyn_A$dynamic.index.w['expansion']
                              range_ecospat_s$range_unfilling_A <- EBBA12_range_dyn_A$dynamic.index.w['unfilling']
                              
                              # save plot of range dynamics:
                              jpeg(filename = file.path(plots_dir, paste0(sel_species[s], ".jpg")),
                                   quality = 100, height = 920, width = 920)
                              
                              ecospat.plot.niche.dyn(grid.range.EBBA1,
                                                     grid.range.EBBA2,
                                                     quant = 0.05,
                                                     interest = 2, # 1 = EBBA1 density, 2 = EBBA2 density
                                                     title = paste(sel_species[s],":",
                                                                   "stability:", round(range_ecospat_s$range_stability_NA,3),
                                                                   "expansion:", round(range_ecospat_s$range_expansion_NA,3),
                                                                   "unfiling:", round(range_ecospat_s$range_unfilling_NA,3)),
                                                     name.axis1 = "X",
                                                     name.axis2 = "Y")
                              dev.off()
                              # blue = stability
                              # red = expansion
                              # green = unfilling
                              
                              range_ecospat_s
                            }


## standardise range dynamic metrics: ------------------------------------------

# non-analogue area (= accessible area in either EBBA1 or EBBA2):
e <- range_ecospat_df$range_expansion_NA # expansion related to whole accessible EBBA area (EBBA1 or EBBA2)
s <- range_ecospat_df$range_stability_NA # stability related to whole accessible EBBA area (EBBA1 or EBBA2)
u <- range_ecospat_df$range_unfilling_NA # unfilling related to whole accessible EBBA area (EBBA1 or EBBA2)

# stability, non-analogue conditions, related to EBBA1 accessible area:
s_x <- 1 - u

# proportion stability related to whole accessible area and stability related to EBBA1 accessible area (= proportion by which one accessible area is larger / smaller than the other)
x <- s / s_x

# total accessible area (non-analogue):
t <- u * x + s + e # unfilling rescaled to relate to EBBA2 accessible area

# unfilling, stability and expansion standardized so that they sum up to one: towards total non-analogue area
u_rel <- u * x / t
s_rel <- s / t
e_rel <- e / t

# target values:
range_ecospat_df$range_stability_std <- s_rel # stability
range_ecospat_df$range_unfilling_std <- u_rel # unfilling
range_ecospat_df$range_expansion_std <- e_rel # expansion

# check: values should sum to 1:
total <-  range_ecospat_df$range_unfilling_std + 
  range_ecospat_df$range_stability_std + 
  range_ecospat_df$range_expansion_std

# merge results:
range_shift_all_df <- cbind(range_shift_df, range_ecospat_df[,-1])

# save resulting data frame:
write.csv(range_shift_all_df, file = results_file, row.names = FALSE)