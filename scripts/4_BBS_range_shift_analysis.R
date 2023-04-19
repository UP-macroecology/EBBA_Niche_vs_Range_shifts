# BBS range overlap between time periods: 

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

# settings:

# environmental background: use presences and absences within 600 km buffer around presences or all true absences within conterminous US:
env_background_contUS <- FALSE # set to FALSE to use presences and absences within 600 km buffer around presences

# which historic time period should be used:
hist_years <- 1980:1983 # maximum gap between historic and recent time period
#hist_years <- 1995:1998 # similar gap between historic and recent time period as in EBBA analysis

# paths to data:

# project data:
#data_dir <- file.path("Data")
data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")

# files to store outputs: 

# folder for range dynamics plots:
if(env_background_contUS){
  plots_dir <- file.path(data_dir, "plots", paste0("BBS_range_dyn_env_contUS_hist", min(hist_years), "_", max(hist_years)))
} else {
  plots_dir <- file.path(data_dir, "plots", paste0("BBS_range_dyn_env_600km_hist", min(hist_years), "_", max(hist_years)))
}
if(!dir.exists(plots_dir)){
  dir.create(plots_dir, recursive = TRUE)
}

# output files:
if(env_background_contUS){
  results_file <- file.path(data_dir, paste0("BBS_range_shift_results_env_contUS_hist", min(hist_years), "_", max(hist_years), ".csv"))
} else {
  results_file <- file.path(data_dir, paste0("BBS_range_shift_results_env_600km_hist", min(hist_years), "_", max(hist_years), ".csv"))
}

# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# species selection: 
sel_species <- read.csv(file = file.path(data_dir, "species_stability_contUS_BL22_120423.csv")) %>% 
  arrange(-stability) %>% 
  filter(stability >= 0.5) %>% 
  pull(species) %>% 
  sort

# BBS data, only selected species:
if(all(hist_years == 1980:1983)){
  
  # which species for 1980:1983:
  load(file = file.path(data_dir, "BBS_prep_steps1-4_V1.RData")) # output of 2_1_BBS_species_filtering_1-4.R
  species_filtered_V1 <- sort(unique(hist_prep_df$species))
  sel_species_final <- species_filtered_V1[which(species_filtered_V1 %in% sel_species)] # 207
  
  BBS_hist <- read_sf(file.path(data_dir, "BBS_historic1_proj_centroids.shp")) %>% 
    filter(species %in% sel_species_final)
  
  BBS_rec <- read_sf(file.path(data_dir, "BBS_recent1_proj_centroids.shp")) %>% 
    filter(species %in% sel_species_final)
}

if(all(hist_years == 1995:1998)){
  
  # which species for V2:
  load(file = file.path(data_dir, "BBS_prep_steps1-4_V2.RData")) # output of 2_1_BBS_species_filtering_1-4.R
  species_filtered_V2 <- sort(unique(hist_prep_df$species))
  sel_species_final <- species_filtered_V2[which(species_filtered_V2 %in% sel_species)] # 296
  
  BBS_hist <- read_sf(file.path(data_dir, "BBS_historic2_proj_centroids.shp")) %>% 
    filter(species %in% sel_species_final)
  
  BBS_rec <- read_sf(file.path(data_dir, "BBS_recent2_proj_centroids.shp")) %>% 
    filter(species %in% sel_species_final)
}


# ------------------------------------------------ #
#       Range shift between time periods:     ####
# ------------------------------------------------ #

# data frame to store results:
range_shift_df <- data.frame("species" = sel_species_final,
                             "centroid_hist_X" = NA,
                             "centroid_hist_Y" = NA,
                             "centroid_rec_X" = NA,
                             "centroid_rec_Y" = NA,
                             "centroid_dist" = NA,
                             "centroid_NS_shift" = NA, # north-south shift
                             "centroid_EW_shift" = NA, # east-west shift
                             "n_routes_hist" = NA,
                             "n_routes_rec" = NA
)

# central range coordinate historic time of each species:
centroids_hist_sf <- BBS_hist %>%
  filter(pres == 1) %>% 
  group_by(species) %>% 
  summarize(geometry = st_union(geometry), .groups = "keep") %>% 
  st_centroid

range_shift_df[, c("centroid_hist_X", "centroid_hist_Y")] <- st_coordinates(centroids_hist_sf)

# check:
#plot(st_geometry(BBS_hist[which(BBS_hist$species == sel_species_final[96] & BBS_hist$pres == 1),]))
#plot(st_geometry(centroids_hist_sf[96,]), col = "red", pch = 19, add = TRUE)

# central range coordinate recent time of each species:
centroids_rec_sf <- BBS_rec %>%
  filter(pres == 1) %>%
  group_by(species) %>% 
  summarize(geometry = st_union(geometry), .groups = "keep") %>% 
  st_centroid

range_shift_df[, c("centroid_rec_X", "centroid_rec_Y")] <- st_coordinates(centroids_rec_sf)

# Euclidian distance between range centroids of historic and recent time period:
range_shift_df$centroid_dist <- st_distance(centroids_hist_sf, centroids_rec_sf, by_element = TRUE) 

# north-south shift:
range_shift_df$centroid_NS_shift <- units::set_units(range_shift_df$centroid_rec_Y - range_shift_df$centroid_hist_Y, "m") # positive values = north shift, negative values = south shift

# east-west shift:
range_shift_df$centroid_EW_shift <- units::set_units(range_shift_df$centroid_rec_X - range_shift_df$centroid_hist_X, "m") # positive values = east shift, negative values = west shift

# number of routes with presences of historic range:
range_shift_df$n_routes_hist <- BBS_hist %>% 
  st_drop_geometry %>% 
  filter(pres == 1) %>%
  group_by(species) %>% 
  summarize(n_routes = n(), .groups = "keep") %>% 
  pull(n_routes)

# number of routes with presences of recent range:
range_shift_df$n_routes_rec <- BBS_rec %>% 
  st_drop_geometry %>% 
  filter(pres == 1) %>%
  group_by(species) %>% 
  summarize(n_routes = n(), .groups = "keep") %>% 
  pull(n_routes)


# ------------------------------------------------ #
#       Range overlap between time periods:     ####
# ------------------------------------------------ #


## mask to restrict occurrence density grids to land area: ----
# conterminous US, Canada, Mexico:
NA_mask <- spData::world %>% 
  filter(name_long %in% c("Canada", "United States", "Mexico")) %>% 
  st_transform(crs = "ESRI:102003") %>% 
  st_union %>% 
  as_Spatial
# using species specific backgrounds, the occurrence density grids of the individual species differ in terms of
# extent and resolution because ecospat only allows to specify the number of grid cells (argument R)


# loop over species:----

# register cores for parallel computation:
registerDoParallel(cores = 2)
#getDoParWorkers() # check registered number of cores

range_ecospat_df <- foreach(s = 1:2,#length(sel_species_final),
                            .combine = rbind,
                            .packages = c("ecospat", "ade4", "dplyr", "sf"),
                            .verbose = TRUE,
                            .errorhandling = "remove",
                            .inorder = FALSE) %dopar% {
                              
                              spec <- sel_species_final[s]
                              
                              # data frame to store results for single species:
                              range_ecospat_s <- data.frame("species" = spec,
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
                              
                              # historic presences:
                              
                              hist_presences <- BBS_hist %>%
                                filter(species == spec & pres == 1) %>% 
                                st_coordinates %>% 
                                as.data.frame
                              
                              if(!env_background_contUS){
                                
                                # environmental background = all presences and true absences within 600 km of presences (Sofaer et al. 2018)
                                
                                ## 600 km buffer around presences:
                                background_bf <- BBS_hist %>%
                                  filter(species == spec & pres == 1) %>% 
                                  st_buffer(dist = 600000) %>% 
                                  st_union
                                
                                ## species records (presences and true absences) within 600 km of presences:
                                background_hist <- BBS_hist %>%
                                  filter(species == spec) %>% 
                                  st_filter(y = background_bf, .predicate = st_within) %>% 
                                  st_coordinates %>% 
                                  as.data.frame
                                
                              } else {
                                
                                # environmental background = all true absences of a species across conterminous US:
                                background_hist <- BBS_hist %>%
                                  filter(species == spec) %>% 
                                  st_coordinates %>% 
                                  as.data.frame
                              }
                              
                              
                              # recent presences:
                              
                              rec_presences <- BBS_rec %>%
                                filter(species == spec & pres == 1) %>% 
                                st_coordinates %>% 
                                as.data.frame
                              
                              if(!env_background_contUS){
                                
                                # environmental background = all presences and true absences within 600 km of presences (Sofaer et al. 2018)
                                
                                ## 600 km buffer around presences:
                                background_bf <- BBS_rec %>%
                                  filter(species == spec & pres == 1) %>% 
                                  st_buffer(dist = 600000) %>% 
                                  st_union
                                
                                ## species records (presences and true absences) within 600 km of presences:
                                background_rec <- BBS_rec %>%
                                  filter(species == spec) %>% 
                                  st_filter(y = background_bf, .predicate = st_within) %>% 
                                  st_coordinates %>% 
                                  as.data.frame
                                
                              } else {
                                
                                # environmental background = all true absences of a species across conterminous US:
                                background_rec <- BBS_rec %>%
                                  filter(species == spec) %>% 
                                  st_coordinates %>% 
                                  as.data.frame
                              }
                              
                              
                              ## occurrence densities grids for historic and recent distribution: ------

                              glob <- rbind(background_hist, background_rec) # coordinates instead of PCA scores
                              glob1_hist <- background_hist
                              glob1_rec <- background_rec

                              grid.range.hist <- ecospat.grid.clim.dyn(glob = glob, # whole study area, both time periods
                                                                       glob1 = glob1_hist, # env. background
                                                                       sp = hist_presences, # species occurrence
                                                                       R = 100,
                                                                       geomask = NA_mask,
                                                                       extend.extent = c(-2e+05, 2e+05, -2e+05, 2e+05) # relative to original extent
                                                                       ) 
                              
                              # check:
                              # ecospat.plot.niche(z = grid.range.hist, name.axis1 = "X", name.axis2 = "Y", title = spec)
                              # plot(NA_mask, add = TRUE, border= "pink")
                              # points(background_hist, col = "blue", cex = 0.3)
                              # points(hist_presences, col = "red", cex = 0.3)
                              
                              # gridding recent range:
                              grid.range.rec <- ecospat.grid.clim.dyn(glob = glob, # whole study area, both time periods
                                                                      glob1 = glob1_rec, # env. background
                                                                      sp = rec_presences, # species occurrence
                                                                      R = 100, # grid resolution, results together with extended extent in spatial resolution of 53 x 49 km when using whole EBBA area as background
                                                                      geomask = NA_mask,
                                                                      extend.extent = c(-2e+05, 2e+05, -2e+05, 2e+05)) # relative to original extent
                              
                              # check:
                              # ecospat.plot.niche(z = grid.range.rec, name.axis1 = "X", name.axis2 = "Y", title = spec)
                              # plot(NA_mask, add = TRUE, border= "pink")
                              # points(background_rec, col = "blue", cex = 0.3)
                              # points(rec_presences, col = "red", cex = 0.3)
                              
                              # range overlap: Schoener's D:
                              range_ecospat_s$D <- ecospat.niche.overlap(grid.range.rec, grid.range.hist, cor = TRUE)$D # correct occurrence densities by prevalence of the environments
                              
                              ## range similarity tests: -----------------------
                              
                              ### tests for range shifting: ----
                              
                              # non-analogue conditions:
                              sim_test_shift_NA <- ecospat.niche.similarity.test(grid.range.hist,
                                                                                 grid.range.rec,
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
                              sim_test_shift_A <- ecospat.niche.similarity.test(grid.range.hist,
                                                                                grid.range.rec,
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
                              sim_test_cons_NA <- ecospat.niche.similarity.test(grid.range.hist, 
                                                                                grid.range.rec, 
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
                              sim_test_cons_A <- ecospat.niche.similarity.test(grid.range.hist, 
                                                                               grid.range.rec, 
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
                              
                              # between historic (as z1) and recent time period (as z2)
                              
                              # non-analogue conditions:
                              range_dyn_NA <- ecospat.niche.dyn.index(grid.range.hist, 
                                                                      grid.range.rec,
                                                                      intersection = NA) # analysis performed on both ranges of EBBA1 and EBBA2
                              
                              range_ecospat_s$range_stability_NA <- range_dyn_NA$dynamic.index.w['stability']
                              range_ecospat_s$range_expansion_NA <- range_dyn_NA$dynamic.index.w['expansion']
                              range_ecospat_s$range_unfilling_NA <- range_dyn_NA$dynamic.index.w['unfilling']
                              
                              # analogue conditions:
                              range_dyn_A <- ecospat.niche.dyn.index(grid.range.hist, 
                                                                     grid.range.rec,
                                                                     intersection = 0) # analysis performed on intersection of ranges of EBBA1 and EBBA2
                              
                              range_ecospat_s$range_stability_A <- range_dyn_A$dynamic.index.w['stability']
                              range_ecospat_s$range_expansion_A <- range_dyn_A$dynamic.index.w['expansion']
                              range_ecospat_s$range_unfilling_A <- range_dyn_A$dynamic.index.w['unfilling']
                              
                              # save plot of range dynamics:
                              jpeg(filename = file.path(plots_dir, paste0(spec, ".jpg")),
                                   quality = 100, height = 920, width = 920)
                              
                              ecospat.plot.niche.dyn(grid.range.hist,
                                                     grid.range.rec,
                                                     quant = 0.05,
                                                     interest = 2, # 1 = historic density, 2 = recent density
                                                     title = paste(spec,":",
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