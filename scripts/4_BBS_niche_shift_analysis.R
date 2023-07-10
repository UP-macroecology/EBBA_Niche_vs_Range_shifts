# Niche overlap between time periods for BBS data:

# note: 
# - ecospat dependencies ‘biomod2’, ‘randomForest’ depend on R >= 4.1.0
# -> before running the script on the cluster run "module load R/4.1.0-foss-2021a" on ecoc9 to ensure that R version 4.1.0 (2021-05-18) is used

# outputs:
# for each species:
# - niche overlap (Schoener's D)
# - results of niche shift tests
# - results of niche conservatism tests
# - niche dynamic metrics: niche stability, expansion and unfilling between historic and recent time period
# - plots of niche dynamics
# - standardised niche dynamic metrics: stability, expansion, unfilling, abandonment, pioneering
# - niche sizes, centroids and distance between centroids?

library(dplyr)
library(sf)
library(terra)
library(doParallel)
library(ecospat)
library(ade4)


# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# environmental background: presences and absences within 600 km buffer around presences (TRUE) or all true absences within conterminous US (FALSE):
bg_spec <- TRUE

# which historic time period should be used:
hist_years <- 1980:1983 # maximum gap between historic and recent time period
#hist_years <- 1987:1990 # similar gap between historic and recent time period as in EBBA analysis

# paths to data:

# project data:
#data_dir <- file.path("data", "BBS_analysis")
data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")

# bioclimatic variables:
#bioclim_data_dir <- file.path("data", "BBS_analysis)
bioclim_data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA")

#plots_dir <- "plots"
plots_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts", "plots")

# folder for niche dynamics plots:
plots_dir <- file.path(plots_dir, "niche_dynamics_species", paste0("BBS_niche_dyn_bg_", ifelse(bg_spec, "spec", "US"), "_hist", 
                                                                   ifelse(all(hist_years == 1980:1983), "81-83", "88-90")))
if(!dir.exists(plots_dir)){dir.create(plots_dir, recursive = TRUE)}

# results table:
results_file <- file.path(data_dir, paste0("BBS_niche_shift_results_bg_", ifelse(bg_spec, "spec", "US"),  "_hist",
                                           ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), "_070723.csv"))


# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# species selection: 
sel_species <- read.csv(file = file.path(data_dir, "BBS_stability_PCA_contUS_BL22_060723.csv")) %>%
  filter(stability >= 0.5) %>% 
  pull(species) %>% 
  sort # 237

# BBS data, only selected species:
load(file = file.path(data_dir, paste0("BBS_prep_steps1-4_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".RData"))) # output of 2_1_BBS_species_filtering_1-4.R
species_filtered <- sort(unique(hist_prep_df$species))
sel_species_final <- species_filtered[which(species_filtered %in% sel_species)] # 196 (hist81-83) / 236 (hist96-98)

BBS_hist <- read_sf(file.path(data_dir, paste0("BBS_historic_centr_proj_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".shp"))) %>% # output of 1_BBS_prep_data.R
  filter(species %in% sel_species_final)

BBS_rec <- read_sf(file.path(data_dir, paste0("BBS_recent_centr_proj_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".shp"))) %>% # output of 1_BBS_prep_data.R
  filter(species %in% sel_species_final)


# ------------------------------------------------ #
#       Niche overlap between time periods:     ####
# ------------------------------------------------ #

# loop over species:

# register cores for parallel computation:
registerDoParallel(cores = 15)
#getDoParWorkers() # check registered number of cores

niche_shift_df <- foreach(s = 1:length(sel_species_final),
                          .combine = rbind,
                          .packages = c("ecospat", "ade4", "dplyr", "sf", "terra"), 
                          .verbose = TRUE,
                          .errorhandling = "remove",
                          .inorder = FALSE) %dopar% {
                            
                            spec <- sel_species_final[s]

                            # climate data:
                            # tifs with bioclim variables for historic and recent time period
                            # need to be loaded within foreach since SpatRasters and SpatVectors are non-exportable objects
                            
                            biovars_hist_rast <- rast(list.files(file.path(bioclim_data_dir, paste0("Bioclim_", min(hist_years), "_", max(hist_years))),
                                                                 pattern = "1km.tif$",
                                                                 full.names = TRUE))
                            
                            biovars_rec_rast <- rast(list.files(file.path(bioclim_data_dir, "Bioclim_2015_2018"),
                                                                pattern = "1km.tif$",
                                                                full.names = TRUE))
                            
                            # data frame to store results for single species:
                            niche_shift_spec_df <- data.frame("species" = spec,
                                                              # niche overlap:  
                                                              "D" = NA, # Schoener's D
                                                              # results of niche similarity test when testing for niche divergence / niche shift:
                                                              "shift_p_D_NA" = NA, # NA = non-analogue conditions
                                                              "shift_p_E_NA" = NA,
                                                              "shift_p_S_NA" = NA,
                                                              "shift_p_U_NA" = NA,
                                                              "shift_p_D_A" = NA, # A = analogue conditions
                                                              "shift_p_E_A" = NA,
                                                              "shift_p_S_A" = NA,
                                                              "shift_p_U_A" = NA,
                                                              # results of niche similarity test when testing for niche conservatism:
                                                              "cons_p_D_NA" = NA,
                                                              "cons_p_E_NA" = NA,
                                                              "cons_p_S_NA" = NA,
                                                              "cons_p_U_NA" = NA,
                                                              "cons_p_D_A" = NA,
                                                              "cons_p_E_A" = NA,
                                                              "cons_p_S_A" = NA,
                                                              "cons_p_U_A" = NA,
                                                              # niche dynamic metrics:
                                                              "nS_NA" = NA,
                                                              "nE_NA" = NA,
                                                              "nU_NA" = NA,
                                                              "nS_A" = NA,
                                                              "nE_A" = NA,
                                                              "nU_A" = NA,
                                                              # niche characteristics:
                                                              "n_centroid_hist_X" = NA, # useful?
                                                              "n_centroid_hist_Y" = NA, # useful?
                                                              "n_centroid_rec_X" = NA, # useful?
                                                              "n_centroid_rec_Y" = NA, # useful?
                                                              "n_centroids_dist" = NA) # useful?
                            
                            
                            ## PCA scores: ---------------------------------------------
                            
                            # historic time period:
                            
                            # species records (presences and true absences):
                            spec_hist <- BBS_hist %>%
                              filter(species == spec)

                            if(bg_spec){
                              
                              # environmental background = all presences and true absences within 500 km of presences (Sofaer et al. 2018)
                              # -> PCA scores for species presence and absence points
                              
                              # 500 km buffer around presences:
                              background_bf <- BBS_hist %>%
                                filter(species == spec & pres == 1) %>% 
                                st_buffer(dist = 500000) %>% 
                                st_union
                              
                              # species records (presences and true absences) within 500 km of presences:
                              background <- spec_hist %>% 
                                st_filter(y = background_bf, .predicate = st_within)
                              
                              # check:
                              #plot(st_geometry(background_bf))
                              #plot(st_geometry(spec_hist), add = TRUE, col = spec_hist$pres+2)
                              #plot(st_geometry(background), add = TRUE, col = "blue")
                            
                              } else {
                                # environmental background = all true absences of a species across conterminous US
                                background <- spec_hist
                            }
                            
                            # data set for niche analysis:
                            BBS_hist_niche_data <- background %>% 
                              st_buffer(dist = 21000) %>% # 21 km as in Sofaer et al. 2018 
                              # add bioclim data:
                              cbind(terra::extract(biovars_hist_rast, 
                                                   y = .,
                                                   fun = mean,
                                                   na.rm = TRUE,
                                                   touches = TRUE, # all cells in the buffer area considered
                                                   weights = TRUE, # weighted mean: fraction of each cell in the buffer is used as weight to calculate mean)
                                                   exact = TRUE)) %>%
                              dplyr::select(-c(ID, AOU, start_X, start_Y)) %>% 
                              st_drop_geometry
                              

                            # recent time period:
                            
                            # species records (presences and true absences):
                            spec_rec <- BBS_rec %>%
                              filter(species == spec)
                            
                            if(bg_spec){
                              
                              # environmental background = all presences and true absences within 500 km of presences (Sofaer et al. 2018)
                              # -> PCA scores for species presence and absence points
                              
                              # 500 km buffer around presences:
                              background_bf <- BBS_rec %>%
                                filter(species == spec & pres == 1) %>% 
                                st_buffer(dist = 500000) %>% 
                                st_union
                              
                              # species records (presences and true absences) within 500 km of presences:
                              background <- spec_rec %>% 
                                st_filter(y = background_bf, .predicate = st_within)
                              
                              # check:
                              #plot(st_geometry(background_bf))
                              #plot(st_geometry(spec_hist), add = TRUE)
                              #plot(st_geometry(background), add = TRUE, col = "green")
                              
                            } else {
                              # environmental background = all true absences of a species across conterminous US
                              background <- spec_rec
                            }
                            
                            # data set for niche analysis:
                            BBS_rec_niche_data <- background %>% 
                              st_buffer(dist = 21000) %>% # 21 km as in Sofaer et al. 2018 
                              # add bioclim data:
                              cbind(terra::extract(biovars_rec_rast, 
                                                   y = .,
                                                   fun = mean,
                                                   na.rm = TRUE,
                                                   touches = TRUE, # all cells in the buffer area considered
                                                   weights = TRUE, # weighted mean: fraction of each cell in the buffer is used as weight to calculate mean)
                                                   exact = TRUE)) %>%
                              dplyr::select(-c(ID, AOU, start_X, start_Y)) %>% 
                              st_drop_geometry
                              
                              
                            # assess climate niche by using the first 2 PCA axes:
                            # calibrating the PCA in the whole study area, including both ranges:
                            pca.env <- dudi.pca(rbind(BBS_hist_niche_data, BBS_rec_niche_data)[, paste0("bio", 1:19)],
                                                scannf = FALSE,
                                                nf = 2) # number of axes
                            
                            # plot variable contribution to PCs:
                            #ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
                            
                            # predict the scores on the PCA axes:
                            # historic used as z1 (corresponds to native distribution in tutorials)
                            # recent used as z2 (corresponds to invasive distribution in tutorials)

                            # PCA scores for whole study area (historic and recent):
                            scores.globclim <- pca.env$li 
                            # PCA scores for historic species distribution:
                            scores.sp.hist <- suprow(pca.env, BBS_hist_niche_data[which(BBS_hist_niche_data$pres == 1), paste0("bio", 1:19)])$li 
                            # PCA scores for whole relevant area:
                            scores.clim.hist <- suprow(pca.env, BBS_hist_niche_data[, paste0("bio", 1:19)])$li
                            # PCA scores for recent species distribution:
                            scores.sp.rec <- suprow(pca.env, BBS_rec_niche_data[which(BBS_rec_niche_data$pres == 1), paste0("bio", 1:19)])$li 
                            # PCA scores for whole relevant area:
                            scores.clim.rec <- suprow(pca.env, BBS_rec_niche_data[, paste0("bio", 1:19)])$li
                            

                            ## occurrence densities grids for historic and recent distribution: ------

                            # gridding historic niche:
                            grid.clim.hist <- ecospat.grid.clim.dyn(glob = scores.globclim, # whole study area, both time periods
                                                                    glob1 = scores.clim.hist, # env. background
                                                                    sp = scores.sp.hist, # species occurrence
                                                                    R = 100, # grid resolution
                                                                    th.sp = 0)
                            # gridding recent niche:
                            grid.clim.rec <- ecospat.grid.clim.dyn(glob = scores.globclim, # whole study area, both time periods
                                                                   glob1 = scores.clim.rec, # env. background
                                                                   sp = scores.sp.rec, # species occurrence
                                                                   R = 100, 
                                                                   th.sp = 0)
                            
                            # niche overlap -  Schoener's D:
                            niche_shift_spec_df$D <- ecospat.niche.overlap(grid.clim.hist, grid.clim.rec, cor = TRUE)$D # occurrence densities corrected with prevalence of the environments
                            
                            
                            ## niche similarity tests: -------------------------------------------
                            
                            ### tests for niche shifting: ----
                            
                            # non-analogue conditions:
                            sim_test_shift_NA <- ecospat.niche.similarity.test(grid.clim.hist,
                                                                               grid.clim.rec,
                                                                               rep = 1000,
                                                                               overlap.alternative = "lower", # test for niche shifting / divergence: niche overlap is lower than random
                                                                               expansion.alternative = "higher", # test for niche shifting: expansion higher than random
                                                                               stability.alternative = "lower", # test for niche shifting: stability lower than random
                                                                               unfilling.alternative = "higher", # test for niche shifting: unfilling higher than random
                                                                               intersection = NA, 
                                                                               rand.type = 2)
                            
                            niche_shift_spec_df$shift_p_D_NA <- sim_test_shift_NA$p.D
                            niche_shift_spec_df$shift_p_E_NA <- sim_test_shift_NA$p.expansion
                            niche_shift_spec_df$shift_p_S_NA <- sim_test_shift_NA$p.stability
                            niche_shift_spec_df$shift_p_U_NA <- sim_test_shift_NA$p.unfilling
                            
                            # plot results of similarity test:
                            #ecospat.plot.overlap.test(sim_test_shift_NA, "D", "Similarity")
                            
                            # analogue conditions:
                            sim_test_shift_A <- ecospat.niche.similarity.test(grid.clim.hist,
                                                                              grid.clim.rec,
                                                                              rep = 1000,
                                                                              overlap.alternative = "lower", # test for niche shifting / divergence: niche overlap is lower than random
                                                                              expansion.alternative = "higher", # test for niche shifting: expansion higher than random
                                                                              stability.alternative = "lower", # test for niche shifting: stability lower than random
                                                                              unfilling.alternative = "higher", # test for niche shifting: unfilling higher than random
                                                                              intersection = 0, 
                                                                              rand.type = 2)
                            
                            niche_shift_spec_df$shift_p_D_A <- sim_test_shift_A$p.D
                            niche_shift_spec_df$shift_p_E_A <- sim_test_shift_A$p.expansion
                            niche_shift_spec_df$shift_p_S_A <- sim_test_shift_A$p.stability
                            niche_shift_spec_df$shift_p_U_A <- sim_test_shift_A$p.unfilling
                            
                            
                            ### tests for niche conservatism: ----
                            
                            # non-analogue conditions:
                            sim_test_cons_NA <- ecospat.niche.similarity.test(grid.clim.hist, 
                                                                              grid.clim.rec, 
                                                                              rep = 1000, 
                                                                              overlap.alternative = "higher", # test for niche conservatism: niche overlap is higher than random
                                                                              expansion.alternative = "lower", # test for niche conservatism: expansion lower than random
                                                                              stability.alternative = "higher", # test for niche conservatism: stability higher than random
                                                                              unfilling.alternative = "lower", # test for niche conservatism: unfilling lower than random
                                                                              intersection = NA, 
                                                                              rand.type = 2)
                            
                            niche_shift_spec_df$cons_p_D_NA <- sim_test_cons_NA$p.D
                            niche_shift_spec_df$cons_p_E_NA <- sim_test_cons_NA$p.expansion
                            niche_shift_spec_df$cons_p_S_NA <- sim_test_cons_NA$p.stability
                            niche_shift_spec_df$cons_p_U_NA <- sim_test_cons_NA$p.unfilling
                            
                            # plot results of similarity test:
                            #ecospat.plot.overlap.test(sim_test_cons_NA, "D", "Similarity")
                            
                            # analogue conditions:
                            sim_test_cons_A <- ecospat.niche.similarity.test(grid.clim.hist,
                                                                             grid.clim.rec, 
                                                                             rep = 1000, 
                                                                             overlap.alternative = "higher", # test for niche conservatism: niche overlap is higher than random
                                                                             expansion.alternative = "lower", # test for niche conservatism: expansion lower than random
                                                                             stability.alternative = "higher", # test for niche conservatism: stability higher than random
                                                                             unfilling.alternative = "lower", # test for niche conservatism: unfilling lower than random
                                                                             intersection = 0, 
                                                                             rand.type = 2)
                            
                            niche_shift_spec_df$cons_p_D_A <- sim_test_cons_A$p.D
                            niche_shift_spec_df$cons_p_E_A <- sim_test_cons_A$p.expansion
                            niche_shift_spec_df$cons_p_S_A <- sim_test_cons_A$p.stability
                            niche_shift_spec_df$cons_p_U_A <- sim_test_cons_A$p.unfilling
                            
                            
                            ## indices for niche dynamics: ---------------------------------------
                            
                            # between historic (as z1) and recent period (as z2)
                            
                            # non-analogue conditions:
                            niche_dyn_NA <- ecospat.niche.dyn.index(grid.clim.hist, 
                                                                    grid.clim.rec,
                                                                    intersection = NA) # analysis performed on niche space of histoirc and recent period
                            
                            niche_shift_spec_df$nS_NA <- niche_dyn_NA$dynamic.index.w['stability']
                            niche_shift_spec_df$nE_NA <- niche_dyn_NA$dynamic.index.w['expansion']
                            niche_shift_spec_df$nU_NA <- niche_dyn_NA$dynamic.index.w['unfilling']
                            
                            # analogue conditions:
                            niche_dyn_A <- ecospat.niche.dyn.index(grid.clim.hist, 
                                                                   grid.clim.rec,
                                                                   intersection = 0) # analysis performed on intersection of niche space of BBS_hist and BBS_rec
                            
                            niche_shift_spec_df$nS_A <- niche_dyn_A$dynamic.index.w['stability']
                            niche_shift_spec_df$nE_A <- niche_dyn_A$dynamic.index.w['expansion']
                            niche_shift_spec_df$nU_A <- niche_dyn_A$dynamic.index.w['unfilling']
                            
                            
                            # save plot of niche dynamics:
                            jpeg(filename = file.path(plots_dir, paste0(spec, ".jpg")),
                                 quality = 100, height = 920, width = 920)
                            
                            ecospat.plot.niche.dyn(grid.clim.hist,
                                                   grid.clim.rec,
                                                   quant = 0.05,
                                                   interest = 2, # 1 = historic density, 2 = recent density
                                                   title = paste(spec,":",
                                                                 "stability:", round(niche_shift_spec_df$nS_NA,3),
                                                                 "expansion:", round(niche_shift_spec_df$nE_NA,3),
                                                                 "unfiling:", round(niche_shift_spec_df$nU_NA,3)),
                                                   name.axis1 = "PC1",
                                                   name.axis2 = "PC2")
                            dev.off()
                            # blue = stability
                            # red = expansion
                            # green = unfilling
                            
                            
                            ## niche characteristics: --------------------------------------------
                            
                            # niche centroids:
                            
                            # as median of distribution scores along the first 2 PCA axes:
                            niche_shift_spec_df$n_centroid_hist_X <- median(scores.sp.hist[,1]) 
                            niche_shift_spec_df$n_centroid_hist_Y <- median(scores.sp.hist[,2])
                            niche_shift_spec_df$n_centroid_rec_X <- median(scores.sp.rec[,1])
                            niche_shift_spec_df$n_centroid_rec_Y <- median(scores.sp.rec[,2])
                            
                            # niche shift:
                            niche_shift_spec_df$n_centroids_dist <- sqrt(((niche_shift_spec_df$n_centroid_rec_X - niche_shift_spec_df$n_centroid_hist_X)^2) + 
                                                                               ((niche_shift_spec_df$n_centroid_rec_Y - niche_shift_spec_df$n_centroid_hist_Y)^2))
                            
                            niche_shift_spec_df
                          }


## standardise niche dynamic metrics: ------------------------------------------

# non-analogue climatic conditions:
e <- niche_shift_df$nE_NA # expansion related to whole niche space (EBBA1 and EBBA2) 
s <- niche_shift_df$nS_NA # stability related to whole niche space (EBBA1 and EBBA2)
u <- niche_shift_df$nU_NA # unfilling related to whole niche space (EBBA1 and EBBA2)

# stability, non-analogue conditions, related to EBBA1 niche:
s_x <- 1 - u

# proportion stability related to whole niche space and stability related to EBBA1 niche (= proportion by which one niche is larger / smaller than the other)
x <- s / s_x

# total niche space (non-analogue):
t <- u * x + s + e # unfilling rescaled to relate to EBBA2 niche space

# unfilling, stability and expansion standardized so that they sum up to one: towards total non-analogue niche
u_rel <- u * x / t
s_rel <- s / t
e_rel <- e / t

# analogue climatic conditions:
e_a <- niche_shift_df$nE_A # expansion related to analogue EBBA2 niche space
s_a <- niche_shift_df$nS_A # stability related to analogue EBBA2 niche space
u_a <- niche_shift_df$nU_A # unfilling related to analogue EBBA1 niche space

# stability, analogue conditions, related to EBBA1 niche space:
s_y <- 1 - u_a

# proportion stability related to EBBA2 analogue niche space and stability related to EBBA2 niche space (= proportion by which one analogue niche space is larger / smaller than the other)
y <- s_a / s_y

# total niche space (analogue)
t_a <- u_a * y + s_a + e_a 

# unfilling, stability and expansion standardized so that they sum up to one, relativ to total analogue niche
s_rel_a1 <- s_a / t_a
e_rel_a1 <- e_a / t_a
u_rel_a1 <- u_a * y / t_a

# size of t_a in relation to t:
i <- s_rel / s_rel_a1 # proportion non-analogue to analogue stability 

# target values:
niche_shift_df$niche_stability_std <- s_rel_a1 * i # stability
niche_shift_df$niche_unfilling_std <- u_rel_a1 * i # unfilling
niche_shift_df$niche_expansion_std <- e_rel_a1 * i # expansion
niche_shift_df$niche_abandonment_std <- u_rel - niche_shift_df$niche_unfilling_std # abandonment
niche_shift_df$niche_abandonment_std[niche_shift_df$niche_abandonment_std < 0] <- 0 # slightly negative values result from rounding errors of ecospat package
niche_shift_df$niche_pioneering_std <- e_rel - niche_shift_df$niche_expansion_std # pioneering
niche_shift_df$niche_pioneering_std[niche_shift_df$niche_pioneering_std < 0] <- 0 # slightly negative values result from rounding errors of ecospat package

# check: values should sum to 1:
total <- niche_shift_df$niche_abandonment_std + 
  niche_shift_df$niche_unfilling_std + 
  niche_shift_df$niche_stability_std + 
  niche_shift_df$niche_expansion_std + 
  niche_shift_df$niche_pioneering_std

# save resulting data frame:
write.csv(niche_shift_df, file = results_file, row.names = FALSE)
