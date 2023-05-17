# Niche overlap between time periods: 

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

# species-specific background (500 km buffer around presences) or whole EBBA area as background:
bg_spec <- FALSE # set to FALSE for whole EBBA area as environmental background

# paths to data:

#data_dir <- file.path("data", "EBBA_analysis")
data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")

#plots_dir <- "plots"
plots_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts", "plots")

#bioclim_data_dir <- file.path("data", "EBBA_analysis")
bioclim_data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA")

# folder for niche dynamics plots:
plots_dir <- file.path(plots_dir, "niche_dynamics_species", paste0("EBBA_niche_dyn_bg_", ifelse(bg_spec, "spec", "EBBA")))
if(!dir.exists(plots_dir)){dir.create(plots_dir, recursive = TRUE)}

# results table:
results_file <- file.path(data_dir, paste0("EBBA_niche_shift_results_bg_", ifelse(bg_spec, "spec", "EBBA"), ".csv"))


# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# comparable EBBA cells:
EBBA_cells <- read_sf(file.path(data_dir, "EBBA_change.shp")) %>% # output of 1_EBBA_prep_data.R
  select(-species, -Change) %>% 
  distinct(cell50x50, .keep_all = TRUE) # keep geometry; would be the same when using EBBA2

# final species selection:
load(file = file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_final.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
sel_species <- sort(unique(EBBA1_prep$species))

# EBBA 1, only selected species: 
EBBA1_sel_spec_df <- read_sf(file.path(data_dir, "EBBA1_change.shp")) %>% # output of 1_EBBA_prep_data.R
  st_drop_geometry %>% 
  filter(species %in% sel_species)
# add all comparable cells to later extract absences:
EBBA1 <- EBBA_cells %>% 
  full_join(y = EBBA1_sel_spec_df, by = "cell50x50")

# EBBA 2, only selected species:
EBBA2_sel_spec_df <- read_sf(file.path(data_dir, "EBBA2_change.shp")) %>% # output of 1_EBBA_prep_data.R
  st_drop_geometry %>% 
  filter(species %in% sel_species)
# add all comparable cells to later extract absences: # xx not necessary?
EBBA2 <- EBBA_cells %>% 
  full_join(y = EBBA2_sel_spec_df, by = "cell50x50")


# ------------------------------------------------ #
#       Niche overlap between time periods:     ####
# ------------------------------------------------ #

## PCA environment (if environmental background = whole EBBA area for all species) ----

if(bg_spec == FALSE){
  
  # environmental background = all comparable EBBA cells:
  
  # EBBA 1:
  # climate data historic period:
  biovars_hist_rast <- rast(list.files(file.path(bioclim_data_dir, "Bioclim_1981_1990"),
                                       pattern = "50km.tif$",
                                       full.names = TRUE))
  
  # climate data corresponding to EBBA 1 cells:
  EBBA1_climate_data <- EBBA_cells %>% 
    vect %>% # convert to terra object
    cbind(terra::extract(biovars_hist_rast, y = .)) %>%
    as.data.frame %>% 
    select(-ID)
  
  # EBBA 2:
  # climate data recent period:
  biovars_rec_rast <- rast(list.files(file.path(bioclim_data_dir, "Bioclim_2009_2018"),
                                      pattern = "50km.tif$",
                                      full.names = TRUE))
  
  # climate data corresponding to EBBA 2 cells:
  EBBA2_climate_data <- EBBA_cells %>% # same as EBBA2 cells
    vect %>% # convert to terra object
    cbind(terra::extract(biovars_rec_rast, y = .)) %>%
    as.data.frame %>% 
    select(-ID)
  
  # assess climate niche by using the first 2 PCA axes:
  # calibrating the PCA in the whole study area:
  pca.env <- dudi.pca(rbind(EBBA1_climate_data, EBBA2_climate_data)[, paste0("bio", 1:19)],
                      scannf = FALSE,
                      nf = 2) # number of axes
  
  # plot variable contribution to PCs:
  #ecospat.plot.contrib(contrib = pca.env$co, eigen = pca.env$eig)
  
  # predict the scores on the PCA axes:
  # EBBA1 used as z1 (corresponds to native distribution in tutorials)
  # EBBA2 used as z2 (corresponds to invasive distribution in tutorials)
  
  # PCA scores for whole study area (EBBA1 and EBBA2):
  scores.globclim <- pca.env$li 
  # PCA scores for whole EBBA1 area:
  scores.clim.EBBA1 <- suprow(pca.env, EBBA1_climate_data[, paste0("bio", 1:19)])$li
  # PCA scores for whole EBBA2 area:
  scores.clim.EBBA2 <- suprow(pca.env, EBBA2_climate_data[, paste0("bio", 1:19)])$li
  
}


# loop over species:

# register cores for parallel computation:
registerDoParallel(cores = 4)
#getDoParWorkers() # check registered number of cores

niche_shift_df <- foreach(s = 1:length(sel_species),
        .combine = rbind,
        .packages = c("ecospat", "ade4", "dplyr", "sf", "terra"), 
        .verbose = TRUE,
        .errorhandling = "remove") %dopar% {
          
          #print(paste(s,":", sel_species[s]))
  
          # climate data:
          # tifs with bioclim variables for historic and recent time period
          # need to be loaded within foreach since SpatRasters and SpatVectors are non-exportable objects
          biovars_hist_rast <- rast(list.files(file.path(bioclim_data_dir, "Bioclim_1981_1990"),
                                               pattern = "50km.tif$",
                                               full.names = TRUE))
          biovars_rec_rast <- rast(list.files(file.path(bioclim_data_dir, "Bioclim_2009_2018"),
                                              pattern = "50km.tif$",
                                              full.names = TRUE))
          
          # data frame to store results for single species:
          niche_shift_spec_df <- data.frame("species" = sel_species[s],
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
                                            "nS_NA" = NA, # niche stability, non-analogue
                                            "nE_NA" = NA, # niche expansion, non-analogue
                                            "nU_NA" = NA, # niche unfilling, non-analogue
                                            "nS_A" = NA, # niche stability, analogue
                                            "nE_A" = NA, # niche expansion, analogue
                                            "nU_A" = NA, # niche unfilling, analogue
                                            # niche characteristics:
                                            "n_centroid_hist_X" = NA, # useful?
                                            "n_centroid_hist_Y" = NA, # useful?
                                            "n_centroid_rec_X" = NA, # useful?
                                            "n_centroid_rec_Y" = NA, # useful?
                                            "n_centroids_dist" = NA) # useful?

          
          ## (missing) PCA scores: ---------------------------------------------
          
          if(bg_spec){
            
            # environmental background = all comparable EBBA1 cells that are within 500 km of the recorded presences of a species
            # -> PCA scores for species presence and absence points
            
            # EBBA1 data for niche analysis:
            
            # presences:
            EBBA1_pres_sf <- EBBA1 %>%
              filter(species == sel_species[s])
            
            # buffer to determine cells used as absences:
            EBBA1_spec_buffer <- EBBA1_pres_sf %>% 
              st_buffer(dist = 500000) %>% 
              st_union
            
            # background environment, presences and absences:
            EBBA1_env_cells_df <- EBBA_cells %>% 
              st_filter(EBBA1_spec_buffer, .pred = st_within) %>% 
              st_drop_geometry
            
            # data set for niche analysis:
            EBBA1_spec_niche_data <- EBBA_cells %>% 
              mutate(species_occ = ifelse(cell50x50 %in% EBBA1_pres_sf$cell50x50, 1, 
                                          ifelse(cell50x50 %in% EBBA1_env_cells_df$cell50x50, 0, NA))) %>% # species occurrence (1 = presence, 0 = absence):
              filter(!is.na(species_occ)) %>% # remove cells outside of buffer
              vect %>% # convert to terra object
              cbind(terra::extract(biovars_hist_rast, y = .)) %>% # add bioclim data
              as.data.frame %>% 
              select(-ID)
            
            # check:
            #plot(st_geometry(EBBA_cells %>% st_filter(EBBA1_spec_buffer, .pred = st_within)))
            #plot(st_geometry(EBBA1_pres_sf), col = "red", add = TRUE)
            #plot(st_geometry(EBBA1_spec_buffer), add = TRUE)
            
            # EBBA2 data for niche analysis:

            # presences:
            EBBA2_pres_sf <- EBBA2 %>%
              filter(species == sel_species[s])
            
            # buffer to determine cells used as absences:
            EBBA2_spec_buffer <- EBBA2_pres_sf %>% 
              st_buffer(dist = 500000) %>% 
              st_union
            
            # background environment, presences and absences:
            EBBA2_env_cells_df <- EBBA_cells %>% 
              st_filter(EBBA2_spec_buffer, .pred = st_within) %>% 
              st_drop_geometry
            
            # data set for niche analysis:
            EBBA2_spec_niche_data <- EBBA_cells %>% 
              mutate(species_occ = ifelse(cell50x50 %in% EBBA2_pres_sf$cell50x50, 1, 
                                          ifelse(cell50x50 %in% EBBA2_env_cells_df$cell50x50, 0, NA))) %>% # species occurrence (1 = presence, 0 = absence):
              filter(!is.na(species_occ)) %>% # remove cells outside of buffer
              vect %>% # convert to terra object
              cbind(terra::extract(biovars_hist_rast, y = .)) %>% # add bioclim data
              as.data.frame %>% 
              select(-ID)
            
            # check:
            #plot(st_geometry(EBBA_cells %>% st_filter(EBBA2_spec_buffer, .pred = st_within)))
            #plot(st_geometry(EBBA2_pres_sf), col = "red", add = TRUE)
            #plot(st_geometry(EBBA2_spec_buffer), add = TRUE)
            
            # assess climate niche by using the first 2 PCA axes:
            # calibrating the PCA in the whole study area, including both ranges:
            pca.env <- dudi.pca(rbind(EBBA1_spec_niche_data, EBBA2_spec_niche_data)[, paste0("bio", 1:19)],
                                scannf = FALSE,
                                nf = 2) # number of axes
            
            # plot variable contribution to PCs:
            #ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)
            
            # predict the scores on the PCA axes:
            # EBBA1 used as z1 (corresponds to native distribution in tutorials)
            # EBBA2 used as z2 (corresponds to invasive distribution in tutorials)
            
            # PCA scores for whole study area (EBBA1 and EBBA2):
            scores.globclim <- pca.env$li 
            # PCA scores for EBBA1 species distribution:
            scores.sp.EBBA1 <- suprow(pca.env, EBBA1_spec_niche_data[which(EBBA1_spec_niche_data$species_occ == 1), paste0("bio", 1:19)])$li 
            # PCA scores for whole relevant EBBA1 area:
            scores.clim.EBBA1 <- suprow(pca.env, EBBA1_spec_niche_data[, paste0("bio", 1:19)])$li
            # PCA scores for EBBA2 species distribution:
            scores.sp.EBBA2 <- suprow(pca.env, EBBA2_spec_niche_data[which(EBBA2_spec_niche_data$species_occ == 1), paste0("bio", 1:19)])$li 
            # PCA scores for whole relevant EBBA2 area:
            scores.clim.EBBA2 <- suprow(pca.env, EBBA2_spec_niche_data[, paste0("bio", 1:19)])$li
          
            
            } else {
            
            # environmental background = whole EBBA area for all species
            # -> PCA scores of species presence points
            # (PCA scores for whole study area and EBBA1 and EBBA2 area seperately calculated above) 
            
            # EBBA 1 presences and matching climate data:
            EBBA1_spec_niche_data <- EBBA1 %>%
              filter(species == sel_species[s]) %>% 
              vect %>% # convert to terra object
              cbind(terra::extract(biovars_hist_rast, y = .)) %>%
              as.data.frame %>% 
              select(-ID)
            
            # EBBA 2 presences and matching climate data:
            EBBA2_spec_niche_data <- EBBA2 %>%
              filter(species == sel_species[s]) %>% 
              vect %>% # convert to terra object
              cbind(terra::extract(biovars_rec_rast, y = .)) %>%
              as.data.frame %>% 
              select(-ID)
            
            # PCA scores for EBBA1 species distribution:
            scores.sp.EBBA1 <- suprow(pca.env, EBBA1_spec_niche_data[, paste0("bio", 1:19)])$li 
            
            # PCA scores for EBBA2 species distribution:
            scores.sp.EBBA2 <- suprow(pca.env, EBBA2_spec_niche_data[, paste0("bio", 1:19)])$li 
            
          }
          
          
          ## occurrence densities grids for EBBA1 and EBBA2 distribution: ------
          
          # gridding EBBA1 niche:
          grid.clim.EBBA1 <- ecospat.grid.clim.dyn(glob = scores.globclim, # whole study area, both time periods
                                                   glob1 = scores.clim.EBBA1, # env. background
                                                   sp = scores.sp.EBBA1, # species occurrence
                                                   R = 100, # grid resolution
                                                   th.sp = 0)
          # gridding EBBA2 niche:
          grid.clim.EBBA2 <- ecospat.grid.clim.dyn(glob = scores.globclim, # whole study area, both time periods
                                                   glob1 = scores.clim.EBBA2, # env. background
                                                   sp = scores.sp.EBBA2, # species occurrence
                                                   R = 100, 
                                                   th.sp = 0)
          
          # niche overlap -  Schoener's D:
          niche_shift_spec_df$D <- ecospat.niche.overlap(grid.clim.EBBA1, grid.clim.EBBA2, cor = TRUE)$D # occurrence densities corrected with prevalence of the environments
        
          ## niche similarity tests: -------------------------------------------
          
          ### tests for niche shifting: ----
          
          # non-analogue conditions:
          sim_test_shift_NA <- ecospat.niche.similarity.test(grid.clim.EBBA1,
                                                             grid.clim.EBBA2,
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
          sim_test_shift_A <- ecospat.niche.similarity.test(grid.clim.EBBA1,
                                                            grid.clim.EBBA2,
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
          sim_test_cons_NA <- ecospat.niche.similarity.test(grid.clim.EBBA1, 
                                                            grid.clim.EBBA2, 
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
          sim_test_cons_A <- ecospat.niche.similarity.test(grid.clim.EBBA1,
                                                           grid.clim.EBBA2, 
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
          
          # between EBBA1 (as z1) and EBBA2 (as z2)
          
          # non-analogue conditions:
          EBBA12_niche_dyn_NA <- ecospat.niche.dyn.index(grid.clim.EBBA1, 
                                                         grid.clim.EBBA2,
                                                         intersection = NA) # analysis performed on niche space of EBBA1 and EBBA2
          
          niche_shift_spec_df$nS_NA <- EBBA12_niche_dyn_NA$dynamic.index.w['stability']
          niche_shift_spec_df$nE_NA <- EBBA12_niche_dyn_NA$dynamic.index.w['expansion']
          niche_shift_spec_df$nU_NA <- EBBA12_niche_dyn_NA$dynamic.index.w['unfilling']
        
          # analogue conditions:
          EBBA12_niche_dyn_A <- ecospat.niche.dyn.index(grid.clim.EBBA1, 
                                                        grid.clim.EBBA2,
                                                        intersection = 0) # analysis performed on intersection of niche space of EBBA1 and EBBA2
          
          niche_shift_spec_df$nS_A <- EBBA12_niche_dyn_A$dynamic.index.w['stability']
          niche_shift_spec_df$nE_A <- EBBA12_niche_dyn_A$dynamic.index.w['expansion']
          niche_shift_spec_df$nU_A <- EBBA12_niche_dyn_A$dynamic.index.w['unfilling']
          
          
          # save plot of niche dynamics:
          jpeg(filename = file.path(plots_dir, paste0(sel_species[s], ".jpg")),
               quality = 100, height = 920, width = 920)
          
          ecospat.plot.niche.dyn(grid.clim.EBBA1,
                                 grid.clim.EBBA2,
                                 quant = 0.05,
                                 interest = 2, # 1 = EBBA1 density, 2 = EBBA2 density
                                 title = paste(sel_species[s],":",
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
          
          # as median of species EBBA distribution scores along the first 2 PCA axes:
          niche_shift_spec_df$n_centroid_hist_X <- median(scores.sp.EBBA1[,1]) 
          niche_shift_spec_df$n_centroid_hist_Y <- median(scores.sp.EBBA1[,2])
          niche_shift_spec_df$n_centroid_rec_X <- median(scores.sp.EBBA2[,1])
          niche_shift_spec_df$n_centroid_rec_Y <- median(scores.sp.EBBA2[,2])
          
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