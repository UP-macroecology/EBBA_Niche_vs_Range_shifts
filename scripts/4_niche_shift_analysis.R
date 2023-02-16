# Niche overlap between time periods: 

# ecospat dependencies ‘biomod2’, ‘randomForest’ depend on R >= 4.1.0
# -> before running the script run "module load R/4.1.0-foss-2021a" on ecoc9 to ensure that R version 4.1.0 (2021-05-18) is used

library(dplyr)
library(sf)
library(terra)
library(doParallel)
library(ecospat)
library(ade4)

# directories:

#EBBA_data_dir <- file.path("Data")
EBBA_data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")

#bioclim_data_dir <- file.path("Data")
bioclim_data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA")

# load data:

# comparable EBBA cells:
EBBA1_cells <- read_sf(file.path(EBBA_data_dir, "EBBA1_comparable.shp")) %>% # output of 1_prep_EBBA_data.R
  select(-species) %>% 
  distinct(cell50x50, .keep_all = TRUE) # keep geometry; would be the same when using EBBA2

# species selection: 
sel_species <- read.csv(file.path(EBBA_data_dir, "EBBA_niche_range_shifts_species_selection.csv"))$species %>% # output of 2_3_species_filtering_5_climatic_niche_analysis.R
  sort # alphabetically sorted (to avoid confusion with indices later)

# EBBA 1, only selected species: 
EBBA1_sel_spec_df <- read_sf(file.path(EBBA_data_dir, "EBBA1_comparable_harmonized.shp")) %>% # output of 1_prep_EBBA_data.R
  st_drop_geometry %>% 
  filter(species %in% sel_species)
# add all comparable cells to later extract absences:
EBBA1 <- EBBA1_cells %>% 
  full_join(y = EBBA1_sel_spec_df, by = "cell50x50")

# EBBA 2, only selected species:
EBBA2_sel_spec_df <- read_sf(file.path(EBBA_data_dir, "EBBA2_comparable_harmonized.shp")) %>% # output of 1_prep_EBBA_data.R
  st_drop_geometry %>% 
  filter(species %in% sel_species)
# add all comparable cells to later extract absences:
EBBA2 <- EBBA1_cells %>% 
  full_join(y = EBBA2_sel_spec_df, by = "cell50x50")

# loop over species:

# register cores for parallel computation:
registerDoParallel(cores = 8)
#getDoParWorkers() # check registered number of cores

niche_shift_df <- foreach(s = 1:length(sel_species),
        .combine = rbind,
        .packages = c("ecospat", "ade4", "dplyr", "sf", "terra"), 
        .verbose = TRUE) %dopar% {
          
          options(warn = 1) # default warn = 0; 1 = warnings are printed as they occur
          
          print(paste(s,":", sel_species[s]))
  
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
                                       "D" = NA,
                                       "shift_simtest_NA_p_D" = NA, # NA = non-analogue conditions
                                       "shift_simtest_NA_p_expansion" = NA,
                                       "shift_simtest_NA_p_stability" = NA,
                                       "shift_simtest_NA_p_unfilling" = NA,
                                       "shift_simtest_A_p_D" = NA, # A = analogue conditions
                                       "shift_simtest_A_p_expansion" = NA,
                                       "shift_simtest_A_p_stability" = NA,
                                       "shift_simtest_A_p_unfilling" = NA,
                                       "cons_simtest_NA_p_D" = NA,
                                       "cons_simtest_NA_p_expansion" = NA,
                                       "cons_simtest_NA_p_stability" = NA,
                                       "cons_simtest_NA_p_unfilling" = NA,
                                       "cons_simtest_A_p_D" = NA,
                                       "cons_simtest_A_p_expansion" = NA,
                                       "cons_simtest_A_p_stability" = NA,
                                       "cons_simtest_A_p_unfilling" = NA,
                                       "niche_stability_NA" = NA,
                                       "niche_expansion_NA" = NA,
                                       "niche_unfilling_NA" = NA,
                                       "niche_stability_A" = NA,
                                       "niche_expansion_A" = NA,
                                       "niche_unfilling_A" = NA,
                                       "niche_centroid_hist_X" = NA, # not sure whether this is useful xx
                                       "niche_centroid_hist_Y" = NA, # not sure whether this is useful xx
                                       "niche_centroid_rec_X" = NA, # not sure whether this is useful xx
                                       "niche_centroid_rec_Y" = NA, # not sure whether this is useful xx
                                       "niche_centroids_dist" = NA) # direction with regard to how env. variables load on PCs? xx
          # information on niche size? xx
          
          # EBBA1 data for niche analysis:
          # environment = all comparable EBBA1 cells that are within 500 km of the recorded presences of a species
          
          # presences:
          EBBA1_pres_sf <- EBBA1 %>%
            filter(species == sel_species[s])
          
          # buffer to determine cells used as absences:
          EBBA1_spec_buffer <- EBBA1_pres_sf %>% 
            st_buffer(dist = 500000) %>% 
            st_union
          
          # background environment, presences and absences:
          EBBA1_env_cells_df <- EBBA1_cells %>% 
              st_filter(EBBA1_spec_buffer, .pred = st_within) %>% 
              st_drop_geometry
        
          # data set for niche analysis:
          EBBA1_spec_niche_data <- EBBA1_cells %>% 
            mutate(species_occ = ifelse(cell50x50 %in% EBBA1_pres_sf$cell50x50, 1, 
                                        ifelse(cell50x50 %in% EBBA1_env_cells_df$cell50x50, 0, NA))) %>% # species occurrence (1 = presence, 0 = absence):
            filter(!is.na(species_occ)) %>% # remove cells outside of buffer
            vect %>% # convert to terra object
            cbind(terra::extract(biovars_hist_rast, y = .)) %>% # add bioclim data
            as.data.frame %>% 
            select(-ID)
        
          # check:
          #plot(st_geometry(EBBA1_cells %>% st_filter(EBBA1_spec_buffer, .pred = st_within)))
          #plot(st_geometry(EBBA1_pres_sf), col = "red", add = TRUE)
          #plot(st_geometry(EBBA1_spec_buffer), add = TRUE)
          
          
          # EBBA2 data for niche analysis:
          # environment = all comparable EBBA2 cells that are within 500 km of the recorded presences of a species
          
          # presences:
          EBBA2_pres_sf <- EBBA2 %>%
            filter(species == sel_species[s])
          
          # buffer to determine cells used as absences:
          EBBA2_spec_buffer <- EBBA2_pres_sf %>% 
            st_buffer(dist = 500000) %>% 
            st_union
          
          # background environment, presences and absences:
          EBBA2_env_cells_df <- EBBA1_cells %>% 
            st_filter(EBBA2_spec_buffer, .pred = st_within) %>% 
            st_drop_geometry
          
          # data set for niche analysis:
          EBBA2_spec_niche_data <- EBBA1_cells %>% 
            mutate(species_occ = ifelse(cell50x50 %in% EBBA2_pres_sf$cell50x50, 1, 
                                        ifelse(cell50x50 %in% EBBA2_env_cells_df$cell50x50, 0, NA))) %>% # species occurrence (1 = presence, 0 = absence):
            filter(!is.na(species_occ)) %>% # remove cells outside of buffer
            vect %>% # convert to terra object
            cbind(terra::extract(biovars_hist_rast, y = .)) %>% # add bioclim data
            as.data.frame %>% 
            select(-ID)
          
          # check:
          #plot(st_geometry(EBBA1_cells %>% st_filter(EBBA2_spec_buffer, .pred = st_within)))
          #plot(st_geometry(EBBA2_pres_sf), col = "red", add = TRUE)
          #plot(st_geometry(EBBA2_spec_buffer), add = TRUE)
          
          
          # niche analysis:
          
          # assess climate niche by using the first 2 PCA axes:
          # calibrating the PCA in the whole study area, including both ranges:
          pca.env <- dudi.pca(rbind(EBBA1_spec_niche_data, EBBA2_spec_niche_data)[,3:21],
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
          scores.sp.EBBA1 <- suprow(pca.env, EBBA1_spec_niche_data[which(EBBA1_spec_niche_data$species_occ == 1),3:21])$li 
          # PCA scores for whole relevant EBBA1 area:
          scores.clim.EBBA1 <- suprow(pca.env, EBBA1_spec_niche_data[,3:21])$li
          # PCA scores for EBBA2 species distribution:
          scores.sp.EBBA2 <- suprow(pca.env, EBBA2_spec_niche_data[which(EBBA2_spec_niche_data$species_occ == 1),3:21])$li 
          # PCA scores for whole relevant EBBA2 area:
          scores.clim.EBBA2 <- suprow(pca.env, EBBA2_spec_niche_data[,3:21])$li
          
          # calculate the Occurrence Densities Grid for EBBA1 and EBBA2 distribution:
          
          # gridding EBBA1 niche:
          grid.clim.EBBA1 <- ecospat.grid.clim.dyn(glob = scores.globclim, # whole study area, both time periods
                                                   glob1 = scores.clim.EBBA1, # env. background
                                                   sp = scores.sp.EBBA1, # species occurrence
                                                   R = 100, # grid resolution
                                                   th.sp = 0)
          # gridding EBBA1 niche:
          grid.clim.EBBA2 <- ecospat.grid.clim.dyn(glob = scores.globclim, # whole study area, both time periods
                                                   glob1 = scores.clim.EBBA2, # env. background
                                                   sp = scores.sp.EBBA2, # species occurrence
                                                   R = 100, 
                                                   th.sp = 0)
          
          # niche overlap: Schoener's D:
          niche_shift_spec_df$D <- ecospat.niche.overlap (grid.clim.EBBA1, grid.clim.EBBA2, cor = TRUE)$D # correct occurrence densities by prevalence of the environments
          
          # niche similarity test - test for niche shifting:
          
          # non-analogue conditions:
          sim_test_niche_shift_NA <- ecospat.niche.similarity.test(grid.clim.EBBA1,
                                                                grid.clim.EBBA2,
                                                                rep = 1000,
                                                                overlap.alternative = "lower", # test for niche shifting / divergence: niche overlap is lower than random
                                                                expansion.alternative = "higher", # test for niche shifting: expansion higher than random
                                                                stability.alternative = "lower", # test for niche shifting: stability lower than random
                                                                unfilling.alternative = "higher", # test for niche shifting: unfilling higher than random
                                                                intersection = NA, 
                                                                rand.type = 2)
          
          niche_shift_spec_df$shift_simtest_NA_p_D <- sim_test_niche_shift_NA$p.D
          niche_shift_spec_df$shift_simtest_NA_p_expansion <- sim_test_niche_shift_NA$p.expansion
          niche_shift_spec_df$shift_simtest_NA_p_stability <- sim_test_niche_shift_NA$p.stability
          niche_shift_spec_df$shift_simtest_NA_p_unfilling <- sim_test_niche_shift_NA$p.unfilling
          
          # plot results of similarity test:
          #ecospat.plot.overlap.test(sim_test_niche_shift_NA, "D", "Similarity")
          
          # analogue conditions:
          sim_test_niche_shift_A <- ecospat.niche.similarity.test(grid.clim.EBBA1,
                                                                   grid.clim.EBBA2,
                                                                   rep = 1000,
                                                                   overlap.alternative = "lower", # test for niche shifting / divergence: niche overlap is lower than random
                                                                   expansion.alternative = "higher", # test for niche shifting: expansion higher than random
                                                                   stability.alternative = "lower", # test for niche shifting: stability lower than random
                                                                   unfilling.alternative = "higher", # test for niche shifting: unfilling higher than random
                                                                   intersection = 0, 
                                                                   rand.type = 2)
          
          niche_shift_spec_df$shift_simtest_A_p_D <- sim_test_niche_shift_A$p.D
          niche_shift_spec_df$shift_simtest_A_p_expansion <- sim_test_niche_shift_A$p.expansion
          niche_shift_spec_df$shift_simtest_A_p_stability <- sim_test_niche_shift_A$p.stability
          niche_shift_spec_df$shift_simtest_A_p_unfilling <- sim_test_niche_shift_A$p.unfilling
          
          # niche similarity test -  test for niche conservatism:
          
          # non-analogue conditions:
          sim_test_niche_cons_NA <- ecospat.niche.similarity.test(grid.clim.EBBA1, 
                                                               grid.clim.EBBA2, 
                                                               rep = 1000, 
                                                               overlap.alternative = "higher", # test for niche conservatism: niche overlap is higher than random
                                                               expansion.alternative = "lower", # test for niche conservatism: expansion lower than random
                                                               stability.alternative = "higher", # test for niche conservatism: stability higher than random
                                                               unfilling.alternative = "lower", # test for niche conservatism: unfilling lower than random
                                                               intersection = NA, 
                                                               rand.type = 2)
          
          niche_shift_spec_df$cons_simtest_NA_p_D <- sim_test_niche_cons_NA$p.D
          niche_shift_spec_df$cons_simtest_NA_p_expansion <- sim_test_niche_cons_NA$p.expansion
          niche_shift_spec_df$cons_simtest_NA_p_stability <- sim_test_niche_cons_NA$p.stability
          niche_shift_spec_df$cons_simtest_NA_p_unfilling <- sim_test_niche_cons_NA$p.unfilling
          
          # plot results of similarity test:
          #ecospat.plot.overlap.test(sim_test_niche_cons_NA, "D", "Similarity")
          
          # analogue conditions:
          sim_test_niche_cons_A <- ecospat.niche.similarity.test(grid.clim.EBBA1, 
                                                                  grid.clim.EBBA2, 
                                                                  rep = 1000, 
                                                                  overlap.alternative = "higher", # test for niche conservatism: niche overlap is higher than random
                                                                  expansion.alternative = "lower", # test for niche conservatism: expansion lower than random
                                                                  stability.alternative = "higher", # test for niche conservatism: stability higher than random
                                                                  unfilling.alternative = "lower", # test for niche conservatism: unfilling lower than random
                                                                  intersection = 0, 
                                                                  rand.type = 2)
          
          niche_shift_spec_df$cons_simtest_A_p_D <- sim_test_niche_cons_A$p.D
          niche_shift_spec_df$cons_simtest_A_p_expansion <- sim_test_niche_cons_A$p.expansion
          niche_shift_spec_df$cons_simtest_A_p_stability <- sim_test_niche_cons_A$p.stability
          niche_shift_spec_df$cons_simtest_A_p_unfilling <- sim_test_niche_cons_A$p.unfilling
          
          
          # indices for niche dynamics between EBBA1 (as z1) and EBBA2 (as z2):
          
          # non-analogue conditions:
          EBBA12_niche_dyn_NA <- ecospat.niche.dyn.index(grid.clim.EBBA1, 
                                                      grid.clim.EBBA2,
                                                      intersection = NA) # analysis performed on niche space of EBBA1 and EBBA2
          
          niche_shift_spec_df$niche_stability_NA <- EBBA12_niche_dyn_NA$dynamic.index.w['stability']
          niche_shift_spec_df$niche_expansion_NA <- EBBA12_niche_dyn_NA$dynamic.index.w['expansion']
          niche_shift_spec_df$niche_unfilling_NA <- EBBA12_niche_dyn_NA$dynamic.index.w['unfilling']
        
          # analogue conditions:
          EBBA12_niche_dyn_A <- ecospat.niche.dyn.index(grid.clim.EBBA1, 
                                                        grid.clim.EBBA2,
                                                        intersection = 0) # analysis performed on intersection of niche space of EBBA1 and EBBA2
          
          niche_shift_spec_df$niche_stability_A <- EBBA12_niche_dyn_A$dynamic.index.w['stability']
          niche_shift_spec_df$niche_expansion_A <- EBBA12_niche_dyn_A$dynamic.index.w['expansion']
          niche_shift_spec_df$niche_unfilling_A <- EBBA12_niche_dyn_A$dynamic.index.w['unfilling']
          
          # niche centroids:
          
          # as median of species EBBA distribution scores along the first 2 PCA axes:
          niche_shift_spec_df$niche_centroid_hist_X <- median(scores.sp.EBBA1[,1]) # range differs for each species due to different PCAs, should it be scaled to be comparable? xx
          niche_shift_spec_df$niche_centroid_hist_Y <- median(scores.sp.EBBA1[,2])
          niche_shift_spec_df$niche_centroid_rec_X <- median(scores.sp.EBBA2[,1])
          niche_shift_spec_df$niche_centroid_rec_Y <- median(scores.sp.EBBA2[,2])
          
          # niche shift:
          niche_shift_spec_df$niche_centroids_dist <- sqrt(((niche_shift_spec_df$niche_centroid_rec_X - niche_shift_spec_df$niche_centroid_hist_X)^2) + 
                                                           ((niche_shift_spec_df$niche_centroid_rec_Y - niche_shift_spec_df$niche_centroid_hist_Y)^2))
          
          # save plot of niche dynamics:
          jpeg(filename = file.path(EBBA_data_dir, "plots", "EBBA_1_2_niche_dyn", paste0(sel_species[s], ".jpg")),
               quality = 100, height = 920, width = 920)

          ecospat.plot.niche.dyn(grid.clim.EBBA1,
                                 grid.clim.EBBA2,
                                 quant = 0.05,
                                 interest = 2, # 1 = EBBA1 density, 2 = EBBA2 density
                                 title = paste(sel_species[s],":",
                                               "stability:", round(niche_shift_spec_df$niche_stability_NA,3),
                                               "expansion:", round(niche_shift_spec_df$niche_expansion_NA,3),
                                               "unfiling:", round(niche_shift_spec_df$niche_unfilling_NA,3)),
                                 name.axis1 = "PC1",
                                 name.axis2 = "PC2")
          dev.off()
          # blue = stability
          # red = expansion
          # green = unfilling
          
          niche_shift_spec_df
          }

# standardise niche dynamic metrics:

# non-analogue climatic conditions:
e <- niche_shift_df$niche_expansion_NA # expansion related to whole niche space (EBBA1 and EBBA2) 
s <- niche_shift_df$niche_stability_NA # stability related to whole niche space (EBBA1 and EBBA2)
u <- niche_shift_df$niche_unfilling_NA # unfilling related to whole niche space (EBBA1 and EBBA2)

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
e_a <- niche_shift_df$niche_expansion_A # expansion related to analogue EBBA2 niche space
s_a <- niche_shift_df$niche_stability_A # stability related to analogue EBBA2 niche space
u_a <- niche_shift_df$niche_unfilling_A # unfilling related to analogue EBBA1 niche space

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
niche_shift_df$std_niche_stability <- s_rel_a1 * i # stability
niche_shift_df$std_niche_unfilling <- u_rel_a1 * i # unfilling
niche_shift_df$std_niche_expansion <- e_rel_a1 * i # expansion
niche_shift_df$std_niche_abandonment <- u_rel - u_rel_a # abandonment
niche_shift_df$std_niche_abandonment[niche_shift_df$std_niche_abandonment < 0] <- 0 # slightly negative values result from rounding errors of ecospat package
niche_shift_df$std_niche_pioneering <- e_rel - e_rel_a # pioneering
niche_shift_df$std_niche_pioneering[niche_shift_df$std_niche_pioneering < 0] <- 0 # slightly negative values result from rounding errors of ecospat package

# check: values should sum to 1:
total <- niche_shift_df$std_niche_abandonment + 
  niche_shift_df$std_niche_unfilling + 
  niche_shift_df$std_niche_stability + 
  niche_shift_df$std_niche_expansion + 
  niche_shift_df$std_niche_pioneering

# save resulting data frame:
write.csv(niche_shift_df, 
          file = file.path(EBBA_data_dir, "niche_shift_results.csv"),
          row.names = FALSE)