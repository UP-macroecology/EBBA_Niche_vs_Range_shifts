# explore the results of the niche and range analysis:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(cowplot)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

datashare_EBCC <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "EBCC")
EBBA_data_dir <- file.path("Data")
plots_dir <- file.path("plots")

# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# species selection: 
sel_species <- read.csv(file.path(EBBA_data_dir, "EBBA_niche_range_shifts_species_selection_change.csv"))$species %>% # output of 2_3_species_filtering_5_climatic_niche_analysis.R
  sort # alphabetically sorted (to avoid confusion with indices later)

# EBBA 1, only selected species:
EBBA1 <- read_sf(file.path(EBBA_data_dir, "EBBA1_change.shp")) %>% # output of 1_prep_EBBA_data.R
  filter(species %in% sel_species)

# EBBA 2, only selected species:
EBBA2 <- read_sf(file.path(EBBA_data_dir, "EBBA2_change.shp")) %>% # output of 1_prep_EBBA_data.R
  filter(species %in% sel_species)

# comparable EBBA cells:
EBBA_cells <- read_sf(file.path(EBBA_data_dir, "EBBA_change.shp")) %>% # output of 1_prep_EBBA_data.R
  dplyr::select(-species) %>% 
  distinct(cell50x50, .keep_all = TRUE) # keep geometry; would be the same when using EBBA2


# ---------------------------- #
#          Explore:       ####
# ---------------------------- #

## species richness maps of selected species: ----------------------------------

# EBBA 1:
# number of species recorded in each EBBA cell:
EBBA1_richness <- EBBA1 %>% 
  group_by(cell50x50) %>% 
  summarise(EBBA1_richness = n(), .groups = "keep") %>% 
  ungroup

# EBBA 2:
# number of species recorded in each EBBA cell:
EBBA2_richness <- EBBA2 %>% 
  group_by(cell50x50) %>% 
  summarise(EBBA2_richness = n(), .groups = "keep") %>% 
  ungroup

# plot based on EBBA grid:

# EBBA grid comparable cells:
# EBBA_change_grid <- read_sf(file.path(datashare_EBCC, "EBBA_change", "ebba2_grid50x50_change_v1", "ebba2_grid50x50_change_v1.shp")) %>% 
#   rename(cell50x50 = cell50x50_)
EBBA_change_grid <- read_sf(file.path("Data", "ebba2_grid50x50_change_v1", "ebba2_grid50x50_change_v1.shp")) %>% # EPSG:3035   
  rename(cell50x50 = cell50x50_)
  
# join richness by cell ID:
EBBA_richness_grid <- EBBA_change_grid %>% 
  left_join(st_drop_geometry(EBBA1_richness)) %>% 
  left_join(st_drop_geometry(EBBA2_richness), by = "cell50x50", suffix = c("_EBBA1", "_EBBA2")) %>% 
  # difference EBBA 1 and EBBA 2 richness:
  mutate(diff = EBBA2_richness - EBBA1_richness) %>% 
  mutate(diff2 = ifelse(diff > quantile(diff, 0.9999, na.rm = TRUE), quantile(diff, 0.9999, na.rm = TRUE), diff)) %>% # cap values below and above 99.99 % quantile to avoid colour scale being dominated by outliers
  mutate(diff2 = ifelse(diff < -quantile(diff, 0.9999, na.rm = TRUE), -quantile(diff, 0.9999, na.rm = TRUE), diff2)) 

# plots:

# EBBA 1 richness:
# world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
#   st_transform("EPSG:3035")
world <- rnaturalearth::ne_coastline(returnclass = "sf") %>% 
  st_transform("EPSG:3035")
ggplot() +
  geom_sf(data = EBBA_richness_grid, aes(fill = EBBA1_richness), colour = NA) +
  scale_fill_viridis_c("richness EBBA 1", na.value = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000))

# EBBA 2 richness:
ggplot() +
  geom_sf(data = EBBA_richness_grid, aes(fill = EBBA2_richness), colour = NA) +
  scale_fill_viridis_c("richness EBBA 2", na.value = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000))

# make colour scale centered at zero:
#limit <- max(abs(EBBA_richness_grid$diff), na.rm = TRUE) * c(-1, 1)
#limit <- c(-23, 23)

# difference:
diffplot <- ggplot() +
  geom_sf(data = EBBA_richness_grid, aes(fill = diff2), colour = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA, linewidth = 0.01) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000)) +
  theme_dark() +
  scale_fill_distiller("EBBA 2 - EBBA 1", type = "div", 
                       palette = "RdBu",#, limit = limit,
                       na.value = NA)

diffplot

# species richness as facets:

# reformat data:
EBBA_richness_grid_lf <- EBBA_richness_grid %>% 
  pivot_longer(cols = ends_with("_richness"), names_to = "atlas", values_to = "richness") %>% 
  mutate(atlas = ifelse(atlas == "EBBA1_richness", "EBBA 1", "EBBA 2"))
  
# plot:
sr <- ggplot(data = EBBA_richness_grid_lf) +
  facet_wrap(~atlas) + 
  geom_sf(aes(fill = richness), colour = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA, linewidth = 0.01) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000)) +
  scale_fill_viridis_c("species richness", na.value = NA)

# save plot:
pdf(file = file.path(plots_dir, "species_richness_EBBAs.pdf"),
    height = 5, width = 10)
sr
dev.off()

# combine species richness and differences in one plot:
sr_diff <- plot_grid(sr, diffplot, labels = "AUTO",
                     align = "h",
                     axis = "tblr",
                     rel_widths = c(2, 1.3),
                     vjust = 7)

# save plot:
pdf(file = file.path(plots_dir, "species_richness_EBBAs_diff.pdf"),
    height = 5, width = 12)
sr_diff
dev.off()


## dynamics boxplot: -----------------------------------------------------------

## species-specific background or whole EBBA area as background:
env_background_species_specific <- TRUE # set to FALSE to use whole EBBA area as environmental background

# load analyses results (output of 4_niche_shift_analysis.R and 4_range_shift_analysis.R)
if(env_background_species_specific){
  niche_results <- read.csv(file = file.path(EBBA_data_dir, "niche_shift_results_env_species_spec_change.csv"))
  range_results <- read.csv(file = file.path(EBBA_data_dir, "range_shift_results_env_species_spec_change.csv"))
  } else {
  niche_results <- read.csv(file = file.path(EBBA_data_dir, "niche_shift_results_env_EBBA_change.csv"))
  range_results <- read.csv(file = file.path(EBBA_data_dir, "range_shift_results_env_EBBA_change.csv"))
  }

# species (%) with significant dynamic values:
## niche:
niche_test_sign <- niche_results %>% 
  select(c(species, matches("_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species)*100, 2)))
## range:
range_test_sign <- range_results %>% 
  select(c(species, matches("_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species)*100, 2)))

# species (%) with significantly higher dynamics than by chance:
spec_sign_higher <- c("NA", # niche abandonment
                      niche_test_sign$shift_p_U_NA_n_sig, # niche unfilling
                      niche_test_sign$cons_p_S_NA_n_sig, # niche stability
                      niche_test_sign$shift_p_E_NA_n_sig, # niche expansion
                      "NA", # niche pioneering,
                      range_test_sign$shift_p_U_NA_n_sig, # range unfilling
                      range_test_sign$cons_p_S_NA_n_sig, # range stabiliy
                      range_test_sign$shift_p_E_NA_n_sig # range expansion
)

# species (%) with significantly lower dynamics than by chance:
spec_sign_lower <- c("NA", # niche abandonment
                      niche_test_sign$cons_p_U_NA_n_sig, # niche unfilling
                      niche_test_sign$shift_p_S_NA_n_sig, # niche stability
                      niche_test_sign$cons_p_E_NA_n_sig, # niche expansion
                      "NA", # niche pioneering,
                      range_test_sign$cons_p_U_NA_n_sig, # range unfilling
                      range_test_sign$shift_p_S_NA_n_sig, # range stabiliy
                      range_test_sign$cons_p_E_NA_n_sig # range expansion
)

# join niche and range shift results, convert to long format:
niche_range_df <- niche_results %>% 
  left_join(range_results, by = "species", suffix = c("_niche", "_range")) %>% 
  pivot_longer(cols = ends_with("_std"),
               names_to = c("category", "metric"), names_pattern = "(.*)_(.*)_",
               values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("abandonment", "unfilling", "stability", "expansion", "pioneering")))

# data set for labels:
labelsdat <- tibble(category = c(rep("niche", 5), rep("range", 3)),
                    metric = factor(c("abandonment", "unfilling", "stability", "expansion", "pioneering", "unfilling", "stability", "expansion"), 
                                    levels = c("abandonment", "unfilling", "stability", "expansion", "pioneering")),
                    sign_higher = spec_sign_higher,
                    sign_lower = spec_sign_lower,
                    ypos_higher = 110,
                    ypos_lower = -15)

# plot:
p <- ggplot(niche_range_df, aes(x = category, y = value*100, fill = metric)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("dynamics [%]") +
  xlab("") +
  scale_fill_viridis_d("dynamics") +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(add = 5)) +
  annotate(geom = "text", label = "Species (%) with significantly higher dynamics than by chance:", 
           x = 0.5, y = 115, hjust = 0, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_sign_higher, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) +
  geom_hline(yintercept = -5) +
  annotate(geom = "text", label = "Species (%) with significantly lower dynamics than by chance:", 
           x = 0.5, y = -10, hjust = 0, vjust = 0, size = 3) +
  geom_text(data = labelsdat, aes(label = spec_sign_lower, y = ypos_lower),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3)

p

# save plot:
# png(filename = file.path(plots_dir, paste0("Boxplot_dynamics_env_species_spec_", env_background_species_specific, ".png")),
#     height = 452, width = 750, res = 100)
# p
# dev.off()

pdf(file = file.path(plots_dir, paste0("Boxplot_dynamics_env_species_spec_", env_background_species_specific, ".pdf")),
    height = 5, width = 8)
p
dev.off()



# for which species are results significant: ----

spec_range_sign_stab <- range_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
  filter(cons_p_S_NA <= 0.05)  %>% # higher than by chance
  pull(species)

spec_range_sign_stab <- range_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
  filter(cons_p_S_A <= 0.05)  %>% # higher than by chance
  pull(species)

spec_niche_sign_stab <- niche_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
  filter(cons_p_S_NA <= 0.05)  %>% # higher than by chance
  pull(species)

spec_niche_sign_low_unf <- niche_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
  filter(cons_p_U_NA <= 0.05)  %>% # higher than by chance
  pull(species)

## differences in bioclim variable values between historic and recent time period: ----

bioclim1_hist <- rast(file.path(EBBA_data_dir, "Bioclim_1981_1990", "CHELSA_bio1_1981_1990_50km.tif"))
bioclim1_rec <- rast(file.path(EBBA_data_dir, "Bioclim_2009_2018", "CHELSA_bio1_2009_2018_50km.tif"))
change_rast <- (bioclim1_hist - bioclim1_rec)/bioclim1_hist * 100 # change in percentage compared to historic time period

biovars_hist_rast <- rast(list.files(file.path(EBBA_data_dir, "Bioclim_1981_1990"),
                                     pattern = "50km.tif$",
                                     full.names = TRUE))
biovars_rec_rast <- rast(list.files(file.path(bioclim_data_dir, "Bioclim_2009_2018"),
                                    pattern = "50km.tif$",
                                    full.names = TRUE))
change_rast <- (biovars_hist_rast - biovars_rec_rast)/biovars_hist_rast * 100 
# mask with EBBA comparable area:
EBBA_mask_sp <- as.polygons(rast(file.path(EBBA_data_dir, "EBBA_mask_ecospat.tif")))
change_rast_masked <- mask(change_rast, EBBA_mask_sp)

# colors:
for(i in 1:nlyr(change_rast_masked)){
  ceil <- values(change_rast_masked[[i]], na.rm=TRUE) |> abs() |> max() |> ceiling() 
  pal <- leaflet::colorNumeric(palette = "RdBu", domain=c(-ceil, ceil), reverse = T)
  # plot:
  plot(change_rast_masked[[i]], range=c(-ceil, ceil), col=pal(seq(-ceil,ceil,.1)), main = change_rast[[i]]@ptr$names)
}


## plot EBBA 1 and EBBA 2 presences of single species: --------------------------

world <- rnaturalearth::ne_coastline(returnclass = "sf") %>% 
  st_transform("EPSG:3035")

for(s in 1:2){#length(sel_species)){
  
  EBBA1_spec <- EBBA1 %>%
    filter(species == sel_species[s])
  
  EBBA2_spec <- EBBA2 %>%
    filter(species == sel_species[s])
  
  map <- ggplot(world) +
    geom_sf(color = "gray50") +
    lims(x = c(1100000, 6500000), y = c(1500000, 6400000)) +
    geom_sf(data = EBBA_change_grid, colour = "black", size = 1) +
    geom_sf(data = EBBA1_spec, colour = "yellow", size = 1) +
    geom_sf(data = EBBA2_spec, colour = "blue", size = 1) +
    ggtitle(sel_species[s])
  
  print(map)
  
  readline(prompt = "Press [enter] to continue") # to go to next species; press [esc] to stop the loop
}


## misc. range / niche analysis explorations: ----------------------------------

res_df <- read.csv(file = file.path(EBBA_data_dir, "range_shift_results_env_EBBA.csv"))
#res_df <- read.csv(file = file.path(EBBA_data_dir, "range_shift_results_env_species_spec.csv"))
#res_df <- read.csv(file = file.path(EBBA_data_dir, "niche_shift_results_env_EBBA.csv"))
#res_df <- read.csv(file = file.path(EBBA_data_dir, "niche_shift_results_env_species_spec.csv"))

# for range analysis:
summary(res_df$centroid_dist)
summary(res_df$centroid_NS_shift)
summary(res_df$centroid_EW_shift)

# species for which similarity test with regard to niche / range shifting yields significant results:

# Schoener's D:
# non-analogue conditions:
res_df %>% 
  filter(shift_p_D_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
res_df %>% 
  filter(shift_p_D_A < 0.05) %>% 
  pull(species)
# -> none

# expansion:
# non-analogue conditions:
res_df %>% 
  filter(shift_p_E_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
res_df %>% 
  filter(shift_p_E_A < 0.05) %>% 
  pull(species)
# -> none

# stability:
res_df %>% 
  filter(shift_p_S_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
res_df %>% 
  filter(shift_p_S_A < 0.05) %>% 
  pull(species)
# -> none

# unfilling:
res_df %>% 
  filter(shift_p_U_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
res_df %>% 
  filter(shift_p_U_A < 0.05) %>% 
  pull(species)
# -> none

# species for which similarity test with regard to niche / range conservatism yields significant results:

# Schoener's D:
# non-analogue conditions:
D_NA_cons_sign <- res_df %>% 
  filter(cons_p_D_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
D_A_cons_sign <- res_df %>% 
  filter(cons_p_D_A < 0.05) %>% 
  pull(species)
identical(D_NA_cons_sign, D_A_cons_sign) # TRUE
length(D_NA_cons_sign)
# ->except for:
res_df %>% 
  filter(cons_p_D_A >= 0.05) %>% 
  pull(species)

# expansion:
# non-analogue conditions:
exp_NA_cons_sign <- res_df %>% 
  filter(cons_p_E_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
exp_A_cons_sign <- res_df %>% 
  filter(cons_p_E_A < 0.05) %>% 
  pull(species)
identical(exp_NA_cons_sign, exp_A_cons_sign)
exp_A_cons_sign[which(!exp_A_cons_sign %in% exp_NA_cons_sign)]
exp_NA_cons_sign[which(!exp_NA_cons_sign %in% exp_A_cons_sign)]

# stability:
# non-analogue conditions:
st_NA_cons_sign <- res_df %>% 
  filter(cons_p_S_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
st_A_cons_sign <- res_df %>% 
  filter(cons_p_S_A < 0.05) %>% 
  pull(species)
identical(st_NA_cons_sign, st_A_cons_sign)
st_A_cons_sign[which(!st_A_cons_sign %in% st_NA_cons_sign)]
st_NA_cons_sign[which(!st_NA_cons_sign %in% st_A_cons_sign)]
identical(exp_NA_cons_sign, st_NA_cons_sign)
identical(exp_A_cons_sign, st_A_cons_sign)
# for which not:
res_df %>% 
  filter(cons_p_S_NA >= 0.05) %>% 
  pull(species) 
res_df %>% 
  filter(cons_p_S_A >= 0.05) %>% 
  pull(species) 

# unfilling:
# non-analogue conditions:
unf_NA_cons_sign <- res_df %>% 
  filter(cons_p_U_NA < 0.05) %>% 
  pull(species)
# analogue conditions:
unf_A_cons_sign <- res_df %>% 
  filter(cons_p_U_A < 0.05) %>% 
  pull(species)
identical(unf_NA_cons_sign, unf_A_cons_sign) # TRUE
unf_NA_cons_sign[which(!unf_NA_cons_sign %in% unf_A_cons_sign)]
unf_A_cons_sign[which(!unf_A_cons_sign %in% unf_NA_cons_sign)]
# except for:
res_df %>% 
  filter(cons_p_U_NA >= 0.05) %>% 
  pull(species)

sel_species[which(sel_species %in% c(D_NA_cons_sign) &
                    sel_species %in% c(D_A_cons_sign) &
                    sel_species %in% c(exp_NA_cons_sign) &
                    sel_species %in% c(exp_A_cons_sign) &
                    sel_species %in% c(st_NA_cons_sign) &
                    sel_species %in% c(st_A_cons_sign) &
                    sel_species %in% c(unf_NA_cons_sign) &
                    sel_species %in% c(unf_A_cons_sign)
                    )]