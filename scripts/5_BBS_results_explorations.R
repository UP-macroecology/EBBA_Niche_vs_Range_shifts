# explore the results of the BBS niche and range analysis:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(cowplot)

data_dir <- file.path("Data")

# which historic time period should be used:
#hist_years <- 1980:1983 # maximum gap between historic and recent time period
hist_years <- 1995:1998 # similar gap between historic and recent time period as in EBBA analysis

# environmental background: use presences and absences within 600 km buffer around presences or all true absences within conterminous US:
env_background_contUS <- TRUE # set to FALSE to use presences and absences within 600 km buffer around presences


# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# species selection (V1 and V2): 
sel_species <- read.csv(file = file.path(data_dir, "species_stability_contUS_BL22_120423.csv")) %>% 
  arrange(-stability) %>% 
  filter(stability >= 0.5) %>% 
  pull(species) %>% 
  sort

if(all(hist_years == 1980:1983)){
  
  # which species for V1:
  load(file = file.path(data_dir, "BBS_prep_steps1-4_V1.RData")) # output of 2_1_BBS_species_filtering_1-4.R
  species_filtered_V1 <- sort(unique(hist_prep_df$species))
  sel_species_final <- species_filtered_V1[which(species_filtered_V1 %in% sel_species)] # 207
  
  # BBS data:
  
  ## historic period, only selected species:
  BBS_hist_sf <- read_sf(file.path(data_dir, "BBS_historic1_proj_centroids.shp")) %>% # output of 1_BBS_prep_data.R
    filter(species %in% sel_species_final)
  
  ## recent period, only selected species:
  BBS_rec_sf <- read_sf(file.path(data_dir, "BBS_recent1_proj_centroids.shp")) %>% # output of 1_BBS_prep_data.R
    filter(species %in% sel_species_final)
  
  # load analyses results (output of 4_BBS_niche_shift_analysis.R and 4_BBS_range_shift_analysis.R)
  if(env_background_contUS){
    niche_results <- read.csv(file = file.path(data_dir, "BBS_niche_shift_results_env_contUS_hist1980_1983.csv"))
    range_results <- read.csv(file = file.path(data_dir, "BBS_range_shift_results_env_contUS_hist1980_1983.csv"))
    # results boxplot name:
    res_boxplot <- "BBS_boxplot_dynamics_env_contUS_hist_1980_1983"
  } else {
    niche_results <- read.csv(file = file.path(data_dir, "BBS_niche_shift_results_env_600km_hist1980_1983.csv"))
    range_results <- read.csv(file = file.path(data_dir, "BBS_range_shift_results_env_600km_hist1980_1983.csv"))
    # results boxplot name:
    res_boxplot <- "BBS_boxplot_dynamics_env_600km_hist_1980_1983"
  }
  
} else {
  
  # which species for V2:
  load(file = file.path(data_dir, "BBS_prep_steps1-4_V2.RData")) # output of 2_1_BBS_species_filtering_1-4.R
  species_filtered_V2 <- sort(unique(hist_prep_df$species))
  sel_species_final <- species_filtered_V2[which(species_filtered_V2 %in% sel_species)] # 296
  
  # BBS data:
  
  ## historic period, only selected species:
  BBS_hist_sf <- read_sf(file.path(data_dir, "BBS_historic2_proj_centroids.shp")) %>% # output of 1_BBS_prep_data.R
    filter(species %in% sel_species_final)
  
  ## recent period, only selected species:
  BBS_rec_sf <- read_sf(file.path(data_dir, "BBS_recent2_proj_centroids.shp")) %>% # output of 1_BBS_prep_data.R
    filter(species %in% sel_species_final)
  
  # load analyses results (output of 4_BBS_niche_shift_analysis.R and 4_BBS_range_shift_analysis.R)
  if(env_background_contUS){
    niche_results <- read.csv(file = file.path(data_dir, "BBS_niche_shift_results_env_contUS_hist1995_1998.csv"))
    range_results <- read.csv(file = file.path(data_dir, "BBS_range_shift_results_env_contUS_hist1995_1998.csv"))
    # results boxplot name:
    res_boxplot <- "BBS_boxplot_dynamics_env_contUS_hist_1995_1998"
  } else {
    niche_results <- read.csv(file = file.path(data_dir, "BBS_niche_shift_results_env_600km_hist1995_1998.csv"))
    range_results <- read.csv(file = file.path(data_dir, "BBS_range_shift_results_env_600km_hist1995_1998.csv"))
    # results boxplot name:
    res_boxplot <- "BBS_boxplot_dynamics_env_600km_hist_1995_1998"
  }
}


# ---------------------------- #
#          Explore:       ####
# ---------------------------- #

## species richness maps of selected species: ----------------------------------

# historic period:
# number of species recorded on each route:
BBS_hist_richness <- BBS_hist_sf %>% 
  filter(pres == 1) %>% 
  group_by(RTENO) %>% 
  summarise(hist_richness = n(), .groups = "keep") %>% 
  ungroup # 522 routes for version 1980-83

# recent period:
# number of species recorded on each route:
BBS_rec_richness <- BBS_rec_sf %>% 
  filter(pres == 1) %>% 
  group_by(RTENO) %>% 
  summarise(rec_richness = n(), .groups = "keep") %>% 
  ungroup 


# plots:

# historic period richness:
 world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
   st_transform("ESRI:102003")
# world <- rnaturalearth::ne_coastline(returnclass = "sf") %>% 
#   st_transform("ESRI:102003")
ggplot() +
  geom_sf(data = BBS_hist_richness, aes(colour = hist_richness)) +
  scale_colour_viridis_c("richness historic period", na.value = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000))

# recent period richness:
ggplot() +
  geom_sf(data = BBS_rec_richness, aes(colour = rec_richness)) +
  scale_colour_viridis_c("richness recent period", na.value = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000))

# make colour scale centered at zero:
#limit <- max(abs(EBBA_richness_grid$diff), na.rm = TRUE) * c(-1, 1)
#limit <- c(-23, 23)
BBS_hist_richness

# join richness by cell ID:
BSS_richness_diff <- BBS_hist_richness %>% 
  left_join(st_drop_geometry(BBS_rec_richness)) %>% 
  # difference EBBA 1 and EBBA 2 richness:
  mutate(diff = rec_richness - hist_richness) %>% 
  mutate(diff2 = ifelse(diff > quantile(diff, 0.99, na.rm = TRUE), quantile(diff, 0.99, na.rm = TRUE), diff)) %>% # cap values below and above 99 % quantile to avoid colour scale being dominated by outliers
  mutate(diff2 = ifelse(diff < -quantile(diff, 0.99, na.rm = TRUE), -quantile(diff, 0.99, na.rm = TRUE), diff2)) 


# difference:
diffplot <- ggplot() +
  geom_sf(data = BSS_richness_diff, aes(colour = diff)) + # diff2
  geom_sf(data = world, colour = "gray20", fill = NA, linewidth = 0.01) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000)) +
  theme_dark() +
  scale_colour_distiller("rec. period - hist. period", type = "div", 
                       palette = "RdBu",#, limit = limit,
                       na.value = NA)

diffplot

# species richness as facets:

# reformat data:
BBS_richness_lf <- BSS_richness_diff %>% 
  pivot_longer(cols = ends_with("_richness"), names_to = "period", values_to = "richness") %>% 
  mutate(period = ifelse(period == "hist_richness", "historic", "recent"))

# plot:
sr <- ggplot(data = BBS_richness_lf) +
  facet_wrap(~period) + 
  geom_sf(aes(colour = richness)) +
  geom_sf(data = world, colour = "gray20", fill = NA, linewidth = 0.01) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000)) +
  scale_colour_viridis_c("species richness", na.value = NA)

# save plot:
pdf(file = file.path("plots", paste0("species_richness_BBS_", min(hist_years), "_", max(hist_years), ".pdf")),
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
pdf(file = file.path("plots", paste0("species_richness_BBS_diff_", min(hist_years), "_", max(hist_years), ".pdf")),
    height = 5, width = 12)
sr_diff
dev.off()



## dynamics boxplot: -----------------------------------------------------------


# species (%) with significant dynamic values:

## niche:
niche_test_sign <- niche_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species)*100, 2)))

## range:
range_test_sign <- range_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
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
# png(filename = file.path("plots", paste0(res_boxplot, ".png")),
#     height = 452, width = 750, res = 100)
# p
# dev.off()

pdf(file = file.path("plots", paste0(res_boxplot, ".pdf")),
    height = 5, width = 8)
p
dev.off()


## for which species are results significant: ----

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
