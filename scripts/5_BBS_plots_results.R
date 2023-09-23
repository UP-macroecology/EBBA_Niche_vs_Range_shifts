# explore the results of the BBS niche and range analysis:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(cowplot)
library(ade4)
library(quantreg)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
data_dir <- file.path("data", "BBS_analysis")
plots_dir <- file.path("plots")

# which historic time period should be used:
hist_years <- 1980:1983 # maximum gap between historic and recent time period
# hist_years <- 1987:1990 # similar gap between historic and recent time period as in EBBA analysis

# environmental background: presences and absences within 500 km buffer around presences (TRUE) or all true absences within conterminous US (FALSE):
bg_spec <- TRUE


# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# species selection: 
sel_species <- read.csv(file = file.path(data_dir, "BBS_stability_PCA_contUS_BL22.csv")) %>% 
  filter(stability >= 0.5) %>% 
  pull(species) %>% 
  sort

# BBS data, only selected species:
load(file = file.path(data_dir, paste0("BBS_prep_steps1-4_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".RData"))) # output of 2_1_BBS_species_filtering_1-4.R
species_filtered <- sort(unique(hist_prep_df$species))
sel_species_final <- species_filtered[which(species_filtered %in% sel_species)]

BBS_hist_sf <- read_sf(file.path(data_dir, paste0("BBS_historic_centr_proj_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".shp"))) %>% # output of 1_BBS_prep_data.R
  filter(species %in% sel_species_final)

BBS_rec_sf <- read_sf(file.path(data_dir, paste0("BBS_recent_centr_proj_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".shp"))) %>% # output of 1_BBS_prep_data.R
  filter(species %in% sel_species_final)

# load analyses results (output of 4_BBS_niche_shift_analysis.R and 4_BBS_range_shift_analysis.R)
niche_results <- read.csv(file.path(data_dir, paste0("BBS_niche_shift_results_bg_",
                                                     ifelse(bg_spec, "spec", "US"), "_hist",
                                                     ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                     ".csv"))) # 195 (81-83), 233 (88-90)

range_results <- read.csv(file.path(data_dir, paste0("BBS_range_shift_results_bg_",
                                                     ifelse(bg_spec, "spec", "US"), "_hist",
                                                     ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                     ".csv")))


# ---------------------------- #
#          Explore:       ####
# ---------------------------- #

## climate change as PCA based maps: ----

# climate data historic period:
biovars_hist_rast <- rast(list.files(file.path(data_dir, paste0("Bioclim_", min(hist_years), "_", max(hist_years))),
                                     pattern = "1km.tif$",
                                     full.names = TRUE))

# aggregate to 50 km resolution (as EBBA):
biovars_hist_rast_50km <- biovars_hist_rast %>% 
  aggregate(fact = 50, fun = "mean", na.rm = TRUE)

# extract cell values:
hist_climate_data <- values(biovars_hist_rast_50km, dataframe = TRUE, na.rm = TRUE)

# climate data recent period:
biovars_rec_rast <- rast(list.files(file.path(data_dir, "Bioclim_2015_2018"),
                                    pattern = "1km.tif$",
                                    full.names = TRUE))

# aggregate to 50 km resolution (as EBBA):
biovars_rec_rast_50km <- biovars_rec_rast %>% 
  aggregate(fact = 50, fun = "mean", na.rm = TRUE)

# extract cell values:
rec_climate_data <- values(biovars_rec_rast_50km, dataframe = TRUE, na.rm = TRUE)


# assess climate niche by using the first 2 PCA axes:
# calibrating the PCA in the whole study area:
pca.env <- dudi.pca(rbind(hist_climate_data, rec_climate_data)[, paste0("bio", 1:19)],
                    scannf = FALSE,
                    nf = 2) # number of axes

# plot(x = pca.env$li$Axis1, y = pca.env$li$Axis2)

p_biplot2 <- factoextra::fviz_pca_biplot(pca.env, repel = TRUE,
                                         col.var = "black", # Variables color
                                         col.ind = c(rep("historic", nrow(hist_climate_data)),
                                                     rep("recent", nrow(rec_climate_data))),  # Individuals color
                                         label = "var",
                                         palette = c( "red", "blue"),
                                         pointsize = 1.5,
                                         alpha.ind = 0.1,
                                         legend.title = "",
                                         title = "",
                                         legend = c(0.9, 0.9))


# PCA scores historic points:
scores.clim.hist <- suprow(pca.env, hist_climate_data[, paste0("bio", 1:19)])$li
# PCA scores recent points:
scores.clim.rec <- suprow(pca.env, rec_climate_data[, paste0("bio", 1:19)])$li

# difference in first 2 PCA axes per cell:
map_df <- data.frame("cell_num" = 1:ncell(biovars_hist_rast_50km[[1]]),
                     "diff_PC1" = NA,
                     "diff_PC2" = NA)
map_df$diff_PC1[as.numeric(row.names(scores.clim.hist))] <- scores.clim.rec$Axis1 - scores.clim.hist$Axis1 # change       
map_df$diff_PC2[as.numeric(row.names(scores.clim.hist))] <- scores.clim.rec$Axis2 - scores.clim.hist$Axis2 # change  

# PC difference rasters:
diff_PC1_rast <- rast(x = biovars_hist_rast_50km[[1]],
                      vals = map_df$diff_PC1,
                      names = "diff_PC1")

diff_PC2_rast <- rast(x = biovars_hist_rast_50km[[1]],
                      vals = map_df$diff_PC2,
                      names = "diff_PC2")

# change to sf object to use same code for plots:
map_sf <- diff_PC1_rast %>% 
  as.polygons(trunc = FALSE) %>% 
  st_as_sf

diff_PC2_sf <- diff_PC2_rast %>% 
  as.polygons(trunc = FALSE) %>% 
  st_as_sf

map_sf$diff_PC2 <- diff_PC2_sf$diff_PC2

# maps:
# plot(map_sf["diff_PC1"])
# plot(map_sf["diff_PC2"])
# plot(map_sf["diff_PC1"], axes = TRUE, graticule = TRUE)

limit <- max(abs(map_sf$diff_PC1)) * c(-1, 1) # to center colour scale

world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
  st_transform("EPSG:3035")

p_diff_pc1 <- ggplot() +
  geom_sf(data = map_sf, aes(fill = diff_PC1), colour = NA) +
  scale_fill_distiller("PC1\n(recent - historic)", type = "div", 
                       palette = "RdBu", limit = limit) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000))

p_diff_pc2 <- ggplot() +
  geom_sf(data = map_sf, aes(fill = diff_PC2), colour = NA) +
  scale_fill_distiller("PC2\n(recent - historic)", type = "div", 
                       palette = "RdBu", limit = limit) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000))

# combine biplot and difference maps in one plot:
p_climate_change <- plot_grid(p_biplot2, p_diff_pc1, p_diff_pc2,
                              labels = "AUTO",
                              align = "h",
                              nrow = 1,
                              axis = "tblr",
                              rel_widths = c(1.5, 1.5, 1.5),
                              rel_heights = c(0.25, 5, 5)
)

# save plot:
pdf(file = file.path(plots_dir, "climate_change_PCAs", paste0("BBS_climate_change_PCA_plots_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".pdf")),
    height = 5, width = 20)
p_climate_change
dev.off()


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

world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
  st_transform("ESRI:102003")

# join richness by cell ID:
BSS_richness_diff <- BBS_hist_richness %>% 
  left_join(st_drop_geometry(BBS_rec_richness)) %>% 
  mutate(diff = rec_richness - hist_richness)

limit <- max(abs(BSS_richness_diff$diff), na.rm = TRUE) * c(-1, 1) # to center colour scale

# difference:
diffplot <- ggplot() +
  geom_sf(data = BSS_richness_diff, aes(colour = diff)) +
  geom_sf(data = world, colour = "gray20", fill = NA, linewidth = 0.01) +
  lims(x = c(-2400000, 2200000), y = c(-1300000, 1600000)) +
  theme_dark() +
  scale_colour_distiller("rec. period - hist. period", type = "div", 
                       palette = "RdBu", limit = limit,
                       na.value = NA)

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
pdf(file = file.path("plots", "species_richness", paste0("BBS_species_richness_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".pdf")),
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
pdf(file = file.path("plots", "species_richness", paste0("BBS_species_richness_hist", ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), "_incl_diff.pdf")),
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
  mutate(across(everything(), ~ round(.x/length(sel_species_final)*100, 2)))

## range:
range_test_sign <- range_results %>% 
  dplyr::select(c(species, matches("_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_final)*100, 2)))

# species (%) with significantly higher dynamics than by chance:
spec_sign_higher <- c("", # niche abandonment
                      niche_test_sign$shift_p_U_A_n_sig, # niche unfilling
                      niche_test_sign$cons_p_S_A_n_sig, # niche stability
                      niche_test_sign$shift_p_E_A_n_sig, # niche expansion
                      "", # niche pioneering,
                      range_test_sign$shift_p_U_A_n_sig, # range unfilling
                      range_test_sign$cons_p_S_A_n_sig, # range stabiliy
                      range_test_sign$shift_p_E_A_n_sig # range expansion
)

# species (%) with significantly lower dynamics than by chance:
spec_sign_lower <- c("", # niche abandonment
                     niche_test_sign$cons_p_U_A_n_sig, # niche unfilling
                     niche_test_sign$shift_p_S_A_n_sig, # niche stability
                     niche_test_sign$cons_p_E_A_n_sig, # niche expansion
                     "", # niche pioneering,
                     range_test_sign$cons_p_U_A_n_sig, # range unfilling
                     range_test_sign$shift_p_S_A_n_sig, # range stabiliy
                     range_test_sign$cons_p_E_A_n_sig # range expansion
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
                    ypos_lower = -22)

# plot:
p <- ggplot(niche_range_df, aes(x = category, y = value*100, fill = metric)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Metrics [%]") +
  xlab("") +
  scale_fill_viridis_d("Metrics") +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(add = 10)) +
  annotate(geom = "text", label = "Species (%) with significantly higher metrics than expected by chance:", 
           x = 0.5, y = 118, hjust = 0, vjust = 0, size = 3) +
  geom_text(data = labelsdat, aes(label = spec_sign_higher, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) +
  geom_hline(yintercept = -5) +
  annotate(geom = "text", label = "Species (%) with significantly lower metrics than expected by chance:", 
           x = 0.5, y = -14, hjust = 0, vjust = 0, size = 3) +
  geom_text(data = labelsdat, aes(label = spec_sign_lower, y = ypos_lower),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3)

p

# save plot:
pdf(file = file.path("plots", "dynamics_boxplots", paste0("BBS_boxplot_dynamics_bg_", 
                                                          ifelse(bg_spec, "spec", "US"), "_hist",
                                                          ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                          ".pdf")),
    height = 3, width = 8)
p
dev.off()


## dynamics correlation plots: ----

metrics_df <- range_results %>% 
  left_join(niche_results, by = c("species" = "species")) %>% 
  dplyr::select(range_stability_std, range_unfilling_std, range_expansion_std,
         niche_stability_std, niche_unfilling_std, niche_expansion_std)

cor1 <- with(metrics_df,cor.test(range_expansion_std, niche_unfilling_std))
xrange <- range(metrics_df$niche_unfilling_std)
yrange <- range(metrics_df$range_expansion_std)
p1 <- ggplot(data = metrics_df, aes(y = range_expansion_std, x = niche_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche U") + ylab("range E") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor1$estimate,2),ifelse(cor1$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")
cor2 <- with(metrics_df,cor.test(range_expansion_std, niche_stability_std))
xrange <- range(metrics_df$niche_stability_std)
yrange <- range(metrics_df$range_expansion_std)
p2 <- ggplot(data = metrics_df, aes(y = range_expansion_std, x = niche_stability_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche S") + ylab("range E") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor2$estimate,2),ifelse(cor2$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")
cor3 <- with(metrics_df,cor.test(range_expansion_std, niche_expansion_std))
xrange <- range(metrics_df$niche_expansion_std)
yrange <- range(metrics_df$range_expansion_std)
p3 <- ggplot(data = metrics_df, aes(y = range_expansion_std, x = niche_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche E") + ylab("range E") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor3$estimate,2),ifelse(cor3$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")

cor4 <- with(metrics_df,cor.test(range_stability_std, niche_unfilling_std))
xrange <- range(metrics_df$niche_unfilling_std)
yrange <- range(metrics_df$range_stability_std)
p4 <- ggplot(data = metrics_df, aes(y = range_stability_std, x = niche_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche U") + ylab("range S") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor4$estimate,2),ifelse(cor4$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")
cor5 <- with(metrics_df,cor.test(range_stability_std, niche_stability_std))
xrange <- range(metrics_df$niche_stability_std)
yrange <- range(metrics_df$range_stability_std)
p5 <- ggplot(data = metrics_df, aes(y = range_stability_std, x = niche_stability_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche S") + ylab("range S") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.95), label = paste0("r=",round(cor5$estimate,2),ifelse(cor5$p.value<0.05,"*","")), vjust=0, hjust=1, color="black")
cor6 <- with(metrics_df,cor.test(range_stability_std, niche_expansion_std))
xrange <- range(metrics_df$niche_expansion_std)
yrange <- range(metrics_df$range_stability_std)
p6 <- ggplot(data = metrics_df, aes(y = range_stability_std, x = niche_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche E") + ylab("range S") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor6$estimate,2),ifelse(cor6$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")

cor7 <- with(metrics_df,cor.test(range_unfilling_std, niche_unfilling_std))
xrange <- range(metrics_df$niche_unfilling_std)
yrange <- range(metrics_df$range_unfilling_std)
p7 <- ggplot(data = metrics_df, aes(y = range_unfilling_std, x = niche_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche U") + ylab("range U") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor7$estimate,2),ifelse(cor7$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")
cor8 <- with(metrics_df,cor.test(range_unfilling_std, niche_stability_std))
xrange <- range(metrics_df$niche_stability_std)
yrange <- range(metrics_df$range_unfilling_std)
p8 <- ggplot(data = metrics_df, aes(y = range_unfilling_std, x = niche_stability_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche S") + ylab("range U") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor8$estimate,2),ifelse(cor8$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")
cor9 <- with(metrics_df,cor.test(range_unfilling_std, niche_expansion_std))
xrange <- range(metrics_df$niche_expansion_std)
yrange <- range(metrics_df$range_unfilling_std)
p9 <- ggplot(data = metrics_df, aes(y = range_unfilling_std, x = niche_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche E") + ylab("range U") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor9$estimate,2),ifelse(cor9$p.value<0.05,"*","")), vjust=1, hjust=1, color="black")


pdf(file = file.path(plots_dir, "dynamics_correlations", paste0("BBS_correlation_dynamics_bg_",
                                                                ifelse(bg_spec, "spec", "US"), "_hist",
                                                                ifelse(all(hist_years == 1980:1983), "81-83", "88-90"), ".pdf")),
    height = 5, width = 8)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9,
          align = "h",
          axis = "tblr",
          nrow = 3,
          ncol = 3)
dev.off()



## direction of range shifts: ----

# same for using species-specific background and whole cont. US as background

# one line for each species:
dir_shift_p <- ggplot(range_results, aes(x = shift_direction, y = Eucl_distance/1000)) +
  geom_point() + 
  geom_segment(aes(x = shift_direction, xend = shift_direction, y=0, yend = Eucl_distance/1000)) +
  ylab("Range shift [km]") +
  xlab("") +
  coord_polar(start = pi, direction = 1) +
  scale_x_continuous(limits = c(-180,180),
                     breaks = seq(-180, 180, by = 30),
                     minor_breaks = seq(-180, 180, by = 15)) +
  theme_minimal()

# one line for each species, log scale
dir_shift_p2 <- ggplot(range_results, aes(x = shift_direction, y = Eucl_distance/1000)) +
  geom_point() + 
  geom_segment(aes(x = shift_direction, xend = shift_direction, y=0, yend = Eucl_distance/1000)) +
  ylab("Range shift [km]") +
  xlab("") +
  coord_polar(start = pi, direction = 1) +
  scale_x_continuous(limits = c(-180,180),
                     breaks = seq(-180, 180, by = 30),
                     minor_breaks = seq(-180, 180, by = 15)) +
  scale_y_continuous(trans = 'log',
                     breaks = c(0, 25, 50,100, 200,400)) +
  theme_minimal()


# histogram:
dir_shift_hist <- ggplot(range_results, aes(x = shift_direction)) +
  geom_histogram(breaks = seq(-180, 180, by = 10)) +
  coord_polar(start = pi, direction = 1) +
  scale_x_continuous(limits = c(-180, 180),
                     breaks = seq(-180, 180, by = 30),
                     minor_breaks = seq(-180, 180, by = 15)) +
  ylab("Number of species") + 
  xlab("") +
  theme_minimal()


# save plots:
pdf(file = file.path(plots_dir, "range_shift_direction", paste0("BBS_rshift_direction", 
                                                                ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                                ".pdf")),
    height = 5, width = 6)
dir_shift_p
dev.off()

pdf(file = file.path(plots_dir, "range_shift_direction", paste0("BBS_rshift_direction_log", 
                                                                ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                                ".pdf")),
    height = 5, width = 6)
dir_shift_p2
dev.off()

pdf(file = file.path(plots_dir, "range_shift_direction", paste0("BBS_rshift_direction_histogram_hist", 
                                                                ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                                ".pdf")),
    height = 5, width = 6)
dir_shift_hist
dev.off()

# combine range shift and direction:
p_climate_change <- plot_grid(dir_shift_p2, dir_shift_hist,
                              labels = "AUTO",
                              align = "h",
                              nrow = 1,
                              axis = "tblr",
                              rel_widths = c(1.5, 1.5, 1.5),
                              rel_heights = c(0.25, 5, 5))

pdf(file = file.path(plots_dir, "range_shift_direction", paste0("BBS_rshift_direction_hist", 
                                                                ifelse(all(hist_years == 1980:1983), "81-83", "88-90"),
                                                                ".pdf")),
    height = 3, width = 8)
p_climate_change
dev.off()
