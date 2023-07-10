# explore the results of the niche and range analysis:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(cowplot)
library(ade4)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

datashare_EBCC <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "EBCC")
data_dir <- file.path("data", "EBBA_analysis")
plots_dir <- file.path("plots")

## species-specific background or whole EBBA area as background:
bg_spec <- TRUE # set to FALSE for whole EBBA area as environmental background

# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# species selection:
load(file = file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_final.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
sel_species <- sort(unique(EBBA1_prep$species))

# EBBA 1, only selected species:
EBBA1 <- read_sf(file.path(data_dir, "EBBA1_change.shp")) %>% # output of 1_EBBA_prep_data.R
  filter(species %in% sel_species)

# EBBA 2, only selected species:
EBBA2 <- read_sf(file.path(data_dir, "EBBA2_change.shp")) %>% # output of 1_EBBA_prep_data.R
  filter(species %in% sel_species)

# comparable EBBA cells (centroids):
EBBA_cells <- read_sf(file.path(data_dir, "EBBA_change.shp")) %>% # output of 1_EBBA_prep_data.R
  dplyr::select(-species) %>% 
  distinct(cell50x50, .keep_all = TRUE) # keep geometry; same when using EBBA2

# comparable EBBA cells (grid):
# EBBA_change_grid <- read_sf(file.path(datashare_EBCC, "EBBA_change", "ebba2_grid50x50_change_v1.shp")) %>% 
#   rename(cell50x50 = cell50x50_)
EBBA_change_grid <- read_sf(file.path("data", "ebba2_grid50x50_change_v1", "ebba2_grid50x50_change_v1.shp")) %>% # EPSG:3035   
  rename(cell50x50 = cell50x50_)

# results of niche and range shift analyses:

# load analyses results (output of 4_EBBA_niche_shift_analysis.R and 4_EBBA_range_shift_analysis.R)
niche_results <- read.csv(file.path(data_dir, paste0("EBBA_niche_shift_results_bg_", ifelse(bg_spec, "spec", "EBBA"), "_070723.csv")))
range_results <- read.csv(file.path(data_dir, paste0("EBBA_range_shift_results_bg_", ifelse(bg_spec, "spec", "EBBA"), "_070723.csv")))

# -------------------------------- #
#     Plots and explorations:   ####
# -------------------------------- #

## climate change as PCA based biplot and maps: --------------------------------

# environmental background = all comparable EBBA cells

# EBBA 1:

# climate data historic period:
biovars_hist_rast <- rast(list.files(file.path(data_dir, "Bioclim_1984_1988"),
                                     pattern = "50km.tif$",
                                     full.names = TRUE))
# join EBBA cells and corresponding climate data:
EBBA1_climate_data <- EBBA_cells %>% 
  vect %>% # convert to terra object
  cbind(terra::extract(biovars_hist_rast, y = .)) %>%
  as.data.frame %>% 
  select(-ID)

# EBBA 2:

# climate data recent period:
biovars_rec_rast <- rast(list.files(file.path(data_dir, "Bioclim_2012_2017"),
                                    pattern = "50km.tif$",
                                    full.names = TRUE))
# join EBBA cells and corresponding climate data:
EBBA2_climate_data <- EBBA_cells %>% # same as EBBA2 cells
  vect %>% # convert to terra object
  cbind(terra::extract(biovars_rec_rast, y = .)) %>%
  as.data.frame %>% 
  select(-ID)


# assess climate niche by using the first 2 PCA axes:
# calibrating the PCA in the whole study area:
pca.env <- dudi.pca(rbind(EBBA1_climate_data, EBBA2_climate_data)[, paste0("bio", 1:19)],
                    scannf = FALSE, nf = 2) # number of axes

#ecospat::ecospat.plot.contrib(contrib=pca.env$co, eigen=pca.env$eig)

# biplot:

# p_biplot1 <- factoextra::fviz_pca_var(pca.env,
#                                       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                                       repel = TRUE,
#                                       title = "EBBA PCA Biplot") # avoid text overlapping

p_biplot2 <- factoextra::fviz_pca_biplot(pca.env, repel = TRUE,
                                         col.var = "black", # variables color
                                         col.ind = c(rep("historic", nrow(EBBA1_climate_data)),
                                                     rep("recent", nrow(EBBA2_climate_data))),  # individuals colour
                                         label = "var",
                                         palette = c( "red", "blue"),
                                         pointsize = 1.5,
                                         alpha.ind = 0.1,
                                         legend.title = "",
                                         title = "",
                                         legend = c(0.9, 0.9))

# maps:

# PCA scores for EBBA1:
scores.clim.EBBA1 <- suprow(pca.env, EBBA1_climate_data[, paste0("bio", 1:19)])$li
# PCA scores for EBBA2:
scores.clim.EBBA2 <- suprow(pca.env, EBBA2_climate_data[, paste0("bio", 1:19)])$li

# difference in first 2 PCA axes per cell:
map_df <- data.frame("cell50x50" = EBBA1_climate_data$cell50x50, 
                     "diff_PC1" = scores.clim.EBBA2$Axis1 - scores.clim.EBBA1$Axis1, # change
                     "diff_PC2" = scores.clim.EBBA2$Axis2 - scores.clim.EBBA1$Axis2) # change

# join to spatial data:
map_sf <- EBBA_change_grid %>% 
  left_join(map_df, by = c("cell50x50" = "cell50x50"))

# plots:
plot(map_sf["diff_PC1"])
plot(map_sf["diff_PC2"])
plot(map_sf["diff_PC1"], axes = TRUE, graticule = TRUE)

world <- rnaturalearth::ne_coastline(returnclass = "sf") %>% 
  st_transform("EPSG:3035")

limit <- max(abs(map_sf$diff_PC1)) * c(-1, 1) # to center colour scale

p_diff_pc1 <- ggplot() +
  geom_sf(data = map_sf, aes(fill = diff_PC1), colour = NA) +
  scale_fill_distiller("PC1\n(recent - historic)", type = "div", 
                       palette = "RdBu", limit = limit) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000))

p_diff_pc2 <- ggplot() +
  geom_sf(data = map_sf, aes(fill = diff_PC2), colour = NA) +
  scale_fill_distiller("PC2\n(recent - historic)", type = "div", 
                       palette = "RdBu", limit = limit) +
  geom_sf(data = world, colour = "gray20", fill = NA) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000))

# combine biplot and difference maps in one plot:
p_climate_change <- plot_grid(p_biplot2, p_diff_pc1, p_diff_pc2,
                              labels = "AUTO",
                              align = "h",
                              nrow = 1,
                              axis = "tblr",
                              rel_widths = c(1.5, 1.5, 1.5),
                              rel_heights = c(0.25, 5, 5)
)
p_climate_change

# save plot:
pdf(file = file.path(plots_dir, "climate_change_PCAs", "EBBA_climate_change_PCA_plots.pdf"), height = 5, width = 20)
p_climate_change
dev.off()


## species richness maps: ------------------------------------------------------

# number of species recorded in each EBBA cell:

## EBBA 1:
EBBA1_richness <- EBBA1 %>% 
  group_by(cell50x50) %>% 
  summarise(EBBA1_richness = n(), .groups = "keep") %>% 
  ungroup

## EBBA 2:
EBBA2_richness <- EBBA2 %>% 
  group_by(cell50x50) %>% 
  summarise(EBBA2_richness = n(), .groups = "keep") %>% 
  ungroup

# join EBBA change grid and richness by cell ID:
EBBA_richness_grid <- EBBA_change_grid %>% # EBBA change grid
  left_join(st_drop_geometry(EBBA1_richness)) %>% 
  left_join(st_drop_geometry(EBBA2_richness), by = "cell50x50", suffix = c("_EBBA1", "_EBBA2")) %>% 
  # difference EBBA 1 and EBBA 2 richness:
  mutate(diff = EBBA2_richness - EBBA1_richness)

# plots:

# world <- rnaturalearth::ne_countries(returnclass = "sf") %>% 
#   st_transform("EPSG:3035")
world <- rnaturalearth::ne_coastline(returnclass = "sf") %>% 
  st_transform("EPSG:3035")

# EBBA 1 richness:
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

# center colour scale at zero:
limit <- max(abs(EBBA_richness_grid$diff), na.rm = TRUE) * c(-1, 1)

# difference species richness EBBA2-EBBA1:
diffplot <- ggplot() +
  geom_sf(data = EBBA_richness_grid, aes(fill = diff), colour = NA) +
  geom_sf(data = world, colour = "gray20", fill = NA, linewidth = 0.01) +
  lims(x = c(1100000, 6500000), y = c(1500000, 6400000)) +
  theme_dark() +
  scale_fill_distiller("EBBA 2 - EBBA 1", type = "div", 
                       palette = "RdBu", limit = limit,
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
pdf(file = file.path(plots_dir, "species_richness", "EBBA_species_richness.pdf"), height = 5, width = 10)
sr
dev.off()

# combine species richness and differences:
sr_diff <- plot_grid(sr, diffplot, labels = "AUTO",
                     align = "h",
                     axis = "tblr",
                     rel_widths = c(2, 1.3),
                     vjust = 7)

# save plot:
pdf(file = file.path(plots_dir, "species_richness", "EBBA_species_richness_incl_diff.pdf"), height = 5, width = 12)
sr_diff
dev.off()


## dynamics boxplot: -----------------------------------------------------------

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
) # rather use non-analogue conditions?

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
           x = 0.5, y = 118, hjust = 0, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_sign_higher, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) +
  geom_hline(yintercept = -5) +
  annotate(geom = "text", label = "Species (%) with significantly lower metrics than expected by chance:", 
           x = 0.5, y = -14, hjust = 0, vjust = 0, size = 3) +
  geom_text(data = labelsdat, aes(label = spec_sign_lower, y = ypos_lower),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3)

# save plot:
pdf(file = file.path(plots_dir, "dynamics_boxplots", paste0("EBBA_boxplot_dynamics_bg_", ifelse(bg_spec, "spec", "EBBA"), ".pdf")),
    height = 3, width = 8)
p
dev.off()


## dynamics correlation plots: -------------------------------------------------

# merge niche and range shift results:
metrics_df <- range_results %>% 
  left_join(niche_results, by = c("species" = "species")) %>% 
  select(range_stability_std, range_unfilling_std, range_expansion_std,
         niche_stability_std, niche_unfilling_std, niche_expansion_std)

cor1 <- with(metrics_df,cor.test(range_expansion_std, niche_unfilling_std))
xrange <- range(metrics_df$niche_unfilling_std)
yrange <- range(metrics_df$range_expansion_std)
p1 <- ggplot(data = metrics_df, aes(y = range_expansion_std, x = niche_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche U") + ylab("range E") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor1$estimate,2),ifelse(cor1$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")
cor2 <- with(metrics_df,cor.test(range_expansion_std, niche_stability_std))
xrange <- range(metrics_df$niche_stability_std)
yrange <- range(metrics_df$range_expansion_std)
p2 <- ggplot(data = metrics_df, aes(y = range_expansion_std, x = niche_stability_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche S") + ylab("range E") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor2$estimate,2),ifelse(cor2$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")
cor3 <- with(metrics_df,cor.test(range_expansion_std, niche_expansion_std))
xrange <- range(metrics_df$niche_expansion_std)
yrange <- range(metrics_df$range_expansion_std)
p3 <- ggplot(data = metrics_df, aes(y = range_expansion_std, x = niche_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche E") + ylab("range E") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor3$estimate,2),ifelse(cor3$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")

cor4 <- with(metrics_df,cor.test(range_stability_std, niche_unfilling_std))
xrange <- range(metrics_df$niche_unfilling_std)
yrange <- range(metrics_df$range_stability_std)
p4 <- ggplot(data = metrics_df, aes(y = range_stability_std, x = niche_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche U") + ylab("range S") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor4$estimate,2),ifelse(cor4$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")
cor5 <- with(metrics_df,cor.test(range_stability_std, niche_stability_std))
xrange <- range(metrics_df$niche_stability_std)
yrange <- range(metrics_df$range_stability_std)
p5 <- ggplot(data = metrics_df, aes(y = range_stability_std, x = niche_stability_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche S") + ylab("range S") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.95), label = paste0("r=",round(cor5$estimate,2),ifelse(cor5$p.value<0.05,"*","")), vjust=0, hjust=1, color="red")
cor6 <- with(metrics_df,cor.test(range_stability_std, niche_expansion_std))
xrange <- range(metrics_df$niche_expansion_std)
yrange <- range(metrics_df$range_stability_std)
p6 <- ggplot(data = metrics_df, aes(y = range_stability_std, x = niche_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche E") + ylab("range S") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor6$estimate,2),ifelse(cor6$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")

cor7 <- with(metrics_df,cor.test(range_unfilling_std, niche_unfilling_std))
xrange <- range(metrics_df$niche_unfilling_std)
yrange <- range(metrics_df$range_unfilling_std)
p7 <- ggplot(data = metrics_df, aes(y = range_unfilling_std, x = niche_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche U") + ylab("range U") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor7$estimate,2),ifelse(cor7$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")
cor8 <- with(metrics_df,cor.test(range_unfilling_std, niche_stability_std))
xrange <- range(metrics_df$niche_stability_std)
yrange <- range(metrics_df$range_unfilling_std)
p8 <- ggplot(data = metrics_df, aes(y = range_unfilling_std, x = niche_stability_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche S") + ylab("range U") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor8$estimate,2),ifelse(cor8$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")
cor9 <- with(metrics_df,cor.test(range_unfilling_std, niche_expansion_std))
xrange <- range(metrics_df$niche_expansion_std)
yrange <- range(metrics_df$range_unfilling_std)
p9 <- ggplot(data = metrics_df, aes(y = range_unfilling_std, x = niche_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche E") + ylab("range U") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor9$estimate,2),ifelse(cor9$p.value<0.05,"*","")), vjust=1, hjust=1, color="red")



pdf(file = file.path(plots_dir, "dynamics_correlations", paste0("EBBA_correlation_dynamics_bg_", 
                                                                ifelse(bg_spec, "spec", "EBBA"), ".pdf")),
    height = 5, width = 8)
plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9,
          align = "h", axis = "tblr", nrow = 3, ncol = 3)
dev.off()


## direction of range shifts: --------------------------------------------------

# same for using species-specific background and whole EBBA area as background

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
pdf(file = file.path(plots_dir, "range_shift_direction", "EBBA_rshift_direction.pdf"),
    height = 5, width = 6)
dir_shift_p
dev.off()

pdf(file = file.path(plots_dir, "range_shift_direction", "EBBA_rshift_direction_log.pdf"),
    height = 5, width = 6)
dir_shift_p2
dev.off()

pdf(file = file.path(plots_dir, "range_shift_direction", "EBBA_rshift_direction_histogram.pdf"),
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

pdf(file = file.path(plots_dir, "range_shift_direction", "EBBA_rshift_direction.pdf"),
    height = 3, width = 8)
p_climate_change
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

## differences in bioclim variable values between historic and recent time period: ----

bioclim1_hist <- rast(file.path(data_dir, "Bioclim_1984_1988", "CHELSA_bio1_1984_1988_50km.tif"))
bioclim1_rec <- rast(file.path(data_dir, "Bioclim_2012_2017", "CHELSA_bio1_2012_2017_50km.tif"))
change_rast <- (bioclim1_hist - bioclim1_rec)/bioclim1_hist * 100 # change in percentage compared to historic time period

biovars_hist_rast <- rast(list.files(file.path(data_dir, "Bioclim_1984_1988"),
                                     pattern = "50km.tif$",
                                     full.names = TRUE))
biovars_rec_rast <- rast(list.files(file.path(data_dir, "Bioclim_2012_2017"),
                                    pattern = "50km.tif$",
                                    full.names = TRUE))
change_rast <- (biovars_hist_rast - biovars_rec_rast)/biovars_hist_rast * 100  # change in percent
# mask with EBBA comparable area:
EBBA_mask_sp <- as.polygons(rast(file.path(data_dir, "EBBA_mask_ecospat.tif")))
change_rast_masked <- mask(change_rast, EBBA_mask_sp)

# colors:
for(i in 1:nlyr(change_rast_masked)){
  ceil <- values(change_rast_masked[[i]], na.rm=TRUE) |> abs() |> max() |> ceiling() 
  pal <- leaflet::colorNumeric(palette = "RdBu", domain=c(-ceil, ceil), reverse = T)
  # plot:
  plot(change_rast_masked[[i]], range=c(-ceil, ceil), col=pal(seq(-ceil,ceil,.1)), 
       main = change_rast[[i]]@ptr$names)
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

res_df <- range_results
#res_df <- niche_results

# for range analysis:
summary(res_df$Eucl_dist)
summary(res_df$NS_shift)
summary(res_df$EW_shift)

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