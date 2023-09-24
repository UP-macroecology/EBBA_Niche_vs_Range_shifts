# explore the results of the EBBA and BBS niche and range analysis:

library(dplyr)
library(sf)
library(terra)
library(ggplot2)
library(tidyr)
library(cowplot)
library(ade4)
library(ggprism)
library(dotwhisker)
library(tidyverse)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
data_dir <- file.path("data")
data_dir_BBS <- file.path("data", "BBS_analysis")
data_dir_EBBA <- file.path("data", "EBBA_analysis")
plots_dir <- file.path("plots")

# which historic time period should be used for BBS:
hist_years <- 1980:1983 # maximum gap between historic and recent time period
# hist_years <- 1987:1990 # similar gap between historic and recent time period as in EBBA analysis

# environmental background: presences and absences within 500 km buffer around presences (TRUE) or all true absences within conterminous US (FALSE):
bg_spec <- TRUE

niche_results_BBS <- read.csv(file.path(data_dir_BBS, "BBS_niche_shift_results_bg_spec_hist81-83.csv"))
range_results_BBS <- read.csv(file.path(data_dir_BBS, "BBS_range_shift_results_bg_spec_hist81-83.csv"))
sel_species_BBS <- niche_results_BBS$species

niche_results_EBBA <- read.csv(file.path(data_dir_EBBA, "EBBA_niche_shift_results_bg_spec.csv"))
range_results_EBBA <- read.csv(file.path(data_dir_EBBA, "EBBA_range_shift_results_bg_spec.csv"))
sel_species_EBBA <- niche_results_EBBA$species


#################################
# Climatic niche coverage

stability_BBS <- read.csv(file.path(data_dir_BBS,"BBS_stability_PCA_contUS_BL22.csv"))
stability_EBBA <- read.csv(file.path(data_dir_EBBA,"species_stability_PCA_EBBA2_BL22.csv"))

stab_df <- data.frame(stability=c(stability_EBBA$stability, stability_BBS$stability),
                      Region=c(rep("Europe",nrow(stability_EBBA)),rep("US",nrow(stability_BBS))))

p <- ggplot(stab_df, aes(x = Region, y = stability*100, fill = Region)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  theme_bw() +
  # scale_fill_viridis_d("Region") +
  theme(panel.grid = element_blank()) +
  ylab("Climatic niche coverage [%]") +
  xlab("") 
  
pdf(file = file.path(plots_dir, "coverage_climatic_niche", "EBBA_BBS_boxplot_niche_coverage.pdf"),
    height = 4, width = 4.5)
p
dev.off()



#--------------------------------

## Niche overlap boxplot: -----------------------------------------------------------

# species (%) with significant dynamic values:

## niche BBS:
niche_test_sign_BBS <- niche_results_BBS %>% 
  dplyr::select(c(species, matches("_p_D_"))) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_BBS)*100, 2)))

## range BBS:
range_test_sign_BBS <- range_results_BBS %>% 
  dplyr::select(c(species, matches("_p_D_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_BBS)*100, 2)))

## Niche EBBA:
niche_test_sign_EBBA <- niche_results_EBBA %>% 
  dplyr::select(c(species, matches("_p_D_"))) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_EBBA)*100, 2)))

## range EBBA:
range_test_sign_EBBA <- range_results_EBBA %>% 
  dplyr::select(c(species, matches("_p_D_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_EBBA)*100, 2)))

# species (%) with significantly lower overlap than by chance:
spec_sign_lower <- c(range_test_sign_EBBA$shift_p_D_A_n_sig, # range overlap
                      range_test_sign_BBS$shift_p_D_A_n_sig, # range overlap
                      niche_test_sign_EBBA$shift_p_D_A_n_sig, # niche overlap
                      niche_test_sign_BBS$shift_p_D_A_n_sig) # niche overlap

# species (%) with significantly higher overlap than by chance:
spec_sign_higher <- c(range_test_sign_EBBA$cons_p_D_A_n_sig, # range overlap
                     range_test_sign_BBS$cons_p_D_A_n_sig, # range overlap
                     niche_test_sign_EBBA$cons_p_D_A_n_sig, # niche overlap
                     niche_test_sign_BBS$cons_p_D_A_n_sig) # niche overlap
 
# join niche and range shift results, in long format:
niche_range_df <- data.frame(
  value=c(range_results_EBBA$D,range_results_BBS$D, niche_results_EBBA$D, range_results_BBS$D), 
  metric=c(rep("Range overlap",length(sel_species_BBS)+length(sel_species_EBBA)), rep ("Niche overlap", length(sel_species_BBS)+length(sel_species_EBBA))),
  Region=rep(c( rep("Europe",length(sel_species_EBBA)),rep("US",length(sel_species_BBS))),2)
  )

# data set for labels:
labelsdat <- tibble(metric = c(rep("Range overlap", 2), rep("Niche overlap", 2)),
                    Region = factor(rep(c("Europe","US"),2),levels = c("Europe","US")),
                    sign_higher = spec_sign_higher,
                    sign_lower = spec_sign_lower,
                    ypos_higher = 1.10,
                    ypos_lower = .15)


wilcox_niche <- wilcox.test(niche_results_EBBA$D,niche_results_BBS$D)$p.value
wilcox_range <- wilcox.test(range_results_EBBA$D,niche_results_BBS$D)$p.value
two.means.grouped1 <- tibble::tribble(
  ~group1, ~group2, ~p.signif,  ~y.position, ~metric,
  "Europe", "US", ifelse(wilcox_niche<0.001,"***",ifelse(wilcox_niche<0.01,"**",ifelse(wilcox_niche<0.05,"*",""))),  1, "Niche overlap",
  "Europe", "US", ifelse(wilcox_range<0.001,"***",ifelse(wilcox_range<0.01,"**",ifelse(wilcox_range<0.05,"*",""))), 1, "Range overlap"
)

# plot:
p <- ggplot(niche_range_df, aes(x = metric, y = value)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(aes(fill = Region),lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Schoener's D") +
  xlab("") +
  # scale_fill_viridis_d("Region") +
  scale_y_continuous(breaks = seq(.4, 1, .2), expand = expansion(add = .05)) +
  annotate(geom = "text", label = "Species (%) conserving niche or range:", 
           x = 0.5, y = 1.15, hjust = 0, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_sign_higher, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 1.05) +
  geom_hline(yintercept = .25) +
  annotate(geom = "text", label = "Species (%) switching niche or range:", 
           x = 0.5, y = .20, hjust = 0, vjust = 0, size = 3) +
  geom_text(data = labelsdat, aes(label = spec_sign_lower, y = ypos_lower),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) + 
  add_pvalue(two.means.grouped1, x="metric")


# save plot:
pdf(file = file.path(plots_dir, "dynamics_boxplots", paste0("EBBA_BBS_boxplot_overlap_bg_", ifelse(bg_spec, "spec", "EBBA"), ".pdf")),
    height = 4, width = 4.5)
p
dev.off()


#####################################################
#------------- Niche dynamic boxplots ---------------


## niche BBS:
niche_test_sign_BBS <- niche_results_BBS %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_BBS)*100, 2)))

## range BBS:
range_test_sign_BBS <- range_results_BBS %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_BBS)*100, 2)))

## Niche EBBA:
niche_test_sign_EBBA <- niche_results_EBBA %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_EBBA)*100, 2)))

## range EBBA:
range_test_sign_EBBA <- range_results_EBBA %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_EBBA)*100, 2)))


# BBS species (%) with significant niche tracking and range lagging
spec_cons_BBS <- c("", # niche abandonment
               niche_test_sign_BBS$cons_p_U_A_n_sig, # niche unfilling
               niche_test_sign_BBS$cons_p_S_A_n_sig, # niche stability
               niche_test_sign_BBS$cons_p_E_A_n_sig, # niche expansion
                      "", # niche pioneering,
               range_test_sign_BBS$cons_p_U_A_n_sig, # range unfilling
               range_test_sign_BBS$cons_p_S_A_n_sig, # range stabiliy
               range_test_sign_BBS$cons_p_E_A_n_sig # range expansion
)

# EBBA species (%) with significant niche tracking and range lagging
spec_cons_EBBA <- c("", # niche abandonment
                   niche_test_sign_EBBA$cons_p_U_A_n_sig, # niche unfilling
                   niche_test_sign_EBBA$cons_p_S_A_n_sig, # niche stability
                   niche_test_sign_EBBA$cons_p_E_A_n_sig, # niche expansion
                   "", # niche pioneering,
                   range_test_sign_EBBA$cons_p_U_A_n_sig, # range unfilling
                   range_test_sign_EBBA$cons_p_S_A_n_sig, # range stabiliy
                   range_test_sign_EBBA$cons_p_E_A_n_sig # range expansion
)

# join niche and range shift results, convert to long format:
niche_range_df_BBS <- niche_results_BBS %>% 
  left_join(range_results_BBS, by = "species", suffix = c("_niche", "_range")) %>% 
  pivot_longer(cols = ends_with("_std"),
               names_to = c("category", "metric"), names_pattern = "(.*)_(.*)_",
               values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("abandonment", "unfilling", "stability", "expansion", "pioneering")))

# join niche and range shift results, convert to long format:
niche_range_df_EBBA <- niche_results_EBBA %>% 
  left_join(range_results_EBBA, by = "species", suffix = c("_niche", "_range")) %>% 
  pivot_longer(cols = ends_with("_std"),
               names_to = c("category", "metric"), names_pattern = "(.*)_(.*)_",
               values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("abandonment", "unfilling", "stability", "expansion", "pioneering")))



# data set for labels:
labelsdat <- tibble(category = c(rep("niche", 5), rep("range", 3)),
                    metric = factor(c("abandonment", "unfilling", "stability", "expansion", "pioneering", "unfilling", "stability", "expansion"), 
                                    levels = c("abandonment", "unfilling", "stability", "expansion", "pioneering")),
                    spec_cons_BBS = spec_cons_BBS,
                    spec_cons_EBBA = spec_cons_EBBA,
                    ypos_higher = 110)

# plot:
p_EBBA <- ggplot(niche_range_df_EBBA, aes(x = category, y = value*100, fill = metric)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Metrics [%]") +
  xlab("") +
  scale_fill_viridis_d("Metrics") +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(add = 10)) +
  annotate(geom = "text", label = "% species niche tracking:", 
           x = 1, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  annotate(geom = "text", label = "% species range lagging:", 
           x = 2, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_cons_EBBA, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) 

p_BBS <- ggplot(niche_range_df_BBS, aes(x = category, y = value*100, fill = metric)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Metrics [%]") +
  xlab("") +
  scale_fill_viridis_d("Metrics") +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(add = 10)) +
  annotate(geom = "text", label = "% species niche tracking:", 
           x = 1, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  annotate(geom = "text", label = "% species range lagging:", 
           x = 2, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_cons_BBS, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) 


# save plot:
pdf(file = file.path("plots", "dynamics_boxplots", "EBBA_boxplot_dynamics_bg_spec.pdf"),
    height = 3, width = 8)
p_EBBA
dev.off()

pdf(file = file.path("plots", "dynamics_boxplots", "BBS_boxplot_dynamics_bg_spec_hist81-83.pdf"),
    height = 3, width = 8)
p_BBS
dev.off()


#------------- Niche dynamic boxplots ---------------
#------------ in analogue climate space -------------


## niche BBS:
niche_test_sign_BBS <- niche_results_BBS %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_BBS)*100, 2)))

## range BBS:
range_test_sign_BBS <- range_results_BBS %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_BBS)*100, 2)))

## Niche EBBA:
niche_test_sign_EBBA <- niche_results_EBBA %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>%
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_EBBA)*100, 2)))

## range EBBA:
range_test_sign_EBBA <- range_results_EBBA %>% 
  dplyr::select(c(species, matches("cons_p_"))) %>% 
  summarise(across(.cols = where(is.numeric), 
                   .fns = ~length(which(.x <= 0.05)), 
                   .names = "{.col}_n_sig")) %>% 
  mutate(across(everything(), ~ round(.x/length(sel_species_EBBA)*100, 2)))


# BBS species (%) with significant niche tracking and range lagging
spec_cons_BBS <- c(#"", # niche abandonment
                   niche_test_sign_BBS$cons_p_U_A_n_sig, # niche unfilling
                   niche_test_sign_BBS$cons_p_S_A_n_sig, # niche stability
                   niche_test_sign_BBS$cons_p_E_A_n_sig, # niche expansion
                   #"", # niche pioneering,
                   range_test_sign_BBS$cons_p_U_A_n_sig, # range unfilling
                   range_test_sign_BBS$cons_p_S_A_n_sig, # range stabiliy
                   range_test_sign_BBS$cons_p_E_A_n_sig # range expansion
)

# EBBA species (%) with significant niche tracking and range lagging
spec_cons_EBBA <- c(#"", # niche abandonment
                    niche_test_sign_EBBA$cons_p_U_A_n_sig, # niche unfilling
                    niche_test_sign_EBBA$cons_p_S_A_n_sig, # niche stability
                    niche_test_sign_EBBA$cons_p_E_A_n_sig, # niche expansion
                    #"", # niche pioneering,
                    range_test_sign_EBBA$cons_p_U_A_n_sig, # range unfilling
                    range_test_sign_EBBA$cons_p_S_A_n_sig, # range stabiliy
                    range_test_sign_EBBA$cons_p_E_A_n_sig # range expansion
)

# rescale
niche_results_BBS_an <- niche_results_BBS %>% 
  select(-c("niche_abandonment_std","niche_pioneering_std"))
sum_analog_BBS <- niche_results_BBS_an %>% 
  select(starts_with("niche")) %>%
  rowSums()
niche_results_BBS_an$niche_expansion_std <- niche_results_BBS_an$niche_expansion_std/sum_analog_BBS
niche_results_BBS_an$niche_unfilling_std <- niche_results_BBS_an$niche_unfilling_std/sum_analog_BBS
niche_results_BBS_an$niche_stability_std <- niche_results_BBS_an$niche_stability_std/sum_analog_BBS

niche_results_EBBA_an <- niche_results_EBBA %>% 
  select(-c("niche_abandonment_std","niche_pioneering_std"))
sum_analog_EBBA <- niche_results_EBBA_an %>% 
  select(starts_with("niche")) %>%
  rowSums()
niche_results_EBBA_an$niche_expansion_std <- niche_results_EBBA_an$niche_expansion_std/sum_analog_EBBA
niche_results_EBBA_an$niche_unfilling_std <- niche_results_EBBA_an$niche_unfilling_std/sum_analog_EBBA
niche_results_EBBA_an$niche_stability_std <- niche_results_EBBA_an$niche_stability_std/sum_analog_EBBA

# join niche and range shift results, convert to long format:
niche_range_df_BBS <- niche_results_BBS_an %>% 
  left_join(range_results_BBS, by = "species", suffix = c("_niche", "_range")) %>% 
  pivot_longer(cols = ends_with("_std"),
               names_to = c("category", "metric"), names_pattern = "(.*)_(.*)_",
               values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("unfilling", "stability", "expansion")))

# join niche and range shift results, convert to long format:
niche_range_df_EBBA <- niche_results_EBBA_an %>% 
  left_join(range_results_EBBA, by = "species", suffix = c("_niche", "_range")) %>% 
  pivot_longer(cols = ends_with("_std"),
               names_to = c("category", "metric"), names_pattern = "(.*)_(.*)_",
               values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("unfilling", "stability", "expansion")))


# data set for labels:
labelsdat <- tibble(category = c(rep("niche", 3), rep("range", 3)),
                    metric = factor(rep(c("unfilling", "stability", "expansion"),2), 
                                    levels = c("abandonment", "unfilling", "stability", "expansion", "pioneering")),
                    spec_cons_BBS = spec_cons_BBS,
                    spec_cons_EBBA = spec_cons_EBBA,
                    ypos_higher = 110)

# plot:
p_EBBA <- ggplot(niche_range_df_EBBA, aes(x = category, y = value*100, fill = metric)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Metrics [%]") +
  xlab("") +
  scale_fill_viridis_d("Metrics") +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(add = 10)) +
  annotate(geom = "text", label = "% species niche tracking:", 
           x = 1, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  annotate(geom = "text", label = "% species range lagging:", 
           x = 2, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_cons_EBBA, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) 

p_BBS <- ggplot(niche_range_df_BBS, aes(x = category, y = value*100, fill = metric)) + # aes given inside ggplot() necessary for geom_text
  geom_boxplot(lwd = 0.1, outlier.size = 0.7, outlier.colour = "grey30", width = 0.9,
               position = position_dodge2(width = 1, preserve = "single")) +
  #facet_grid(~ category, scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  ylab("Metrics [%]") +
  xlab("") +
  scale_fill_viridis_d("Metrics") +
  scale_y_continuous(breaks = seq(0, 100, 20), expand = expansion(add = 10)) +
  annotate(geom = "text", label = "% species niche tracking:", 
           x = 1, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  annotate(geom = "text", label = "% species range lagging:", 
           x = 2, y = 118, hjust = 0.5, vjust = 0, size = 3) + # xx
  geom_text(data = labelsdat, aes(label = spec_cons_BBS, y = ypos_higher),
            position = position_dodge2(width = 0.9, preserve = "single"), size = 3) +
  geom_hline(yintercept = 105) 


# save plot:
pdf(file = file.path("plots", "dynamics_boxplots", "EBBA_boxplot_dynamics_bg_spec_analogClim.pdf"),
    height = 3, width = 5.5)
p_EBBA
dev.off()

pdf(file = file.path("plots", "dynamics_boxplots", "BBS_boxplot_dynamics_bg_spec_hist81-83_analogClim.pdf"),
    height = 3, width = 5.5)
p_BBS
dev.off()




#################################
# Correlation between niche and range metrics

metrics_BBS_df <- range_results_BBS %>% 
  rename(r_D=D) %>%
  left_join(niche_results_BBS, by = c("species" = "species")) %>% 
  rename(n_D=D) %>%
  dplyr::select(r_D, n_D, range_stability_std, range_unfilling_std, range_expansion_std,
         niche_stability_std, niche_unfilling_std, niche_expansion_std)

metrics_EBBA_df <- range_results_EBBA %>% 
  rename(r_D=D) %>%
  left_join(niche_results_EBBA, by = c("species" = "species")) %>% 
  rename(n_D=D) %>%
  dplyr::select(r_D, n_D, range_stability_std, range_unfilling_std, range_expansion_std,
         niche_stability_std, niche_unfilling_std, niche_expansion_std)

cor1 <- with(metrics_BBS_df,cor.test(r_D, n_D))
xrange <- range(metrics_BBS_df$n_D)
yrange <- range(metrics_BBS_df$r_D)
p1 <- ggplot(data = metrics_BBS_df, aes(y = r_D, x = n_D)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche D") + ylab("range D") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.95), label = paste0("r=",round(cor1$estimate,2),ifelse(cor1$p.value<0.05,"*","")), vjust=0, hjust=1)

cor2 <- with(metrics_EBBA_df,cor.test(r_D, n_D))
xrange <- range(metrics_EBBA_df$n_D)
yrange <- range(metrics_EBBA_df$r_D)
p2 <- ggplot(data = metrics_EBBA_df, aes(y = r_D, x = n_D)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("niche D") + ylab("range D") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.95), label = paste0("r=",round(cor2$estimate,2),ifelse(cor2$p.value<0.05,"*","")), vjust=0, hjust=1)

pdf(file = file.path(plots_dir, "dynamics_correlations", "dynamics_correlations_BBS_EBBA_SchoenersD_spec.pdf"),
    height = 3, width = 5)
plot_grid(p2, p1,
          align = "h",
          labels = "AUTO",
          axis = "tblr")
dev.off()



####

# correlation between range metrics

cor1 <- with(metrics_EBBA_df,cor.test(range_stability_std, range_unfilling_std))
xrange <- range(metrics_EBBA_df$range_unfilling_std)
yrange <- range(metrics_EBBA_df$range_stability_std)
p1 <- ggplot(data = metrics_EBBA_df, aes(y = range_stability_std, x = range_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("range unfilling") + ylab("range stability") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor1$estimate,2),ifelse(cor1$p.value<0.05,"*","")), vjust=1, hjust=1)

cor2 <- with(metrics_EBBA_df,cor.test(range_stability_std, range_expansion_std))
xrange <- range(metrics_EBBA_df$range_expansion_std)
yrange <- range(metrics_EBBA_df$range_stability_std)
p2 <- ggplot(data = metrics_EBBA_df, aes(y = range_stability_std, x = range_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("range expansion") + ylab("range stability") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor2$estimate,2),ifelse(cor2$p.value<0.05,"*","")), vjust=1, hjust=1)

cor3 <- with(metrics_EBBA_df,cor.test(range_unfilling_std, range_expansion_std))
xrange <- range(metrics_EBBA_df$range_expansion_std)
yrange <- range(metrics_EBBA_df$range_unfilling_std)
p3 <- ggplot(data = metrics_EBBA_df, aes(y = range_unfilling_std, x = range_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("range expansion") + ylab("range unfilling") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor3$estimate,2),ifelse(cor3$p.value<0.05,"*","")), vjust=1, hjust=1)

cor4 <- with(metrics_BBS_df,cor.test(range_stability_std, range_unfilling_std))
xrange <- range(metrics_BBS_df$range_unfilling_std)
yrange <- range(metrics_BBS_df$range_stability_std)
p4 <- ggplot(data = metrics_BBS_df, aes(y = range_stability_std, x = range_unfilling_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("range unfilling") + ylab("range stability") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor4$estimate,2),ifelse(cor4$p.value<0.05,"*","")), vjust=1, hjust=1)

cor5 <- with(metrics_BBS_df,cor.test(range_stability_std, range_expansion_std))
xrange <- range(metrics_BBS_df$range_expansion_std)
yrange <- range(metrics_BBS_df$range_stability_std)
p5 <- ggplot(data = metrics_BBS_df, aes(y = range_stability_std, x = range_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("range expansion") + ylab("range stability") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor5$estimate,2),ifelse(cor5$p.value<0.05,"*","")), vjust=1, hjust=1)

cor6 <- with(metrics_BBS_df,cor.test(range_unfilling_std, range_expansion_std))
xrange <- range(metrics_BBS_df$range_expansion_std)
yrange <- range(metrics_BBS_df$range_unfilling_std)
p6 <- ggplot(data = metrics_BBS_df, aes(y = range_unfilling_std, x = range_expansion_std)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_bw() + xlab("range expansion") + ylab("range unfilling") +
  geom_text(x=xrange[2]-(diff(xrange)*0.05), y = yrange[2]-(diff(yrange)*0.05), label = paste0("r=",round(cor6$estimate,2),ifelse(cor6$p.value<0.05,"*","")), vjust=1, hjust=1)


pdf(file = file.path(plots_dir, "dynamics_correlations", "dynamics_correlations_BBS_EBBA_rangemetrics_spec.pdf"),
    height = 3, width = 7)
plot_grid(p1, p2, p3, p4, p5, p6,
          align = "h",
          axis = "tblr",
          nrow=2)
dev.off()


#--------------------------------------------
# trait analyses - plot effect sizes:

trait_res_BBS <- read.csv(file.path(data_dir,"results_TraitAnal_df_US.csv"))
trait_res_EBBA <- read.csv(file.path(data_dir,"results_TraitAnal_df_Europe.csv"))
traits_res <- rbind(data.frame(trait_res_BBS, region='US'),data.frame(trait_res_EBBA, region='Europe'))

covariates <- c("Mass", "Hand.Wing.Index", "Trophic.Level", "Habitat.Density", "Range.Size", "Migration", "Centroid.Latitude", "niche_breadth_zcor")


traits_res_tb <- traits_res %>%
  as_tibble() %>% 
  # mutate(term=Trait, estimate=rangeD_coef, std.error=rangeD_stderr, statistic=rangeD_coef/rangeD_stderr, p.value=rangeD_p, model=region, varimp=rangeD_varimp) %>%
  # mutate(term=Trait, estimate=range_stab_coef, std.error=range_stab_stderr, statistic=range_stab_coef/range_stab_stderr, p.value=range_stab_p, model=region, varimp=range_stab_varimp) %>%
  # mutate(term=Trait, estimate=range_unf_coef, std.error=range_unf_stderr, statistic=range_unf_coef/range_unf_stderr, p.value=range_unf_p, model=region, varimp=range_unf_varimp) %>%
  # mutate(term=Trait, estimate=range_exp_coef, std.error=range_exp_stderr, statistic=range_exp_coef/range_exp_stderr, p.value=range_exp_p, model=region, varimp=range_exp_varimp) %>%
  mutate(term=Trait, estimate=nicheD_coef, std.error=nicheD_stderr, statistic=nicheD_coef/nicheD_stderr, p.value=nicheD_p, model=region, varimp=nicheD_varimp) %>%
  # mutate(term=Trait, estimate=niche_stab_coef, std.error=niche_stab_stderr, statistic=niche_stab_coef/niche_stab_stderr, p.value=niche_stab_p, model=region, varimp=niche_stab_varimp) %>%
  # mutate(term=Trait, estimate=niche_unf_coef, std.error=niche_unf_stderr, statistic=niche_unf_coef/niche_unf_stderr, p.value=niche_unf_p, model=region, varimp=niche_unf_varimp) %>%
  # mutate(term=Trait, estimate=niche_exp_coef, std.error=niche_exp_stderr, statistic=niche_exp_coef/niche_exp_stderr, p.value=niche_exp_p, model=region, varimp=niche_exp_varimp) %>%
  select(term, estimate, std.error, statistic, p.value, model, varimp) %>% 
  filter(term %in% covariates)
varimp_US <- traits_res_tb %>% filter(model=="US") %>% select(varimp)
varimp_Europe <- traits_res_tb %>% filter(model=="Europe") %>% select(varimp)
  
quartz(w=5,h=3)
dwplot(traits_res_tb, 
       vline = geom_vline(xintercept = 0, colour = "grey60", linetype = 2)
       ) %>% 
  relabel_predictors(
    c(
      Mass = "Body mass",
      Hand.Wing.Index = "HWI",
      Trophic.Level = "Trophic level",
      Habitat.Density = "Habitat openness",
      Range.Size = "Range size",
      Migration = "Migration",
      Centroid.Latitude = "Breeding latitude",
      niche_breadth_zcor = "Niche breadth"
    )) + 
  # ggtitle("Range overlap - Schoener's D") +
  # ggtitle("Range stability") +
  # ggtitle("Range unfilling") +
  # ggtitle("Range expansion") +
  ggtitle("Niche overlap - Schoener's D") +
  # ggtitle("Niche stability") +
  # ggtitle("Niche unfilling") +
  # ggtitle("Niche expansion") +
  theme_bw(base_size = 12) + 
  xlab("Coefficient Estimate") + 
  ylab("") + 
  scale_colour_hue(
    name = "Region") +
  # annotate("text", y = 8:1-.2,
  #          x = max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error), na.rm=T) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
  #          label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), ""),
  #          hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[8]-.2, 
           x = traits_res_tb$estimate[8] + (1.96*traits_res_tb$std.error[8]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[8], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[7]-.2, 
           x = traits_res_tb$estimate[7] + (1.96*traits_res_tb$std.error[7]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[7], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[6]-.2, 
           x = traits_res_tb$estimate[6] + (1.96*traits_res_tb$std.error[6]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[6], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[5]-.2, 
           x = traits_res_tb$estimate[5] + (1.96*traits_res_tb$std.error[5]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[5], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[4]-.2, 
           x = traits_res_tb$estimate[4] + (1.96*traits_res_tb$std.error[4]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[4], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[3]-.2, 
           x = traits_res_tb$estimate[3] + (1.96*traits_res_tb$std.error[3]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[3], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[2]-.2, 
           x = traits_res_tb$estimate[2] + (1.96*traits_res_tb$std.error[2]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[2], hjust=1, size=3, colour="#F8766D") +
  annotate("text", y = c(8:1)[1]-.2, 
           x = traits_res_tb$estimate[1] + (1.96*traits_res_tb$std.error[1]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_US[[1]])), paste0(round(as.numeric(varimp_US[[1]])*100,0),"%"), "")[1], hjust=1, size=3, colour="#F8766D") +
  # annotate("text", y = 8:1+.2, 
  #          x = max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error), na.rm=T) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
  #          label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), ""), hjust=1, size=3, colour="#00BFC4")
  annotate("text", y = c(8:1)[8]+.2, 
           x = traits_res_tb$estimate[8+8] + (1.96*traits_res_tb$std.error[8+8]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[8], hjust=1, size=3, colour="#00BFC4") +
annotate("text", y = c(8:1)[7]+.2, 
         x = traits_res_tb$estimate[8+7] + (1.96*traits_res_tb$std.error[8+7]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
         label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[7], hjust=1, size=3, colour="#00BFC4") +
  annotate("text", y = c(8:1)[6]+.2, 
           x = traits_res_tb$estimate[8+6] + (1.96*traits_res_tb$std.error[8+6]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[6], hjust=1, size=3, colour="#00BFC4") +
  annotate("text", y = c(8:1)[5]+.2, 
           x = traits_res_tb$estimate[8+5] + (1.96*traits_res_tb$std.error[8+5]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[5], hjust=1, size=3, colour="#00BFC4") +
  annotate("text", y = c(8:1)[4]+.2, 
           x = traits_res_tb$estimate[8+4] + (1.96*traits_res_tb$std.error[8+4]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[4], hjust=1, size=3, colour="#00BFC4") +
  annotate("text", y = c(8:1)[3]+.2, 
           x = traits_res_tb$estimate[8+3] + (1.96*traits_res_tb$std.error[8+3]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[3], hjust=1, size=3, colour="#00BFC4") +
  annotate("text", y = c(8:1)[2]+.2, 
           x = traits_res_tb$estimate[8+2] + (1.96*traits_res_tb$std.error[8+2]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[2], hjust=1, size=3, colour="#00BFC4") +
  annotate("text", y = c(8:1)[1]+.2, 
           x = traits_res_tb$estimate[8+1] + (1.96*traits_res_tb$std.error[8+1]) + diff(c(min(traits_res_tb$estimate - (1.96*traits_res_tb$std.error),na.rm=T), max(traits_res_tb$estimate + (1.96*traits_res_tb$std.error),na.rm=T)))*.15, 
           label = ifelse(!is.na(as.numeric(varimp_Europe[[1]])), paste0(round(as.numeric(varimp_Europe[[1]])*100,0),"%"), "")[1], hjust=1, size=3, colour="#00BFC4")
  
