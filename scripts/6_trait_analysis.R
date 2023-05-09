# explorative plots for trait analysis:

library(dplyr)
library(ggplot2)

# project data:
data_dir <- file.path("Data")

# use EBBA or BBS results:
dataset <- "BBS"
#dataset <- "EBBA" 

# use results of background = buffer around presences (TRUE) (or background = whole EBBA-area / conterminous US):
env_background_species_specific <- TRUE

# folder to store plots:
plots_dir <- file.path("plots", "traits_explorations", dataset)
if(!dir.exists(plots_dir)){dir.create(plots_dir)}


# functions: -------------------------------------------------------------------

# function to get number of observations:
give.n <- function(x){
  return(c(y = 1, label = length(x))) 
}

# function to print and save boxplots:
boxplot_fun <- function(data, x, y, ylab, title, ylim, filename, height = 5, width = 5){
  p <- ggplot(data = data, aes(x = {{x}}, y = {{y}})) + 
    geom_boxplot() +
    theme_bw() +
    ylab(ylab) +
    ylim(ylim) +
    stat_summary(fun.data = give.n, geom = "text") +
    ggtitle(title)
  print(p)
  # save plot:
  pdf(file = file.path(plots_dir, filename),
      height = height, width = width)
  print(p)
  dev.off()
}

# function to print and save scatterplots:
scatterplot_fun <- function(data, x, y, ylab, xlab, title,
                            filename, height = 5, width = 8, ylim = c(0,1)){
  p <- ggplot(data = data, aes(x = {{x}}, y = {{y}})) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_bw() +
    ylim(ylim) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title)
  print(p)
  # save plot:
  pdf(file = file.path(plots_dir, filename),
      height = height, width = width)
  print(p)
  dev.off()
}


# load data: -------------------------------------------------------------------

if(dataset == "EBBA"){
  
  avonet <- read.csv(file = file.path(data_dir, "AVONET_EBBA_species.csv"))
  
  # niche breadth:
  nb <- read.csv(file = file.path(data_dir, "EBBA_niche_breadth.csv"))
  
  
  # load analyses results (output of 4_EBBA_niche_shift_analysis.R and 4_EBBA_range_shift_analysis.R)
  if(env_background_species_specific){
    niche_results <- read.csv(file = file.path(data_dir, "niche_shift_results_env_species_spec_change.csv"))
    range_results <- read.csv(file = file.path(data_dir, "range_shift_results_env_species_spec_change.csv"))
  } else {
    niche_results <- read.csv(file = file.path(data_dir, "niche_shift_results_env_EBBA_change.csv"))
    range_results <- read.csv(file = file.path(data_dir, "range_shift_results_env_EBBA_change.csv"))
  }
  
  # join analysis results and trait data:
  res_traits_df <- avonet %>% 
    left_join(niche_results, by = c("Species1" = "species")) %>% 
    left_join(range_results, by = c("Species1" = "species"), suffix = c("_niche", "_range")) %>% 
    left_join(nb, by = c("Species1" = "species")) %>% 
    mutate(diff_range_niche_D = D_range - D_niche)
  # missing data for 6 species!
  
} else {
  
  avonet <- read.csv(file = file.path(data_dir, "AVONET_BBS_species.csv"))
  
  # niche breadth:
  nb <- read.csv(file = file.path(data_dir, "BBS_niche_breadth.csv"))
  
  # load analyses results (output of 4_BBS_niche_shift_analysis.R and 4_BBS_range_shift_analysis.R)
  if(env_background_species_specific){
    niche_results <- read.csv(file = file.path(data_dir, "BBS_niche_shift_results_env_600km_hist1995_1998.csv"))
    range_results <- read.csv(file = file.path(data_dir, "BBS_range_shift_results_env_600km_hist1995_1998.csv"))
  } else {
    niche_results <- read.csv(file = file.path(data_dir, "BBS_niche_shift_results_env_contUS_hist1995_1998.csv"))
    range_results <- read.csv(file = file.path(data_dir, "BBS_range_shift_results_env_contUS_hist1995_1998.csv"))
  } # xx other time period!
  
  # join analysis results and trait data:
  res_traits_df <- avonet %>% 
    left_join(niche_results, by = c("BBS_species" = "species")) %>% 
    left_join(range_results, by = c("BBS_species" = "species"), suffix = c("_niche", "_range")) %>% 
    left_join(nb, by = c("BBS_species" = "species")) %>% 
    mutate(diff_range_niche_D = D_range - D_niche)
}


# explorative plots:------------------------------------------------------------

## Niche overlap ~ range overlap:

scatterplot_fun(data = res_traits_df, x = D_niche, y = D_range,
                ylab = "range overlap D", xlab = "niche overlap D",
                title = paste0(dataset, ", background = buffer"),
                filename = "range_D~niche_D.pdf",
                width = 5, height = 5)
cor.test(res_traits_df$D_niche, res_traits_df$D_range, alternative = "two.sided",  method = "s")
# niche D and range D correlated!


## Niche breadth: ----

# range overlap:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = D_range,
                ylab = "range overlap D", xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), 
                filename = "range_D~Niche breadth.pdf",
                width = 8, height = 5)
# niche overlap:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = D_niche,
                ylab = "niche overlap D", xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), 
                filename = "niche_D~Niche breadth.pdf",
                width = 8, height = 5)

# niche expansion:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = niche_expansion_std,
                ylab = "niche E", ylim = c(0, 0.2), xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), filename = "niche_E~Niche breadth.pdf",
                width = 8)
# range expansion:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = range_expansion_std,
                ylab = "range E", ylim = c(0, 0.2), xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), filename = "range_E~Niche breadth.pdf")

# niche unfilling:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = niche_unfilling_std,
                ylab = "niche U", ylim = c(0, 0.2), xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), filename = "niche_U~Niche breadth.pdf",
                width = 8)
# range unfilling:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = range_unfilling_std,
                ylab = "range U", ylim = c(0, 0.2), xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), filename = "range_U~Niche breadth.pdf")

# niche stability:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = niche_stability_std,
                ylab = "niche S", ylim = c(0, 1), xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), filename = "niche_S~Niche breadth.pdf",
                width = 8)
# range stability:
scatterplot_fun(data = res_traits_df, x = niche_breadth_zcor, y = range_stability_std,
                ylab = "range S", ylim = c(0, 1), xlab = "Niche breadth",
                title = paste0(dataset, ", background = buffer"), filename = "range_S~Niche breadth.pdf")


## Trophic level:----

boxplot_fun(data = res_traits_df, x = Trophic.Level, y = D_range,
            ylab = "range overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_D~trophic_level.pdf")

boxplot_fun(data = res_traits_df, x = Trophic.Level, y = D_niche,
            ylab = "niche overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_D~trophic_level.pdf")

# difference range D - niche D:
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = diff_range_niche_D,
            ylab = "range D - niche D", ylim = range(res_traits_df$diff_range_niche_D), 
            title = paste0(dataset, ", background = buffer"), filename = "diff_rangeD_niche_D~trophic_level.pdf",
            width = 8)

### stability ~ trophic level:
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = range_stability_std,
            ylab = "range S", ylim = c(0.7, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_S~trophic_level.pdf")
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = niche_stability_std,
            ylab = "niche S", ylim = c(0.7, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_S~trophic_level.pdf")

### expansion ~ trophic level:
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = range_expansion_std,
            ylab = "range E", ylim = c(0, 0.5), 
            title = paste0(dataset, ", background = buffer"), filename = "range_E~trophic_level.pdf")
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = niche_expansion_std,
            ylab = "niche E", ylim = c(0, 0.5), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_E~trophic_level.pdf")

### unfilling ~ trophic level:
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = range_unfilling_std,
            ylab = "range U", ylim = c(0, 0.5), 
            title = paste0(dataset, ", background = buffer"), filename = "range_U~trophic_level.pdf")
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = niche_unfilling_std,
            ylab = "niche U", ylim = c(0, 0.5), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_U~trophic_level.pdf")

### range conservatism ~ trophic level:
# S:
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = cons_p_S_A_range,
            ylab = "range conservatism p-value S", ylim = c(0, 1),
            title = paste0(dataset, ", background = buffer"), filename = "range_cons_p_S~trophic_level.pdf")
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = cons_p_S_A_niche,
            ylab = "niche conservatism p-value S", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_cons_p_S~trophic_level.pdf")

# range:
res_traits_df %>% 
  group_by(Trophic.Level) %>% 
  mutate(n_tl = n()) %>% 
  filter(cons_p_S_A_range <= 0.05) %>% 
  group_by(Trophic.Level) %>% 
  summarise(n = n(), n_tl = unique(n_tl)) %>% 
  mutate(perc = n / n_tl * 100)

# niche:
res_traits_df %>% 
  group_by(Trophic.Level) %>% 
  mutate(n_tl = n()) %>% 
  filter(cons_p_S_A_niche <= 0.05) %>% 
  group_by(Trophic.Level) %>% 
  summarise(n = n(), n_tl = unique(n_tl)) %>% 
  mutate(perc = n / n_tl * 100)

### range shift ~ trophic level:
# S:
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = shift_p_S_A_range,
            ylab = "range shift p-value S", ylim = c(0, 1),
            title = paste0(dataset, ", background = buffer"), filename = "range_shift_p_S~trophic_level.pdf")
boxplot_fun(data = res_traits_df, x = Trophic.Level, y = shift_p_S_A_niche,
            ylab = "niche shift p-value S", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_shift_p_S~trophic_level.pdf")


## Habitat: ----

boxplot_fun(data = res_traits_df, x = Habitat, y = D_range,
            ylab = "range overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_D~Habitat.pdf",
            width = 10)
boxplot_fun(data = res_traits_df, x = Habitat, y = D_niche,
            ylab = "niche overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_D~Habitat.pdf",
            width = 10)

# difference range D - niche D:
res_traits_df %>% 
  mutate(diff_range_niche_D = D_range - D_niche) %>% 
  group_by(Habitat) %>% 
  summarise(mean_diff_range_niche_D = mean(diff_range_niche_D, na.rm = TRUE)) %>% 
  arrange(-mean_diff_range_niche_D)

boxplot_fun(data = res_traits_df, x = Habitat, y = diff_range_niche_D,
            ylab = "range D - niche D", ylim = range(res_traits_df$diff_range_niche_D), 
            title = paste0(dataset, ", background = buffer"), filename = "diff_rangeD_niche_D~Habitat.pdf",
            width = 8)


### stability ~ habitat:
boxplot_fun(data = res_traits_df, x = Habitat, y = range_stability_std,
            ylab = "range S", ylim = c(0.5, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_S~Habitat.pdf",
            width = 8)
boxplot_fun(data = res_traits_df, x = Habitat, y = niche_stability_std,
            ylab = "niche S", ylim = c(0.5, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_S~Habitat.pdf",
            width = 8)


## Body weight:----

scatterplot_fun(data = res_traits_df, x = Mass, y = D_range,
                ylab = "range overlap D", xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), 
                filename = "range_D~body_weight.pdf",
                width = 8, height = 5)

res_traits_df_ss <- res_traits_df %>% 
  filter(Mass <= quantile(res_traits_df$Mass, 0.95))

scatterplot_fun(data = res_traits_df_ss, x = Mass, y = D_range,
                ylab = "range overlap D", xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), 
                filename = "range_D~body_weight_subset.pdf",
                width = 8, height = 5)

scatterplot_fun(data = res_traits_df, x = Mass, y = D_niche,
                ylab = "niche overlap D", xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), 
                filename = "niche_D~body_weight.pdf",
                width = 8, height = 5)

scatterplot_fun(data = res_traits_df_ss, x = Mass, y = D_niche,
                ylab = "niche overlap D", xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), 
                filename = "niche_D~body_weight_subset.pdf",
                width = 8, height = 5)

# difference range D - niche D:
scatterplot_fun(data = res_traits_df, x = Mass, y = diff_range_niche_D,
                ylab = "range D - niche D", xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), 
                filename = "diff_rangeD_niche_D~body_weight.pdf",
                width = 8, height = 5, ylim = range(res_traits_df$diff_range_niche_D))

### stability ~ body weight:
scatterplot_fun(data = res_traits_df, x = Mass, y = range_stability_std,
            ylab = "range S", ylim = c(0.5, 1),xlab = "Body weight",
            title = paste0(dataset, ", background = buffer"), filename = "range_S~body_weight.pdf")
scatterplot_fun(data = res_traits_df, x = Mass, y = niche_stability_std,
            ylab = "niche S", ylim = c(0.5, 1), xlab = "Body weight",
            title = paste0(dataset, ", background = buffer"), filename = "niche_S~body_weight.pdf",
            width = 8)

### expansion ~ body weight:
scatterplot_fun(data = res_traits_df, x = Mass, y = range_expansion_std,
                ylab = "range E", ylim = c(0, 0.5),xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), filename = "range_E~body_weight.pdf")
scatterplot_fun(data = res_traits_df, x = Mass, y = niche_expansion_std,
                ylab = "niche E", ylim = c(0, 0.5), xlab = "Body weight",
                title = paste0(dataset, ", background = buffer"), filename = "niche_E~body_weight.pdf",
                width = 8)



## Latitudinal range centroid: ----

scatterplot_fun(data = res_traits_df, x = Centroid.Latitude, y = D_range,
                ylab = "range overlap D", xlab = "Lat. range centroid",
                title = paste0(dataset, ", background = buffer"), 
                filename = "range_D~lat_range_centroid.pdf",
                width = 8, height = 5)

scatterplot_fun(data = res_traits_df, x = Centroid.Latitude, y = D_niche,
                ylab = "niche overlap D", xlab = "Lat. range centroid",
                title = paste0(dataset, ", background = buffer"), 
                filename = "niche_D~lat_range_centroid.pdf",
                width = 8, height = 5)

# difference range D - niche D:
scatterplot_fun(data = res_traits_df, x = Centroid.Latitude, y = diff_range_niche_D,
                ylab = "range D - niche D", xlab = "Lat. range centroid",
                title = paste0(dataset, ", background = buffer"), 
                filename = "diff_rangeD_niche_D~lat_range_centroid.pdf",
                width = 8, height = 5, ylim = range(res_traits_df$diff_range_niche_D))

## Range size:

scatterplot_fun(data = res_traits_df, x = Range.Size, y = D_range,
                ylab = "range overlap D", xlab = "Range size",
                title = paste0(dataset, ", background = buffer"), 
                filename = "range_D~range_size.pdf",
                width = 8, height = 5)

scatterplot_fun(data = res_traits_df, x = Range.Size, y = D_niche,
                ylab = "niche overlap D", xlab = "Range size",
                title = paste0(dataset, ", background = buffer"), 
                filename = "niche_D~range_size.pdf",
                width = 8, height = 5)

# difference range D - niche D:

scatterplot_fun(data = res_traits_df, x = Range.Size, y = diff_range_niche_D,
                ylab = "range D - niche D", xlab = "Range size",
                title = paste0(dataset, ", background = buffer"), 
                filename = "diff_rangeD_niche_D~range_size.pdf",
                width = 8, height = 5, ylim = range(res_traits_df$diff_range_niche_D))

# S:
scatterplot_fun(data = res_traits_df, x = Range.Size, y = range_stability_std,
                ylab = "range S", xlab = "Range size",
                title = paste0(dataset, ", background = buffer"), 
                filename = "range_S~range_size.pdf",
                width = 8, height = 5)

## Others: ----
## Trophic niche:

boxplot_fun(data = res_traits_df, x = Trophic.Niche, y = D_range,
            ylab = "range overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_D~trophic_niche.pdf",
            width = 10)
boxplot_fun(data = res_traits_df, x = Trophic.Niche, y = D_niche,
            ylab = "niche overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_D~trophic_niche.pdf",
            width = 10)
# difference range D - niche D:
boxplot_fun(data = res_traits_df, x = Trophic.Niche, y = diff_range_niche_D,
            ylab = "range D - niche D", ylim = range(res_traits_df$diff_range_niche_D), 
            title = paste0(dataset, ", background = buffer"), filename = "diff_rangeD_niche_D~trophic_niche.pdf",
            width = 8)

## Primary lifestyle:

boxplot_fun(data = res_traits_df, x = Primary.Lifestyle, y = D_range,
            ylab = "range overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_D~primary_lifestyle.pdf",
            width = 10)
boxplot_fun(data = res_traits_df, x = Primary.Lifestyle, y = D_niche,
            ylab = "niche overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_D~primary_lifestyle.pdf",
            width = 10)

# difference range D - niche D:
boxplot_fun(data = res_traits_df, x = Primary.Lifestyle, y = diff_range_niche_D,
            ylab = "range D - niche D", ylim = range(res_traits_df$diff_range_niche_D), 
            title = paste0(dataset, ", background = buffer"), filename = "diff_rangeD_niche_D~primary_lifestyle.pdf",
            width = 8)

## Migration:

boxplot_fun(data = res_traits_df, x = as.factor(Migration), y = D_range,
            ylab = "range overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "range_D~Migration.pdf")
boxplot_fun(data = res_traits_df, x = as.factor(Migration), y = D_niche,
            ylab = "niche overlap D", ylim = c(0, 1), 
            title = paste0(dataset, ", background = buffer"), filename = "niche_D~Migration.pdf")

# difference range D - niche D:
boxplot_fun(data = res_traits_df, x = as.factor(Migration), y = diff_range_niche_D,
            ylab = "range D - niche D", ylim = range(res_traits_df$diff_range_niche_D), 
            title = paste0(dataset, ", background = buffer"), filename = "diff_rangeD_niche_D~migration.pdf",
            width = 8)