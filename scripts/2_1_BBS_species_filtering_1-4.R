# Filter BBS species to use in analysis:

# - Exclude pelagic specialists (according to Wilman et al. 2014) 
# - Exclude rare species with n<20 occurrences in any of the two periods 
# - Exclude very common species with >90% prevalence in at least one of the two periods 
# - Only those species that occur in both periods 
# - Only records with information on species level

library(dplyr)
library(sf)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
data_dir <- file.path("data", "BBS_analysis")

# pelagic specialists (according to Wilman et al. 2014):
SeaBirds <- read.csv2(file.path("data", "BirdFuncDat.txt"), header = TRUE, sep = "\t") # BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
SeaBirds <- subset(SeaBirds, PelagicSpecialist==1)
SeaBirds <- SeaBirds$Scientific


# ----------------------------- #
#       Filter species:      ####
# ----------------------------- #

historic_periods <- list(1981:1983, 1988:1990) 

for(i in 1:length(historic_periods)){
  
  # read prepared BBS data (output of 1_BBS_prep_data.R):
  
  if(i == 1){
    hist_prep_sf <- st_read(file.path(data_dir, "BBS_historic_centr_proj_hist81-83.shp")) # output of 1_BBS_prep_data.R
    rec_prep_sf <- st_read(file.path(data_dir, "BBS_recent_centr_proj_hist81-83.shp")) # output of 1_BBS_prep_data.R
    } else {
      hist_prep_sf <- st_read(file.path(data_dir, "BBS_historic_centr_proj_hist88-90.shp")) # output of 1_BBS_prep_data.R
      rec_prep_sf <- st_read(file.path(data_dir, "BBS_recent_centr_proj_hist88-90.shp")) # output of 1_BBS_prep_data.R
      }
  
  # number of validation routes:
  n_valid_routes <- length(unique(hist_prep_sf$RTENO)) # V1: 522; V2: 721
  
  # convert to data frame and join species names:
  
  # historic time period:
  hist_prep_df <- hist_prep_sf %>% 
    st_drop_geometry() %>% 
    filter(pres == 1)# only keep species with presence records
  
  length(unique(hist_prep_df$species)) # 459 species for historic time period 1981-1983 (V1); 497 for historic time period 1988-1990 (V2) 
  
  # recent time period:
  rec_prep_df <- rec_prep_sf %>% 
    st_drop_geometry() %>% 
    filter(pres == 1) # only keep species with presence records
  
  length(unique(rec_prep_df$species)) # 480 species in recent time period (V1), 522 (V2)
  
  
  # 1. exclude pelagic specialists (according to Wilman et al. 2014): ----
  
  hist_prep_df <- hist_prep_df[which(!hist_prep_df$species %in% SeaBirds),]
  rec_prep_df <- rec_prep_df[which(!rec_prep_df$species %in% SeaBirds),]
  
  length(unique(hist_prep_df$species)) # 452 species left (V1: 7 species removed); 488 (V2)
  length(unique(rec_prep_df$species)) # 477 species left (V1: 3 species removed); 512 (V2)
  
  
  # 2. exclude rare species with n<20 occurrences in any of the two periods: ----
  
  hist_prep_df <- hist_prep_df %>% 
    group_by(species) %>% 
    mutate(n_occurrences = n()) %>%
    filter(n_occurrences >= 20) %>% 
    ungroup # 34,468
  
  length(unique(hist_prep_df$species)) # 230 species left (V1); 285 (V2)
  
  rec_prep_df <- rec_prep_df %>% 
    group_by(species) %>% 
    mutate(n_occurrences = n()) %>%
    filter(n_occurrences >= 20) %>% 
    ungroup # 36,494
  
  length(unique(rec_prep_df$species)) # 235 species left (V1); 288 (V2)
  
  
  # 3. exclude very common species with >90% prevalence in both periods: ----
  # = occur in >90% of validation routes:
  
  # which species are excluded:
  hist_prep_df %>% 
    filter(n_occurrences > 0.9 * n_valid_routes) %>% 
    distinct(species)
  rec_prep_df %>% 
    filter(n_occurrences > 0.9 * n_valid_routes) %>% 
    distinct(species)
  
  # exclude species:
  hist_prep_df <- hist_prep_df %>% 
    filter(n_occurrences <= 0.9 * n_valid_routes)
  
  length(unique(hist_prep_df$species)) # 221 species left (V1); 280 (V2)
  
  rec_prep_df <- rec_prep_df %>% 
    filter(n_occurrences <= 0.9 * n_valid_routes)
  
  length(unique(rec_prep_df$species)) # 228 species left (V1); 284 (V2)
  
  
  # 4. keep only those species that occur in both periods: ----
  
  hist_prep_df <- hist_prep_df %>% 
    filter(species %in% rec_prep_df$species)
  
  length(unique(hist_prep_df$species)) # 212 species left (V1); 264 (V2)
  
  rec_prep_df <- rec_prep_df %>% 
    filter(species %in% hist_prep_df$species)
  
  length(unique(rec_prep_df$species)) # 212 species left (V1); 264 (V2)
  
  # 5. remove records where information is not on species level: ----
  
  hist_prep_df <- hist_prep_df %>% 
    filter(!grepl(pattern = "sp\\.|/", x = species))
  
  length(unique(hist_prep_df$species)) # 212 (V1); 264 (V2)
  
  rec_prep_df <- rec_prep_df %>% 
    filter(!grepl(pattern = "sp\\.|/", x = species))
  
  length(unique(rec_prep_df$species)) # 212 (V1); 264 (V2)
  
  # how often does each species occur:
  hist_prep_df <- hist_prep_df %>%
    arrange(-n_occurrences) %>%
    dplyr::select(species, n_occurrences) %>% 
    distinct(species, n_occurrences)
  
  rec_prep_df <- rec_prep_df %>%
    arrange(-n_occurrences) %>%
    dplyr::select(species, n_occurrences) %>% 
    distinct(species, n_occurrences)
  
  if(i == 1){
    save(hist_prep_df, rec_prep_df, file = file.path(data_dir, "BBS_prep_steps1-4_hist81-83.RData"))
    } else {
      save(hist_prep_df, rec_prep_df, file = file.path(data_dir, "BBS_prep_steps1-4_hist88-90.RData"))
      }
}