# Filter species to use in analysis, steps 1 to 4:

# - Exclude pelagic specialists (according to Wilman et al. 2014) 
# - Exclude rare species with n<20 occurrences in any of the two atlas periods 
# - Exclude very common species with >90% prevalence in at least one of the two atlas periods 
# - Only those species that occur in both atlas periods 

library(sf)
library(dplyr)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# EBBA data folder:
data_dir <- file.path("data", "EBBA_analysis")

# final or preliminary filtering:
final_filtering <- TRUE # TRUE: filtering of EBBA change data, FALSE = preliminary filtering to request EBBA change data

# results file:
if(final_filtering){
  res_file <- file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_final.RData")
} else {
  res_file <- file.path(data_dir, "EBBA1_EBBA2_prep_steps1-4_prelim.RData")
}


# ---------------------------- #
#          Load data:       ####
# ---------------------------- #

# read prepared EBBA data:

if(final_filtering){
  
  # species with stability >= 50 %:
  sel_species_stab <- read.csv(file = file.path(data_dir, "species_stability_EBBA2_BL22_030423.csv")) %>% # output of 2_3_EBBA_species_filtering_5_climatic_niche_analysis.R
    filter(stability >= 0.5) %>% 
    pull(species)
  
  EBBA1_prep_sf <- st_read(file.path(data_dir, "EBBA1_change.shp")) %>% 
    filter(species %in% sel_species_stab) # output of 1_EBBA_prep_data.R
  
  EBBA2_prep_sf <- st_read(file.path(data_dir, "EBBA2_change.shp")) %>% 
    filter(species %in% sel_species_stab) # output of 1_EBBA_prep_data.R
  
  EBBA_cells_sf <- st_read(file.path(data_dir, "EBBA_change.shp")) # output of 1_EBBA_prep_data.R
  
} else {
  
  EBBA1_prep_sf <- st_read(file.path(data_dir, "EBBA1_prelim_comparable_harmonized.shp")) # output of 1_EBBA_prep_data.R
  
  EBBA2_prep_sf <- st_read(file.path(data_dir, "EBBA2_prelim_comparable_harmonized.shp")) # output of 1_EBBA_prep_data.R
  
  # comparable EBBA cells:
  EBBA_cells_sf <- st_read(file.path(data_dir, "EBBA1_prelim_comparable_cells.shp")) # output of 1_EBBA_prep_data.R
}


# drop geometry (not needed to create species list)
EBBA1_prep <- EBBA1_prep_sf %>% 
  st_drop_geometry()
EBBA2_prep <- EBBA2_prep_sf %>% 
  st_drop_geometry()

length(unique(EBBA1_prep$species))
length(unique(EBBA2_prep$species))


# ---------------------------- #
#      Filter species:      ####
# ---------------------------- #

## 1. exclude pelagic specialists (according to Wilman et al. 2014): ----

SeaBirds <- read.csv2(file.path("data", "BirdFuncDat.txt"), header = TRUE, sep = "\t") # BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
SeaBirds <- subset(SeaBirds, PelagicSpecialist==1)
SeaBirds <- SeaBirds$Scientific

EBBA1_prep <- EBBA1_prep[which(!EBBA1_prep$species %in% SeaBirds),]
EBBA2_prep <- EBBA2_prep[which(!EBBA2_prep$species %in% SeaBirds),]

length(unique(EBBA1_prep$species))
length(unique(EBBA2_prep$species))


## 2. exclude rare species with n<20 occurrences in any of the two atlas periods: ----

EBBA1_prep <- EBBA1_prep %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup

length(unique(EBBA1_prep$species))

EBBA2_prep <- EBBA2_prep %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup

length(unique(EBBA2_prep$species))


## 3. exclude very common species with >90% prevalence in both atlas periods: ----
# = occur in >90% EBBA cells:

# number of comparable EBBA cells:
nEBBAcells <- length(unique(EBBA_cells_sf$cell50x50))

# which species are excluded:
EBBA1_prep %>% 
  filter(n_occurrences > 0.9 * nEBBAcells) %>% 
  distinct(species)
EBBA2_prep %>% 
  filter(n_occurrences > 0.9 * nEBBAcells) %>% 
  distinct(species)

# exclude species:
EBBA1_prep <- EBBA1_prep %>% 
  filter(n_occurrences <= 0.9 * nEBBAcells)

length(unique(EBBA1_prep$species))

EBBA2_prep <- EBBA2_prep %>% 
  filter(n_occurrences <= 0.9 * nEBBAcells)

length(unique(EBBA2_prep$species))


## 4. keep only those species that occur in both atlas periods: ----

EBBA1_prep <- EBBA1_prep %>% 
  filter(species %in% EBBA2_prep$species)

length(unique(EBBA1_prep$species))

EBBA2_prep <- EBBA2_prep %>% 
  filter(species %in% EBBA1_prep$species)

length(unique(EBBA2_prep$species))

## how often does each species occur in each of the EBBAs:
EBBA1_prep %>%
  arrange(-n_occurrences) %>%
  distinct(species, n_occurrences)

EBBA2_prep %>%
  arrange(-n_occurrences) %>%
  distinct(species, n_occurrences)

save(EBBA1_prep, EBBA2_prep, file = res_file)
