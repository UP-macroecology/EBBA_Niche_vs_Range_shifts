# Filter species to use in analysis, steps 1 to 4:

# - Exclude pelagic specialists (according to Wilman et al. 2014) 
# - Exclude rare species with n<20 occurrences in any of the two atlas periods 
# - Exclude very common species with >90% prevalence in at least one of the two atlas periods 
# - Only those species that occur in both atlas periods 

library(sf)
library(dplyr)

# read prepared EBBA data (output of 1_EBBA_prep_data.R):

EBBA1_prep_sf <- st_read(file.path("Data", "EBBA1_comparable_harmonized.shp"))
EBBA2_prep_sf <- st_read(file.path("Data", "EBBA2_comparable_harmonized.shp"))

#EBBA1_prep_sf <- st_read(file.path("Data", "EBBA1_change.shp"))
#EBBA2_prep_sf <- st_read(file.path("Data", "EBBA2_change.shp"))

# drop geometry (not needed to create species list)
EBBA1_prep <- EBBA1_prep_sf %>% 
  st_drop_geometry()
EBBA2_prep <- EBBA2_prep_sf %>% 
  st_drop_geometry()

length(unique(EBBA1_prep$species)) # 445 species
length(unique(EBBA2_prep$species)) # 517 species


## 1. exclude pelagic specialists (according to Wilman et al. 2014): ----

SeaBirds <- read.csv2(file.path("Data", "BirdFuncDat.txt"), header = TRUE, sep = "\t") # BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
SeaBirds <- subset(SeaBirds, PelagicSpecialist==1)
SeaBirds <- SeaBirds$Scientific

EBBA1_prep <- EBBA1_prep[which(!EBBA1_prep$species %in% SeaBirds),]
EBBA2_prep <- EBBA2_prep[which(!EBBA2_prep$species %in% SeaBirds),]

length(unique(EBBA1_prep$species)) # 417 species left
length(unique(EBBA2_prep$species)) # 491 species left


## 2. exclude rare species with n<20 occurrences in any of the two atlas periods: ----

EBBA1_prep <- EBBA1_prep %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup

length(unique(EBBA1_prep$species)) # 371 species left

EBBA2_prep <- EBBA2_prep %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup

length(unique(EBBA2_prep$species)) # 390 species left


## 3. exclude very common species with >90% prevalence in both atlas periods: ----
# = occur in >90% EBBA cells:

# read prepared EBBA data (output of 1_prep_EBBA_data.R, before taxonomic harmonization, to get number of comparable EBBA cells):
EBBA1_comp_sf <- st_read(file.path("Data", "EBBA1_comparable.shp"))
EBBA2_comp_sf <- st_read(file.path("Data", "EBBA2_comparable.shp"))

#EBBA_change_sf <- st_read(file.path("Data", "EBBA_change.shp"))

nEBBAcells <- length(unique(EBBA1_comp_sf$cell50x50)) # 2845 (EBBA2 methods chapter: "change map thus encompassed a total of 2,972 50-km squares"; here less because I excluded cells from EBBA2 not in EBBA1, the 2,972 squares contain cells not included in change analysis due to insufficient coverage?!
nEBBAcells <- length(unique(EBBA_change_sf$cell50x50)) # 2913 xx

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

length(unique(EBBA1_prep$species)) # 365 species left

EBBA2_prep <- EBBA2_prep %>% 
  filter(n_occurrences <= 0.9 * nEBBAcells)

length(unique(EBBA2_prep$species)) # 379 species left


## 4. keep only those species that occur in both atlas periods: ----

EBBA1_prep <- EBBA1_prep %>% 
  filter(species %in% EBBA2_prep$species)

length(unique(EBBA1_prep$species)) # 325 species left

EBBA2_prep <- EBBA2_prep %>% 
  filter(species %in% EBBA1_prep$species)

length(unique(EBBA2_prep$species)) # 325 species left

## how often does each species occur in each of the EBBAs:
EBBA1_prep %>%
  arrange(-n_occurrences) %>%
  distinct(species, n_occurrences)

EBBA2_prep %>%
  arrange(-n_occurrences) %>%
  distinct(species, n_occurrences) %>% View

save(EBBA1_prep, EBBA2_prep, file = file.path("Data", "EBBA1_EBBA2_prep_steps1-4_prelim.RData"))
