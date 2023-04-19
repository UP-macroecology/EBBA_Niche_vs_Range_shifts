# Filter BBS species to use in analysis, steps 1 to 4:

# - Exclude pelagic specialists (according to Wilman et al. 2014) 
# - Exclude rare species with n<20 occurrences in any of the two periods 
# - Exclude very common species with >90% prevalence in at least one of the two periods 
# - Only those species that occur in both periods 

# run once for each of the two historic time periods used

library(dplyr)
library(sf)

#data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")
data_dir <- file.path("Data")

# which historic time period should be used (V1: 1981-1983, V2: 1996-1998)
historic <- 1996:1998 # 1981:1983 or 1996:1998

# read prepared BBS data (output of 1_BBS_prep_data.R):

if(all(historic == 1981:1983)){
  hist_prep_sf <- st_read(file.path(data_dir, "BBS_historic1_proj.shp")) # output of 1_BBS_prep_data.R
  rec_prep_sf <- st_read(file.path(data_dir, "BBS_recent1_proj.shp")) # output of 1_BBS_prep_data.R
}
if(all(historic == 1996:1998)){
  hist_prep_sf <- st_read(file.path(data_dir, "BBS_historic2_proj.shp")) # output of 1_BBS_prep_data.R
  rec_prep_sf <- st_read(file.path(data_dir, "BBS_recent2_proj.shp")) # output of 1_BBS_prep_data.R
}

# number of validation routes:
n_valid_routes <- length(unique(hist_prep_sf$RTENO)) # V1: 522; V2: 933

# convert to data frame and join species names:

# historic time period:
hist_prep_df <- hist_prep_sf %>% 
  st_drop_geometry() %>% 
  filter(pres == 1)# only keep species with presence records

length(unique(hist_prep_df$species)) # 459 species for historic time period 1981-1983 (V1); 513 for historic time period 1996-1998 (V2) 

# recent time period:
rec_prep_df <- rec_prep_sf %>% 
  st_drop_geometry() %>% 
  filter(pres == 1) # only keep species with presence records

length(unique(rec_prep_df$species)) # 480 species in recent time period (V1), 525 (V2)


## 1. exclude pelagic specialists (according to Wilman et al. 2014): ----

SeaBirds <- read.csv2(file.path(data_dir, "BirdFuncDat.txt"), header = TRUE, sep = "\t") # BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
SeaBirds <- subset(SeaBirds, PelagicSpecialist==1)
SeaBirds <- SeaBirds$Scientific # $English

hist_prep_df <- hist_prep_df[which(!hist_prep_df$species %in% SeaBirds),]
rec_prep_df <- rec_prep_df[which(!rec_prep_df$species %in% SeaBirds),]

length(unique(hist_prep_df$species)) # 452 species left (V1: 7 species removed); 505 (V2, 8 species removed)
length(unique(rec_prep_df$species)) # 477 species left (V1: 3 species removed); 518 (V2, 7 species removed)


## 2. exclude rare species with n<20 occurrences in any of the two periods: ----
hist_prep_df <- hist_prep_df %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup # 34,468

length(unique(hist_prep_df$species)) # 230 species left (V1); 326 (V2)

rec_prep_df <- rec_prep_df %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup # 36,494

length(unique(rec_prep_df$species)) # 235 species left (V1); 337 (V2)


## 3. exclude very common species with >90% prevalence in both periods: ----
# = occur in >90% of validation routes: xx

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

length(unique(hist_prep_df$species)) # 221 species left (V1); 322 (V2)

rec_prep_df <- rec_prep_df %>% 
  filter(n_occurrences <= 0.9 * n_valid_routes)

length(unique(rec_prep_df$species)) # 228 species left (V1); 334 (V2)


## 4. keep only those species that occur in both atlas periods: ----

hist_prep_df <- hist_prep_df %>% 
  filter(species %in% rec_prep_df$species)

length(unique(hist_prep_df$species)) # 212 species left (V1); 314 (V2)

rec_prep_df <- rec_prep_df %>% 
  filter(species %in% hist_prep_df$species)

length(unique(rec_prep_df$species)) # 212 species left (V1); 314 (V2)

## 5. remove records where information is not on species level:
hist_prep_df <- hist_prep_df %>% 
  filter(!grepl(pattern = "sp\\.|/", x = species))

length(unique(hist_prep_df$species)) # 212 (V1); 312 (V2)

rec_prep_df <- rec_prep_df %>% 
  filter(!grepl(pattern = "sp\\.|/", x = species))

length(unique(rec_prep_df$species)) # 212 (V1); 312 (V2)

## how often does each species occur:
hist_prep_df %>%
  arrange(-n_occurrences) %>%
  distinct(species, n_occurrences)

rec_prep_df %>%
  arrange(-n_occurrences) %>%
  distinct(species, n_occurrences)

if(all(historic == 1981:1983)){
  save(hist_prep_df, rec_prep_df, file = file.path(data_dir, "BBS_prep_steps1-4_V1.RData"))
}
if(all(historic == 1996:1998)){
  save(hist_prep_df, rec_prep_df, file = file.path(data_dir, "BBS_prep_steps1-4_V2.RData"))
}
