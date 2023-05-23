# Script to explore the BioShifts dataset (Lenoir et al.) with regard to selected species of EBBA and BBS analysis:

# (https://figshare.com/articles/dataset/BioShifts_a_global_geodatabase_of_climate-induced_species_redistribution_over_land_and_sea/7413365/1)
# Lenoir J., Bertrand, Comte L., R., Murienne J., Bourgeaud L., Hattab T., Grenouillet G.  (2020) Species better track climate warming in the oceans than on land. TO BE COMPLETED
# Comte L., Grenouillet G., Bertrand, R., Murienne J., Bourgeaud L., Hattab T., Lenoir J. (2020) BioShifts: a global geodatabase of climate-induced species redistribution over land and sea. Figshare, 10.6084/m9.figshare.7413365.

library(dplyr)
library(sf)

# read and subset BioShifts data: ----

BS <- read.csv(file.path("data", "BioShifts", "BioShifts", "BioShifts.csv"))

BS_subset <- BS %>% 
  as_tibble() %>% 
  filter(Class == "Aves") %>% 
  filter(Ecosystem == "Terrestrial") %>% 
  filter(Gradient == "Latitudinal") %>% 
  filter(Hemisphere == "North")

length(unique(BS_subset$Species)) # 687 bird species

# study areas of the studies included in the BioShifts dataset: 
study_areas_dir <- file.path("data", "BioShifts", "Bioshifts", "Study_Areas.gdb")
layer_list <- st_layers(study_areas_dir)

# BioShifts references:
BioShifts_ref <- read.csv(file.path("data", "BioShifts", "Bioshifts", "References.csv"), header = FALSE)


# EBBA analysis: ----

# BioShifts data for EBBA species names:

# EBBA selected species:
load(file = file.path("data", "EBBA_analysis", "EBBA1_EBBA2_prep_steps1-4_final.RData")) # output of 2_1_EBBA_species_filtering_1-4.R
sel_species <- sort(unique(EBBA1_prep$species))

BS_subset_spec <- BS_subset %>% 
  filter(Species %in% sel_species)

sel_species[which(sel_species %in% BS_subset_spec$Species)] 
sel_species[which(!sel_species %in% BS_subset_spec$Species)] # 39 of selected 117 EBBA species not directly found in BioShifts (without name harmonization effort)

# study areas of the studies containing EBBA species included in the BioShifts dataset: 
rel_areas_EBBA <- layer_list$name[which(layer_list$name %in% BS_subset_spec$Source)]

study_areas_EBBA_sf <- list()
for(i in 1:length(rel_areas_EBBA)) {
  print(i)
  study_areas_EBBA_sf[[rel_areas_EBBA[i]]] = read_sf(study_areas_dir, layer = rel_areas_EBBA[i])
  plot(st_geometry(study_areas_EBBA_sf[[i]]), main = i)
}

EBBA_polys <- names(study_areas_EBBA_sf[c(17,26)]) # polygons covering most of EBBA area

# corresponding references:
ref_ids <- BS_subset_spec %>% 
  filter(Source %in% EBBA_polys) %>% 
  pull(Reference) %>% 
  unique

BioShifts_ref[ref_ids, 2]


# BBS analysis: ----

# BioShifts data for BBS species:

# species selection: 
BBS_sel_species <- read.csv(file = file.path("data", "BBS_analysis", "species_stability_contUS_BL22.csv")) %>% 
  filter(stability >= 0.5) %>% 
  pull(species) %>% 
  sort # 297

BBS_sel_species[which(!BBS_sel_species %in% BS_subset$Species)] # 41 of selected 297 BBS species not directly found in BioShifts

BS_subset_spec_BBS <- BS_subset %>% 
  filter(Species %in% BBS_sel_species)

# study areas corresponding to selected BBS species:

rel_areas_BBS <- layer_list$name[which(layer_list$name %in% BS_subset_spec_BBS$Source)]

study_areas_BBS_sf <- list()
for(i in 1:length(rel_areas_BBS)) {
  print(i)
  study_areas_BBS_sf[[rel_areas_BBS[i]]] = read_sf(study_areas_dir, layer = rel_areas_BBS[i])
  plot(st_geometry(study_areas_BBS_sf[[i]]), main = i)
}

contUS_polys <- names(study_areas_BBS_sf[c(1,11,23,25,26,27,28,30)]) # polygons covering most of conterminous US

# corresponding references:
ref_ids <- BS_subset_spec_BBS %>% 
  filter(Source %in% contUS_polys) %>% 
  pull(Reference) %>% 
  unique

BioShifts_ref[ref_ids, 2]