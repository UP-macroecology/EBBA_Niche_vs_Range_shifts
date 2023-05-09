# 1) Assemble trait data for European and North American bird species left after 5 filtering steps
# (2_1_species_filtering_1-4.R and 2_3_species_filtering_5_climatic_niche_analysis.R)
# using AVONET data

library(dplyr)

# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
data_dir <- file.path("Data")

# AVONET data:
# Tobias et al. 2022: AVONET: morphological, ecological and geographical data for all birds.
# I downloaded the AVONET dataset from https://figshare.com/s/b990722d72a26b5bfead -> AVONET Supplementary dataset 1
# (link in publication)
# (I didn't use the package "traitdata" to access AVONET because the package seems to use 
# the data provided to reproduce the analyses of Tobias at al. 2022 instead of the full AVONET data set
# (e.g. species names only available in one unspecified format while main AVONET data set provides 3 formats))


# ----------------------------------------------- #
#   Match AVONET trait data to species names:  ####
# ----------------------------------------------- #

# read species names:

EBBA_species_final <- read.csv(file = file.path(data_dir, "species_stability_EBBA2_BL22_030423.csv")) %>% 
  arrange(-stability) %>% 
  filter(stability >= 0.5) %>% 
  pull(species) %>% 
  sort

BBS_species_final <- read.csv(file = file.path(data_dir, "species_stability_contUS_BL22_060423.csv")) %>% 
  arrange(-stability) %>% 
  filter(stability >= 0.5) %>% 
  pull(species)


# read AVONET data:

# sheet with taxonomy according to Birdlife International
avonet_BL <- readxl::read_xlsx(path = file.path(data_dir, "AVONET Supplementary dataset 1.xlsx"),
                 sheet = "AVONET1_BirdLife")

# AVONET data for selected EBBA species: ---------------------------------------

EBBA_avonet <- subset(avonet_BL, avonet_BL$Species1 %in% EBBA_species_final)
# all selected EBBA species are found in AVONET
write.csv(EBBA_avonet, file = file.path(data_dir, "AVONET_EBBA_species.csv"),
          row.names = FALSE)


# AVONET data for selected BBS species: ----------------------------------------

BBS_species_final %in% avonet_BL$Species1
# not for all BBS names matching AVONET data are found directly

# for some species BBS uses different names than BL:
# (used also in 2_3_BBS_species_filtering_5_climatic_niche_analysis.R)
spec_name_change_df <- data.frame("BBS_name" = "Antigone canadensis", "BL_name" = "Grus canadensis")
spec_name_change_df[2,] <- c("Aphelocoma woodhouseii", "Aphelocoma californica")
spec_name_change_df[3,] <- c("Centronyx bairdii", "Passerculus bairdii")
spec_name_change_df[4,] <- c("Centronyx henslowii", "Passerculus henslowii")
spec_name_change_df[5,] <- c("Chroicocephalus philadelphia", "Larus philadelphia")
spec_name_change_df[6,] <- c("Coccothraustes vespertinus", "Hesperiphona vespertina")
spec_name_change_df[7,] <- c("Colaptes auratus cafer", "Colaptes cafer")
spec_name_change_df[8,] <- c("Dryobates albolarvatus", "Leuconotopicus albolarvatus")
spec_name_change_df[9,] <- c("Himantopus mexicanus", "Himantopus himantopus")
spec_name_change_df[10,] <- c("Icterus bullockii", "Icterus bullockiorum")
spec_name_change_df[11,] <- c("Junco hyemalis caniceps", "Junco hyemalis")
spec_name_change_df[12,] <- c("Junco hyemalis hyemalis", "Junco hyemalis")
spec_name_change_df[13,] <- c("Junco hyemalis oreganus", "Junco hyemalis")
spec_name_change_df[14,] <- c("Leucophaeus atricilla", "Larus atricilla")
spec_name_change_df[15,] <- c("Leucophaeus pipixcan", "Larus pipixcan")
spec_name_change_df[16,] <- c("Phalacrocorax auritus", "Nannopterum auritus")
spec_name_change_df[17,] <- c("Phalaropus tricolor", "Steganopus tricolor")
spec_name_change_df[18,] <- c("Pica nuttalli", "Pica nutalli")
spec_name_change_df[19,] <- c("Regulus calendula", "Corthylio calendula")
spec_name_change_df[20,] <- c("Setophaga coronata audoboni", "Setophaga auduboni")
spec_name_change_df[21,] <- c("Setophaga coronata coronata", "Setophaga coronata")
spec_name_change_df[22,] <- c("Butorides virescens", "Butorides striata")
spec_name_change_df[23,] <- c("Dryobates villosus", "Leuconotopicus villosus")
spec_name_change_df[24,] <- c("Dryocopus pileatus", "Hylatomus pileatus")
spec_name_change_df[25,] <- c("Colaptes auratus auratus", "Colaptes auratus")

# match AVONET data either on BBS species name or on BL species name:

BBS_sel_species_df <- bind_rows( # join on either BBS name or BL name and bind rows
  
  # join on Birdlife species name:
  data.frame("BBS_species" = BBS_species_final) %>% 
  left_join(spec_name_change_df, by = c("BBS_species" = "BBS_name")) %>% 
  inner_join(avonet_BL, by = c("BL_name" = "Species1")),
  
  # join on BBS species name:
  data.frame("BBS_species" = BBS_species_final) %>% 
  left_join(spec_name_change_df, by = c("BBS_species" = "BBS_name")) %>% 
  inner_join(avonet_BL, by = c("BBS_species" = "Species1"))
  ) %>% 
  mutate(BL_name = ifelse(is.na(BL_name), BBS_species, BL_name)) %>% 
  arrange(BBS_species)

write.csv(BBS_sel_species_df, file = file.path(data_dir, "AVONET_BBS_species.csv"),
          row.names = FALSE)
