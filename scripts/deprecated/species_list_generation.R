#### Generating species list and preliminary species inputs for requesting aligned data from raw EBBA1 and EBBA2

# datashare:
datashare_EBCC <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "EBCC")
datashare_Chelsa <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "envidat", "biophysical", "CHELSA_V2", "global")
transfer <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Transfer", "EBBA_niche_vs_range_shift", "Data")

library(dplyr)
library(sf)

# EBBA1: 3950 cells and 495 species according to https://ebba2.info/about/ebba1/
# (3912 locations in GBIF data ?)
# EBBA2: 5303 grid cells (5079 with data?), 596 species (data in csv-file, locations in shapefile)

# EBBA1 data: ----
EBBA1 <- read.csv(file.path(datashare_EBCC, "EBBA1", "EBBA1.csv"), sep="\t", header=T)
nrow(EBBA1)# 1'339'711
# EBBA1.csv contains for each species presence AND absence records
# this information is contained in column "issue"
# use only EBBA1 presences:
EBBA1 <- EBBA1 %>% 
  filter(issue %in% c("COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;INDIVIDUAL_COUNT_INVALID",
                      "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_FUZZY;INDIVIDUAL_COUNT_INVALID"))
nrow(EBBA1)# 355'488

# explorations regarding EBBA 1: ----

library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)

unique(EBBA1$issue) # 7 issues
table(EBBA1$issue)
# https://gbif.github.io/gbif-api/apidocs/org/gbif/api/vocabulary/OccurrenceIssue.html:
# COORDINATE ROUNDED = Original coordinate modified by rounding to 5 decimals.
# COUNTRY_DERIVED_FROM_COORDINATES = The interpreted country is based on the coordinates, not the verbatim string information.
# RECORDED_DATE_INVALID = A (partial) invalid date is given, such as a non existing date, zero month, etc.
# TAXON_MATCH_HIGHERRANK = Matching to the taxonomic backbone can only be done on a higher rank and not the scientific name.
# TAXON_MATCH_FUZZY = Matching to the taxonomic backbone can only be done using a fuzzy, non exact match.
# INDIVIDUAL_COUNT_INVALID = The individual count value is not a positive integer.
# TAXON_MATCH_NONE = Matching to the taxonomic backbone cannot be done because there was no match at all, or several matches with too little information to keep them apart (potentially homonyms).
# BASIS_OF_RECORD_INVALID = The given basis of record is impossible to interpret or significantly different from the recommended vocabulary.

# check for some species:
unique(EBBA1$species) # 496 species

world <- ne_coastline(returnclass = "sf")

Razorbill <- EBBA1 %>% 
  filter(species == "Alca torda")
Razorbill_sf <- st_as_sf(Razorbill, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
ggplot(world) +
  geom_sf(color = "gray50") +
  geom_sf(data = Razorbill_sf, aes(color = issue), size = 1) +
  lims(x = c(-30, 50), y = c(30, 80)) +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))

Eagle_owl <- EBBA1 %>% 
  filter(species == "Bubo bubo")
Eagle_owl_sf <- st_as_sf(Eagle_owl, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
ggplot(world) +
  geom_sf(color = "gray50") +
  geom_sf(data = Eagle_owl_sf, aes(color = issue), size = 1) +
  lims(x = c(-30, 50), y = c(30, 80)) +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))

Capercaillie <- EBBA1 %>% 
  filter(species == "Tetrao urogallus")
Capercaillie_sf <- st_as_sf(Capercaillie, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
ggplot(world) +
  geom_sf(color = "gray50") +
  geom_sf(data = Capercaillie_sf, aes(color = issue), size = 1) +
  lims(x = c(-30, 50), y = c(30, 80)) +
  theme(legend.position = "bottom") +
  guides(colour=guide_legend(nrow=2,byrow=TRUE))
# species absence = issue: "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID"
# species presence = issue: "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;INDIVIDUAL_COUNT_INVALID"

# species with other than those 2 issues:
issues_spec <- EBBA1 %>% 
  filter(issue %in% c("COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_FUZZY",
                      "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_FUZZY;INDIVIDUAL_COUNT_INVALID",
                      "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_HIGHERRANK",
                      "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_HIGHERRANK;INDIVIDUAL_COUNT_INVALID",
                      "TAXON_MATCH_NONE;BASIS_OF_RECORD_INVALID")) %>% 
  distinct(species) %>% 
  pull(species)# 15 species

# maps:
for(spec in issues_spec[-1]){
  print(spec)
  spec_data <- EBBA1 %>% 
    filter(species == !!spec)
  spec_sf <- st_as_sf(spec_data, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
  print(ggplot(world) +
          geom_sf(color = "gray50") +
          geom_sf(data = spec_sf, aes(color = issue), size = 1) +
          ggtitle(paste(spec)) +
          lims(x = c(-30, 50), y = c(30, 80)) +
          theme(legend.position = "bottom") +
          guides(colour=guide_legend(nrow=2, byrow=TRUE)))
}

# -> similar as before:
# presence = "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_FUZZY;INDIVIDUAL_COUNT_INVALID"
# absence = "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_FUZZY"
# (Vanellus leucurus: only absences)

# other issues don't yield species level information:

# "TAXON_MATCH_HIGHERRANK" issue:
EBBA1 %>% 
  filter(issue == "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_HIGHERRANK" |
           issue == "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_HIGHERRANK;INDIVIDUAL_COUNT_INVALID") %>% 
  distinct(genus) # 2 genera: Stictocarbo, Picoides
# are there information for these genera on species level at all:
EBBA1 %>% 
  filter(genus == "Stictocarbo") %>%
  distinct(taxonRank) %>% 
  pull(taxonRank) # no
EBBA1 %>% 
  filter(genus == "Picoides") %>%
  distinct(taxonRank) %>% 
  pull(taxonRank) # yes
# which: 
EBBA1 %>% 
  filter(genus == "Picoides") %>%
  distinct(species)

# "TAXON_MATCH_NONE;BASIS_OF_RECORD_INVALID" issue:
EBBA1 %>% 
  filter(issue == "TAXON_MATCH_NONE;BASIS_OF_RECORD_INVALID") # 1 entry


# EBBA2 data: ----
EBBA2_dt <- read.csv(file.path(datashare_EBCC, "EBBA2", "ebba2_data_occurrence_50km.csv"), sep=";", header=T)
nrow(EBBA2_dt) # 579'263
# each row indicates the presence of a given species in a given square
# cell50x50 = square code from the 50x50 km EBBA2 grid
EBBA2_sf <- read_sf(file.path(datashare_EBCC, "EBBA2", "ebba2_grid50x50_v1", "ebba2_grid50x50_v1.shp"))
nrow(EBBA2_sf) # 5303

# spatial processing of EBBA data: ----

# CRS of EBBA2 grid is ETRS89-extended / LAEA Europe (EPSG:3035)
# EBBA1 data have been downloaded from GBIF, which only allows data in WGS84 longitude/latitude (EPSG:4326)
# thus, EBBA1 data must be transformed from EPSG:4326 to EPSG:3035

# I use EPSG:3035 throughout the project, so also transform Chelsa data to EPSG:3035
# European Environment Agency recommends usage of Lambert azimuthal equal-area projection for 
# European mapping for statistical analysis and display (https://en.wikipedia.org/wiki/European_grid)

EBBA1_sf <- st_as_sf(EBBA1, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
EBBA1_sf_EPSG3035 <- st_transform(EBBA1_sf, crs = 3035)
rm(EBBA1_sf)

# subset to save storage space:
EBBA1_sf_EPSG3035_ss <- EBBA1_sf_EPSG3035 %>% 
  select(occurrenceID, countryCode, species, scientificName, issue)
#st_write(EBBA1_sf_EPSG3035_ss, file.path(transfer, "input_data", "EBCC", "EBBA1_EPSG3035_KS.shp"))
#st_write(EBBA1_sf_EPSG3035_ss, "EBBA1_EPSG3035_KS.shp")

# convert EBBA 2 grid cells to point data (as EBBA1) and add EBBA 2 data:
EBBA2_sf_points <- st_centroid(EBBA2_sf)
length(unique(EBBA2_sf_points$cell50x50)) # 5,303

EBBA2_sf_pts_dt <- full_join(EBBA2_sf_points, EBBA2_dt, by = "cell50x50") %>% 
  rename(species = birdlife_scientific_name) %>% 
  select(-birdlife_code)

nrow(EBBA2_sf_pts_dt) # 579,487 (more because there are 224 cells for which no presence is recorded)
#st_write(EBBA2_sf_pts_dt, file.path(transfer, "input_data", "EBCC", "EBBA2", "EBBA2_EPSG3035_KS.shp"))
#st_write(EBBA2_sf_pts_dt, "EBBA2_EPSG3035_KS.shp")


# drop non-comparable cells according to EBBA2 methods chapter: ----

# "We decided to restrict the change map to a well defined, continuous geographical area based on coverage in EBBA1. 
# Cyprus, the Asian part of Turkey and the Canary Islands were not included at all in EBBA1. 
# Moreover, a very large area including most of the European parts of Russia and Kazakhstan, Georgia, Armenia and Azerbaijan was only partially surveyd at the time.
# Consequently, this vast area in eastern Europe was excluded from the geographical area used in the change maps.
# More specifically, 50-km squares that lie totally or mostly (>70% of their area) within the above mentioned areas
# were not included."

# (GBIF countryCode follows ISO 3166-1-alpha-2 standard)

# add cell-IDs to EBBA1 (spatially join EBBA1 observations (with country code) and EBBA2 grid cells (with cellID)):
EBBA1_sf_EPSG3035_cellID <- st_join(EBBA1_sf_EPSG3035_ss, EBBA2_sf, join = st_intersects)

# add country codes to EBBA2:
EBBA2_sf_pts_dt_country_code <- left_join(EBBA2_sf_pts_dt, 
                                          distinct(st_drop_geometry(EBBA1_sf_EPSG3035_cellID[, c("cell50x50", "countryCode")])), 
                                          by = c("cell50x50" = "cell50x50"))

# countries / areas to exclude because of low sampling effort:
low_survey_effort <- c("CY", "TR", "RU", "KZ", "GE", "AM", "AZ", "IR") # IR = Iran, some observations are located just across the border

# exclude cell-IDs Canary Islands: 
Canaries <- c("28RDU2", "27RYL3", "27RYM3", "28RBR1", "28RCS2", "28RCR1", "28RBS4", "28RBS1", 
              "28RDR3", "28RDR1", "28RCS4", "28RCS3", "28RES3", "28RES2", "28RDS4", "28RDS2", 
              "28RFT1", "28RFS2", "28RFS1", "28RES4", "28RFT4", "28RFT2")

# include cell-IDs Russia between PL and LT:
Russia_Europe <- c("34UDG4", "34UDF1", "34UDF3", "34UEF1", "34UEF3", "34UEF4")

# include cell-IDs European Turkey: 
Turkey_Europe <- c("35TMG4", "35TNG2", "35TNG4", "35TMF3", "35TNF1", "35TNF3", "35TPF1", 
                   "35TPF3", "35TMF4", "35TNF2", "35TNF4", "35TPF2", "35TPF4", "35TME3")
# (cellIDs that lie totally or mostly (>70% of their area) within the above mentioned areas xx?)

# clip EBBA1 data to cells to keep:
EBBA1_clipped_sf <- EBBA1_sf_EPSG3035_cellID %>%
  mutate(include = case_when(cell50x50 %in% Canaries ~ FALSE,
                             cell50x50 %in% Turkey_Europe ~ TRUE,
                             cell50x50 %in% Russia_Europe ~ TRUE,
                             countryCode %in% low_survey_effort ~ FALSE,
                             TRUE ~ TRUE)) %>%
  filter(include == TRUE & !is.na(cell50x50)) %>% # also exclude observations from EBBA1 that are not inside an EBBA2 cell (151 observations, all along the coast)
  select(-include)

nrow(EBBA1_sf_EPSG3035_cellID) # 355'488
nrow(EBBA1_clipped_sf) # 330'337
#st_write(EBBA1_clipped_sf, "EBBA1_clipped_KS.shp")

# clip EBBA2 data to cells to keep:
EBBA2_clipped_sf <- EBBA2_sf_pts_dt_country_code %>%
  mutate(include = case_when(cell50x50 %in% Canaries ~ FALSE,
                             cell50x50 %in% Turkey_Europe ~ TRUE,
                             cell50x50 %in% Russia_Europe ~ TRUE,
                             countryCode %in% low_survey_effort ~ FALSE,
                             is.na(countryCode) ~ FALSE, # cell not included in EBBA1 -> drop #xx
                             TRUE ~ TRUE)) %>% 
  filter(include == TRUE) %>% 
  select(-include)

nrow(EBBA2_sf_pts_dt_country_code) # 579'487
nrow(EBBA2_clipped_sf) # 365'389
#st_write(EBBA2_clipped_sf, "EBBA2_clipped_KS.shp")


# unify taxonomy between EBBAs: ----

# keep only species and cell-ID:
EBBA1_clipped <- EBBA1_clipped_sf %>% 
  st_drop_geometry() %>%
  select(cell50x50, species)

EBBA2_clipped <- EBBA2_clipped_sf %>% 
  st_drop_geometry() %>%
  select(cell50x50, species)

# based on EBBA2 methods chapter:
# "EBBA1 data were attributed, where possible, to the species recognised in EBBA2 (Table 6)"

EBBA1_taxunif <- EBBA1_clipped
EBBA2_taxunif <- EBBA2_clipped

# account for name change from EBBA1 to EBBA2:
EBBA1_taxunif[which(EBBA1_taxunif$species=="Acanthis hornemanni"),'species'] <- "Acanthis flammea"
EBBA1_taxunif[which(EBBA1_taxunif$species=="Oceanodroma castro"),'species'] <- "Hydrobates castro"
EBBA1_taxunif[which(EBBA1_taxunif$species=="Parus lugubris"),'species'] <- "Poecile lugubris"
EBBA1_taxunif[which(EBBA1_taxunif$species=="Regulus ignicapillus"),'species'] <- "Regulus ignicapilla"
EBBA1_taxunif[which(EBBA1_taxunif$species=="Serinus citrinella"),'species'] <- "Carduelis citrinella"

# remove all non-comparable species according to Table 6 in EBBA2 methods chapter:
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Branta hutchinsii"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Calonectris diomedea"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Calonectris borealis"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Calonectris borealis"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Hydrobates monteiroi"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Larus cachinnas"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Larus michahellis"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Larus michahellis"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Picus viridis"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Picus sharpei"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Picus viridis"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Lanius excubitor"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Lanius meridionalis"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Lanius excubitor"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Poecile hyrcanus"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Oenanthe xanthoprymna"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Oenanthe chrysopygia"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Oenanthe xanthoprymna"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Phylloscopus collybita"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Phylloscopus tristis"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Phylloscopus ibericus"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Phylloscopus collybita"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Phylloscopus nitidus"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Sylvia sarda"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Sylvia balearica"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Sylvia sarda"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Sylvia cantillans"),]
EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Sylvia subalpina"),]
EBBA1_taxunif <- EBBA1_taxunif[which(EBBA1_taxunif$species!="Sylvia cantillans"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Regulus  madeirensis"),]

EBBA2_taxunif <- EBBA2_taxunif[which(EBBA2_taxunif$species!="Carduelis corsicana"),]

length(unique(EBBA1_taxunif$species)) # 445 species left
length(unique(EBBA2_taxunif$species)) # 517 species left


# Species filtering: ----

EBBA1_filtered <- EBBA1_taxunif
EBBA2_filtered <- EBBA2_taxunif

## 1. exclude pelagic specialists (according to Wilman et al. 2014): ----

SeaBirds <- read.csv2(file.path(transfer, "input_data", "BirdFuncDat.txt"), header=TRUE, sep="\t") #BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
SeaBirds <- subset(SeaBirds, PelagicSpecialist==1)
SeaBirds <- SeaBirds$Scientific

EBBA1_filtered <- EBBA1_filtered[which(!EBBA1_filtered$species %in% SeaBirds),]
EBBA2_filtered <- EBBA2_filtered[which(!EBBA2_filtered$species %in% SeaBirds),]

length(unique(EBBA1_filtered$species)) # 417 species left
length(unique(EBBA2_filtered$species)) # 491 species left


## 2. exclude rare species with n<20 occurrences in any of the two atlas periods: ----

EBBA1_filtered <- EBBA1_filtered %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup

length(unique(EBBA1_filtered$species)) # 371 species left

EBBA2_filtered <- EBBA2_filtered %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20) %>% 
  ungroup

length(unique(EBBA2_filtered$species)) # 390 species left


## 3. exclude very common species with >90% prevalence in both atlas periods: ----
# = occur in >90% EBBA cells:

nEBBAcells <- length(unique(EBBA1_clipped$cell50x50)) # 2845 (EBBA2 methods chapter: "change map thus encompassed a total of 2,972 50-km squares" -> more because I excluded cells from EBBA2 not in EBBA1, the 2,972 squares contain cells not included in change analysis due to insufficient coverage?!
length(unique(EBBA2_clipped$cell50x50)) # 2845

# which species are excluded:
EBBA1_filtered %>% 
  filter(n_occurrences > 0.9 * nEBBAcells) %>% 
  distinct(species)
EBBA2_filtered %>% 
  filter(n_occurrences > 0.9 * nEBBAcells) %>% 
  distinct(species)

# exclude species:
EBBA1_filtered <- EBBA1_filtered %>% 
  filter(n_occurrences <= 0.9 * nEBBAcells)

length(unique(EBBA1_filtered$species)) # 365 species left

EBBA2_filtered <- EBBA2_filtered %>% 
  filter(n_occurrences <= 0.9 * nEBBAcells)

length(unique(EBBA2_filtered$species)) # 379 species left


## 4. keep only those species that occur in both atlas periods: ----

EBBA1_filtered <- EBBA1_filtered %>% 
  filter(species %in% EBBA2_filtered$species)

length(unique(EBBA1_filtered$species)) # 325 species left

EBBA2_filtered <- EBBA2_filtered %>% 
  filter(species %in% EBBA1_filtered$species)

length(unique(EBBA2_filtered$species)) # 325 species left

## how often does each species occur in each of the EBBAs:
# EBBA1_filtered %>%
#   arrange(-n_occurrences) %>% 
#   distinct(species, n_occurrences)
# 
# EBBA2_filtered %>%
#   arrange(-n_occurrences) %>%
#   distinct(species, n_occurrences)

#save(EBBA1_filtered, EBBA2_filtered, file = file.path("Data", "EBBA1_EBBA2_filtered_steps1-4.RData"))



## 5. exclude species that have <70% ? of climatic niche covered in Europe: ----

# packages:

library(archive) # to extract zipped Birdlife ranges
library(gdalUtilities)
library(terra)
library(dismo) # requires spatial data formats of raster package -> no transition to terra possible
library(raster)
library(stringr)
library(ecospat)
library(ade4)

# setup:

datashare_Birdlife <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "Birdlife")
load(file = file.path("Data", "EBBA1_EBBA2_filtered_steps1-4.RData"))

# species left after filtering step 4 (pelagic specialists, rare and very common species 
# and species occurring only in one atlas period excluded)
species_filtered <- sort(sub(" ", "_", unique(EBBA2_filtered$species)))


### 1. Project Chelsa data of EBBA1 period (1981 - 1990): ----

# note: this code section was run from the cluster
# note: auf ecoc9 module load R/4.1.0-foss-2021a  machst und dann die neuere R version nutzt. Du musst dann alle notwendigen Pakete nochmal neu installieren


library(gdalUtilities)
#xx
#datashare <- file.path("/mnt","ibb_share","zurell","envidat","biophysical","CHELSA_V2","global") 
datashare <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "envidat","biophysical","CHELSA_V2","global")

chelsa_tifs <- list.files(datashare, full.names = FALSE, pattern = paste0("(", paste(1981:1990, collapse = "|"), ")_V.2.1.tif")) # 480

# folder to store reprojected CHELSA data:
#chelsa_birdlife_path <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Chelsa_projected")
chelsa_birdlife_path <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "users$", "schifferle1", "Documents", "EBBA_Niche_vs_Range_shifts", "Data", "Chelsa_projected_GDAL3_4")


# create folder if it doesn't exist yet:
if(!dir.exists(chelsa_birdlife_path)){
  dir.create(chelsa_birdlife_path, recursive = TRUE)
}

# create list of file paths for the reprojected data:
names <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_50km.tif")
names <- file.path(chelsa_birdlife_path, names)

# global equal area projection used by SoilGrids (Homolosine projection applied to the WGS84 datum)
# from: https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection
homolosine <- 'PROJCS["Homolosine", 
                     GEOGCS["WGS 84", 
                            DATUM["WGS_1984", 
                                  SPHEROID["WGS 84",6378137,298.257223563, 
                                           AUTHORITY["EPSG","7030"]], 
                                  AUTHORITY["EPSG","6326"]], 
                            PRIMEM["Greenwich",0, 
                                   AUTHORITY["EPSG","8901"]], 
                            UNIT["degree",0.0174532925199433, 
                                 AUTHORITY["EPSG","9122"]], 
                            AUTHORITY["EPSG","4326"]], 
                     PROJECTION["Interrupted_Goode_Homolosine"], 
                     UNIT["Meter",1]]'

# reproject downloaded CHELSA layers:
options(warn = 1) # default warn = 0; 1 = warnings are printed as they occur

for(i in 1:length(chelsa_tifs)){
  
  print(paste(i, "of", length(chelsa_tifs)))
  
  # reproject Chelsa data:
  gdalUtilities::gdalwarp(srcfile = file.path(datashare, 
                                              chelsa_tifs[i]),
                          dstfile = names[i],
                          overwrite = TRUE,
                          tr = c(50000, 50000), # target resolution
                          r = "bilinear", # resampling method
                          t_srs = homolosine)
}


# In CPL_gdalwarp(source, destination, options, oo, doo,  ... :
# GDAL Error 1: Too many points (526 out of 529) failed to transform, unable to compute output bounds.
# 2: In CPL_gdalwarp(source, destination, options, oo, doo,  ... :
#    GDAL Message 1: Unable to compute source region for output window 588,0,13,10, skipping.

#xx

#---
# whether I get Error messages depends on GDAL / PROJ version:
chelsa_test_tif <- "CHELSA_pr_01_1981_V.2.1.tif"
names <- paste0(strsplit(chelsa_test_tif, "\\.tif"), "_50km.tif")
names <- file.path("Data", "Chelsa_projected_GDAL_3_4", names)

# reproject Chelsa data:
gdalUtilities::gdalwarp(srcfile = file.path("Data", chelsa_test_tif),
                        dstfile = names,
                        overwrite = TRUE,
                        tr = c(50000, 50000), # target resolution
                        r = "bilinear", # resampling method
                        t_srs = homolosine)

# Cluster:
#Linking to GEOS 3.9.0, GDAL 3.2.2, PROJ 7.2.1

# my laptop:
#Linking to GEOS 3.10.2, GDAL 3.4.1, PROJ 7.2.1; sf_use_s2() is TRUE
# Jettes other cluster option:
# Linking to GEOS 3.9.1, GDAL 3.3.0, PROJ 8.0.1; sf_use_s2() is TRUE

# results differ!

# GDAL Error 1: PROJ: igh: Invalid latitude
# may be due to proj version difference

# testen, ob Ergebnis gleich ist:
gdal32 <- rast(file.path("Data", "Chelsa_Birdlife_ranges_50km","CHELSA_pr_01_1981_V.2.1_50km.tif"))
gdal33 <- rast(file.path("Data", "Chelsa_projected_GDAL_3_3", "CHELSA_pr_01_1981_V.2.1_50km.tif"))
gdal34 <- rast(file.path("Data", "Chelsa_projected_GDAL_3_4", "CHELSA_pr_01_1981_V.2.1_50km.tif"))

plot(gdal32)
plot(gdal33)
plot(gdal34)

all.equal(gdal33, gdal34)
compareGeom(gdal33, gdal34)
m <- minmax(gdal33 - gdal34)
all(abs(m) < 1e-7)

#---



### 2. Create mask to extract relevant Chelsa data: ----

# mask covers Birdlife range polygons of all species left after filtering step 4 and EBBA grid
# (resolution = 50 km, as EBBAs)

# extract zipped Birdlife ranges:

# all files in zipped archive:
#Birdlife_files <- archive(file.path(datashare_Birdlife, "All_Shapefiles.7z"))
Birdlife_files <- archive(file.path("Data", "Birdlife_shapes", "All_Shapefiles.7z"))

# files of relevant bird species:
Birdlife_rel_files <- lapply(X = species_filtered, 
                             FUN = function(x) Birdlife_files$path[which(grepl(pattern = paste0(x, "_"), 
                                                                               x = Birdlife_files$path))]) # shapefiles with 4 sub-files per species

# some species have other names in Birdlife data than in EBBA data:
spec_synonyms_df <- data.frame("EBBA_name" = species_filtered[which(lengths(Birdlife_rel_files) == 0)])
spec_synonyms_df$BL_name <- NA
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Acanthis_flammea")] <- "Carduelis_flammea"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Cecropis_daurica")] <- "Hirundo_daurica"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Cyanistes_caeruleus")] <- "Parus_caeruleus"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Emberiza_calandra")] <- "Miliaria_calandra"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Iduna_pallida")] <- "Hippolais_pallida"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Linaria_cannabina")] <- "Carduelis_cannabina"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Linaria_flavirostris")] <- "Carduelis_flavirostris"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Lophophanes_cristatus")] <- "Parus_cristatus"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Periparus_ater")] <- "Parus_ater"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Poecile_cinctus")] <- "Parus_cinctus"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Poecile_lugubris")] <- "Parus_lugubris"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Poecile_montanus")] <- "Parus_montanus"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Poecile_palustris")] <- "Parus_palustris"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Ptyonoprogne_rupestris")] <- "Hirundo_rupestris"
spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == "Spinus_spinus")] <- "Carduelis_spinus"

# data for "Passer_italiae" and "Sylvia_ruppeli" not included in zipped archive on datashare
# stored data on Github / transfer xx

# add files of species with synonym name:
Birdlife_rel_files[which(lengths(Birdlife_rel_files) == 0)] <- lapply(X = spec_synonyms_df$BL_name, 
       FUN = function(x) Birdlife_files$path[which(grepl(pattern = x, x = Birdlife_files$path))])
Birdlife_rel_files <- unlist(Birdlife_rel_files)

# extract files:
archive_extract(#archive = file.path(datashare_Birdlife, "All_Shapefiles.7z"),
                archive = file.path("Data", "Birdlife_shapes", "All_Shapefiles.7z"),
                dir = file.path("Data", "Birdlife_shapes"),
                #dir = file.path(transfer, "IUCN_range_polys"),
                files = Birdlife_rel_files)


# for each species: convert Birdlife range shapefile to raster, 
# project raster and change resolution to 50 km:

# input shapefiles:
Birdlife_rel_files_shp <- list.files(file.path("Data", "Birdlife_shapes", "All_shapefiles"),
                                     full.names = FALSE, pattern = ".shp")

# output tifs:
Birdlife_tifs <- sub(".shp", ".tif", Birdlife_rel_files_shp)

# global equal area projection used by SoilGrids (Homolosine projection applied to the WGS84 datum)
# from: https://www.isric.org/explore/soilgrids/faq-soilgrids#How_can_I_use_the_Homolosine_projection
homolosine <- 'PROJCS["Homolosine", 
                     GEOGCS["WGS 84", 
                            DATUM["WGS_1984", 
                                  SPHEROID["WGS 84",6378137,298.257223563, 
                                           AUTHORITY["EPSG","7030"]], 
                                  AUTHORITY["EPSG","6326"]], 
                            PRIMEM["Greenwich",0, 
                                   AUTHORITY["EPSG","8901"]], 
                            UNIT["degree",0.0174532925199433, 
                                 AUTHORITY["EPSG","9122"]], 
                            AUTHORITY["EPSG","4326"]], 
                     PROJECTION["Interrupted_Goode_Homolosine"], 
                     UNIT["Meter",1]]'

options(warn = 1) # default warn = 0; 1 = warnings are printed as they occur

for(i in 1:length(Birdlife_rel_files_shp)){# 325
  
  print(paste(i, Birdlife_rel_files_shp[i]))
  
  # rasterize shapefile:
  gdalUtilities::gdal_rasterize(src_datasource = file.path("Data", "Birdlife_shapes", "All_shapefiles", Birdlife_rel_files_shp[i]),
                                where = "SEASONAL = '1' OR SEASONAL = '2'", # SEASONAL = 1: resident throughout the year, SEASONAL = 2: breeding season
                                dst_filename = file.path("Data", "Birdlife_shapes", "Raster", Birdlife_tifs[i]),
                                burn = 1,
                                tr = c(0.1, 0.1), # target resolution in degrees (same unit as src_datasource) (0.1° enough since final resolution will be 50 km)
                                a_nodata = -99999) # value for cells with missing data
  
  # project and change resolution to 50 km (as EBBA):
  gdalUtilities::gdalwarp(srcfile = file.path("Data", "Birdlife_shapes", "Raster", Birdlife_tifs[i]), 
                          dstfile = file.path("Data", "Birdlife_shapes", "Raster", sub("\\.tif", "_50km.tif", Birdlife_tifs[i])),
                          overwrite = TRUE,
                          tr = c(50000, 50000), # target resolution in meters (same as unit of target srs)
                          r = "max", # "near", resampling method
                          t_srs = homolosine # target spatial reference
  )
} 

# same for EBBA grid:

# rasterize shapefile:
gdalUtilities::gdal_rasterize(src_datasource = file.path(datashare_EBCC, "EBBA2", "ebba2_grid50x50_v1", "ebba2_grid50x50_v1.shp"),
                              dst_filename = file.path("Data", "EBBA2.tif"),
                              burn = 1,
                              tr = c(50000, 50000), # target resolution in degrees (same unit as src_datasource)
                              a_nodata = -99999) # value for cells with missing data

# project and change resolution to 50 km:
gdalUtilities::gdalwarp(srcfile = file.path("Data", "EBBA2.tif"), 
                        dstfile = file.path("Data", "EBBA2_hom.tif"), 
                        overwrite = TRUE,
                        tr = c(50000, 50000), # target resolution in meters (same as unit of target srs)
                        r = "max", # "near", resampling method
                        t_srs = homolosine # target spatial reference
                        )

# overlay Birdlife ranges of all species and EBBA grid:

# Birdlife ranges of all species:
ranges <- lapply(file.path("Data", "Birdlife_shapes", "Raster", sub("\\.tif", "_50km.tif", Birdlife_tifs)),
                 rast)
# add EBBA grid:
ranges <- append(ranges, rast(file.path("Data", "EBBA2_hom.tif")))

# overlay rasters, add 50 km buffer and resample to match projected Chelsa rasters:
# (buffering makes sure that despite of resampling there is in the end Chelsa data for each grid point available)
# (resampling aligns raster origins, necessary for masking)

ranges_combined_bf <- do.call(terra::merge, ranges) %>% 
  buffer(width = 100000) %>% # Warning: "[merge] rasters did not align and were resampled" -> fine #xx
  #resample(y = rast(file.path("Data", "Chelsa_Birdlife_ranges_50km", "CHELSA_pr_01_1981_V.2.1_50km.tif")),
  #         method = "near") # projected Chelsa tif as template
  resample(y = rast(file.path("//ibb-fs01.ibb.uni-potsdam.de", "users$",
                              "schifferle1", "Documents", "EBBA_Niche_vs_Range_shifts", "Data", "Chelsa_projected_GDAL3_3", "CHELSA_pr_01_1981_V.2.1_50km.tif")),
                    method = "near") %>% 
  crop(y = rast(file.path("//ibb-fs01.ibb.uni-potsdam.de", "users$",
                         "schifferle1", "Documents", "EBBA_Niche_vs_Range_shifts", "Data", "Chelsa_projected_GDAL3_3", "CHELSA_pr_01_1981_V.2.1_50km.tif")))

# save mask:
writeRaster(ranges_combined_bf, 
            filename = file.path("Data", "Birdlife_ranges_mask_310123.tif"),
            overwrite = TRUE,
            NAflag = FALSE)

ext(ranges_combined_bf)
ext(rast(file.path("//ibb-fs01.ibb.uni-potsdam.de", "users$",
                   "schifferle1", "Documents", "EBBA_Niche_vs_Range_shifts", "Data", "Chelsa_projected_GDAL3_3", "CHELSA_pr_01_1981_V.2.1_50km.tif")))
### exploration: compare EBBA1 range with range maps from Birdlife: ----

library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)

# Birdlife range map:
BL_range_sf <- read_sf(file.path("Data", "Birdlife_shapes", Birdlife_rel_files[87]))

# EBBA distribution:
EBBA1_range_sf <- EBBA1_filtered %>%
  filter(species == sub("_", " ", species_filtered[2])) %>%
  left_join(EBBA2_sf, by = c("cell50x50" = "cell50x50")) %>%
  st_as_sf()

world <- ne_countries(scale = "small", returnclass = "sf")

ggplot(world) +
  geom_sf(color = "gray50") +
  geom_sf(data = st_transform(BL_range_sf, crs = st_crs(world)), fill = "red", alpha = 0.5) +
  geom_sf(data = st_transform(EBBA1_range_sf, crs = st_crs(world)), fill = "blue", colour = NA)


### 3. Mask Chelsa data and calculate bioclim variables: ----

# note: this code section was run from the cluster

library(terra)
library(dismo) # requires spatial data formats of raster package -> no transition to terra possible
library(raster)
library(stringr)

datashare <- file.path("/mnt","ibb_share","zurell","envidat","biophysical","CHELSA_V2","global") 
chelsa_tifs <- list.files(datashare, full.names = FALSE, pattern = paste0("(", paste(1981:1990, collapse = "|"), ")_V.2.1.tif")) # 480

# folder where projected CHELSA data are stored:
chelsa_birdlife_path <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Chelsa_projected")

# list of file paths of projected data:
names <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_50km.tif")
names <- file.path(chelsa_birdlife_path, names)

# mask:
ranges_mask <- rast(file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Birdlife_ranges_mask2.tif"))

# folder to store masked CHELSA data:
chelsa_masked_path <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Chelsa_masked")
# create folder if it doesn't exist yet:
if(!dir.exists(chelsa_masked_path)){
  dir.create(chelsa_masked_path, recursive = TRUE)
}

# create list of file paths for the reprojected and masked data:
names_masked <- paste0(unlist(lapply(chelsa_tifs, FUN = function(x) {strsplit(x, "\\.tif")})), "_50km_masked.tif")
names_masked <- file.path(chelsa_masked_path, names_masked)

for(i in 1:length(chelsa_tifs)){
  
  print(i)
  
  # mask reprojected raster so all non-terrestrial areas become NA:
  names[i] %>%
    rast %>%
    mask(mask = ranges_mask) %>%
    writeRaster(names_masked[i], overwrite = TRUE)
}

# calculate bioclim variables from Chelsa data:

vars <- c("pr", "tas", "tasmin", "tasmax")
months <- str_pad(1:12, width = 2, pad = "0")
years <- 1981:1990

# Chelsa files:
chelsa_masked_files <- list.files(chelsa_masked_path,
                                  pattern = "tif$",
                                  full.names = TRUE)

# load mask (using raster package since this is required by dismo::biovars())
mask_ranges_combined <- raster(file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Birdlife_ranges_mask2.tif"))

# create a spatial brick template using the mask and 
# create bricks to store the month-wise means across the selected years for each bioclim variable
out <- brick(mask_ranges_combined, values = FALSE)
pr_mean <- tasmin_mean <- tasmax_mean <- out

# loop over months:
for(j in as.numeric(months)){
  
  print(paste("month", j))
  
  # calculate precipitation mean for the current month (across all available years)
  pr_mean[[j]] <- chelsa_masked_files[which(grepl(paste0("pr_", months[j]),
                                                  chelsa_masked_files))] %>% 
    stack %>% 
    mean
  
  # calculate tasmin mean for the current month (across all available years)
  tasmin_mean[[j]] <- chelsa_masked_files[which(grepl(paste0("tasmin_", months[j]), 
                                                      chelsa_masked_files))] %>% 
    stack %>% 
    mean
  
  # calculate tasmax mean for the current month (across all available years)
  tasmax_mean[[j]] <- chelsa_masked_files[which(grepl(paste0("tasmax_", months[j]), chelsa_masked_files))] %>% 
    stack %>% 
    mean
}

# calculate rasters of bioclim variables:
biovars <- dismo::biovars(prec = pr_mean,  
                          tmin = tasmin_mean, 
                          tmax = tasmax_mean)

biovars_rast <- rast(biovars) # convert to terra object

# save tifs:
chelsa_bioclim_path <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "Chelsa_for_EBBA", "Bioclim_raster")
writeRaster(biovars_rast,
            filename = file.path(chelsa_bioclim_path, paste0("CHELSA_", names(biovars), "_", 
                                                             min(years), "_", max(years), "_", 
                                                             "50km.tif")), 
            overwrite = TRUE)


### 4. Climatic niche analyses: ---- 

# calculate stability index: how much of the species global climatic niche is covered in Europe

biovars_rast <- rast(file.path("Data", "Bioclim_raster_EBBA1_310123", paste0("CHELSA_", names(biovars), "_", 
                                           min(years), "_", max(years), "_", 
                                           "50km.tif")))

stability_df <- data.frame("species" = sub("_", " ", species_filtered),
                           "stability" = NA)

for(i in 1:length(species_filtered)){
  
  print(i)
  print(species_filtered[i])
  
  # EBBA1 occurrence point of species i:
  
  EBBA1_spec_sf <- EBBA1_clipped_sf %>% 
    st_transform(crs = homolosine) %>% # transform to Homolosine projection to match bioclim variables
    filter(species == sub("_", " ", species_filtered[i])) %>% 
    vect # convert to terra object
  
  # add bioclim variables to occurrences:
  EBBA1_spec_vars_df <- cbind(values(EBBA1_spec_sf[, c("species", "cell50x50")]),
                              terra::extract(biovars_rast, EBBA1_spec_sf))
  
  # Birdlife range of species i:
  
  BL_file <- grep(pattern = species_filtered[i], ranges_proj_files, value = TRUE)
  if(length(BL_file) == 0){ # Birdlife used another species name:
    BL_file <- grep(pattern = spec_synonyms_df$BL_name[which(spec_synonyms_df$EBBA_name == species_filtered[i])], 
                    ranges_proj_files, 
                    value = TRUE)
    }
  
  BL_pts_spec <- BL_file %>% 
    rast %>% # load raster
    as.points # convert raster to grid points
  
  # add bioclim variables to occurrences:
  BL_spec_vars_df <- cbind(data.frame(species = sub("_", " ", species_filtered[i]),
                                      cell50x50 = "BL"),
                           terra::extract(biovars_rast, BL_pts_spec)) %>% 
    filter(complete.cases(.)) # for few occurrence points "on the edge of the globe" no Chelsa data are available; because points are raster grid points, not cell centroids, and resampling is done, few points end up outside of global Chelsa range)
  
  # assess climate niche by using the first 2 PCA axes:
  # calibrating the PCA in the whole study area, including both ranges:
  pca.env <- dudi.pca(rbind(EBBA1_spec_vars_df, BL_spec_vars_df)[,4:22],
                      scannf=FALSE,
                      nf=2) # number of axes
  
  # predict the scores on the PCA axes:
  # EBBA1 used as z1 (corresponds to native distribution in tutorials)
  # Birdlife range used as z2 (corresponds to invasive distribution in tutorials)
  scores.globclim <- pca.env$li # PCA scores for EBBA + Birdlife distribution
  scores.clim.EBBA <- suprow(pca.env,EBBA1_spec_vars_df[,4:22])$li # PCA scores for EBBA distribution
  scores.clim.BL <- suprow(pca.env,BL_spec_vars_df[,4:22])$li # PCA scores for Birdlife distribution
  
  # calculate the Occurrence Densities Grid for EBBA and Birdlife ditribution:
  
  # EBBA1 distribution area:
  grid.clim.EBBA <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                          glob1 = scores.clim.EBBA, 
                                          sp = scores.clim.EBBA, # same for background an occurrences -> fine? xx
                                          R = 100, # grid resolution
                                          th.sp = 0)
  # Birdlife range area:
  grid.clim.BL <- ecospat.grid.clim.dyn(glob = scores.globclim,
                                        glob1 = scores.clim.BL, 
                                        sp = scores.clim.BL, # same for background an occurrences -> fine? xx
                                        R = 100, 
                                        th.sp = 0)
  
  # assess niche dynamics between EBBA1 (as z1) and Birdlife range map (as z2) 
  # => the stability index shows how much of the species global climatic niche 
  # is covered in Europe (should be e.g. >0.7):
  
  EBBA_BL_niche_dyn <- ecospat.niche.dyn.index(grid.clim.EBBA, 
                                               grid.clim.BL,
                                               intersection=NA) # analysis performed on EBBA + Birdlife extent
  stability_df$stability[i] <- EBBA_BL_niche_dyn$dynamic.index.w['stability']
  
  print(paste("Stability =", round(stability_df$stability[i],2)))
  
  # save plot of niche dynamics:
  jpeg(filename = file.path("plots", "EBBA_Birdlife_niche_dyn2", paste0(species_filtered[i], ".jpg")),
       quality = 100, height = 920, width = 920
       )

  ecospat.plot.niche.dyn(grid.clim.EBBA, 
                         grid.clim.BL, 
                         quant = 0.1,
                         interest = 2, # 1 = EBBA density, 2 = BL density
                         title = paste(species_filtered[i], "stability:", round(stability_df$stability[i],2)), 
                         name.axis1 = "PC1",
                         name.axis2 = "PC2")
  dev.off()
  # blue = stability
  # red = expansion
  # green = unfilling
}

write.csv(stability_df, 
          file = file.path("Data", "EBBA_niche_range_shifts_selected_species_310123.csv"),
          row.names = FALSE)


# Plot regarding stability threshold:
stability_df %>% 
  arrange(-stability) %>%  View

plot(sort(stability_df$stability, decreasing = TRUE), 
     ylab = "stability", xlab = "number of species", las = 1)
abline(h = 0.7, col = "red")
abline(h = 0.5, col = "blue")

# test:
alt <- read.csv(file.path("Data", "EBBA_niche_range_shifts_selected_species.csv"))
neu <- read.csv(file.path("Data", "EBBA_niche_range_shifts_selected_species_310123.csv"))

alt$species[!alt$species %in% neu$species]
neu$species[!neu$species %in% alt$species]
