# create comparable spatial data sets based on EBBA1 and EBBA2

library(dplyr)
library(sf)

# -------------------------------------------------- #
#            Preliminary data preparation:        ####
# -------------------------------------------------- #

# preliminary = analysis to identify suitable species for which we request EBBA1-EBBA2 change dataset
# (EBBA1-EBBA2 change dataset used for main niche vs. range shift analysis,
# but can also be run using publicly available EBBA1 and EBBA 2 data)


## load and explore EBBA data: -------------------------------------------------

# datashare:
datashare_EBCC <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "EBCC")

### EBBA1 data: --------------------------------------

EBBA1 <- read.csv(file.path(datashare_EBCC, "EBBA1", "EBBA1.csv"), sep = "\t", header = TRUE)
nrow(EBBA1)# 1'339'711


#### explorations regarding EBBA 1: ------------------

library("rnaturalearth")
library("rnaturalearthdata")
library(ggplot2)

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

# use only EBBA1 presences:

# EBBA1.csv contains presence AND absence records of each species
# this information is contained in column "issue"

EBBA1 <- EBBA1 %>% 
  filter(issue %in% c("COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;INDIVIDUAL_COUNT_INVALID",
                      "COORDINATE_ROUNDED;COUNTRY_DERIVED_FROM_COORDINATES;RECORDED_DATE_INVALID;TAXON_MATCH_FUZZY;INDIVIDUAL_COUNT_INVALID"))
nrow(EBBA1)# 355'488


### EBBA2 data: ---------------------------------------

EBBA2_dt <- read.csv(file.path(datashare_EBCC, "EBBA2", "ebba2_data_occurrence_50km.csv"), sep=";", header=T)
nrow(EBBA2_dt) # 579'263

# each row indicates the presence of a given species in a given square
# cell50x50 = square code from the 50x50 km EBBA2 grid

EBBA2_grid <- read_sf(file.path(datashare_EBCC, "EBBA2", "ebba2_grid50x50_v1", "ebba2_grid50x50_v1.shp"))
nrow(EBBA2_grid) # 5303


## spatial processing of EBBA data:---------------------------------------------

# CRS of EBBA2 grid is ETRS89-extended / LAEA Europe (EPSG:3035)
# EPSG:3035 is used throughout the project
# (European Environment Agency recommends usage of Lambert azimuthal equal-area projection for 
# European mapping for statistical analysis and display (https://en.wikipedia.org/wiki/European_grid))

# EBBA1 data have been downloaded from GBIF, which only allows data in WGS84 longitude/latitude (EPSG:4326)
# thus, EBBA1 data must be transformed from EPSG:4326 to EPSG:3035

EBBA1_sf <- st_as_sf(EBBA1, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% 
  st_transform(crs = 3035) %>% 
  select(occurrenceID, countryCode, species, scientificName, issue) # keep relevant attributes only

# convert EBBA 2 grid cells to point data (as EBBA1) and add EBBA 2 data:
EBBA2_sf <- st_centroid(EBBA2_grid) %>% 
  full_join(EBBA2_dt, by = "cell50x50") %>% 
  rename(species = birdlife_scientific_name) %>% 
  select(-birdlife_code)

# drop non-comparable cells according to EBBA2 methods chapter:

# EBBA2 methods chapter:
# "We decided to restrict the change map to a well defined, continuous geographical area based on coverage in EBBA1. 
# Cyprus, the Asian part of Turkey and the Canary Islands were not included at all in EBBA1. 
# Moreover, a very large area including most of the European parts of Russia and Kazakhstan, Georgia, Armenia and Azerbaijan was only partially surveyd at the time.
# Consequently, this vast area in eastern Europe was excluded from the geographical area used in the change maps.
# More specifically, 50-km squares that lie totally or mostly (>70% of their area) within the above mentioned areas
# were not included."

# drop cells based on country code or cell ID:
# (EBBA1 contains GBIF country codes, follows ISO 3166-1-alpha-2 standard; EBBA2 contains cell IDs)

# add cell-IDs to EBBA1 (spatially join EBBA1 observations (with country code) and EBBA2 grid cells (with cellID)):
EBBA1_sf_cellID <- st_join(EBBA1_sf, EBBA2_grid, join = st_intersects)

# add country codes to EBBA2:
EBBA2_sf_country_code <- left_join(EBBA2_sf,
                                   distinct(st_drop_geometry(EBBA1_sf_cellID[, c("cell50x50", "countryCode")])), 
                                   by = c("cell50x50" = "cell50x50"))

# countries / areas to exclude because of low sampling effort:
low_survey_effort <- c("CY", "TR", "RU", "KZ", "GE", "AM", "AZ", "IR") # IR = Iran, included since some observations are located just across the border in Iran

# exclude cell-IDs Canary Islands: 
Canaries <- c("28RDU2", "27RYL3", "27RYM3", "28RBR1", "28RCS2", "28RCR1", "28RBS4", "28RBS1", 
              "28RDR3", "28RDR1", "28RCS4", "28RCS3", "28RES3", "28RES2", "28RDS4", "28RDS2", 
              "28RFT1", "28RFS2", "28RFS1", "28RES4", "28RFT4", "28RFT2")

# include cell-IDs Russia between PL and LT:
Russia_Europe <- c("34UDG4", "34UDF1", "34UDF3", "34UEF1", "34UEF3", "34UEF4")

# include cell-IDs European Turkey: 
Turkey_Europe <- c("35TMG4", "35TNG2", "35TNG4", "35TMF3", "35TNF1", "35TNF3", "35TPF1", 
                   "35TPF3", "35TMF4", "35TNF2", "35TNF4", "35TPF2", "35TPF4", "35TME3")

# clip EBBA1 data to cells to keep:
EBBA1_clipped_sf <- EBBA1_sf_cellID %>%
  mutate(include = case_when(cell50x50 %in% Canaries ~ FALSE,
                             cell50x50 %in% Turkey_Europe ~ TRUE,
                             cell50x50 %in% Russia_Europe ~ TRUE,
                             countryCode %in% low_survey_effort ~ FALSE,
                             TRUE ~ TRUE)) %>%
  filter(include == TRUE & !is.na(cell50x50)) %>% # also exclude observations from EBBA1 that are not inside an EBBA2 cell (151 observations, all along the coast)
  select(cell50x50, species)

nrow(EBBA1_sf_cellID) # 355'488
nrow(EBBA1_clipped_sf) # 330'337

# clip EBBA2 data to cells to keep:
EBBA2_clipped_sf <- EBBA2_sf_country_code %>%
  mutate(include = case_when(cell50x50 %in% Canaries ~ FALSE,
                             cell50x50 %in% Turkey_Europe ~ TRUE,
                             cell50x50 %in% Russia_Europe ~ TRUE,
                             countryCode %in% low_survey_effort ~ FALSE,
                             is.na(countryCode) ~ FALSE, # drop cell not included in EBBA1
                             TRUE ~ TRUE)) %>% 
  filter(include == TRUE) %>% 
  select(cell50x50, species)

nrow(EBBA2_sf_country_code) # 579'487
nrow(EBBA2_clipped_sf) # 365'389

# save:
st_write(EBBA1_clipped_sf, file.path("data", "EBBA_analysis", "EBBA1_prelim_comparable_cells.shp")) # save before taxonomic harmonization because otherwise 5 EBBA2 cells are dropped that don't contain any record of comparable species
st_write(EBBA2_clipped_sf, file.path("data", "EBBA_analysis", "EBBA2_prelim_comparable_cells.shp"))


## taxonomic harmonization of EBBAs:--------------------------------------------

# based on EBBA2 methods chapter:
# "EBBA1 data were attributed, where possible, to the species recognised in EBBA2 (Table 6)"

EBBA1_taxunif_sf <- EBBA1_clipped_sf
EBBA2_taxunif_sf <- EBBA2_clipped_sf

# account for name change from EBBA1 to EBBA2:
EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species=="Acanthis hornemanni"),'species'] <- "Acanthis flammea"
EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species=="Oceanodroma castro"),'species'] <- "Hydrobates castro"
EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species=="Parus lugubris"),'species'] <- "Poecile lugubris"
EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species=="Regulus ignicapillus"),'species'] <- "Regulus ignicapilla"
EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species=="Serinus citrinella"),'species'] <- "Carduelis citrinella"

# remove all non-comparable species according to Table 6 in EBBA2 methods chapter:
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Branta hutchinsii"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Calonectris diomedea"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Calonectris borealis"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Calonectris borealis"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Hydrobates monteiroi"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Larus cachinnas"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Larus michahellis"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Larus michahellis"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Picus viridis"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Picus sharpei"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Picus viridis"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Lanius excubitor"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Lanius meridionalis"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Lanius excubitor"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Poecile hyrcanus"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Oenanthe xanthoprymna"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Oenanthe chrysopygia"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Oenanthe xanthoprymna"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Phylloscopus collybita"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Phylloscopus tristis"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Phylloscopus ibericus"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Phylloscopus collybita"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Phylloscopus nitidus"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Sylvia sarda"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Sylvia balearica"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Sylvia sarda"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Sylvia cantillans"),]
EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Sylvia subalpina"),]
EBBA1_taxunif_sf <- EBBA1_taxunif_sf[which(EBBA1_taxunif_sf$species!="Sylvia cantillans"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Regulus  madeirensis"),]

EBBA2_taxunif_sf <- EBBA2_taxunif_sf[which(EBBA2_taxunif_sf$species!="Carduelis corsicana"),]

length(unique(EBBA1_taxunif_sf$species)) # 445 species left
length(unique(EBBA2_taxunif_sf$species)) # 517 species left

# save:
st_write(EBBA1_taxunif_sf, file.path("data", "EBBA_analysis", "EBBA1_prelim_comparable_harmonized.shp"))
st_write(EBBA2_taxunif_sf, file.path("data", "EBBA_analysis", "EBBA2_prelim_comparable_harmonized.shp"))


# --------------------------------------------- #
#            EBBA change dataset:            ####
# --------------------------------------------- #

# requested for species selected with 2_3_EBBA_species_filtering_5_climatic_niche_analysis.R

# Metadata: 
# Change: A = Loss, B = Apparent Loss, C = Stable, D = Apparent Stable, E = Gain, F = Apparent Gain, G = Present in EBBA1 or EBBA2 but taxonomic assignment uncertain or wrong
# the word "Apparent" is preceded when the coverage in that square was considered insufficient for a proper comparison in any of the two atlases

EBBA_change <- read.csv2(file.path(datashare_EBCC, "EBBA_change", "ebba2_data_change_50km.csv"))
# (Sylvia ruppeli missing because "the species info did not pass the steps to have enough squares good enough in terms of comparability of efforts of coverage between the two atlases. Observed change for this species in the few squares in S Greece are thus considered not reliable at all (just effort differences)"

# include all cell-species combinations with change codes A, C or E
# sort to EBBA1 and EBBA2:
## A = Loss = species occupies cell in EBBA1 but not in EBBA2
## E = Gain = species occupies cell in EBBA2 but not in EBBA1
## C = Stable = species occupies cell in both EBBA1 and EBBA2

# EBBA comparable cells according to change data set:
EBBA_change_grid <- read_sf(file.path(datashare_EBCC, "EBBA_change", "ebba2_grid50x50_change_v1", "ebba2_grid50x50_change_v1.shp"))

EBBA_change_sf <- st_centroid(EBBA_change_grid) %>% 
  full_join(EBBA_change, by = c("cell50x50_" = "cell50x50_change")) %>% 
  rename(species = birdlife_scientific_name, cell50x50 = cell50x50_) %>% 
  select(-birdlife_code)

# EBBA1 change data:
EBBA1_ch_sf <- EBBA_change_sf %>% 
  filter(Change %in% c("A", "C"))

# EBBA2 change data:
EBBA2_ch_sf <- EBBA_change_sf %>% 
  filter(Change %in% c("E", "C"))
# not necessary to harmonize taxonomy or to remove non comparable cells

# save:
st_write(EBBA_change_sf, file.path("data", "EBBA_analysis", "EBBA_change.shp"), append = FALSE)
st_write(EBBA1_ch_sf, file.path("data", "EBBA_analysis", "EBBA1_change.shp"), append = FALSE)
st_write(EBBA2_ch_sf, file.path("data", "EBBA_analysis", "EBBA2_change.shp"), append = FALSE)