#### Generating species list and preliminary species inputs for requesting aligned data from raw EBBA1 and EBBA2

# datashare:
datashare_EBCC <- file.path("//ibb-fs01.ibb.uni-potsdam.de", "daten$", "AG26", "Arbeit", "datashare", "data", "biodat", "distribution", "EBCC")
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
# -> looks like:
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
length(unique(EBBA2_sf_points$cell50x50)) # 5303

EBBA2_sf_pts_dt <- full_join(EBBA2_sf_points, EBBA2_dt, by = "cell50x50") %>% 
  rename(species = birdlife_scientific_name) %>% 
  select(-birdlife_code)

nrow(EBBA2_sf_pts_dt) # 579487 (more because there are 224 cells for which no presence is recorded)
#st_write(EBBA2_sf_pts_dt, file.path(transfer, "input_data", "EBCC", "EBBA2", "EBBA2_EPSG3035_KS.shp"))
#st_write(EBBA2_sf_pts_dt, "EBBA2_EPSG3035_KS.shp")

# explorations regarding EBBA1 & convert_UTMcoord.R: ----------

library(sp)
library(rgdal)
library(tibble)

# originally, EBBA1 data were in UTM coordinates, script convert_UTMcoord.R contains examples how
# to convert lon-lat to UTM coordinates, but was not extensively tested

# Helper functions

find_UTM_zone <- function(longitude, latitude) {
  
  # Special zones for Svalbard and Norway
  if (latitude >= 72.0 && latitude < 84.0 ) 
    if (longitude >= 0.0  && longitude <  9.0) 
      return(31);
  if (longitude >= 9.0  && longitude < 21.0)
    return(33)
  if (longitude >= 21.0 && longitude < 33.0)
    return(35)
  if (longitude >= 33.0 && longitude < 42.0) 
    return(37)
  
  (floor((longitude + 180) / 6) %% 60) + 1 # UTM: 60 zones, each 6° of longitude width, zone 1 = 180°- 174°W
}

find_UTM_hemisphere <- function(latitude) {
  
  ifelse(latitude > 0, "north", "south")
}

# returns a DF containing the UTM values, the zone and the hemisphere
longlat_to_UTM <- function(long, lat, units = 'm') {
  
  df <- data.frame(
    id = seq_along(long), 
    x = long, 
    y = lat
  )
  sp::coordinates(df) <- c("x", "y")
  
  hemisphere <- find_UTM_hemisphere(lat)
  zone <- find_UTM_zone(long, lat) # not vectorized
  
  sp::proj4string(df) <- sp::CRS("+init=epsg:4326") 
  CRSstring <- paste0(
    "+proj=utm +zone=", zone,
    " +ellps=WGS84",
    " +", hemisphere,
    " +units=", units)
  
  if (dplyr::n_distinct(CRSstring) > 1L) 
    stop("multiple zone/hemisphere detected")
  
  res <- sp::spTransform(df, sp::CRS(CRSstring[1L])) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      zone = zone,
      hemisphere = hemisphere
    )
  
  res
}

# conversion:
EBBA1 <- EBBA1[!is.na(EBBA1$decimalLongitude),] # none

EBBA1_coords <- longlat_to_UTM(EBBA1$decimalLongitude, EBBA1$decimalLatitude) # only returns zone 35 (but see: https://de.wikipedia.org/wiki/UTM-Koordinatensystem#/media/Datei:Utm-zones.jpg)
# add UTM coordinates (WGS84, zone 35 north) to EBBA1:
EBBA1 <- cbind(EBBA1, EBBA1_coords[,c("x", "y")])

write.csv(EBBA1, file = file.path(transfer, "input_data", "EBCC", "EBBA1", "EBBA1_conv_coords_KS.csv"))

EBBA1_coords %>% 
  distinct(zone) %>% 
  pull(zone) # all zone 35?

EBBA1_coords %>% 
  distinct(x, y) # 3912 unique raster locations

EBBA1_coords_zone <- find_UTM_zone(EBBA1$decimalLongitude[1:2], EBBA1$decimalLatitude[1:2])
# not vectorized?!

find_UTM_zone2 <- function(longitude, latitude) {
  
  zones <- vector(mode = "numeric", length = length(longitude))
  
  for(i in 1:length(longitude)){
    
    # Special zones for Svalbard and Norway
    if (latitude[i] >= 72.0 && latitude[i] < 84.0) 
      if (longitude[i] >= 0.0  && longitude[i] <  9.0) 
        zones[i] <- 31

    if (longitude[i] >= 9.0  && longitude[i] < 21.0)
      zones[i] <- 33

    if (longitude[i] >= 21.0 && longitude[i] < 33.0)
      zones[i] <- 35

    if (longitude[i] >= 33.0 && longitude[i] < 42.0) 
      zones[i] <- 37

    zones[i] <- (floor((longitude[i] + 180) / 6) %% 60) + 1 # UTM: 60 zones, each 6° of longitude width, zone 1 = 180°- 174°W
  }
return(zones)
}

longlat_to_UTM2 <- function(long, lat, units = 'm') { # doesn't work, uses only first zone when transforming
  
  df <- data.frame(
    id = seq_along(long), 
    x = long, 
    y = lat
  )
  sp::coordinates(df) <- c("x", "y")
  
  hemisphere <- find_UTM_hemisphere(lat)
  zone <- find_UTM_zone2(long, lat)
  
  sp::proj4string(df) <- sp::CRS("+init=epsg:4326") 
  CRSstring <- paste0(
    "+proj=utm +zone=", zone,
    " +ellps=WGS84",
    " +", hemisphere,
    " +units=", units)
  
  # if (dplyr::n_distinct(CRSstring) > 1L) 
  #   stop("multiple zone/hemisphere detected")
  
  # loop hier:
  #length(unique(CRSstring)) # 18 zones, use 18 transformations, one for each zone
  
  res <- sp::spTransform(df, sp::CRS(CRSstring)) %>%  # not vectorized
    tibble::as_tibble() %>%
    dplyr::mutate(
      zone = zone,
      hemisphere = hemisphere
    )
  
  res
}

EBBA1_coords2 <- longlat_to_UTM2(EBBA1$decimalLongitude, EBBA1$decimalLatitude)

EBBA1_with_coords2 <- cbind(EBBA1[,c("occurrenceID", "species", "taxonRank", 
                                     "scientificName", "decimalLatitude", "decimalLongitude", "issue")], 
                            EBBA1_coords2[,c("x", "y")])
write.csv(EBBA1_with_coords2, file = file.path(transfer, "input_data", "EBCC", "EBBA1", "EBBA1_conv_coords_KS2.csv"))

# mehrere UTM Zonen verwendet, in QGIS kann man nur eine auswählen?
# spTransform benutzt nur ersten Wert 
# alle Koordinaten werden zu erster Zone im df transformiert, 
# Koordinaten beziehen sich auf den entsprechenden Mittelmeridian
# nehmen daher eigentlich unzulässige (?) UTM-Werte an (außerhalb von 100'000 bis 899'999 (wikipedia), aber QGIS kommt damit klar)
# kann CRS einer einzige Shapefile überhaupt mehrere UTM-Zonen beinhalten?


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

# cell-IDs Canary Islands: 
Canaries <- c("28RDU2", "27RYL3", "27RYM3", "28RBR1", "28RCS2", "28RCR1", "28RBS4", "28RBS1", 
              "28RDR3", "28RDR1", "28RCS4", "28RCS3", "28RES3", "28RES2", "28RDS4", "28RDS2", 
              "28RFT1", "28RFS2", "28RFS1", "28RES4", "28RFT4", "28RFT2")

# cell-IDs Russia between PL and LT:
Russia_Europe <- c("34UDG4", "34UDF1", "34UDF3", "34UEF1", "34UEF3", "34UEF4")

# cell-IDs European Turkey: 
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
  filter(include == TRUE) %>% 
  select(-include)

nrow(EBBA1_sf_EPSG3035_cellID) # 355'488
nrow(EBBA1_clipped_sf) # 330'486
#st_write(EBBA1_clipped_sf, "EBBA1_clipped_KS.shp")

# clip EBBA2 data to cells to keep:
EBBA2_clipped_sf <- EBBA2_sf_pts_dt_country_code %>%
  mutate(include = case_when(cell50x50 %in% Canaries ~ FALSE,
                             cell50x50 %in% Turkey_Europe ~ TRUE,
                             cell50x50 %in% Russia_Europe ~ TRUE,
                             countryCode %in% low_survey_effort ~ FALSE,
                             is.na(countryCode) ~ FALSE, # cell not included in EBBA1 -> drop
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

# hier weiter: xx
# based on EBBA methods chapter - removing all species from Table 6 unless necessary taxon changes clear!
EBBA1[which(EBBA1$species=="Acanthis hornemanni"),'species'] <- "Acanthis flammea"

EBBA1[which(EBBA1$species=="Oceanodroma castro"),'species'] <- "Hydrobates castro"

EBBA1[which(EBBA1$species=="Parus lugubris"),'species'] <- "Poecile lugubris"

EBBA1[which(EBBA1$species=="Regulus ignicapillus"),'species'] <- "Regulus ignicapilla"

EBBA1[which(EBBA1$species=="Serinus citrinella"),'species'] <- "Carduelis citrinella"

# remove non-comparable species:
EBBA2 <- EBBA2[which(EBBA2$species!="Branta hutchinsii"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Calonectris diomedea"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Calonectris borealis"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Calonectris borealis"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Hydrobates monteiroi"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Larus cachinnas"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Larus michahellis"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Larus michahellis"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Picus viridis"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Picus sharpei"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Picus viridis"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Lanius excubitor"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Lanius meridionalis"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Lanius excubitor"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Poecile hyrcanus"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Oenanthe xanthoprymna"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Oenanthe chrysopygia"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Oenanthe xanthoprymna"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Oenanthe xanthoprymna"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Oenanthe chrysopygia"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Oenanthe xanthoprymna"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Phylloscopus collybita"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Phylloscopus tristis"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Phylloscopus ibericus"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Phylloscopus collybita"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Phylloscopus bonelli"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Phylloscopus orientalis"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Phylloscopus bonelli"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Phylloscopus nitidus"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Sylvia hortensis"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Sylvia crassirostris"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Sylvia hortensis"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Sylvia sarda"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Sylvia balearica"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Sylvia sarda"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Sylvia cantillans"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Sylvia subalpina"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Sylvia cantillans"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Iduna pallida"),]
EBBA2 <- EBBA2[which(EBBA2$species!="Iduna opaca"),]

EBBA1 <- EBBA1[which(EBBA1$species!="Hippolais pallida"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Regulus  madeirensis"),]

EBBA2 <- EBBA2[which(EBBA2$species!="Carduelis corsicana"),]


# Species filtering:

## 1. exclude pelagic specialists (according to Wilman et al. 2014):

SeaBirds <- read.csv2(file.path(transfer, "input_data", "BirdFuncDat.txt"), header=TRUE, sep="\t") #BirdFuncDat from the EltonTraits database: https://figshare.com/articles/dataset/Data_Paper_Data_Paper/3559887?backTo=/collections/EltonTraits_1_0_Species-level_foraging_attributes_of_the_world_s_birds_and_mammals/3306933
SeaBirds <- subset(SeaBirds, PelagicSpecialist==1)
SeaBirds <- SeaBirds$Scientific

EBBA1 <- EBBA1[which(!EBBA1$species %in% SeaBirds),]
EBBA2 <- EBBA2[which(!EBBA2$species %in% SeaBirds),]


## 2. exclude rare species with n<20 occurrences in any of the two atlas periods:
EBBA1_1 <- EBBA1 %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20)

EBBA2_1 <- EBBA2 %>% 
  group_by(species) %>% 
  mutate(n_occurrences = n()) %>%
  filter(n_occurrences >= 20)


## 3. keep only those species that occur in both atlas periods:
EBBA1_2 <- EBBA1_1 %>% 
  filter(species %in% EBBA2_1$species)

EBBA2_2 <- EBBA2_1 %>% 
  filter(species %in% EBBA1_1$species)


## 4. exclude very common species with >90% prevalence in both atlas periods:
# = occur in >90% EBBA cells:

# number of EBBA cells:
EBBA2_Grid <- read_sf(file.path(transfer, "input_data", "EBCC", "EBBA2", "ebba2_grid50x50_v1", "ebba2_grid50x50_v1.shp"))
# each cell = one row
KeepCountries <- read_sf(file.path(transfer, "Comparable_region_mask.shp"))
RetainedEBBA_Grid <- terra::intersect(EBBA2_Grid,KeepCountries)
EBBAcells <- nrow(RetainedEBBA_Grid)


##Exporting species distributions for examination and identification of remaining issues##
for(i in unique(EBBA1$species)){
  
  #Isolate single species#
  Sp_EBBA2 <- subset(EBBA2,EBBA2$brdlf__==i)
  Sp_EBBA1 <- subset(EBBA1,EBBA1$species==i)
  
  #Plot species distributions from EBBA1 and 2, exporting so I can go through and identify species with unrealistic distributions and/or seabirds with coastal distributions that are not suited to these analyses#
  tiff(filename = paste("./Data/EBBA_distribution_maps/",i,"_EBBA_map.tiff",sep=""),width=1000,height=1000,res=125)
  
  plot(Sp_EBBA1,col=alpha('red',1),pch=1,cex=0.75,main=paste(i,"red=EBBA1", "blue=EBBA2",sep=" "),xlim=c(xmin(EBBA1), xmax(EBBA1)), ylim=c(ymin(EBBA1), ymax(EBBA1)))
  
  plot(Sp_EBBA2,col=alpha('blue',1),pch=16,cex=0.5,add=T)
  
  plot(KeepCountries,add=T)
  
  dev.off()
}


##Filtering out the most common birds##

EBBA2_common<-as.data.frame(table(EBBA2$brdlf__))

EBBA2_common<-as.character(subset(EBBA2_common,EBBA2_common$Freq>(0.9*EBBAcells))$Var1)#common if occur in >90% EBBA cells, removes 20 species

EBBA2<- as.data.frame(EBBA2) %>% group_by(brdlf__) %>% filter(n()<(0.9*EBBAcells)) %>% ungroup()#Use EBBA2 to filter as less confident presences in EBBA1 are accurate

Sp_List<-unique(EBBA2$brdlf__) #Now 326 species

##IUCN range overlap analysis##
##Filter Out species with <5% range in EBBA comparable cells. Range from IUCN Red List global polygons##

Sp_List <- sub(" ", "_", Sp_List)#replace space with underscore to match IUCN file names

#List IUCN species polygons# 

#!!Make sure polygons are unzipped in folder before attempting!!#

target_path <- file.path(".Data//IUCN_range_polys/All_shapefiles")

IUCN_species <- list.files(target_path, full.names = F, pattern = ".shp")

IUCN_species <- substr(IUCN_species,1,nchar(IUCN_species)-13)

#Match IUCN names to EBBA species and identify unmatched EBBA species to harmonize#

Matched_sp<-IUCN_species[which(IUCN_species %in% Sp_List)]

unmatched_sp<-Sp_List[which(!Sp_List %in% IUCN_species)]##Examine these names and find suitable replacement for harmonization

#Harmonization notes (changed IUCN file names to match those used in EBBA)
# - Carduelis_chloris -> Chloris_chloris
# - Carduelis_flammea -> Acanthis_flammea
# - Carduelis_spinus -> Spinus_spinus
# - Carduelis_cannabina -> Linaria_cannabina
# - Milaria_calandra -> Emberiza_calandra
# - Parus_ater -> Periparus_ater
# - Hirundo_rupestris -> Ptyonoprogne_rupestris
# - Hirundo_daurica -> Cecropis_daurica
# - Parus_caeruleus -> Cyanistes_caeruleus
# - Parus_cristatus -> Lophophanes_cristatus
# - parus_palustris -> Poecile_palustris
# - Carduelis_flavirostris -> Linaria_flavirostris
# - Parus_montanus -> Poecile_montanus
# - Parus_lugubris -> Poecile_lugubris
# - Parus_cinctus -> Poecile_cinctus
# - Sturnus_roseus -> Pastor_roseus
# - Hippolais_caligata -> Iduna_caligata
#!! Passer_italiae not mapped by IUCN!! But is European, so should be retained in species list#

#List IUCN files with full names#

IUCN_polys <- list.files(target_path, full.names = T, pattern = ".shp")

#Restrict IUCN polygons to just species included in the species list

IUCN_polys <- IUCN_polys[grep(paste(Sp_List, collapse="|"),IUCN_polys)]

##Read in EBBA polygons and restrict to retained region
EBBA2_poly <- readOGR("./Data/input_data/EBCC/EBBA2/ebba2_grid50x50_v1/ebba2_grid50x50_v1.shp") %>%
  spTransform(CRS('+proj=laea +lat_0=10 +lon_0=-81 +ellps=WGS84 +units=m +no_defs')) %>%
  gBuffer(byid=TRUE, width=0)

KeepCountries <- readOGR(".Data/","Comparable_region_mask") %>%
  spTransform(CRS('+proj=laea +lat_0=10 +lon_0=-81 +ellps=WGS84 +units=m +no_defs')) %>%
  gBuffer(byid=TRUE, width=0)#Put into global CRS and buffer to avoid gIntersect error

EBBA2_poly <- gIntersection(EBBA2_poly,KeepCountries) %>%
  gBuffer(byid=TRUE, width=0)

EBBA_range_cover<-list()

for(i in IUCN_polys[1:306]){
  
  glob_range<-readOGR(i) %>% 
    spTransform(CRS('+proj=laea +lat_0=10 +lon_0=-81 +ellps=WGS84 +units=m +no_defs')) %>%
    
    gBuffer(byid=TRUE, width=0)
  
  EBBA_range<-gIntersection(glob_range,EBBA2_poly)
  
  if (is.null(EBBA_range)){
    EBBA_range_cover[i]<-0
  } else {
    EBBA_range_cover[i]<-round(area(EBBA_range)/sum(area(glob_range))*100,2)
  }}

#Coerce list into table with species names, then subset by %

RangeCoverdf <- melt(EBBA_range_cover)
RangeCoverdf <- RangeCoverdf[,c(2,1)]
names(RangeCoverdf) <- c("species","range_cover_percent")
RangeCoverdf$species <- substr(RangeCoverdf$species,77,nchar(RangeCoverdf$species)-13)

#Restrict to finalized species list
write.csv(RangeCoverdf,"./Data/IUCN_poly_covered.csv")

hist(RangeCoverdf$range_cover_percent, main="Range within EBBA Region", xlab="Percentage range covered")

Sp_list <- subset(RangeCoverdf,range_cover_percent>5)$species #Likely threshold not defensible, method needs reconsidering

write.csv(Sp_List,"./Data/Change_Data_request_species_list.csv")

writeOGR(EBBA1,"./Data/input_data/EBCC/","Filtered_EBBA1",driver="ESRI Shapefile",overwrite=T)

writeOGR(EBBA2,"./Data/input_data/EBCC/","Filtered_EBBA2",driver="ESRI Shapefile",overwrite=T)