# prepare BBS data based on: https://github.com/UP-macroecology/bbs_dataprep/blob/main/BBS_info.md

library(dplyr)
library(tidyr)
library(sf)
# use the package devtools to install bbsAssistant directly from github
#devtools::install_github("trashbirdecology/bbsAssistant")
library(bbsAssistant)


# ------------------------------ #
#            Set-up:          ####
# ------------------------------ #

# project data:
#data_dir <- file.path("/import", "ecoc9z", "data-zurell", "schifferle", "EBBA_niche_range_shifts")
data_dir <- file.path("Data")

# folder for BBS data:
# create folder if it doesn't exist yet:
if(!dir.exists(file.path(data_dir, "BBS"))){dir.create(file.path(data_dir, "BBS"), recursive = TRUE)}

# function: 
# modified version of bbsAssistant::import_bbs_data (replaced "50-StopData.zip" with "States.zip")
# only needed for data before 1997

import_bbs_data_states <- function(bbs_dir, sb_id) {
  
  ObsN <- RTENO <- Date <- TotalSpp  <- NULL 
  
  # where to save the unzipped files
  tempdir = tempdir()
  
  # create a vector of desired file locations.
  zipF <-
    list.files(path = paste0(bbs_dir),
               pattern = "States.zip",
               full.names = TRUE)
  utils::unzip(zipF, exdir = tempdir) # unzip the top directory
  fns.states <-
    paste0(tempdir,
           "/",
           unzip(
             zipfile = zipF,
             list = TRUE,
             exdir = tempdir()
           )$Name)
  
  fns.states <-
    fns.states[stringr::str_detect(tolower(fns.states), pattern = ".zip")] # to remove the dir that isnt a .zip
  
  fns.routes <- list.files(path = paste0(bbs_dir),
                           pattern = "routes.zip",
                           full.names = TRUE)
  fns.vehicle <- list.files(path = paste0(bbs_dir),
                            pattern = "ehicle",
                            full.names = TRUE)
  fns.weather <- list.files(path = paste0(bbs_dir),
                            pattern = "eather.zip",
                            full.names = TRUE)
  
  # define potential columns and desired types to ensure consistency across data files
  col_types <- readr::cols(
    AOU = readr::col_integer(),
    CountryNum = readr::col_integer(),
    Route = readr::col_character(),
    RouteDataID = readr::col_integer(),
    RPID = readr::col_integer(),
    StateNum = readr::col_integer(),
    Year = readr::col_integer()
  )
  
  # Get observations and routes ---
  
  observations <- list()
  for (i in seq_along(fns.states)) {
    f <- fns.states[i]
    observations[[i]]  <- readr::read_csv(f, col_types = col_types)
  }
  observations <- dplyr::bind_rows(observations)
  
  # Get dataset citation(s) ---
  #citation <- sbtools::item_get_fields(sb_id, "citation")
  
  # Get species list ---
  species_list <- import_species_list(bbs_dir)
  
  # Get route metadata ---
  routes <- suppressWarnings(readr::read_csv(fns.routes, col_types = col_types))
  weather <- suppressWarnings(readr::read_csv(fns.weather, col_types = col_types))
  vehicle_data <- suppressWarnings(readr::read_csv(unzip(zipfile = fns.vehicle, exdir = tempdir), col_types = col_types)) 
  
  observers <- weather %>%
    make.dates() %>% 
    make.rteno() %>%
    dplyr::select(ObsN, RTENO, Date, TotalSpp) %>%
    ##create binary for if observer's first year on the BBS and on the route
    dplyr::group_by(ObsN) %>% #observation identifier (number)
    dplyr::mutate(ObsFirstYearOnBBS = ifelse(Date==min(Date), 1, 0)) %>%
    dplyr::group_by(ObsN, RTENO) %>%
    dplyr::mutate(ObsFirstYearOnRTENO = ifelse(Date==min(Date), 1, 0)) %>%
    dplyr::ungroup() # to be safe
  
  # Create a list of data and information to export or return
  list.elements <-
    list("observations",
         "routes",
         "observers",
         "weather",
         "species_list",
         #"citation",
         "vehicle_data"
    )
  bbs <- lapply(
    list.elements,
    FUN = function(x) {
      eval(parse(text = paste(x))) %>%
        make.rteno()
    }
  )
  names(bbs) <- list.elements
  
  # END FUNCTION ---
  return(bbs)
}


# --------------------------------------------------- #
#      Create presence / absence data from BBS:    ####
# --------------------------------------------------- #

# (In Sofaer, Jarnevich, and Flather (2018), presences correspond to any time a species was 
# observed on a route within a time period. Absences are assigned any time a species was not 
# observed on routes that were sampled in each year of a time period.)

# import BBS counts (aggregated to 5 sections - States.zip):
bbs_agg <- import_bbs_data_states(bbs_dir = file.path(data_dir, "BBS"), 
                                  sb_id = "5ea04e9a82cefae35a129d65") # sb_id = sb_item of the current BBS dataset (can be looked up by typing sb_items)
# write species list (to later merge ids to names):
write.csv(bbs_agg$species_list, file = file.path(data_dir, "BBS_species_list.csv"))

# merge BBS datasets:
bbs_merged <- bbs_agg$observations %>% 
  left_join(bbs_agg$routes) %>% 
  left_join(bbs_agg$weather) %>% 
  # only conterminous United States routes (we have no detailed spatial info for other BBS countries)
  filter(CountryNum == 840) %>% 
  filter(StateNum != 3) # exclude Alaska

rm(bbs_agg)

# time periods to compare:
# 3-year periods as in Sofaer et al. 2018 (Sofaer et al. use 1977-1979 and 2012-2014)
# historic: 2 versions (see following for loop)
recent <- 2016:2018 # 2018: end of available Chelsa data

# route sections from which to use data (all five sections are about 40 km)
sections <- paste0("Count", seq(10, 50, by = 10))

# loop over 2 versions for historic period:
for(i in c(1,2)){ 
  
  if(i == 1){
    historic <- 1980:1982 # 1980: start of available Chelsa data
    } else {
      historic <- 1996:1998 # version that yields similar time gap as between EBBA1 and EBBA2
    }
  
  ## validation route IDs: ----
  
  # routes that were sampled in all three years (to get true absences) of both periods (comparison of historic and recent time period based on same routes)
  
  valid_routes <- bbs_merged %>%
    
    # filter the data by recorded quality index:
    filter(RunType == 1) %>% # RunType 1 means QualityCurrentID == 1 & RouteTypeDetailID == 1 & RPID == 101
    
    # remove water based routes:
    filter(RouteTypeID == 1) %>% 
    
    select(RTENO, Year) %>% 
    filter(Year %in% c(historic, recent)) %>% 
    distinct %>% 
    
    group_by(RTENO) %>% 
    count(RTENO) %>% 
    filter(n == length(c(historic, recent))) %>% # routes sampled in all years of both periods
    pull(RTENO)

  
  # filter data, add presence information: ----
  
  bbs_cleaned <- bbs_merged %>% # 5,726,806
    
    # use only validation routes (= sampled in all years of both time periods):
    filter(RTENO %in% valid_routes) %>% 
    
    # filter the data by recorded quality index:
    filter(RunType == 1) %>% # RunType 1 means QualityCurrentID == 1 & RouteTypeDetailID == 1 & RPID == 101
    
    # remove water based routes:
    filter(RouteTypeID == 1) %>% 
    
    # select only data from two time periods we want to compare (historic and recent):
    filter(Year %in% c(historic, recent)) %>% 
    # assign a flag for each time period
    mutate(period = if_else(Year %in% historic, 1, 2)) %>% # 2,063,001
    
    # add presence information:
    
    # a) assign presence where a species was observed at least once on a route in a year
    mutate(present_year = if_else(rowSums(select(., all_of(sections))) > 0, 1, 0)) %>% 
    
    # b) assign presence where a species was present at least once in a time period on a route (= final presences)
    group_by(AOU, period, RTENO) %>% # species, time period, route ID
    mutate(present_period = if_else(sum(present_year) > 0, 1, 0)) %>% 
    ungroup
  
  
  # split data into historic and recent period, add true absence information, convert to long format: ----
  
  # (true absence = species was not observed in each year of the time period)
  
  # historic time period:
  model_df_hist <- bbs_cleaned %>% 
    select(RTENO, Latitude, Longitude, AOU, period, present_period) %>% 
    distinct() %>% # one row for each route and time period
    pivot_wider(names_from = AOU, # species as columns, to get each route-species combination
                values_from = present_period) %>% # 1 = present, NA = true absence (no presence recorded, although sampled in each year of the time period)
    mutate(across(.cols = c(5:ncol(.)), .fns = ~ ifelse(is.na(.x), 0, .)))  %>% # change NA to 0 for true absence
    filter(period == 1) %>% # historic time period (filter at the end to have column for every species)
    select(-period) %>% 
    pivot_longer(cols = 4:ncol(.), names_to = "AOU", values_to = "pres") %>%  # convert back to long format with presence/absence column
    mutate(AOU = as.numeric(AOU))
  
  # recent time period:
  model_df_rec <- bbs_cleaned %>% 
    select(RTENO, Latitude, Longitude, AOU, period, present_period) %>% 
    distinct() %>% # one row for each route and time period
    pivot_wider(names_from = AOU, # species as columns, to get each route-species combination
                values_from = present_period) %>% #  1 = present, NA = true absence (no presence recorded, although sampled in each year of the time period)
    mutate(across(.cols = c(5:ncol(.)), .fns = ~ ifelse(is.na(.x), 0, .)))  %>% # change NA to 0 for true absence
    filter(period == 2) %>% # recent time period (filter at the end to have column for every species)
    select(-period) %>% 
    pivot_longer(cols = 4:ncol(.), names_to = "AOU", values_to = "pres") %>% # convert back to long format with presence/absence column
    mutate(AOU = as.numeric(AOU))

  # add species names:
  species_names <- read.csv(file.path(data_dir, "BBS_species_list.csv")) %>% 
    select(c(AOU, Scientific_Name))
  
  model_df_hist <- model_df_hist %>% 
    left_join(species_names) %>% 
    rename(species = Scientific_Name)
    
  model_df_rec <- model_df_rec %>% 
    left_join(species_names) %>% 
    rename(species = Scientific_Name)
  
  # write to shapefiles:
  BBS_hist_sf <- st_as_sf(model_df_hist, coords = c("Longitude", "Latitude"), crs = 4269)
  st_write(BBS_hist_sf, file.path(data_dir, paste0("BBS_historic", i, ".shp")), append = FALSE)
  
  BBS_rec_sf <- st_as_sf(model_df_rec, coords = c("Longitude", "Latitude"), crs = 4269)
  st_write(BBS_rec_sf, file.path(data_dir, paste0("BBS_recent", i, ".shp")), append = FALSE) # recent contains different validation routes depending on historic years considered
}