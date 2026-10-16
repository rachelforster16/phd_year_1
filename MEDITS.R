# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data")

# Load packages
library(dplyr)
library(tidyverse)
library(rfishbase)
library(tidyr)
library(ggplot2)
library(ggOceanMaps)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(broom)

# Load taxonomy
taxonomy_lookup <- read.csv("taxonomy_lookup.csv")


#=======================================================================
# MEDITS: Gulf of Lion and eastern Corsica (North-western Mediterranean)
#=======================================================================

# Load data
MEDITS_nw_catch <- read.csv("MEDITS_nw_catch.csv", sep = ";")
MEDITS_nw_operation <- read.csv("MEDITS_nw_operation.csv", sep = ";")
MEDITS_nw_catch_operation <- left_join(MEDITS_nw_catch, MEDITS_nw_operation, by = "haulID") %>%
  rename(accepted_name = scientificName)
MEDITS_nw_full <- left_join(MEDITS_nw_catch_operation, taxonomy_lookup, by = "accepted_name")

MEDITS_nw_full <- MEDITS_nw_full %>%
  ungroup() %>%
  mutate(
    start_lat_cell = floor(startLatDD),
    start_lon_cell = floor(startLongDD),
    end_lat_cell = floor(endLatDD),
    end_lon_cell = floor(endLongDD)
  ) %>%
  mutate(
    matching_lat = start_lat_cell == end_lat_cell,
    matching_lon = start_lon_cell == end_lon_cell
  )


# **CHECK** Each transect_id corresponds with only one year/grid cell combination
MEDITS_nw_full %>%
  distinct(
    haulID,
    year.x,
    start_lat_cell,
    start_lon_cell
  ) %>%
  count(haulID) %>%
  filter(n > 1)

# Rename columns
MEDITS_nw <- MEDITS_nw_full %>%
  mutate(num_cpue = (totalNumber/haulDurMin)) %>%
  mutate(survey = "MEDITS_nw") %>%
  rename(year = year.x, haul_id = haulID, lat_cell = start_lat_cell,
         lon_cell = start_lon_cell) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)

# Find species that have ever been recorded in a grid square
MEDITS_nw_species <- MEDITS_nw %>%
  filter(num_cpue > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class, superclass)

# **CHECK** Only one row per species per grid cell
MEDITS_nw_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
MEDITS_nw_hauls <- MEDITS_nw %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
MEDITS_nw_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
MEDITS_nw_complete <- MEDITS_nw_hauls %>%
  inner_join(
    MEDITS_nw_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    MEDITS_nw %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
MEDITS_nw_mean_num_cpue <- MEDITS_nw_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
MEDITS_nw_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(MEDITS_nw$superclass)

MEDITS_nw_elasmo_ts <- MEDITS_nw_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(MEDITS_nw_catch)
rm(MEDITS_nw_operation)
rm(MEDITS_nw_catch_operation)
rm(MEDITS_nw_full)
rm(MEDITS_nw_species)
rm(MEDITS_nw_hauls)
rm(MEDITS_nw_complete)

write_csv(MEDITS_nw_elasmo_ts, "MEDITS_nw_timeseries.csv")

#==================================
# MEDITS: Spanish Continental Shelf
#==================================

MEDITS_spain_occurrence <- read_tsv("MEDITS_spain_occurrence.txt")
MEDITS_spain_extendedmeasurement<- read_tsv("MEDITS_spain_extendedmeasurementorfact.txt")
MEDITS_spain_event <- read_tsv("MEDITS_spain_event.txt")

MEDITS_spain_measurement_occurrence <- left_join(MEDITS_spain_extendedmeasurement, MEDITS_spain_occurrence, by = "occurrenceID")
MEDITS_spain_all <- left_join(MEDITS_spain_measurement_occurrence, MEDITS_spain_event, by = "eventID")

# Filter data so there are no measurements of e.g. net size etc.
MEDITS_spain_full <- MEDITS_spain_all %>%
  select(collectionCode, id.x, occurrenceID, eventID, measurementType, measurementValue, scientificName, year.y, 
         decimalLatitude.y, decimalLongitude.y, class, order) %>%
  filter(measurementType %in% c("Count (in assayed sample) of biological entity specified elsewhere",
                                "Abundance of biological entity specified elsewhere per unit area of the bed",
                                "Wet weight biomass of biological entity specified elsewhere per unit area of the bed"))


# Rearrange the data so that each occurrence ID has a column for count, abundance, and weight
MEDITS_spain_full_wide <- MEDITS_spain_full %>%
  pivot_wider(
    names_from = measurementType,
    values_from = measurementValue
  ) %>%
  rename(
    count = "Count (in assayed sample) of biological entity specified elsewhere",
    abundance = "Abundance of biological entity specified elsewhere per unit area of the bed",
    weight = "Wet weight biomass of biological entity specified elsewhere per unit area of the bed"
  )

# Make the column names consistent with other databases
MEDITS_spain <- MEDITS_spain_full_wide %>%
  mutate(
    lat_cell = floor(decimalLatitude.y),
    lon_cell = floor(decimalLongitude.y),
    abundance = as.numeric(abundance)
  ) %>%
  rename(survey = collectionCode, haul_id = id.x, year = year.y, num_cpua = abundance, accepted_name = scientificName) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpua, accepted_name, order, class)
  

# **CHECK** Each haul_id corresponds with only one year/grid cell combination
MEDITS_spain %>%
  distinct(
    haul_id,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(haul_id) %>%
  filter(n > 1)


# Find species that have ever been recorded in a grid square
MEDITS_spain_species <- MEDITS_spain %>%
  filter(num_cpua > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class)

# **CHECK** Only one row per species per grid cell
MEDITS_spain_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
MEDITS_spain_hauls <- MEDITS_spain %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
MEDITS_spain_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
MEDITS_spain_complete <- MEDITS_spain_hauls %>%
  inner_join(
    MEDITS_spain_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    MEDITS_spain %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, num_cpua),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class"
    )
  ) %>%
  mutate(num_cpua = replace_na(num_cpua, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
MEDITS_spain_mean_num_cpua <- MEDITS_spain_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpua = mean(num_cpua, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
MEDITS_spain_mean_num_cpua %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpua = max(mean_num_cpua),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpua = sum(max_mean_num_cpua > 0),
    n_with_zero_only = sum(max_mean_num_cpua == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(MEDITS_spain$class)

MEDITS_spain_elasmo_ts <- MEDITS_spain_mean_num_cpua %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

# Clean house
rm(MEDITS_spain_all)
rm(MEDITS_spain_event)
rm(MEDITS_spain_extendedmeasurement)
rm(MEDITS_spain_full)
rm(MEDITS_spain_species)
rm(MEDITS_spain_hauls)
rm(MEDITS_spain_complete)
rm(MEDITS_spain_full_wide)
rm(MEDITS_spain_measurement_occurrence)
rm(MEDITS_spain_occurrence)

write_csv(MEDITS_spain_elasmo_ts, "MEDITS_spain_timeseries.csv")



#==================
# MEDITS: Black Sea
#==================

# Load RoME - contains data on species codes used in Black Sea survey data
# install.packages("RoME")
library(RoME)

# Load the TM (taxonomic list) data
data("TM_list")

# Set up a dataframe with the MEDITS codes and corresponding species names
Sp_Medits <- data.frame(
  MeditsCode = as.character(TM_list$MeditsCode),
  Scientific_name = as.character(
    TM_list[["Scientific.Name.................................................valid"]]
  )
)

# Load the black sea data
MEDITS_blacksea_speciescatch <- read.csv("MEDITS_blacksea_speciescatch.csv")

# Combine the genus and species columns into a MeditsCode column
MEDITS_blacksea_speciescatch <- MEDITS_blacksea_speciescatch %>%
  unite(MeditsCode, genus, species, sep = "")

# Check how many codes match the TM list for each reference list
MEDITS_blacksea_speciescatch %>%
  distinct(lref, MeditsCode) %>%
  mutate(MeditsCode = as.character(MeditsCode)) %>%
  left_join(Sp_Medits, by = "MeditsCode") %>%
  group_by(lref) %>%
  summarise(
    total_codes = n(),
    matched = sum(!is.na(Scientific_name)),
    unmatched = sum(is.na(Scientific_name)),
    .groups = "drop"
  )

# Find which species haven't had codes matched
unmatched <- MEDITS_blacksea_speciescatch %>%
  distinct(lref, MeditsCode) %>%
  mutate(MeditsCode = as.character(MeditsCode)) %>%
  anti_join(tm_lookup, by = "MeditsCode") %>%
  arrange(lref, MeditsCode) %>% 
  as_tibble()

head(unmatched, 8)

# Manually add names for those where the code and species name can be found online (e.g. the MEDITS handbook)
manual_names <- tibble(
  MeditsCode = c("ARGRACU", "ASCDCEA", "BUCCSPP", "FMBONEL", "GASTRDA", "ILLESPP", "PTEDGRI", "SOLEAEG"),
  Scientific_name = c("Argyropelecus aculeatus", "Ascidiacea spp.", "Buccinum spp", "Bonellia viridis", "Gastropoda spp.", 
                      "Illex spp.", "Pteroeides griseum", "Solea aegyptiaca")
)

Sp_Medits <- bind_rows(Sp_Medits, manual_names) %>%
  distinct(MeditsCode, .keep_all = TRUE) %>%
  arrange(MeditsCode)

write_csv(Sp_Medits, "Sp_Medits.csv")


# Add the species name column to the black sea data
MEDITS_blacksea_speciescatch_names <- left_join(MEDITS_blacksea_speciescatch, Sp_Medits, by = "MeditsCode")

# Load haul data
MEDITS_blacksea_haul <- read.csv("MEDITS_blacksea_haul.csv")

# Changing latitude and longitude columns so they are in degrees (not degrees minutes seconds)
dms_to_decimal <- function(x) {
  x <- as.numeric(x)
  degrees <- floor(x/100)
  minutes <- x - degrees * 100
  degrees + minutes / 60
}

MEDITS_blacksea_haul <- MEDITS_blacksea_haul %>%
  mutate(
    latitude = dms_to_decimal(hauling_latitude),
    longitude = dms_to_decimal(hauling_longitude)
  )

# Create unique code for each survey with country, area, vessel, year, month, haul number?
MEDITS_blacksea_speciescatch_names <- MEDITS_blacksea_speciescatch_names %>%
  unite(haul_id, country, area, vessel, year, month, haul_number, sep = ",", remove = FALSE)

MEDITS_blacksea_haul <- MEDITS_blacksea_haul %>%
  unite(haul_id, country, area, vessel, year, month, haul_number, sep = ",", remove = FALSE)

# Join datasets by haul number
MEDITS_blacksea_full <- left_join(MEDITS_blacksea_speciescatch_names, MEDITS_blacksea_haul, by = "haul_id") %>%
  rename(accepted_name = Scientific_name) %>%
  left_join(taxonomy_lookup, by = "accepted_name")

# Make columns consistent with other databases
MEDITS_blacksea <- MEDITS_blacksea_full %>%
  mutate(
    survey = "MEDITS_blacksea",
    lat_cell = floor(latitude),
    lon_cell = floor(longitude),
    num_cpue = (nbtot/hauling_duration)
  ) %>%
  rename(year = year.x) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name, order, class, superclass)


# **CHECK** Each haul_id corresponds with only one year/grid cell combination
MEDITS_blacksea %>%
  distinct(
    haul_id,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(haul_id) %>%
  filter(n > 1)


# Find species that have ever been recorded in a grid square
MEDITS_blacksea_species <- MEDITS_blacksea %>%
  filter(num_cpue > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class, superclass)

# **CHECK** Only one row per species per grid cell
MEDITS_blacksea_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
MEDITS_blacksea_hauls <- MEDITS_blacksea %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
MEDITS_blacksea_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
MEDITS_blacksea_complete <- MEDITS_blacksea_hauls %>%
  inner_join(
    MEDITS_blacksea_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    MEDITS_blacksea %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
MEDITS_blacksea_mean_num_cpue <- MEDITS_blacksea_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
MEDITS_blacksea_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(MEDITS_blacksea$superclass)

MEDITS_blacksea_elasmo_ts <- MEDITS_blacksea_mean_num_cpue %>%
  filter(class %in% c("Chondrichthyes"))

# Clean house
rm(MEDITS_blacksea_haul)
rm(MEDITS_blacksea_speciescatch)
rm(MEDITS_blacksea_speciescatch_names)
rm(MEDITS_blacksea_full)
rm(MEDITS_blacksea_species)
rm(MEDITS_blacksea_hauls)
rm(MEDITS_blacksea_complete)
rm(manual_names)
rm(unmatched)

write_csv(MEDITS_blacksea_elasmo_ts, "MEDITS_blacksea_timeseries.csv")




