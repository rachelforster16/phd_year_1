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



#==================================
# MEDITS: Spanish Continental Shelf
#==================================

MEDITS_spain_occurrence <- read_tsv("MEDITS_spain_occurrence.txt") %>%
  select(id, scientificName, class, order)
MEDITS_spain_extendedmeasurement_all <- read_tsv("MEDITS_spain_extendedmeasurementorfact.txt") %>%
  select(id, measurementType, measurementValue)
MEDITS_spain_event <- read_tsv("MEDITS_spain_event.txt") %>%
  select(id, year, decimalLatitude, decimalLongitude)

MEDITS_spain_extendedmeasurement_all %>%
  count(measurementType, sort = TRUE)

MEDITS_spain_extendedmeasurement <- MEDITS_spain_extendedmeasurement_all %>%
  filter(measurementType == "Count (in assayed sample) of biological entity specified elsewhere")

MEDITS_spain_occurrence_measurement <- left_join(MEDITS_spain_occurrence, MEDITS_spain_extendedmeasurement, by = "id", relationship = "many-to-many")
MEDITS_spain_full <- left_join(MEDITS_spain_occurrence_measurement, MEDITS_spain_event, by = "id", relationship = "many-to-many") %>%
  rename(accepted_name = scientificName)

MEDITS_spain_full <- MEDITS_spain_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(decimalLatitude),
    lon_cell = floor(decimalLongitude)
  )

# **CHECK** Each id corresponds with only one year/grid cell combination
MEDITS_spain_full %>%
  distinct(
    id,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(id) %>%
  filter(n > 1)

# Rename columns
MEDITS_spain <- MEDITS_spain_full %>%
  mutate(survey = "MEDITS_spain",
    num_cpue = as.numeric(measurementValue)) %>%
  rename(haul_id = id) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue,
    accepted_name, order, class)

# Find species that have ever been recorded in a grid square
MEDITS_spain_species <- MEDITS_spain %>%
  filter(num_cpue > 0) %>% 
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
             order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
MEDITS_spain_mean_num_cpue <- MEDITS_spain_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
MEDITS_spain_mean_num_cpue %>%
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
unique(MEDITS_spain$class)

MEDITS_spain_elasmo_ts <- MEDITS_spain_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

# Clean house
rm(MEDITS_spain_event)
rm(MEDITS_spain_extendedmeasurement)
rm(MEDITS_spain_extendedmeasurement_all)
rm(MEDITS_spain_occurrence)
rm(MEDITS_spain_occurrence_measurement)
rm(MEDITS_spain_full)
rm(MEDITS_spain_species)
rm(MEDITS_spain_hauls)
rm(MEDITS_spain_complete)



#==============
# Bind together
#==============

MEDITS_elasmo_ts <- bind_rows(MEDITS_nw_elasmo_ts, MEDITS_spain_elasmo_ts)

write.csv(MEDITS_elasmo_ts, "MEDITS_timeseries.csv")









