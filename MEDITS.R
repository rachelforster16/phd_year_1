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

# Load data
MEDITS_catch <- read.csv("MEDITS_catch.csv", sep = ";")
MEDITS_operation <- read.csv("MEDITS_operation.csv", sep = ";")
MEDITS_catch_operation <- left_join(MEDITS_catch, MEDITS_operation, by = "haulID") %>%
  rename(accepted_name = scientificName)
MEDITS_full <- left_join(MEDITS_catch_operation, taxonomy_lookup, by = "accepted_name")

MEDITS_full <- MEDITS_full %>%
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
MEDITS_full %>%
  distinct(
    haulID,
    year.x,
    start_lat_cell,
    start_lon_cell
  ) %>%
  count(haulID) %>%
  filter(n > 1)

# Rename columns
MEDITS <- MEDITS_full %>%
  mutate(num_cpue = (totalNumber/haulDurMin)) %>%
  rename(survey = serie.x, year = year.x, haul_id = haulID, lat_cell = start_lat_cell,
         lon_cell = start_lon_cell) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)

# Find species that have ever been recorded in a grid square
MEDITS_species <- MEDITS %>%
  filter(num_cpue > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class, superclass)

# **CHECK** Only one row per species per grid cell
MEDITS_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
MEDITS_hauls <- MEDITS %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
MEDITS_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
MEDITS_complete <- MEDITS_hauls %>%
  inner_join(
    MEDITS_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    MEDITS %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
MEDITS_mean_num_cpue <- MEDITS_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
MEDITS_mean_num_cpue %>%
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
unique(MEDITS$superclass)

MEDITS_elasmo_ts <- MEDITS_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(MEDITS_catch)
rm(MEDITS_catch_operation)
rm(MEDITS_full)
rm(MEDITS_species)
rm(MEDITS_hauls)
rm(MEDITS_complete)

write.csv(MEDITS_elasmo_ts, "MEDITS_timeseries.csv")










