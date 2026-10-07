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

# Load data
tf_census <- read.csv("TimeFISH_census_data.csv")
tf_location <- read.csv("TimeFISH_location_information.csv", fileEncoding = "latin1")
tf_taxonomy <- read.csv("TimeFISH_taxonomic_information.csv")
taxonomy_lookup <- read.csv("taxonomy_lookup.csv")

# Join together
tf_census_location <- left_join(tf_census, tf_location, by = "transect_id", relationship = "many-to-many")
tf_census_location_taxonomy <- left_join(tf_census_location, tf_taxonomy, by = "species_name", relationship = "many-to-many")
tf_full <- left_join(tf_census_location_taxonomy, taxonomy_lookup, by = "fishbase_code", relationship = "many-to-many")

tf_check <- tf_full %>%
  select(transect_id, accepted_name, abundance, longitude, latitude, year, order, class, superclass) %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  )

# **CHECK** Each transect_id corresponds with only one year/grid cell combination
tf_check %>%
  distinct(
    transect_id,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(transect_id) %>%
  filter(n > 1)

# Rename columns
tf <- tf_check %>%
  mutate(survey = "TimeFISH") %>%
  rename(haul_id = transect_id, num_cpue = abundance) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name, order, class, superclass)

# Find species that have ever been recorded in a grid square
tf_species <- tf %>%
  filter(num_cpue > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class, superclass)

# **CHECK** Only one row per species per grid cell
tf_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
tf_hauls <- tf %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
tf_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
tf_complete <- tf_hauls %>%
  inner_join(
    tf_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    tf %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
tf_mean_num_cpue <- tf_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
tf_mean_num_cpue %>%
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
unique(tf$superclass)

tf_elasmo_ts <- tf_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(tf_census)
rm(tf_location)
rm(tf_taxonomy)
rm(tf_check)
rm(tf_census_location)
rm(tf_census_location_taxonomy)
rm(tf_full)
rm(tf_species)
rm(tf_hauls)
rm(tf_complete)

























