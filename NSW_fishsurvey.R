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
soviet_trawl_full <- read_tsv("soviet_trawl.txt")

soviet_trawl_clean <- soviet_trawl_full %>%
  select(id, individualCount, occurrenceStatus, year, decimalLatitude, decimalLongitude, scientificName) %>%
  mutate(
    lat_cell = floor(decimalLatitude),
    lon_cell = floor(decimalLongitude)
  )

# **CHECK** Each ID corresponds with only one year/grid cell combination
soviet_trawl_clean %>%
  distinct(
    id,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(id) %>%
  filter(n > 1)  

# Make variable names consistent and add species data
soviet_trawl <- soviet_trawl_clean %>%
  rename(accepted_name = scientificName, num_cpue = individualCount, haul_id = id, presenceabsence = occurrenceStatus) %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
  mutate(survey = "SovietTrawl") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, presenceabsence, accepted_name,
         order, class, superclass)

# Find species that have ever been recorded in a grid square
soviet_trawl_species <- soviet_trawl %>%
  filter(num_cpue > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class, superclass)

# **CHECK** Only one row per species per grid cell
soviet_trawl_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
soviet_trawl_hauls <- soviet_trawl %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
soviet_trawl_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
soviet_trawl_complete <- soviet_trawl_hauls %>%
  inner_join(
    soviet_trawl_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    soviet_trawl %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
soviet_trawl_mean_num_cpue <- soviet_trawl_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
soviet_trawl_mean_num_cpue %>%
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
unique(soviet_trawl$superclass)

soviet_trawl_elasmo_ts <- soviet_trawl_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(soviet_trawl_clean)
rm(soviet_trawl_full)
rm(soviet_trawl_species)
rm(soviet_trawl_hauls)
rm(soviet_trawl_complete)

write.csv(soviet_trawl_elasmo_ts, "SovietTrawl_timeseries.csv")