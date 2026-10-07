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

#=========================================
# IMOS Global Reef Fish Abundance
#=========================================
IMOS_reef_full <- read.csv('IMOS_reef.csv', header = T, sep = ',', skip = 71)

IMOS_reef_full$survey <- "IMOS_reef"

IMOS_reef <- IMOS_reef_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(year = survey_date %>%
           str_split("-", simplify = TRUE) %>%
           .[,1] %>%
           as.numeric())

# **CHECK** Each FID corresponds with only one year/grid cell combination
IMOS_reef %>%
  distinct(
    FID,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(FID) %>%
  filter(n > 1)

# Make variable names consistent and add species data
IMOS_reef <- IMOS_reef %>%
  rename(accepted_name = species_name, num_cpue = total, haul_id = FID) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class)

# ------------------------------------------------
# Calculate zero-filled mean CPUE without creating
# every haul x species combination (too many rows)
# ------------------------------------------------

# Find species that have ever been recorded in each grid square
IMOS_reef_species <- IMOS_reef %>%
  filter(num_cpue > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK** Only one row per species per grid cell
IMOS_reef_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
IMOS_reef_hauls <- IMOS_reef %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
IMOS_reef_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Number of sampled hauls in each grid cell and year
IMOS_reef_hauls_year <- IMOS_reef_hauls %>%
  count(
    survey,
    year,
    lat_cell,
    lon_cell,
    name = "n_hauls"
  )

# Sum of observed CPUE for each species, grid cell and year
IMOS_reef_observed <- IMOS_reef %>%
  group_by(
    survey,
    year,
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  ) %>%
  summarise(
    sum_cpue = sum(num_cpue, na.rm = TRUE),
    .groups = "drop"
  )


# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
IMOS_reef_mean_cpue <- IMOS_reef_species %>%
  inner_join(
    IMOS_reef_hauls_year,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    IMOS_reef_observed,
    by = c(
      "survey",
      "year",
      "lat_cell",
      "lon_cell",
      "accepted_name",
      "order",
      "class"
    )
  ) %>%
  mutate(
    sum_cpue = replace_na(sum_cpue, 0),
    mean_cpue = sum_cpue / n_hauls
  ) %>%
  select(survey, year, lat_cell, lon_cell, accepted_name, order, class, mean_cpue, n_hauls)


# **CHECK** Every species/grid cell combination has >1 positive combination
IMOS_reef_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(IMOS_reef_mean_cpue$class)

IMOS_reef_elasmo_ts <- IMOS_reef_mean_cpue %>%
  filter(class %in% c("Elasmobranchii"))

# Clean house
rm(IMOS_reef_full)
rm(IMOS_reef_species)
rm(IMOS_reef_hauls)
rm(IMOS_reef_hauls_year)
rm(IMOS_reef_observed)



#=========================================
# IMOS Global Cryptobenthic Fish Abundance
#=========================================

IMOS_crypto_full <- read.csv('IMOS_cryptobenthic.csv', header = T, sep = ',', skip = 71)

IMOS_crypto_full$survey <- "IMOS_cryptobenthic"

IMOS_crypto <- IMOS_crypto_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(year = survey_date %>%
           str_split("-", simplify = TRUE) %>%
           .[,1] %>%
           as.numeric())

# **CHECK** Each FID corresponds with only one year/grid cell combination
IMOS_crypto %>%
  distinct(
    FID,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(FID) %>%
  filter(n > 1)


# Make variable names consistent and add species data
IMOS_crypto <- IMOS_crypto %>%
  rename(accepted_name = species_name, num_cpue = total, haul_id = FID) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class)

# Find species that have ever been recorded in each grid square
IMOS_crypto_species <- IMOS_crypto %>%
  filter(num_cpue > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK** Only one row per species per grid cell
IMOS_crypto_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
IMOS_crypto_hauls <- IMOS_crypto %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
IMOS_crypto_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
IMOS_crypto_complete <- IMOS_crypto_hauls %>%
  inner_join(
    IMOS_crypto_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%left_join(
    IMOS_crypto %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
IMOS_crypto_mean_cpue <- IMOS_crypto_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
IMOS_crypto_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(IMOS_crypto$class)

IMOS_crypto_elasmo_ts <- IMOS_crypto_mean_cpue %>%
  filter(class %in% c("Chondrichthyes"))

# Clean house
rm(IMOS_crypto_full)
rm(IMOS_crypto_species)
rm(IMOS_crypto_hauls)
rm(IMOS_crypto_complete)


#==============
# Bind together
#==============

IMOS_ts <- bind_rows(IMOS_reef_elasmo_ts, IMOS_crypto_elasmo_ts)

write.csv(IMOS_ts, "IMOS_timeseries.csv")





