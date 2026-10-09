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
library(lubridate) # For extracting year from the eventDate column

# Load taxonomy
taxonomy_lookup <- read.csv("taxonomy_lookup.csv")


# Load data
NSW_full <- read_tsv("NSW.txt")
print(problems(NSW_full), n = 31)

class(NSW_full$eventDate)
head(NSW_full$eventDate)
sum(is.na(NSW_full$eventDate))

NSW_clean <- NSW_full %>%
  select(id, individualCount, occurrenceStatus, eventDate, decimalLatitude, decimalLongitude, scientificName) %>%
  mutate(
    year = year(eventDate),
    lat_cell = floor(decimalLatitude),
    lon_cell = floor(decimalLongitude)
  )

# **CHECK** Each ID corresponds with only one year/grid cell combination
NSW_clean %>%
  distinct(
    id,
    year,
    lat_cell,
    lon_cell
  ) %>%
  count(id) %>%
  filter(n > 1)  

# Make variable names consistent and add species data
NSW <- NSW_clean %>%
  rename(accepted_name = scientificName, num_cpue = individualCount, haul_id = id, presenceabsence = occurrenceStatus) %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
  mutate(survey = "NSW") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, presenceabsence, accepted_name,
         order, class, superclass)

# Find species that have ever been recorded in a grid square
NSW_species <- NSW %>%
  filter(num_cpue > 0) %>% 
  distinct(lat_cell, lon_cell, accepted_name, order, class, superclass)

# **CHECK** Only one row per species per grid cell
NSW_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
NSW_hauls <- NSW %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
NSW_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
NSW_complete <- NSW_hauls %>%
  inner_join(
    NSW_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    NSW %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
NSW_mean_num_cpue <- NSW_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
NSW_mean_num_cpue %>%
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
unique(NSW$superclass)

NSW_elasmo_ts <- NSW_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(NSW_clean)
rm(NSW_full)
rm(NSW_species)
rm(NSW_hauls)
rm(NSW_complete)

write.csv(NSW_elasmo_ts, "NSW_timeseries.csv")

