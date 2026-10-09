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

#================================
# Taxonomy - create lookup tables
#================================

# Load taxonomy
taxonomy <- load_taxa(server=getOption("fishbase"))

# Create species lookup table
species_lookup <- taxonomy %>%
  select(Species, SpecCode, Order, Class, SuperClass) %>%
  rename(accepted_name = Species) %>%
  distinct()

# Create family lookup table
family_lookup <- taxonomy %>%
  select(Family, SpecCode, Order, Class, SuperClass) %>%
  rename(accepted_name = Family) %>%
  distinct()

# Combine
taxonomy_lookup <- bind_rows(species_lookup, family_lookup) %>%
  distinct(accepted_name, .keep_all = TRUE) %>%
  rename(fishbase_code = SpecCode, order = Order, class = Class, superclass = SuperClass)

# Check each name only occurs once
taxonomy_lookup %>%
  count(accepted_name) %>%
  filter(n > 1)

# Save as CSV
write_csv(taxonomy_lookup, "taxonomy_lookup.csv")

#=======================================================================
# MCR LTER Coral Reef Long-term Population and Community Dynamics Fishes
#=======================================================================

mcr_full <- read.csv("MCR-LTER.csv")

mcr_full$survey <- "MCR-LTER"

mcr <- mcr_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
mcr %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)


# Make variable names consistent and add species data
mcr <- mcr %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)

# Find species that have ever been recorded in each grid square
mcr_species <- mcr %>%
  filter(num_cpue > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
mcr_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
mcr_hauls <- mcr %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
mcr_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
mcr_complete <- mcr_hauls %>%
  inner_join(
    mcr_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  
  # Add actual observations back in
  left_join(
    mcr %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, 
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass"
    )
  ) %>%
  
  # Any missing observations become zero
  mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
mcr_mean_num_cpue <- mcr_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
mcr_mean_num_cpue %>%
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
unique(mcr$superclass)

mcr_elasmo_ts <- mcr_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(mcr_full)
rm(mcr_species)
rm(mcr_hauls)
rm(mcr_complete)


#====================================
# Pelagic Fish Observations 1968-1999
#====================================

pelagic_full <- read.csv("pelagic_fish.csv")

pelagic_full$survey <- "Pelagic"

pelagic <- pelagic_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
pelagic %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

pelagic <- pelagic %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


pelagic_species <- pelagic %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
pelagic_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

pelagic_hauls <- pelagic %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
pelagic_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

pelagic_complete <- pelagic_hauls %>%
  inner_join(
    pelagic_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    pelagic %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

pelagic_mean_num_cpue <- pelagic_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
pelagic_mean_num_cpue %>%
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
unique(pelagic$superclass)

pelagic_elasmo_ts <- pelagic_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(pelagic_full)
rm(pelagic_species)
rm(pelagic_hauls)
rm(pelagic_complete)


#=====================
# Seychelles reef fish
#=====================

seych_full <- read.csv("seychelles.csv")

seych_full$survey <- "Seychelles"

seych <- seych_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
seych %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

seych <- seych %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


seych_species <- seych %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
seych_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

seych_hauls <- seych %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
seych_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

seych_complete <- seych_hauls %>%
  inner_join(
    seych_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    seych %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

seych_mean_num_cpue <- seych_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
seych_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

# No chondrichthyes
unique(seych$superclass)

# Clean house
rm(seych_full)
rm(seych_species)
rm(seych_hauls)
rm(seych_complete)

#================================================================================
# St. John. USVI Fish Assessment and Monitoring Data (2002 - Present) (NOAA-CCMA)
#================================================================================

stj_full <- read.csv("st_john.csv")

stj_full$survey <- "St John USVI"

stj <- stj_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
stj %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

stj <- stj %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


stj_species <- stj %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
stj_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

stj_hauls <- stj %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
stj_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

stj_complete <- stj_hauls %>%
  inner_join(
    stj_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    stj %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

stj_mean_num_cpue <- stj_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
stj_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(stj$superclass)
stj_elasmo_ts <- stj_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(stj_full)
rm(stj_species)
rm(stj_hauls)
rm(stj_complete)


#=================================================================================
# St. Croix. USVI Fish Assessment and Monitoring Data (2002 - Present) (NOAA-CCMA)
#=================================================================================

croix_full <- read.csv("st_croix.csv")

croix_full$survey <- "St Croix USVI"

croix <- croix_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
croix %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

croix <- croix %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


croix_species <- croix %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
croix_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

croix_hauls <- croix %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
croix_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

croix_complete <- croix_hauls %>%
  inner_join(
    croix_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    croix %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

croix_mean_num_cpue <- croix_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
croix_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(croix$superclass)
croix_elasmo_ts <- croix_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(croix_full)
rm(croix_species)
rm(croix_hauls)
rm(croix_complete)



#===========================================================================
# Operation Wallacea marine site in Hoga, Indonesia - Wakatobi National Park
#===========================================================================

hoga_full <- read.csv("hoga_indonesia.csv")

hoga_full$survey <- "Hoga"

hoga <- hoga_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
hoga %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

hoga <- hoga %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


hoga_species <- hoga %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
hoga_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

hoga_hauls <- hoga %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
hoga_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

hoga_complete <- hoga_hauls %>%
  inner_join(
    hoga_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    hoga %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

hoga_mean_num_cpue <- hoga_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
hoga_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(hoga$superclass)
hoga_elasmo_ts <- hoga_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(hoga_full)
rm(hoga_species)
rm(hoga_hauls)
rm(hoga_complete)



#============================================================
# Marine Fish Underwater Surveys in The Israeli Mediterranean
#============================================================

is_med_full <- read.csv("Israeli-Med.csv")

is_med_full$survey <- "Israeli Med"

is_med <- is_med_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
is_med %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

is_med <- is_med %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


is_med_species <- is_med %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
is_med_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

is_med_hauls <- is_med %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

# **CHECK** Only one row per haul
is_med_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

is_med_complete <- is_med_hauls %>%
  inner_join(
    is_med_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    is_med %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

is_med_mean_num_cpue <- is_med_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
is_med_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(is_med$superclass)
is_med_elasmo_ts <- is_med_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(is_med_full)
rm(is_med_species)
rm(is_med_hauls)
rm(is_med_complete)


#====================
# REVIZEE Program 285
#====================

rev_285_full <- read.csv("REVIZEE_285.csv")

rev_285_full$survey <- "OBIS Brazil"

rev_285 <- rev_285_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
rev_285 %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

rev_285 <- rev_285 %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


rev_285_species <- rev_285 %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
rev_285_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

rev_285_hauls <- rev_285 %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

rev_285_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

rev_285_complete <- rev_285_hauls %>%
  inner_join(
    rev_285_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    rev_285 %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

rev_285_mean_num_cpue <- rev_285_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
rev_285_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(rev_285$superclass)
rev_285_elasmo_ts <- rev_285_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(rev_285_full)
rm(rev_285_species)
rm(rev_285_hauls)
rm(rev_285_complete)


#====================
# REVIZEE Program 284
#====================

rev_284_full <- read.csv("REVIZEE_284.csv")

rev_284_full$survey <- "OBIS Brazil"

rev_284 <- rev_284_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
rev_284 %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

rev_284 <- rev_284 %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


rev_284_species <- rev_284 %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
rev_284_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

rev_284_hauls <- rev_284 %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

rev_284_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

rev_284_complete <- rev_284_hauls %>%
  inner_join(
    rev_284_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    rev_284 %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

rev_284_mean_num_cpue <- rev_284_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
rev_284_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(rev_284$superclass)
rev_284_elasmo_ts <- rev_284_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(rev_284_full)
rm(rev_284_species)
rm(rev_284_hauls)
rm(rev_284_complete)


#====================
# REVIZEE Program 135
#====================

rev_135_full <- read.csv("REVIZEE_135.csv")

rev_135_full$survey <- "OBIS Brazil"

rev_135 <- rev_135_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
rev_135 %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

rev_135 <- rev_135 %>%
  rename(accepted_name = valid_name, wgt_cpue = BIOMAS, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, wgt_cpue, accepted_name,
         order, class, superclass)


rev_135_species <- rev_135 %>%
  filter(wgt_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
rev_135_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

rev_135_hauls <- rev_135 %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

rev_135_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

rev_135_complete <- rev_135_hauls %>%
  inner_join(
    rev_135_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    rev_135 %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

rev_135_mean_wgt_cpue <- rev_135_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
rev_135_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )

unique(rev_135$superclass)
rev_135_elasmo_ts <- rev_135_mean_wgt_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(rev_135_full)
rm(rev_135_species)
rm(rev_135_hauls)
rm(rev_135_complete)


#==============================
# Alcatrazes monitoring program
#==============================

alcatrazes_full <- read.csv("alcatrazes.csv")

alcatrazes_full$survey <- "Alcatrazes"

alcatrazes <- alcatrazes_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
alcatrazes %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

alcatrazes <- alcatrazes %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


alcatrazes_species <- alcatrazes %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
alcatrazes_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

alcatrazes_hauls <- alcatrazes %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

alcatrazes_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

alcatrazes_complete <- alcatrazes_hauls %>%
  inner_join(
    alcatrazes_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    alcatrazes %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

alcatrazes_mean_num_cpue <- alcatrazes_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
alcatrazes_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(alcatrazes$superclass)
alcatrazes_elasmo_ts <- alcatrazes_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(alcatrazes_full)
rm(alcatrazes_species)
rm(alcatrazes_hauls)
rm(alcatrazes_complete)



#==================================
# CRED Rapid Ecological Assessments
#==================================

cred_full <- read.csv("CRED.csv")

cred_full$survey <- "CRED"

cred <- cred_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
cred %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

cred <- cred %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


cred_species <- cred %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

# **CHECK** Only one row per species per grid cell
cred_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

cred_hauls <- cred %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

cred_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

cred_complete <- cred_hauls %>%
  inner_join(
    cred_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    cred %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

cred_mean_num_cpue <- cred_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
cred_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(cred$superclass)
cred_elasmo_ts <- cred_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(cred_full)
rm(cred_species)
rm(cred_hauls)
rm(cred_complete)



#============================
# CSIRO Marine Data Warehouse
#============================

csiro_full <- read.csv("CSIRO.csv")

csiro_full$survey <- "CSIRO"

csiro <- csiro_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(LATITUDE),
    lon_cell = floor(LONGITUDE)
  )

# **CHECK** Each SAMPLE_DESC corresponds with only one year/grid cell combination
csiro %>%
  distinct(
    SAMPLE_DESC,
    YEAR,
    lat_cell,
    lon_cell
  ) %>%
  count(SAMPLE_DESC) %>%
  filter(n > 1)

csiro <- csiro %>%
  rename(accepted_name = valid_name, num_cpue = ABUNDANCE, year = YEAR, haul_id = SAMPLE_DESC) %>%
  left_join(taxonomy_lookup, by="accepted_name") %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, accepted_name,
         order, class, superclass)


csiro_species <- csiro %>%
  filter(num_cpue > 0) %>% 
  distinct(                
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class,
    superclass
  )

csiro_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  filter(n > 1)

csiro_hauls <- csiro %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell)

csiro_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

csiro_complete <- csiro_hauls %>%
  inner_join(
    csiro_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    csiro %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name,
             order, class, superclass, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name",
           "order", "class", "superclass")
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

csiro_mean_num_cpue <- csiro_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
csiro_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )

unique(csiro$superclass)
csiro_elasmo_ts <- csiro_mean_num_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(csiro_full)
rm(csiro_species)
rm(csiro_hauls)
rm(csiro_complete)


#==================
# Bind all together
#==================

biotime_ts <- bind_rows(alcatrazes_elasmo_ts, cred_elasmo_ts, csiro_elasmo_ts, croix_elasmo_ts, hoga_elasmo_ts, is_med_elasmo_ts,
                        mcr_elasmo_ts, pelagic_elasmo_ts, rev_135_elasmo_ts, rev_284_elasmo_ts, rev_285_elasmo_ts,
                        stj_elasmo_ts)


write_csv(biotime_ts, "Biotime_timeseries.csv")


