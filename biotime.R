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
  select(Species, Order, Class, SuperClass) %>%
  rename(accepted_name = Species) %>%
  distinct()

# Create family lookup table
family_lookup <- taxonomy %>%
  select(Family, Order, Class, SuperClass) %>%
  rename(accepted_name = Family) %>%
  distinct()

# Combine
taxonomy_lookup <- bind_rows(species_lookup, family_lookup) %>%
  distinct(accepted_name, .keep_all = TRUE)

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
  rename(order = Order, class = Class, superclass = SuperClass) %>%
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
mcr_mean_cpue <- mcr_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
mcr_mean_cpue %>%
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
unique(mcr$superclass)

mcr_elasmo_ts <- mcr_mean_cpue %>%
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
  rename(order = Order, class = Class, superclass = SuperClass) %>%
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

pelagic_mean_cpue <- pelagic_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
pelagic_mean_cpue %>%
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
unique(pelagic$superclass)

pelagic_elasmo_ts <- pelagic_mean_cpue %>%
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
  rename(order = Order, class = Class, superclass = SuperClass) %>%
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

seych_mean_cpue <- seych_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
seych_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
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
  rename(order = Order, class = Class, superclass = SuperClass) %>%
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

stj_mean_cpue <- stj_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
stj_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )

unique(stj$superclass)
stj_elasmo_ts <- stj_mean_cpue %>%
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
  rename(order = Order, class = Class, superclass = SuperClass) %>%
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

croix_mean_cpue <- croix_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
croix_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )

unique(croix$superclass)
croix_elasmo_ts <- croix_mean_cpue %>%
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
  rename(order = Order, class = Class, superclass = SuperClass) %>%
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

hoga_mean_cpue <- hoga_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class, superclass) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
hoga_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )

unique(hoga$superclass)
hoga_elasmo_ts <- hoga_mean_cpue %>%
  filter(superclass %in% c("Chondrichthyes"))

# Clean house
rm(hoga_full)
rm(hoga_species)
rm(hoga_hauls)
rm(hoga_complete)















#==================
# Bind all together
#==================

biotime_ts <- bind_rows(croix_elasmo_ts, hoga_elasmo_ts, mcr_elasmo_ts, pelagic_elasmo_ts,
                        stj_elasmo_ts)


write_csv(biotime_ts, "Biotime_timeseries.csv")


