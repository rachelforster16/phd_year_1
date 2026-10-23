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

#=====================
# Aleutian Island (AI)
#=====================
ai_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/AI_clean.RData",
     envir = ai_env)
ai_full <- ai_env$data

# Modify dataset to get 1x1 grid square for each observation and select necessary variables
ai <- ai_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpua, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(ai_full$haul_id)   # before
n_distinct(ai$haul_id)        # after

# Firstly, find species that have ever been recorded in each grid square
ai_species <- ai %>%
  filter(wgt_cpua > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK** Only one row per species per grid cell
ai_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
ai_hauls <- ai %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK** Only one row per haul
ai_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
ai_complete <- ai_hauls %>%
  inner_join(
    ai_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  
  # Add actual observations back in
  left_join(
    ai %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpua),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name","order", "class"
    )
  ) %>%
  
  # Any missing observations become zero
  mutate(wgt_cpua = replace_na(wgt_cpua, 0))

# Now calculate mean_cpua of all events in a grid square, per species, per year, per survey
ai_mean_wgt_cpua <- ai_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpua = mean(wgt_cpua, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
ai_mean_wgt_cpua %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpua = max(mean_wgt_cpua),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpua = sum(max_mean_wgt_cpua > 0),
    n_with_zero_only = sum(max_mean_wgt_cpua == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(ai_full$class)

ai_elasmo_ts <- ai_mean_wgt_cpua %>%
  filter(class %in% c("Elasmobranchii"))

# Clean house
rm(ai_env)
rm(ai_full)
rm(ai_species)
rm(ai_hauls)
rm(ai_complete)


#==================
# Baltic Sea (BITS)
#==================
bits_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/BITS_clean.RData",
     envir = bits_env)
bits_full <- bits_env$data

# Modify dataset to get 1x1 grid square for each observation and select necessary variables
bits <- bits_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(bits_full$haul_id)   # before
n_distinct(bits$haul_id)        # after

# Firstly, find species that have ever been recorded in each grid square
bits_species <- bits %>%
  filter(num_cpue > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK** Only one row per species per grid cell
bits_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# Using original full dataset as some hauls were lost in the filtering 
bits_hauls <- bits_full %>%
  ungroup() %>%
  mutate(lat_cell = floor(latitude), lon_cell = floor(longitude)) %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK** Only one row per haul
bits_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
bits_complete <- bits_hauls %>%
  inner_join(
    bits_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%

# Add actual observations back in
left_join(
  bits %>%
    select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
  by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name","order", "class"
  )
) %>%

# Any missing observations become zero
mutate(num_cpue = replace_na(num_cpue, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
bits_mean_num_cpue <- bits_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
bits_mean_num_cpue %>%
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
unique(bits_full$class)

bits_elasmo_ts <- bits_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

# Clean house
rm(bits_env)
rm(bits_full)
rm(bits_species)
rm(bits_hauls)
rm(bits_complete)


#=========================
# Eastern Bering Sea (EBS)
#=========================
ebs_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/EBS_clean.RData",
     envir = ebs_env)
ebs_full <- ebs_env$data

# Modify dataset to get 1x1 grid square for each observation and select necessary variables
ebs <- ebs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpua, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(ebs_full$haul_id)   # before
n_distinct(ebs$haul_id)        # after

# Firstly, find species that have ever been recorded in each grid square
ebs_species <- ebs %>%
  filter(wgt_cpua > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK** Only one row per species per grid cell
ebs_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
ebs_hauls <- ebs %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK** Only one row per haul
ebs_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
ebs_complete <- ebs_hauls %>%
  inner_join(
    ebs_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    ebs %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpua),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name","order", "class"
    )
  ) %>%
  mutate(wgt_cpua = replace_na(wgt_cpua, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
ebs_mean_wgt_cpua <- ebs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpua = mean(wgt_cpua, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
ebs_mean_wgt_cpua %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpua = max(mean_wgt_cpua),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpua = sum(max_mean_wgt_cpua > 0),
    n_with_zero_only = sum(max_mean_wgt_cpua == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(ebs_full$class)

ebs_elasmo_ts <- ebs_mean_wgt_cpua %>%
  filter(class %in% c("Elasmobranchii"))

# Clean house
rm(ebs_env)
rm(ebs_full)
rm(ebs_species)
rm(ebs_hauls)
rm(ebs_complete)


#======================
# Bay of Biscay (EVHOE)
#======================
evhoe_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/EVHOE_clean.RData",
     envir = evhoe_env)
evhoe_full <- evhoe_env$data

evhoe <- evhoe_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(evhoe_full$haul_id)   # before
n_distinct(evhoe$haul_id)        # after

evhoe_species <- evhoe %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
evhoe_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

evhoe_hauls <- evhoe %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
evhoe_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

evhoe_complete <- evhoe_hauls %>%
  inner_join(
    evhoe_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    evhoe %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

evhoe_mean_num_cpue <- evhoe_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
evhoe_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(evhoe_full$class)

evhoe_elasmo_ts <- evhoe_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(evhoe_env)
rm(evhoe_full)
rm(evhoe_species)
rm(evhoe_hauls)
rm(evhoe_complete)


#==========================
# English Channel (FR-CGFS)
#==========================
fr_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/FR-CGFS_clean.RData",
     envir = fr_env)
fr_full <- fr_env$data

fr <- fr_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(fr_full$haul_id)   # before
n_distinct(fr$haul_id)        # after

fr_species <- fr %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
fr_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

fr_hauls <- fr %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
fr_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

fr_complete <- fr_hauls %>%
  inner_join(
    fr_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    fr %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

fr_mean_num_cpue <- fr_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
fr_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(fr_full$class)

fr_elasmo_ts <- fr_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(fr_env)
rm(fr_full)
rm(fr_species)
rm(fr_hauls)
rm(fr_complete)


#=======================
# Gulf of Mexico (GMEX)
#=======================
gmex_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/GMEX_clean.RData",
     envir = gmex_env)
gmex_full <- gmex_env$data

gmex <- gmex_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(gmex_full$haul_id)   # before
n_distinct(gmex$haul_id)        # after

gmex_species <- gmex %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
gmex_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

gmex_hauls <- gmex %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
gmex_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

gmex_complete <- gmex_hauls %>%
  inner_join(
    gmex_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    gmex %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

gmex_mean_num_cpue <- gmex_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
gmex_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(gmex_full$class)

gmex_elasmo_ts <- gmex_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(gmex_env)
rm(gmex_full)
rm(gmex_species)
rm(gmex_hauls)
rm(gmex_complete)


#=====================
# Gulf of Alaska (GOA)
#=====================
goa_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/GOA_clean.RData",
     envir = goa_env)
goa_full <- goa_env$data

goa <- goa_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpua, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(goa_full$haul_id)   # before
n_distinct(goa$haul_id)        # after

# Firstly, find species that have ever been recorded in each grid square
goa_species <- goa %>%
  filter(wgt_cpua > 0) %>% # Keeping only positive observations
  distinct(                # Removing duplicates for each grid cell x species combination
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK** Only one row per species per grid cell
goa_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Get a list of all the individual hauls
# As may have many rows per haul
goa_hauls <- goa %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK** Only one row per haul
goa_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

# Create every possible haul x species combination
goa_complete <- goa_hauls %>%
  inner_join(
    goa_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    goa %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpua),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name","order", "class"
    )
  ) %>%
  mutate(wgt_cpua = replace_na(wgt_cpua, 0))

# Now calculate mean_cpue of all events in a grid square, per species, per year, per survey
goa_mean_wgt_cpua <- goa_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpua = mean(wgt_cpua, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
goa_mean_wgt_cpua %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpua = max(mean_wgt_cpua),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpua = sum(max_mean_wgt_cpua > 0),
    n_with_zero_only = sum(max_mean_wgt_cpua == 0)
  )

# Filter for chondrichthyes (elasmobranchii and holocephali/chimaeriformes)
unique(goa_full$class)

goa_elasmo_ts <- goa_mean_wgt_cpua %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

# Clean house
rm(goa_env)
rm(goa_full)
rm(goa_species)
rm(goa_hauls)
rm(goa_complete)

#=====================================
# Northern Gulf of St Lawrence (GSL-N)
#=====================================
gsln_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/GSL-N_clean.RData",
     envir = gsln_env)
gsln_full <- gsln_env$data

# Fixing encoding issue because of é
gsln_full <- gsln_full %>%
  mutate(
    haul_id = iconv(haul_id, from = "latin1", to = "UTF-8")
  )

gsln <- gsln_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(gsln_full$haul_id)   # before
n_distinct(gsln$haul_id)        # after

gsln_species <- gsln %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
gsln_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

gsln_hauls <- gsln %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
gsln_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

gsln_complete <- gsln_hauls %>%
  inner_join(
    gsln_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    gsln %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

gsln_mean_num_cpue <- gsln_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
gsln_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(gsln_full$class)

gsln_elasmo_ts <- gsln_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(gsln_env)
rm(gsln_full)
rm(gsln_species)
rm(gsln_hauls)
rm(gsln_complete)


#=====================================
# Southern Gulf of St Lawrence (GSL-S)
#=====================================
gsls_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/GSL-S_clean.RData",
     envir = gsls_env)
gsls_full <- gsls_env$data

gsls <- gsls_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(gsls_full$haul_id)   # before
n_distinct(gsls$haul_id)        # after

# Some NAs in num_cpue when wgt_cpue has values, so estimate the values using least squares regression

# Settings (See below dataset (HS) for more explanation)
min_n_species <- 3    # min paired rows to trust a species-specific slope
min_n_order   <- 10   # min paired rows to trust an order-level slope
min_n_class   <- 10   # min paired rows to trust a class-level slope

# Training data: rows with BOTH number and weight
train_gsls <- gsls %>%
  filter(!is.na(num_cpue), !is.na(wgt_cpue), wgt_cpue > 0)

# Tier 2: species-specific slopes
sp_slopes_gsls <- train_gsls %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_gsls <- train_gsls %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_gsls <- train_gsls %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
global_slope_gsls <- with(train_gsls, sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2))

# Apply the hierarchy
gsls <- gsls %>%
  left_join(sp_slopes_gsls,    by = "accepted_name") %>%
  left_join(order_slopes_gsls, by = "order") %>%
  left_join(class_slopes_gsls, by = "class") %>%
  mutate(
    num_source = case_when(
      !is.na(num_cpue)                        ~ "observed",
      !is.na(wgt_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(wgt_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(wgt_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(wgt_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    num_cpue = case_when(
      num_source == "observed"      ~ num_cpue,
      num_source == "species_model" ~ sp_slope * wgt_cpue,
      num_source == "order_model"   ~ order_slope * wgt_cpue,
      num_source == "class_model"   ~ class_slope * wgt_cpue,
      num_source == "global_model"  ~ global_slope_gsls * wgt_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# ** CHECK **
count(gsls, num_source)
sum(is.na(gsls$num_cpue))   # should be 0 if all-NA rows were already removed
gsls %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, num_source)

gsls_species <- gsls %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
gsls_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

gsls_hauls <- gsls %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
gsls_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

gsls_complete <- gsls_hauls %>%
  inner_join(
    gsls_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    gsls %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

gsls_mean_num_cpue <- gsls_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
gsls_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(gsls_full$class)

gsls_elasmo_ts <- gsls_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(gsls_env)
rm(gsls_full)
rm(gsls_species)
rm(gsls_hauls)
rm(gsls_complete)
rm(sp_slopes_gsls)
rm(order_slopes_gsls)
rm(class_slopes_gsls)
rm(train_gsls)
rm(global_slope_gsls)


#===========================
# Canada, Hecate Strait (HS)
#===========================
hs_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/HS_clean.RData",
     envir = hs_env)
hs_full <- hs_env$data

# Prepare the dataset
hs <- hs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(hs_full$haul_id)   # before
n_distinct(hs$haul_id)        # after

# Build model to predict missing wgt_cpue based on num_cpue
# hs_model <- lm(wgt_cpue ~ accepted_name * num_cpue, data = hs)
# summary(hs_model)
# anova(hs_model)

# Extract rows where wgt_cpue in NA
# missing_rows_hs <- hs %>%
#  filter(is.na(wgt_cpue)) %>%
#  mutate(.row_id = which(is.na(hs$wgt_cpue))) # Create a temporary row identifier so know where to put predictions back

# Identify species levels known to the model
# known_species_hs <- hs_model$xlevels$accepted_name

# Keep rows with known species levels
# predictable_rows_hs <- missing_rows_hs %>%
#  filter(accepted_name %in% known_species_hs)

# Predict one row at a time, leaving NA if the prediction fails
# predictions_hs <- lapply(seq_len(nrow(predictable_rows_hs)), function(i) {
#  row <- predictable_rows_hs[i,]
#  tryCatch(
#    as.numeric(predict(hs_model, newdata = row)),
#    error = function(e) NA_real_
#  )
# })

# Combine predictions
# predicted_values_hs <- unlist(predictions_hs)

# Fill successful predictions into the original dataset
# valid_predictions_hs <- !is.na(predicted_values_hs)

# hs$wgt_cpue[predictable_rows_hs$.row_id[valid_predictions_hs]] <- predicted_values_hs[valid_predictions_hs]

# Settings
# Setting thresholds, a slope estimated from 1 or 2 points is unreliable, so a group needs at least this many paired rows before its slope is used
min_n_species <- 3    # min paired rows to trust a species-specific slope
min_n_order   <- 10   # min paired rows to trust an order-level slope
min_n_class   <- 10   # min paired rows to trust a class-level slope

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_hs <- hs %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_hs <- train_hs %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_hs <- train_hs %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_hs <- train_hs %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_hs <- with(train_hs, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
hs <- hs %>%
  left_join(sp_slopes_hs,    by = "accepted_name") %>%
  left_join(order_slopes_hs, by = "order") %>%
  left_join(class_slopes_hs, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_hs * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# ** CHECKS **
count(hs, wgt_source)     # How many rows come from each tier
sum(is.na(hs$wgt_cpue))   # Should be 0
hs %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)


# Complete the analysis the same as other datasets
hs <- hs %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

hs_species <- hs %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
hs_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

hs_hauls <- hs %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
hs_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

hs_complete <- hs_hauls %>%
  inner_join(
    hs_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    hs %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

hs_mean_wgt_cpue <- hs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
hs_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )

unique(hs_full$class)

hs_elasmo_ts <- hs_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))


# Clean house
rm(hs_env)
rm(hs_full)
rm(hs_hauls)
rm(hs_complete)
rm(hs_species)
rm(sp_slopes_hs)
rm(order_slopes_hs)
rm(class_slopes_hs)
rm(train_hs)
rm(global_slope_hs)


#====================
# Irish Sea (IE-IGFS)
#====================
ie_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/IE-IGFS_clean.RData",
     envir = ie_env)
ie_full <- ie_env$data

ie <- ie_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(ie_full$haul_id)   # before
n_distinct(ie$haul_id)        # after

ie_species <- ie %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
ie_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

ie_hauls <- ie %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
ie_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

ie_complete <- ie_hauls %>%
  inner_join(
    ie_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    ie %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

ie_mean_num_cpue <- ie_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
ie_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(ie_full$class)

ie_elasmo_ts <- ie_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(ie_env)
rm(ie_full)
rm(ie_species)
rm(ie_hauls)
rm(ie_complete)


#====================
# Northeast US (NEUS)
#====================
neus_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/NEUS_clean.RData",
     envir = neus_env)
neus_full <- neus_env$data

neus <- neus_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(neus_full$haul_id)   # before
n_distinct(neus$haul_id)        # after

# Some NAs in num_cpue when wgt_cpue has values, so estimate the values using least squares regression

# Training data: rows with BOTH number and weight
train_neus <- neus %>%
  filter(!is.na(num_cpue), !is.na(wgt_cpue), wgt_cpue > 0)

# Tier 2: species-specific slopes
sp_slopes_neus <- train_neus %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_neus <- train_neus %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_neus <- train_neus %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
global_slope_neus <- with(train_neus, sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2))

# Apply the hierarchy
neus <- neus %>%
  left_join(sp_slopes_neus,    by = "accepted_name") %>%
  left_join(order_slopes_neus, by = "order") %>%
  left_join(class_slopes_neus, by = "class") %>%
  mutate(
    num_source = case_when(
      !is.na(num_cpue)                        ~ "observed",
      !is.na(wgt_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(wgt_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(wgt_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(wgt_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    num_cpue = case_when(
      num_source == "observed"      ~ num_cpue,
      num_source == "species_model" ~ sp_slope * wgt_cpue,
      num_source == "order_model"   ~ order_slope * wgt_cpue,
      num_source == "class_model"   ~ class_slope * wgt_cpue,
      num_source == "global_model"  ~ global_slope_neus * wgt_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# ** CHECK **
count(neus, num_source)
sum(is.na(neus$num_cpue))   # should be 0 if all-NA rows were already removed
neus %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, num_source)

neus_species <- neus %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
neus_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

neus_hauls <- neus %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
neus_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

neus_complete <- neus_hauls %>%
  inner_join(
    neus_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    neus %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

neus_mean_num_cpue <- neus_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
neus_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(neus_full$class)

neus_elasmo_ts <- neus_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(neus_env)
rm(neus_full)
rm(neus_species)
rm(neus_hauls)
rm(neus_complete)
rm(sp_slopes_neus)
rm(order_slopes_neus)
rm(class_slopes_neus)
rm(train_neus)
rm(global_slope_neus)


#=========================
# Northern Ireland (NIGFS)
#=========================
nigfs_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/NIGFS_clean.RData",
     envir = nigfs_env)
nigfs_full <- nigfs_env$data

nigfs <- nigfs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(nigfs_full$haul_id)   # before
n_distinct(nigfs$haul_id)        # after

nigfs_species <- nigfs %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
nigfs_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

nigfs_hauls <- nigfs %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
nigfs_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

nigfs_complete <- nigfs_hauls %>%
  inner_join(
    nigfs_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    nigfs %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

nigfs_mean_num_cpue <- nigfs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
nigfs_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(nigfs_full$class)

nigfs_elasmo_ts <- nigfs_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(nigfs_env)
rm(nigfs_full)
rm(nigfs_species)
rm(nigfs_hauls)
rm(nigfs_complete)


#=================
# Norway (NOR-BTS)
#=================
nor_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/NOR-BTS_clean.RData",
     envir = nor_env)
nor_full <- nor_env$data

nor <- nor_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(nor_full$haul_id)   # before
n_distinct(nor$haul_id)        # after

nor_species <- nor %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
nor_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

nor_hauls <- nor %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
nor_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

nor_complete <- nor_hauls %>%
  inner_join(
    nor_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    nor %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

nor_mean_num_cpue <- nor_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
nor_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(nor_full$class)

nor_elasmo_ts <- nor_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(nor_env)
rm(nor_full)
rm(nor_species)
rm(nor_hauls)
rm(nor_complete)


#====================
# North Sea (NS-IBTS)
#====================
ns_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/NS-IBTS_clean.RData",
     envir = ns_env)
ns_full <- ns_env$data

ns <- ns_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(ns_full$haul_id)   # before
n_distinct(ns$haul_id)        # after

ns_species <- ns %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
ns_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Using original full dataset as some hauls were lost in the filtering 
ns_hauls <- ns_full %>%
  ungroup() %>%
  mutate(lat_cell = floor(latitude), lon_cell = floor(longitude)) %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
ns_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

ns_complete <- ns_hauls %>%
  inner_join(
    ns_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    ns %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

ns_mean_num_cpue <- ns_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
ns_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(ns_full$class)

ns_elasmo_ts <- ns_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(ns_env)
rm(ns_full)
rm(ns_species)
rm(ns_hauls)
rm(ns_complete)


#===================
# Portugal (PT-IBTS)
#===================
pt_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/PT-IBTS_clean.RData",
     envir = pt_env)
pt_full <- pt_env$data

pt <- pt_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua))) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num_cpue, accepted_name, order, class)

# **CHECK** Not losing any hauls
n_distinct(pt_full$haul_id)   # before
n_distinct(pt$haul_id)        # after

pt_species <- pt %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
pt_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

# Using original full dataset as some hauls were lost in the filtering 
pt_hauls <- pt_full %>%
  ungroup() %>%
  mutate(lat_cell = floor(latitude), lon_cell = floor(longitude)) %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
pt_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

pt_complete <- pt_hauls %>%
  inner_join(
    pt_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    pt %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

pt_mean_num_cpue <- pt_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
pt_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(pt_full$class)

pt_elasmo_ts <- pt_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(pt_env)
rm(pt_full)
rm(pt_species)
rm(pt_hauls)
rm(pt_complete)


#==============================
# Canada, Queen Charlotte (QCS)
#==============================
qcs_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/QCS_clean.RData",
     envir = qcs_env)
qcs_full <- qcs_env$data

# Prepare the dataset
qcs <- qcs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(qcs_full$haul_id)   # before
n_distinct(qcs$haul_id)        # after

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_qcs <- qcs %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_qcs <- train_qcs %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_qcs <- train_qcs %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_qcs <- train_qcs %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_qcs <- with(train_qcs, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
qcs <- qcs %>%
  left_join(sp_slopes_qcs,    by = "accepted_name") %>%
  left_join(order_slopes_qcs, by = "order") %>%
  left_join(class_slopes_qcs, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_qcs * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# Checks
count(qcs, wgt_source)     # How many rows come from each tier
sum(is.na(qcs$wgt_cpue))   # Should be 0
qcs %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)


# Complete the analysis the same as other datasets
qcs <- qcs %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

qcs_species <- qcs %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
qcs_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

qcs_hauls <- qcs %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
qcs_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

qcs_complete <- qcs_hauls %>%
  inner_join(
    qcs_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    qcs %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

qcs_mean_wgt_cpue <- qcs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
qcs_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )

unique(qcs_full$class)

qcs_elasmo_ts <- qcs_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))


# Clean house
rm(qcs_env)
rm(qcs_full)
rm(qcs_hauls)
rm(qcs_complete)
rm(qcs_species)
rm(sp_slopes_qcs)
rm(order_slopes_qcs)
rm(class_slopes_qcs)
rm(train_qcs)
rm(global_slope_qcs)


#==========================
# Rockall Plateau (ROCKALL)
#==========================
rock_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/ROCKALL_clean.RData",
     envir = rock_env)
rock_full <- rock_env$data

# Prepare the dataset
rock <- rock_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(rock_full$haul_id)   # before
n_distinct(rock$haul_id)        # after

rock_species <- rock %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
rock_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

rock_hauls <- rock %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
rock_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

rock_complete <- rock_hauls %>%
  inner_join(
    rock_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    rock %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

rock_mean_num_cpue <- rock_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
rock_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(rock_full$class)

rock_elasmo_ts <- rock_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(rock_env)
rm(rock_full)
rm(rock_species)
rm(rock_hauls)
rm(rock_complete)


#====================
# Scotian Shelf (SCS)
#====================
scs_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SCS_clean.RData",
     envir = scs_env)
scs_full <- scs_env$data

scs <- scs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(scs_full$haul_id)   # before
n_distinct(scs$haul_id)        # after

# Some NAs in num_cpue when wgt_cpue has values, so estimate the values using least squares regression

# Training data: rows with BOTH number and weight
train_scs <- scs %>%
  filter(!is.na(num_cpue), !is.na(wgt_cpue), wgt_cpue > 0)

# Tier 2: species-specific slopes
sp_slopes_scs <- train_scs %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_scs <- train_scs %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_scs <- train_scs %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
global_slope_scs <- with(train_scs, sum(num_cpue * wgt_cpue) / sum(wgt_cpue^2))

# Apply the hierarchy
scs <- scs %>%
  left_join(sp_slopes_scs,    by = "accepted_name") %>%
  left_join(order_slopes_scs, by = "order") %>%
  left_join(class_slopes_scs, by = "class") %>%
  mutate(
    num_source = case_when(
      !is.na(num_cpue)                        ~ "observed",
      !is.na(wgt_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(wgt_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(wgt_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(wgt_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    num_cpue = case_when(
      num_source == "observed"      ~ num_cpue,
      num_source == "species_model" ~ sp_slope * wgt_cpue,
      num_source == "order_model"   ~ order_slope * wgt_cpue,
      num_source == "class_model"   ~ class_slope * wgt_cpue,
      num_source == "global_model"  ~ global_slope_scs * wgt_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# ** CHECK **
count(scs, num_source)
sum(is.na(scs$num_cpue))   # should be 0 if all-NA rows were already removed
scs %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, num_source)

scs_species <- scs %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
scs_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

scs_hauls <- scs %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
scs_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

scs_complete <- scs_hauls %>%
  inner_join(
    scs_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    scs %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

scs_mean_num_cpue <- scs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
scs_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(scs_full$class)

scs_elasmo_ts <- scs_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(scs_env)
rm(scs_full)
rm(scs_species)
rm(scs_hauls)
rm(scs_complete)
rm(sp_slopes_scs)
rm(order_slopes_scs)
rm(class_slopes_scs)
rm(train_scs)
rm(global_slope_scs)


#=====================
# Southeast US (SEUS)
#=====================
seus_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SEUS_clean.RData",
     envir = seus_env)
seus_full <- seus_env$data

seus <- seus_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(seus_full$haul_id)   # before
n_distinct(seus$haul_id)        # after

seus_species <- seus %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
seus_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

seus_hauls <- seus %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
seus_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

seus_complete <- seus_hauls %>%
  inner_join(
    seus_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    seus %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

seus_mean_num_cpue <- seus_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
seus_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(seus_full$class)

seus_elasmo_ts <- seus_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(seus_env)
rm(seus_full)
rm(seus_species)
rm(seus_hauls)
rm(seus_complete)


#=================================
# Canada, Strait of Georgia (SOG)
#=================================
sog_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SOG_clean.RData",
     envir = sog_env)
sog_full <- sog_env$data

sog <- sog_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(sog_full$haul_id)   # before
n_distinct(sog$haul_id)        # after

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_sog <- sog %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_sog <- train_sog %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_sog <- train_sog %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_sog <- train_sog %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_sog <- with(train_sog, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
sog <- sog %>%
  left_join(sp_slopes_sog,    by = "accepted_name") %>%
  left_join(order_slopes_sog, by = "order") %>%
  left_join(class_slopes_sog, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_sog * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# Checks
count(sog, wgt_source)     # How many rows come from each tier
sum(is.na(sog$wgt_cpue))   # Should be 0
sog %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)


# Complete the analysis the same as other datasets
sog <- sog %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

sog_species <- sog %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
sog_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

sog_hauls <- sog %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
sog_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

sog_complete <- sog_hauls %>%
  inner_join(
    sog_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    sog %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

sog_mean_wgt_cpue <- sog_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
sog_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )

unique(sog_full$class)

sog_elasmo_ts <- sog_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))


# Clean house
rm(sog_env)
rm(sog_full)
rm(sog_hauls)
rm(sog_complete)
rm(sog_species)
rm(sp_slopes_sog)
rm(order_slopes_sog)
rm(class_slopes_sog)
rm(train_sog)
rm(global_slope_sog)


#========================
# Gulf of Cadiz (SP-ARSA)
#========================
sparsa_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SP-ARSA_clean.RData",
     envir = sparsa_env)
sparsa_full <- sparsa_env$data

sparsa <- sparsa_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(sparsa_full$haul_id)   # before
n_distinct(sparsa$haul_id)        # after

sparsa_species <- sparsa %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
sparsa_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

sparsa_hauls <- sparsa %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
sparsa_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

sparsa_complete <- sparsa_hauls %>%
  inner_join(
    sparsa_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    sparsa %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

sparsa_mean_num_cpue <- sparsa_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
sparsa_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(sparsa_full$class)

sparsa_elasmo_ts <- sparsa_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(sparsa_env)
rm(sparsa_full)
rm(sparsa_species)
rm(sparsa_hauls)
rm(sparsa_complete)


#=======================
# North Spain (SP-NORTH)
#=======================
spnorth_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SP-NORTH_clean.RData",
     envir = spnorth_env)
spnorth_full <- spnorth_env$data

spnorth <- spnorth_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(spnorth_full$haul_id)   # before
n_distinct(spnorth$haul_id)        # after

spnorth_species <- spnorth %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
spnorth_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

spnorth_hauls <- spnorth %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
spnorth_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

spnorth_complete <- spnorth_hauls %>%
  inner_join(
    spnorth_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    spnorth %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

spnorth_mean_num_cpue <- spnorth_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
spnorth_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(spnorth_full$class)

spnorth_elasmo_ts <- spnorth_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(spnorth_env)
rm(spnorth_full)
rm(spnorth_species)
rm(spnorth_hauls)
rm(spnorth_complete)


#=========================
# Porcupine Bank (SP-PORC)
#=========================
spporc_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SP-PORC_clean.RData",
     envir = spporc_env)
spporc_full <- spporc_env$data

spporc <- spporc_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(spporc_full$haul_id)   # before
n_distinct(spporc$haul_id)        # after

spporc_species <- spporc %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
spporc_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

spporc_hauls <- spporc %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
spporc_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

spporc_complete <- spporc_hauls %>%
  inner_join(
    spporc_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    spporc %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

spporc_mean_num_cpue <- spporc_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
spporc_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(spporc_full$class)

spporc_elasmo_ts <- spporc_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(spporc_env)
rm(spporc_full)
rm(spporc_species)
rm(spporc_hauls)
rm(spporc_complete)


#==============================
# Scotland Shelf Sea (SWC-IBTS)
#==============================
swc_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/SWC-IBTS_clean.RData",
     envir = swc_env)
swc_full <- swc_env$data

swc <- swc_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** Not losing any hauls
n_distinct(swc_full$haul_id)   # before
n_distinct(swc$haul_id)        # after

swc_species <- swc %>%
  filter(num_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
swc_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

swc_hauls <- swc %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
swc_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

swc_complete <- swc_hauls %>%
  inner_join(
    swc_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    swc %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

swc_mean_num_cpue <- swc_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_num_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
swc_mean_num_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_num_cpue = max(mean_num_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_num_cpue > 0),
    n_with_zero_only = sum(max_mean_num_cpue == 0)
  )


unique(swc_full$class)

swc_elasmo_ts <- swc_mean_num_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(swc_env)
rm(swc_full)
rm(swc_species)
rm(swc_hauls)
rm(swc_complete)


#====================================
# California Current (Annual) (WCANN)
#====================================
wcann_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/WCANN_clean.RData",
     envir = wcann_env)
wcann_full <- wcann_env$data

# Prepare the dataset
wcann <- wcann_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(wcann_full$haul_id)   # before
n_distinct(wcann$haul_id)        # after

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_wcann <- wcann %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_wcann <- train_wcann %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_wcann <- train_wcann %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_wcann <- train_wcann %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_wcann <- with(train_wcann, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
wcann <- wcann %>%
  left_join(sp_slopes_wcann,    by = "accepted_name") %>%
  left_join(order_slopes_wcann, by = "order") %>%
  left_join(class_slopes_wcann, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_wcann * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# Checks
count(wcann, wgt_source)     # How many rows come from each tier
sum(is.na(wcann$wgt_cpue))   # Should be 0
wcann %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)

# Analysis as usual
wcann <- wcann %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

wcann_species <- wcann %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
wcann_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

wcann_hauls <- wcann %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
wcann_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

wcann_complete <- wcann_hauls %>%
  inner_join(
    wcann_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    wcann %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

wcann_mean_wgt_cpue <- wcann_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
wcann_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )


unique(wcann_full$class)

wcann_elasmo_ts <- wcann_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(wcann_env)
rm(wcann_full)
rm(wcann_species)
rm(wcann_hauls)
rm(wcann_complete)
rm(sp_slopes_wcann)
rm(order_slopes_wcann)
rm(class_slopes_wcann)
rm(train_wcann)
rm(global_slope_wcann)



#======================================
# Canada, West Coast Haida Gwaii (WCHG)
#======================================
wchg_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/WCHG_clean.RData",
     envir = wchg_env)
wchg_full <- wchg_env$data

# Prepare the dataset
wchg <- wchg_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(wchg_full$haul_id)   # before
n_distinct(wchg$haul_id)        # after

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_wchg <- wchg %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_wchg <- train_wchg %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_wchg <- train_wchg %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_wchg <- train_wchg %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_wchg <- with(train_wchg, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
wchg <- wchg %>%
  left_join(sp_slopes_wchg,    by = "accepted_name") %>%
  left_join(order_slopes_wchg, by = "order") %>%
  left_join(class_slopes_wchg, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_wchg * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# Checks
count(wchg, wgt_source)     # How many rows come from each tier
sum(is.na(wchg$wgt_cpue))   # Should be 0
wchg %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)

# Complete the analysis the same as other datasets
wchg <- wchg %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

wchg_species <- wchg %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
wchg_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

wchg_hauls <- wchg %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
wchg_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

wchg_complete <- wchg_hauls %>%
  inner_join(
    wchg_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    wchg %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

wchg_mean_wgt_cpue <- wchg_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
wchg_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )

unique(wchg_full$class)

wchg_elasmo_ts <- wchg_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))


# Clean house
rm(wchg_env)
rm(wchg_full)
rm(wchg_hauls)
rm(wchg_complete)
rm(wchg_species)
rm(sp_slopes_wchg)
rm(order_slopes_wchg)
rm(class_slopes_wchg)
rm(train_wchg)
rm(global_slope_wchg)



#========================================
# California Current (Trienniall) (WCTRI)
#========================================
wctri_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/WCTRI_clean.RData",
     envir = wctri_env)
wctri_full <- wctri_env$data

# Prepare the dataset
wctri <- wctri_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(wctri_full$haul_id)   # before
n_distinct(wctri$haul_id)        # after

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_wctri <- wctri %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_wctri <- train_wctri %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_wctri <- train_wctri %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_wctri <- train_wctri  %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_wctri <- with(train_wctri, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
wctri <- wctri %>%
  left_join(sp_slopes_wctri,    by = "accepted_name") %>%
  left_join(order_slopes_wctri, by = "order") %>%
  left_join(class_slopes_wctri, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_wctri * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# Checks
count(wctri, wgt_source)     # How many rows come from each tier
sum(is.na(wctri$wgt_cpue))   # Should be 0
wctri %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)

# Analysis as usual
wctri <- wctri %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

wctri_species <- wctri %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
wctri_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

wctri_hauls <- wctri %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
wctri_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

wctri_complete <- wctri_hauls %>%
  inner_join(
    wctri_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    wctri %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

wctri_mean_wgt_cpue <- wctri_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
wctri_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )


unique(wctri_full$class)

wctri_elasmo_ts <- wctri_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(wctri_env)
rm(wctri_full)
rm(wctri_species)
rm(wctri_hauls)
rm(wctri_complete)
rm(sp_slopes_wctri)
rm(order_slopes_wctri)
rm(class_slopes_wctri)
rm(train_wctri)
rm(global_slope_wctri)


#===========================================
# Canada, West Coast Vancouver Island (WCVI)
#===========================================
wcvi_env <- new.env()
load("/Users/mh26992/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data/WCVI_clean.RData",
     envir = wcvi_env)
wcvi_full <- wcvi_env$data

# Prepare the dataset
wcvi <- wcvi_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur, num_cpue, num_cpua,
         wgt_cpue, wgt_cpua, accepted_name, order, class) %>%
  mutate(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua), ~ na_if(.x, 0))) %>% # Turning zeros into NAs
  filter( # Removing any rows which are all NA
    !(is.na(num_cpue) &
        is.na(num_cpua) &
        is.na(wgt_cpue) &
        is.na(wgt_cpua)))

# **CHECK** That not losing any hauls
n_distinct(wcvi_full$haul_id)   # before
n_distinct(wcvi$haul_id)        # after

# Training data: rows with BOTH weight and number
# Also removing rows with zeros as they carry no information
train_wcvi <- wcvi %>%
  filter(!is.na(wgt_cpue), !is.na(num_cpue), num_cpue > 0)

# Tier 2: species-specific slopes
# Take the training data, group by species, summarise the number of rows per species and calculate the 
# slope from the available data. Filter where there aren't enough rows for a reliable estimate 
sp_slopes_wcvi <- train_wcvi %>%
  group_by(accepted_name) %>%
  summarise(
    n_train  = n(),
    sp_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups  = "drop"
  ) %>%
  filter(n_train >= min_n_species) %>%
  select(accepted_name, sp_slope)

# Tier 3: order-level slopes
order_slopes_wcvi <- train_wcvi %>%
  group_by(order) %>%
  summarise(
    n_train     = n(),
    order_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_order) %>%
  select(order, order_slope)

# Tier 4: class-level slopes
class_slopes_wcvi <- train_wcvi  %>%
  group_by(class) %>%
  summarise(
    n_train     = n(),
    class_slope = sum(wgt_cpue * num_cpue) / sum(num_cpue^2),
    .groups     = "drop"
  ) %>%
  filter(n_train >= min_n_class) %>%
  select(class, class_slope)

# Tier 5: global slope
# Calculating a single slope using every paired row in the dataset
global_slope_wcvi <- with(train_wcvi, sum(wgt_cpue * num_cpue) / sum(num_cpue^2))

# Apply the hierarchy
wcvi <- wcvi %>%
  left_join(sp_slopes_wcvi,    by = "accepted_name") %>%
  left_join(order_slopes_wcvi, by = "order") %>%
  left_join(class_slopes_wcvi, by = "class") %>%
  mutate(
    wgt_source = case_when(
      !is.na(wgt_cpue)                        ~ "observed",
      !is.na(num_cpue) & !is.na(sp_slope)     ~ "species_model",
      !is.na(num_cpue) & !is.na(order_slope)  ~ "order_model",
      !is.na(num_cpue) & !is.na(class_slope)  ~ "class_model",
      !is.na(num_cpue)                        ~ "global_model",
      TRUE                                    ~ NA_character_
    ),
    wgt_cpue = case_when(
      wgt_source == "observed"      ~ wgt_cpue,
      wgt_source == "species_model" ~ sp_slope * num_cpue,
      wgt_source == "order_model"   ~ order_slope * num_cpue,
      wgt_source == "class_model"   ~ class_slope * num_cpue,
      wgt_source == "global_model"  ~ global_slope_wcvi * num_cpue
    )
  ) %>%
  select(-sp_slope, -order_slope, -class_slope)

# Checks
count(wcvi, wgt_source)     # How many rows come from each tier
sum(is.na(wcvi$wgt_cpue))   # Should be 0
wcvi %>%                    # Checking how the chondrichthyes have been calculated
  filter(class %in% c("Elasmobranchii", "Holocephali")) %>%
  count(class, wgt_source)

# Analysis as usual
wcvi <- wcvi %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, wgt_cpue, accepted_name, order, class)

wcvi_species <- wcvi %>%
  filter(wgt_cpue > 0) %>%
  distinct(
    lat_cell,
    lon_cell,
    accepted_name,
    order,
    class
  )

# **CHECK**
wcvi_species %>%
  count(lat_cell, lon_cell, accepted_name, order, class) %>%
  filter(n > 1)

wcvi_hauls <- wcvi %>%
  distinct(survey, haul_id, year, lat_cell, lon_cell, depth, haul_dur)

# **CHECK**
wcvi_hauls %>%
  count(survey, haul_id) %>%
  filter(n > 1)

wcvi_complete <- wcvi_hauls %>%
  inner_join(
    wcvi_species,
    by = c("lat_cell", "lon_cell"),
    relationship = "many-to-many"
  ) %>%
  left_join(
    wcvi %>%
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, wgt_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(wgt_cpue = replace_na(wgt_cpue, 0))

wcvi_mean_wgt_cpue <- wcvi_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_wgt_cpue = mean(wgt_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
wcvi_mean_wgt_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_wgt_cpue = max(mean_wgt_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_wgt_cpue > 0),
    n_with_zero_only = sum(max_mean_wgt_cpue == 0)
  )


unique(wcvi_full$class)

wcvi_elasmo_ts <- wcvi_mean_wgt_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(wcvi_env)
rm(wcvi_full)
rm(wcvi_species)
rm(wcvi_hauls)
rm(wcvi_complete)
rm(sp_slopes_wcvi)
rm(order_slopes_wcvi)
rm(class_slopes_wcvi)
rm(train_wcvi)
rm(global_slope_wcvi)



#==================
# Bind all together
#==================

fishglob_ts <- bind_rows(ai_elasmo_ts, bits_elasmo_ts, ebs_elasmo_ts, evhoe_elasmo_ts, 
                         fr_elasmo_ts, gmex_elasmo_ts, goa_elasmo_ts, gsln_elasmo_ts, 
                         gsls_elasmo_ts, hs_elasmo_ts, ie_elasmo_ts, neus_elasmo_ts, 
                         nigfs_elasmo_ts, nor_elasmo_ts, ns_elasmo_ts, pt_elasmo_ts, 
                         qcs_elasmo_ts, rock_elasmo_ts, scs_elasmo_ts, seus_elasmo_ts, 
                         sog_elasmo_ts, sparsa_elasmo_ts, spnorth_elasmo_ts, spporc_elasmo_ts, 
                         swc_elasmo_ts, wcann_elasmo_ts, wchg_elasmo_ts, wctri_elasmo_ts,
                         wcvi_elasmo_ts) %>%
  filter( # Removing rows which are all NA
    !(is.na(mean_num_cpue) & is.na(mean_wgt_cpue) & is.na(mean_wgt_cpua))) %>%
  select(survey, year, lat_cell, lon_cell, accepted_name, order, class, mean_num_cpue,
         mean_wgt_cpue, mean_wgt_cpua, n_hauls)

# **CHECK** That each row has just one metric
fishglob_ts %>%
  mutate(
    n_metrics = rowSums(
      !is.na(select(., mean_num_cpue, mean_wgt_cpue, mean_wgt_cpua)))) %>%
  count(n_metrics)


write_csv(fishglob_ts, "Fishglob_timeseries.csv")


#==========
# Visualise
#==========

# World map
world <- ne_countries(
  scale = "medium",
  returnclass = "sf")

# Rajiformes
rajiformes <- fishglob_ts %>%
  filter(order == "Rajiformes")

ggplot() +
  geom_sf(data = world, fill = "lightblue", colour = "grey40", linewidth = 0.2) +
  geom_tile(data = rajiformes, aes(x = lon_cell, y = lat_cell, fill = mean_cpue)) +
  coord_sf(
    xlim = c(-128, 43),
    ylim = c(24, 81),
    expand = FALSE
  ) +
  facet_wrap(~ year) +
  scale_fill_viridis_c(option = "turbo", trans = "log1p") +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Mean CPUE"
  ) +
  theme_minimal()


# Carcharhiniformes
carcharhiniformes <- fishglob_ts %>%
  filter(order == "Carcharhiniformes")

ggplot() +
  geom_sf(data = world, fill = "lightblue", colour = "grey40", linewidth = 0.2) +
  geom_tile(data = carcharhiniformes, aes(x = lon_cell, y = lat_cell, fill = mean_cpue)) +
  coord_sf(
    xlim = c(-127, 12),
    ylim = c(24, 61),
    expand = FALSE
  ) +
  facet_wrap(~ year) +
  scale_fill_viridis_c(option = "turbo", trans = "log1p") +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = "Mean CPUE"
  ) +
  theme_minimal()



