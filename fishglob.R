# Set working directory
setwd("~/OneDrive - University of Bristol/PhD year 1/R code and data")

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

# No CPUE measurement

ai_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/AI_clean.RData", envir = ai_env)
ai_full <- ai_env$data

ai <- ai_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(ai_env)
rm(ai_full)
rm(ai)


#==================
# Baltic Sea (BITS)
#==================
bits_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/BITS_clean.RData", envir = bits_env)
bits_full <- bits_env$data

# Modify dataset to get 1x1 grid square for each observation and select necessary variables
bits <- bits_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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
# As may have many rows per haul
bits_hauls <- bits %>%
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
bits_mean_cpue <- bits_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK** Every species/grid cell combination has >1 positive combination
bits_mean_cpue %>%
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
unique(bits_full$class)

bits_elasmo_ts <- bits_mean_cpue %>%
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

# No CPUE measurement

ebs_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/EBS_clean.RData", envir = ebs_env)
ebs_full <- ebs_env$data

ebs <- ebs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(ebs_env)
rm(ebs_full)
rm(ebs)


#======================
# Bay of Biscay (EVHOE)
#======================
evhoe_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/EVHOE_clean.RData", envir = evhoe_env)
evhoe_full <- evhoe_env$data

evhoe <- evhoe_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

evhoe_mean_cpue <- evhoe_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
evhoe_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(evhoe_full$class)

evhoe_elasmo_ts <- evhoe_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/FR-CGFS_clean.RData", envir = fr_env)
fr_full <- fr_env$data

fr <- fr_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

fr_mean_cpue <- fr_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
fr_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(fr_full$class)

fr_elasmo_ts <- fr_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/GMEX_clean.RData", envir = gmex_env)
gmex_full <- gmex_env$data

gmex <- gmex_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

gmex_mean_cpue <- gmex_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
gmex_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(gmex_full$class)

gmex_elasmo_ts <- gmex_mean_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(gmex_env)
rm(gmex_full)
rm(gmex_species)
rm(gmex_hauls)
rm(gmex_complete)


#=====================
# Gulf of Alaska (GOA)
#=====================

# No CPUE measurement

goa_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/GOA_clean.RData", envir = goa_env)
goa_full <- goa_env$data

goa <- goa_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(goa_env)
rm(goa_full)
rm(goa)


#=====================================
# Northern Gulf of St Lawrence (GSL-N)
#=====================================
gsln_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/GSL-N_clean.RData", envir = gsln_env)
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
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

gsln_mean_cpue <- gsln_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
gsln_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(gsln_full$class)

gsln_elasmo_ts <- gsln_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/GSL-S_clean.RData", envir = gsls_env)
gsls_full <- gsls_env$data

gsls <- gsls_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

gsls_mean_cpue <- gsls_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
gsls_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(gsls_full$class)

gsls_elasmo_ts <- gsls_mean_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(gsls_env)
rm(gsls_full)
rm(gsls_species)
rm(gsls_hauls)
rm(gsls_complete)


#===========================
# Canada, Hecate Strait (HS)
#===========================

# Some observations do not have CPUE values

hs_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/HS_clean.RData", envir = hs_env)
hs_full <- hs_env$data

hs <- hs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(hs_env)
rm(hs_full)
rm(hs)


#====================
# Irish Sea (IE-IGFS)
#====================
ie_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/IE-IGFS_clean.RData", envir = ie_env)
ie_full <- ie_env$data

ie <- ie_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

ie_mean_cpue <- ie_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
ie_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(ie_full$class)

ie_elasmo_ts <- ie_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/NEUS_clean.RData", envir = neus_env)
neus_full <- neus_env$data

neus <- neus_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

neus_mean_cpue <- neus_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
neus_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(neus_full$class)

neus_elasmo_ts <- neus_mean_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(neus_env)
rm(neus_full)
rm(neus_species)
rm(neus_hauls)
rm(neus_complete)


#=========================
# Northern Ireland (NIGFS)
#=========================
nigfs_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/NIGFS_clean.RData", envir = nigfs_env)
nigfs_full <- nigfs_env$data

nigfs <- nigfs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

nigfs_mean_cpue <- nigfs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
nigfs_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(nigfs_full$class)

nigfs_elasmo_ts <- nigfs_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/NOR-BTS_clean.RData", envir = nor_env)
nor_full <- nor_env$data

nor <- nor_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

nor_mean_cpue <- nor_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
nor_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(nor_full$class)

nor_elasmo_ts <- nor_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/NS-IBTS_clean.RData", envir = ns_env)
ns_full <- ns_env$data

ns <- ns_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

ns_hauls <- ns %>%
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

ns_mean_cpue <- ns_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
ns_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(ns_full$class)

ns_elasmo_ts <- ns_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/PT-IBTS_clean.RData", envir = pt_env)
pt_full <- pt_env$data

pt <- pt_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

pt_hauls <- pt %>%
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

pt_mean_cpue <- pt_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
pt_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(pt_full$class)

pt_elasmo_ts <- pt_mean_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(pt_env)
rm(pt_full)
rm(pt_species)
rm(pt_hauls)
rm(pt_complete)


#==============================
# Canada, Queen Charlotte (QCS)
#==============================

# No CPUE measurement

qcs_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/QCS_clean.RData", envir = qcs_env)
qcs_full <- qcs_env$data

qcs <- qcs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(qcs_env)
rm(qcs_full)
rm(qcs)


#==========================
# Rockall Plateau (ROCKALL)
#==========================
rock_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/ROCKALL_clean.RData", envir = rock_env)
rock_full <- rock_env$data

rock <- rock_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

rock_mean_cpue <- rock_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
rock_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(rock_full$class)

rock_elasmo_ts <- rock_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/SCS_clean.RData", envir = scs_env)
scs_full <- scs_env$data

scs <- scs_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

scs_mean_cpue <- scs_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
scs_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(scs_full$class)

scs_elasmo_ts <- scs_mean_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(scs_env)
rm(scs_full)
rm(scs_species)
rm(scs_hauls)
rm(scs_complete)


#=====================
# Southeast US (SEUS)
#=====================
seus_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/SEUS_clean.RData", envir = seus_env)
seus_full <- seus_env$data

seus <- seus_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

seus_mean_cpue <- seus_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
seus_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(seus_full$class)

seus_elasmo_ts <- seus_mean_cpue %>%
  filter(class %in% c("Elasmobranchii"))

rm(seus_env)
rm(seus_full)
rm(seus_species)
rm(seus_hauls)
rm(seus_complete)


#=================================
# Canada, Strait of Georgia (SOG)
#=================================

# No CPUE measurement

sog_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/SOG_clean.RData", envir = sog_env)
sog_full <- sog_env$data

sog <- sog_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(sog_env)
rm(sog_full)
rm(sog)


#========================
# Gulf of Cadiz (SP-ARSA)
#========================
sparsa_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/SP-ARSA_clean.RData", envir = sparsa_env)
sparsa_full <- sparsa_env$data

sparsa <- sparsa_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

sparsa_mean_cpue <- sparsa_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
sparsa_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(sparsa_full$class)

sparsa_elasmo_ts <- sparsa_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/SP-NORTH_clean.RData", envir = spnorth_env)
spnorth_full <- spnorth_env$data

spnorth <- spnorth_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

spnorth_mean_cpue <- spnorth_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
spnorth_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(spnorth_full$class)

spnorth_elasmo_ts <- spnorth_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/SP-PORC_clean.RData", envir = spporc_env)
spporc_full <- spporc_env$data

spporc <- spporc_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

spporc_mean_cpue <- spporc_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
spporc_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(spporc_full$class)

spporc_elasmo_ts <- spporc_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/SWC-IBTS_clean.RData", envir = swc_env)
swc_full <- swc_env$data

swc <- swc_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

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

swc_mean_cpue <- swc_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
swc_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(swc_full$class)

swc_elasmo_ts <- swc_mean_cpue %>%
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
load("/Users/rachelforster/Documents/R code and data/Raw data/WCANN_clean.RData", envir = wcann_env)
wcann_full <- wcann_env$data

wcann <- wcann_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

wcann_species <- wcann %>%
  filter(num_cpue > 0) %>%
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
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

wcann_mean_cpue <- wcann_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
wcann_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(wcann_full$class)

wcann_elasmo_ts <- wcann_mean_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(wcann_env)
rm(wcann_full)
rm(wcann_species)
rm(wcann_hauls)
rm(wcann_complete)


#======================================
# Canada, West Coast Haida Gwaii (WCHG)
#======================================

# No CPUE measurement

wchg_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/WCHG_clean.RData", envir = wchg_env)
wchg_full <- wchg_env$data

wchg <- wchg_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(wchg_env)
rm(wchg_full)
rm(wchg)


#========================================
# California Current (Trienniall) (WCTRI)
#========================================
wctri_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/WCTRI_clean.RData", envir = wctri_env)
wctri_full <- wctri_env$data

wctri <- wctri_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

wctri_species <- wctri %>%
  filter(num_cpue > 0) %>%
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
      select(survey, haul_id, year, lat_cell, lon_cell, accepted_name, order, class, num_cpue),
    by = c("survey", "haul_id", "year", "lat_cell", "lon_cell", "accepted_name", "order", "class"
    )
  ) %>%
  mutate(num_cpue = replace_na(num_cpue, 0))

wctri_mean_cpue <- wctri_complete %>%
  group_by(survey, year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    mean_cpue = mean(num_cpue, na.rm = TRUE),
    n_hauls = n(),
    .groups = "drop"
  )

# **CHECK**
wctri_mean_cpue %>%
  group_by(lat_cell, lon_cell, accepted_name) %>%
  summarise(
    max_mean_cpue = max(mean_cpue),
    .groups = "drop") %>%
  summarise(
    n_species = n(),
    n_with_positive_cpue = sum(max_mean_cpue > 0),
    n_with_zero_only = sum(max_mean_cpue == 0)
  )


unique(wctri_full$class)

wctri_elasmo_ts <- wctri_mean_cpue %>%
  filter(class %in% c("Elasmobranchii", "Holocephali"))

rm(wctri_env)
rm(wctri_full)
rm(wctri_species)
rm(wctri_hauls)
rm(wctri_complete)


#===========================================
# Canada, West Coast Vancouver Island (WCVI)
#===========================================
wcvi_env <- new.env()
load("/Users/rachelforster/Documents/R code and data/Raw data/WCVI_clean.RData", envir = wcvi_env)
wcvi_full <- wcvi_env$data

wcvi <- wcvi_full %>%
  ungroup() %>%
  mutate(
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, depth,
         haul_dur, num, num_cpue, accepted_name, order, class)

rm(wcvi_env)
rm(wcvi_full)
rm(wcvi)


#==================
# Bind all together
#==================

fishglob_ts <- bind_rows(bits_elasmo_ts, evhoe_elasmo_ts, fr_elasmo_ts, gmex_elasmo_ts,
                         gsln_elasmo_ts, gsls_elasmo_ts, ie_elasmo_ts, neus_elasmo_ts,
                         nigfs_elasmo_ts, nor_elasmo_ts, ns_elasmo_ts, pt_elasmo_ts,
                         rock_elasmo_ts, scs_elasmo_ts, seus_elasmo_ts, sparsa_elasmo_ts,
                         spnorth_elasmo_ts, spporc_elasmo_ts, swc_elasmo_ts, wcann_elasmo_ts,
                         wctri_elasmo_ts)

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



