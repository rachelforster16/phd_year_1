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


#==========================
# IATTC - Purse Seine Catch
#==========================

iattc_ps_flag_wide <- read.csv("IATTC_shark_flag_purse_seine.csv")
iattc_ps_set_wide <- read.csv("IATTC_shark_settype_purse_seine.csv")

# Convert data from wide format to long format
iattc_ps_flag <- iattc_ps_flag_wide %>%
  pivot_longer(
    cols = matches("^[A-Z]+n$"),
    names_to = "species_code",
    names_pattern = "(.*)n$",
    values_to = "count"
  ) %>%
  mutate(
    survey = "IATTC_ps",
    accepted_name = recode(
      species_code,
      BSH = "Prionace glauca",
      CCL = "Carcharhinus limbatus",
      FAL = "Carcharhinus falciformis",
      OCS = "Carcharhinus longimanus",
      SEPN = "Carcharhinus spp.",
      SKH = "Euselachii",
      SPL = "Sphyrna lewini",
      SPN = "Sphyrna spp.",
      SPZ = "Sphyrna zygaena",
      THR = "Alopias spp."
    ),
    lat_cell = floor(LatC1),
    lon_cell = floor(LonC1),
    num_cpue = (count/NumSets)
  ) %>%
  rename(year = Year) %>%
  select(
    survey, year, lat_cell, lon_cell, num_cpue, accepted_name
  )

iattc_ps_flag_elasmo_ts <- iattc_ps_flag %>%
  left_join(taxonomy_lookup, by = "accepted_name") %>%
  select(-fishbase_code)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
