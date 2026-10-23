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
library(stringr)
library(readr)

# Load taxonomy
taxonomy_lookup <- read.csv("taxonomy_lookup.csv")

# Load CSV
sharkipedia <- read.csv("sharkipedia_trends_combined.csv")

print(n = 96, result %>% count(general_data_type, unit, unit_freeform, sort = TRUE))

# Remove NAs, confirmed that these are in years when surveys haven't been done
# VERY important that these are NOT turned into zeros
sharkipedia_clean <- sharkipedia %>%
  filter(!is.na(suppressWarnings(as.numeric(measurement))))

# Into wide format
# Placing measurements into num_cpue etc. columns where appropriate, or into an other column where not possible
rate_types <- c("sCPUE", "nCPUE", "sSPUE", "nSPUE")

sharkipedia_wide <- sharkipedia_clean %>%
  mutate(
    value        = suppressWarnings(as.numeric(measurement)),
    unit_chr     = coalesce(unit, ""),
    is_num       = str_detect(unit_chr, "^individual"),
    is_wgt       = str_detect(unit_chr, "^(kg|ton)(\\s|$)"),
    is_area      = str_detect(unit_chr, "per (km2|m2)$"),
    is_rate_type = general_data_type %in% rate_types,
    standardized = str_starts(general_data_type, "s"),
    effort_unit  = str_remove(unit_chr, "^(individual|kg|ton)( per )?"),
    target = case_when(
      is_num & is_area      ~ "num_cpua",
      is_wgt & is_area      ~ "wgt_cpua",
      is_num & is_rate_type ~ "num_cpue",
      is_wgt & is_rate_type ~ "wgt_cpue",
      TRUE                  ~ "other"
    ),
    num_cpue = if_else(target == "num_cpue", value, NA_real_),
    num_cpua = if_else(target == "num_cpua", value, NA_real_),
    wgt_cpue = if_else(target == "wgt_cpue", value, NA_real_),
    wgt_cpua = if_else(target == "wgt_cpua", value, NA_real_),
    other    = if_else(target == "other",    value, NA_real_)
  ) %>%
  select(dataset_id, year,
         num_cpue, num_cpua, wgt_cpue, wgt_cpua, other,
         standardized, general_data_type, unit, effort_unit, unit_freeform,
         everything(),
         -measurement, -value, -unit_chr,
         -is_num, -is_wgt, -is_area, -is_rate_type, -target)

# Sanity checks
nrow(sharkipedia_wide) == nrow(sharkipedia_clean)   # should be TRUE: no rows lost
sharkipedia_wide %>%
  summarise(across(c(num_cpue, num_cpua, wgt_cpue, wgt_cpua, other),
                   ~ sum(!is.na(.x))))

# Checking what ended up in "other"
sharkipedia_wide %>%
     filter(!is.na(other)) %>%
     count(general_data_type, unit, sort = TRUE)


# Checking that each dataset uses one unit
result %>%
  group_by(dataset_id) %>%
  summarise(n_units = n_distinct(unit), n_freeform = n_distinct(unit_freeform)) %>%
  filter(n_units > 1 | n_freeform > 1)


# Rearranging the dataset for consistency with others
sharkipedia_wide_edited <- sharkipedia_wide %>%
  unite(survey, dataset_id, author_year, sep = "_", remove = FALSE) %>%
  unite(haul_id, dataset_id, author_year, sep = "_") %>%
  mutate(
    survey = if_else(str_starts(survey, "sharkipedia_"),
                     survey,
                     paste0("sharkipedia_", survey)),
    lat_cell = floor(latitude),
    lon_cell = floor(longitude)
  ) %>%
  rename(accepted_name = species, class = subclass) %>%
  mutate(
    accepted_name = case_when(
      haul_id %in% c("3453_baum2010",
                     "3454_baum2003",
                     "3455_cortes2007a",
                     "3456_cortes2007a",
                     "3457_cortes2007a",
                     "3458_cortes2002b",
                     "3459_rice2015") ~ "Alopiidae",
      haul_id %in% c("3496_wueringer2017",
                     "3495_wueringer2017",
                     "4012_everett2015") ~ "Pristidae",
      haul_id %in% c("3548_baum2003",
                     "3554_Steer2020",
                     "3591_Grubert2013",
                     "3642_Roff2018",
                     "3643_Roff2018") ~ "Carcharhinidae",
      haul_id %in% c("3656_walker2007b",
                     "3707_ICES2019",
                     "3708_ICES2019",
                     "4306_kyne2007k",
                     "3657_walker2007b") ~ "Centrophoridae",
      haul_id %in% c("3716_ICES2018",
                     "3717_ICES2018",
                     "3887_Marandel2019") ~ "Rajidae",
      haul_id %in% c("3871_baum2003",
                     "3872_rice2015",
                     "3873_rice2015",
                     "3874_cortes2002b") ~ "Lamnidae",
      haul_id %in% c("3951_Pollackunpubl9",
                     "3952_Pollackunpubl9",
                     "3953_SEDAR2015b") ~ "Triakidae",
      haul_id %in% c("4010_tuck2018a") ~ "Pristiophoridae",
      haul_id %in% c("4153_Roff2018",
                     "4154_Roff2018",
                     "4155_simpfendorfer2011b",
                     "4156_baum2010",
                     "4157_bezerra2016",
                     "4158_baum2003",
                     "4164_reid2011",
                     "4179_McAuley2012",
                     "4180_rice2015") ~ "Sphyrnidae",
      TRUE ~ accepted_name)) %>%
  left_join(taxonomy_lookup %>%
  rename(
  order_lookup = order,
  class_lookup = class), by = "accepted_name") %>%
  mutate(
    order = coalesce(na_if(order, ""), order_lookup),
    class = coalesce(na_if(class, ""), class_lookup)) %>%
  select(-order_lookup, -class_lookup) %>%
  select(survey, haul_id, year, lat_cell, lon_cell, num_cpue, num_cpua, wgt_cpue, 
         wgt_cpua, other, accepted_name, comments, order, class) %>%
  filter(!(haul_id == "4012_everett2015" & year %in% c(1964, 1965)))

# Checking for blank values
sharkipedia_wide_edited %>%
  summarise(
    blank_species = sum(is.na(accepted_name) | accepted_name == ""),
    blank_lat = sum(is.na(lat_cell) | lat_cell == ""),
    blank_lon = sum(is.na(lon_cell) | lon_cell == ""))

# Checking each haul_id corresponds with one species, one lat_cell and one lon_cell
sharkipedia_wide_edited %>%
  group_by(haul_id) %>%
  summarise(
    n_species = n_distinct(accepted_name),
    n_lat = n_distinct(lat_cell),
    n_lon = n_distinct(lon_cell),
    .groups = "drop"
  ) %>%
  filter(n_species != 1 | n_lat != 1 | n_lon != 1)

# Saving the dataset
write_csv(sharkipedia_wide_edited, "sharkipedia_timeseries.csv")


