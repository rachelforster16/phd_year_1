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

# Using the set type is better because can standardise by gear type
iattc_ps_set_wide <- read.csv("IATTC_shark_settype_purse_seine.csv")

# **CHECK** Do the two IATTC files agree on total sets?
flag_file <- read.csv("IATTC_shark_flag_purse_seine.csv")

by_settype <- iattc_ps_set_wide %>%
  group_by(Year, Month, LatC1, LonC1) %>%
  summarise(sets_settype = sum(NumSets, na.rm = TRUE), .groups = "drop")

by_flag <- flag_file %>%
  group_by(Year, Month, LatC1, LonC1) %>%
  summarise(sets_flag = sum(NumSets, na.rm = TRUE), .groups = "drop")

full_join(by_settype, by_flag, by = c("Year", "Month", "LatC1", "LonC1")) %>%
  summarise(
     n_cell_months  = n(),
     totals_match   = sum(sets_settype == sets_flag, na.rm = TRUE),
     missing_in_one = sum(is.na(sets_settype) | is.na(sets_flag))
    )

# **CHECK** Is there one row per year / month / cell / set type? (should return 0 rows)
iattc_ps_set_wide %>%
  count(Year, Month, LatC1, LonC1, SetType) %>%
  filter(n > 1)

#  **CHECK** Cell centres (should be .5 becasue then using floor() will match the other datasets)
sort(unique(iattc_ps_set_wide$LatC1))
sort(unique(iattc_ps_set_wide$LonC1))

# **CHECK** No zeros or NA in NumSets (would create problems when calculating CPUE)
sum(iattc_ps_set_wide$NumSets == 0, na.rm = TRUE)
sum(is.na(iattc_ps_set_wide$NumSets))

# **CHECK** Not recorded values stored as NA or as 0 in the species columns
# Zeros can mean two things: No sharks were caught or sharks were caught but recorded in weight instead
iattc_ps_set_wide %>%
  summarise(across(
    matches("^[A-Z]+(n|mt)$"),
    list(na = ~ sum(is.na(.x)), zero = ~ sum(.x == 0, na.rm = TRUE))
  )) %>%
  pivot_longer(everything(),
               names_to = c("column", ".value"),
               names_pattern = "(.*)_(na|zero)")


# **CHECK** Compare NumSets for the same year/month/cell across types for duplication
sets_by_type <- iattc_ps_set_wide %>%
  select(Year, Month, LatC1, LonC1, SetType, NumSets) %>%
  pivot_wider(names_from = SetType, values_from = NumSets, values_fn = sum)

# If set types are copies of each other, eq count is close to both count
# Looking for equality counts to be small (just chance matches)
sets_by_type %>%
  summarise(
    cell_months   = n(),
    one_type_only = sum(is.na(DEL) + is.na(NOA) + is.na(OBJ) == 2),
    del_eq_noa    = sum(DEL == NOA, na.rm = TRUE),
    del_eq_obj    = sum(DEL == OBJ, na.rm = TRUE),
    noa_eq_obj    = sum(NOA == OBJ, na.rm = TRUE),
    del_noa_both  = sum(!is.na(DEL) & !is.na(NOA)),
    del_obj_both  = sum(!is.na(DEL) & !is.na(OBJ)),
    noa_obj_both  = sum(!is.na(NOA) & !is.na(OBJ))
  )


# Convert from wide to long. Each species has a number and weight (in tonnes) column
iattc_long <- iattc_ps_set_wide %>%
  pivot_longer(
    cols = matches("^[A-Z]+(n|mt)$"),
    names_to = c("species_code", ".value"),
    names_pattern = "^([A-Z]+)(n|mt)$"
  ) %>%
  rename(count = n, wgt_mt = mt)

# **CHECK** Most common measurement
unit_check <- iattc_long %>%
  mutate(
    record_type = case_when(
      count > 0 & wgt_mt > 0 ~ "both",
      count > 0              ~ "number_only",
      wgt_mt > 0             ~ "weight_only",
      TRUE                   ~ "no_catch"
    )
  )

unit_check %>%
  filter(record_type != "no_catch") %>%
  count(record_type) %>%
  mutate(percent = round(100 * n / sum(n), 1))

# Most common measurement by species
unit_check %>%
  filter(record_type != "no_catch") %>%
  count(species_code, record_type) %>%
  pivot_wider(names_from = record_type, values_from = n, values_fill = 0)


# Change zeros to NAs (because no NAs anyway, zero just means "nothing recorded")
# And adjust species names
# Then left_join the taxonomy database
iattc_long <- iattc_long %>%
  filter(species_code != "SKH") %>%            # Not identified to genus level
  mutate(
    count    = na_if(count, 0),
    wgt_mt   = na_if(wgt_mt, 0),
    lat_cell = floor(LatC1),
    lon_cell = floor(LonC1),
    accepted_name = recode(
      species_code,
      BSH  = "Prionace glauca",
      CCL  = "Carcharhinus",               # Silky/blacktip complex (see PDF for explanation)
      FAL  = "Carcharhinus",               # Silky/blacktip complex
      OCS  = "Carcharhinus longimanus",
      SEPN = "Carcharhinus",               # Silky/blacktip complex
      SPL  = "Sphyrna lewini",
      SPN  = "Sphyrna",
      SPZ  = "Sphyrna zygaena",
      THR  = "Alopias"
    )
  ) %>%
  left_join(
    taxonomy_lookup %>% select(accepted_name, order, class) %>% distinct(),
    by = "accepted_name"
  )

# Working out whether weight can be estimated from number
# iattc_long %>%
#  filter(!is.na(count), !is.na(wgt_mt)) %>%
#  mutate(kg_per_individual = wgt_mt * 1000 / count) %>%
#  group_by(species_code) %>%
#  summarise(n = n(), min = min(kg_per_individual),
#            median = median(kg_per_individual), max = max(kg_per_individual))


# The weights are implausible - likely that a row's weight and its count are describing different sharks
# So the problem is the rows where there is either weight only or numnber and weight (these are incomplete)
# Only ~ 0.5% of catch records so will drop them

# Classify, merge, and drop weight recorded rows
iattc_long <- iattc_long %>%
  mutate(
    source = case_when(
      !is.na(wgt_mt) ~ "weight_recorded",
      !is.na(count)  ~ "observed",
      TRUE           ~ "no_catch"
    )
  ) %>%
  # One row per original row (year/month/cell/set type) and name.
  # FAL, CCL and SEPN merge here, with NumSets counted once.
  group_by(Year, Month, SetType, LatC1, LonC1, NumSets,
           lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    catch      = sum(count[source == "observed"]),
    any_weight = any(source == "weight_recorded"),
    .groups = "drop"
  )

# How many rows are dropped, by year?
iattc_long %>%
  filter(catch > 0 | any_weight) %>%
  group_by(Year) %>%
  summarise(pct_dropped = round(100 * mean(any_weight), 1)) %>%
  print(n = Inf)

iattc_long <- iattc_long %>%
  filter(!any_weight) %>%
  select(-any_weight)

# Weights only recorded between 1993 and 2004
# Equal across set type and small
# "Records reported in weight (about 0.5% overall, all before 2005) were excluded"

# Sets should NOT be triple-counted: totals should be almost identical across names
iattc_long %>%
  group_by(accepted_name) %>%
  summarise(total_sets = sum(NumSets))

# Creating num_cpue column
iattc_ps_set <- iattc_long %>%
  group_by(set = SetType, year = Year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    catch    = sum(catch),
    num_sets = sum(NumSets),
    .groups = "drop"
  ) %>%
  mutate(
    num_cpue = catch / num_sets,
    survey   = paste0("IATTC_ps_", set)
  ) %>%
  group_by(survey, accepted_name, lat_cell, lon_cell) %>%
  filter(any(num_cpue > 0)) %>%
  ungroup()

# Zeros should still be there
iattc_ps_set %>%
  summarise(
    n_rows   = n(),
    n_zero   = sum(num_cpue == 0),
    pct_zero = round(100 * mean(num_cpue == 0), 1)
  )

# Series per set type
iattc_ps_set %>%
  distinct(survey, lat_cell, lon_cell, accepted_name) %>%
  count(survey)

# Anything odd in CPUE? (no NA or Inf)
summary(iattc_ps_set$num_cpue)

# Filtering out any time series where less than 5 years were sampled
iattc_ps_elasmo_ts <- iattc_ps_set %>%
  group_by(survey, lat_cell, lon_cell, accepted_name) %>%
  filter(n_distinct(year) >= 5) %>%
  ungroup() %>%
  select(survey, set, year, lat_cell, lon_cell, num_cpue, accepted_name, order, class)

 
#=======================
# IATTC - Longline Catch
#=======================  

# Using the number dataset for consistency the the purse seine data
# There is no way of converting weight to numbers because both datasets have been segregated
iattc_ll_wide <- read.csv("IATTC_longline.csv")

# **CHECK** Is there one row per year / month / cell / flag? (should return 0 rows)
iattc_ll_wide %>%
  count(Year, Month, LatC5, LonC5, Flag) %>%
  filter(n > 1)

#  **CHECK** Cell centres (should be .5 becasue then using floor() will match the other datasets)
sort(unique(iattc_ll_wide$LatC5))
sort(unique(iattc_ll_wide$LonC5))

# **CHECK** No zeros or NA in Hooks (would create problems when calculating CPUE)
sum(iattc_ll_wide$Hooks== 0, na.rm = TRUE)
sum(is.na(iattc_ll_wide$Hooks))

# Remove rows with no effort and before 2009
# Species level reporting only began in 2009. Including data before 2009
# would result in zeros recorded where there has actually been catch
iattc_ll_wide <- iattc_ll_wide %>%
  filter(Hooks > 0, Year >= 2009)

# **CHECK** Not recorded values stored as NA or as 0 in the species columns
# Zeros can mean two things: No sharks were caught or sharks were caught but recorded in weight instead
iattc_ll_wide %>%
  summarise(across(
    matches("^[A-Z]+(n|mt)$"),
    list(na = ~ sum(is.na(.x)), zero = ~ sum(.x == 0, na.rm = TRUE))
  )) %>%
  pivot_longer(everything(),
               names_to = c("column", ".value"),
               names_pattern = "(.*)_(na|zero)")


# Trnasform to long format
iattc_ll_long <- iattc_ll_wide %>%
  pivot_longer(
    cols = matches("^[A-Z]+(n|mt)$"),
    names_to = c("species_code", ".value"),
    names_pattern = "^([A-Z]+)(n|mt)$"
  ) %>%
  rename(count = n, wgt_mt = mt) %>%
  filter(!species_code %in% c("RSK", "SKH")) %>%          # Remove species not identified to genus level
  mutate(
    count    = na_if(count, 0),
    wgt_mt   = na_if(wgt_mt, 0),
    lat_cell = floor(LatC5),
    lon_cell = floor(LonC5),
    accepted_name = recode(
      species_code,
      BSH = "Prionace glauca",
      CCL = "Carcharhinus limbatus",
      FAL = "Carcharhinus falciformis",
      OCS = "Carcharhinus longimanus",
      MAK = "Isurus",
      SMA = "Isurus oxyrinchus",
      SPN = "Sphyrna",
      THR = "Alopias"
    )
  ) %>%
  left_join(
    taxonomy_lookup %>% select(accepted_name, order, class) %>% distinct(),
    by = "accepted_name"
  ) %>%
  mutate(
    source = case_when(
      !is.na(wgt_mt) ~ "weight_recorded",
      !is.na(count)  ~ "observed",
      TRUE           ~ "no_catch"
    )
  )

# Drop weight-recorded rows and create catch
iattc_ll_long <- iattc_ll_long %>%
  filter(source != "weight_recorded") %>%
  mutate(catch = if_else(source == "observed", count, 0))

# Create num_cpue
iattc_ll <- iattc_ll_long %>%
  group_by(Flag, year = Year, lat_cell, lon_cell, accepted_name, order, class) %>%
  summarise(
    catch   = sum(catch),
    hooks   = sum(Hooks),
    .groups = "drop"
  ) %>%
  mutate(
    num_cpue = catch / hooks * 1000,
    survey   = paste0("IATTC_ll_", Flag)
  ) %>%
  group_by(survey, accepted_name, lat_cell, lon_cell) %>%
  filter(any(num_cpue > 0)) %>%
  ungroup()

# Filter for series which are 5 years or more
iattc_ll_elasmo_ts <- iattc_ll %>%
  group_by(survey, lat_cell, lon_cell, accepted_name) %>%
  filter(n_distinct(year) >= 5) %>%
  ungroup() %>%
  select(survey, year, lat_cell, lon_cell, num_cpue, accepted_name, order, class)


# Bind rows

iattc_elasmo_ts <- bind_rows(iattc_ps_elasmo_ts, iattc_ll_elasmo_ts) %>%
  rename(mean_num_cpue = num_cpue)

write_csv(iattc_elasmo_ts, "IATTC_timeseries.csv")

  
