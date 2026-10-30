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
IOTC_catch <- read.csv("IOTC-catch.csv")
IOTC_effort <- read.csv("IOTC-effort.csv")

# For identifiying rows and strata
keys <- c("YEAR", "QUARTER", "MONTH_START", "MONTH_END", "FISHING_GROUND_CODE",
          "FLEET_CODE", "FISHERY_TYPE_CODE", "FISHERY_GROUP_CODE", "FISHERY_CODE",
          "GEAR_CODE", "SCHOOL_TYPE_CODE")

# Getting the location information from the catch database
# Make the code a clean 7-character string (keeps any leading zeros)
IOTC_catch$grid_code <- sprintf("%07d", as.integer(IOTC_catch$FISHING_GROUND_CODE))

# Keep only 1x1 (size digit 5) and 5x5 (size digit 6) cells
# This information was retrieved from:
# https://www.fao.org/cwp-on-fishery-statistics/handbook/general-concepts/main-water-areas/en/#c737533
IOTC_catch$size_digit <- substr(IOTC_catch$grid_code, 1, 1)
IOTC_catch_latlon <- IOTC_catch[IOTC_catch$size_digit %in% c("5", "6"), ]

# Decode the lookup (unique codes only) into cell-centre coordinates
decode_grid <- function(code) {
  type <- substr(code, 1, 1)
  quad <- as.integer(substr(code, 2, 2))
  la   <- as.numeric(substr(code, 3, 4))   # latitude of corner nearest equator
  lo   <- as.numeric(substr(code, 5, 7))   # longitude of corner nearest Greenwich
  
  size <- ifelse(type == "5", 1, ifelse(type == "6", 5, NA))   # square cells
  
  sy <- ifelse(quad %in% c(1, 4),  1, -1)  # 1 NE, 2 SE, 3 SW, 4 NW
  sx <- ifelse(quad %in% c(1, 2),  1, -1)
  
  data.frame(
    grid_code = code,
    grid_size = size,
    grid_lat  = sy * (la + size / 2),
    grid_lon  = sx * (lo + size / 2)
  )
}

lookup <- decode_grid(unique(IOTC_catch_latlon$grid_code))

# 4. Join back to the data
IOTC_catch_latlon <- merge(IOTC_catch_latlon, lookup, by = "grid_code", all.x = TRUE)

# 5. Checks
table(IOTC_catch_latlon$grid_size, useNA = "ifany")          # should be only 1 and 5
sum(is.na(IOTC_catch_latlon$grid_lat))                        # should be 0
range(IOTC_catch_latlon$grid_lat); range(IOTC_catch_latlon$grid_lon)  # should sit in the Indian Ocean
plot(IOTC_catch_latlon$grid_lon, IOTC_catch_latlon$grid_lat, pch = 16, cex = 0.4, col = "#00000030")
maps::map("world", add = TRUE)   # needs the maps package

# Should both return true
# The corner of 1 degree cells can be any integer
# The corner of 5 degree cells should be multiples of 5
d <- IOTC_catch_latlon[!duplicated(IOTC_catch_latlon$grid_code), ]
la <- as.integer(substr(d$grid_code, 3, 4)); lo <- as.integer(substr(d$grid_code, 5, 7))
all(la[d$size_digit == "6"] %% 5 == 0); all(lo[d$size_digit == "6"] %% 5 == 0)


# Filtering for elasmobranchs
IOTC_catch_latlon %>%
  distinct(SPECIES_CATEGORY_CODE, SPECIES_CATEGORY)

IOTC_catch_latlon %>%
  filter(SPECIES_CATEGORY_CODE %in% c("SHARKS", "RAYS")) %>%
  count(SPECIES_CODE, SPECIES, SPECIES_SCIENTIFIC, CATCH_UNIT_CODE)

elasmobranchs <- IOTC_catch_latlon %>%
  filter(SPECIES_CATEGORY_CODE %in% c("SHARKS", "RAYS"))

# Filtering out where sharks and rays were recorded not to genus level
nei_codes <- c("SKH", "UNCL", "SRX", "RAJ", "EAG", "MAN")

# Adding a flag in case needed for sensitivity analyses
nei_strata <- elasmobranchs %>%
  filter(SPECIES_CODE %in% nei_codes) %>%
  distinct(across(all_of(keys))) %>%
  mutate(has_nei_report = TRUE)

elasmo_species <- elasmobranchs %>%
  filter(!SPECIES_CODE %in% nei_codes) %>%
  left_join(nei_strata, by = keys, relationship = 'many-to-one') %>%
  mutate(has_nei_report = replace_na(has_nei_report, FALSE))

nrow(elasmobranchs); nrow(elasmo_species)

# Attach effort
effort_wide <- IOTC_effort %>%
  select(all_of(keys), EFFORT_UNIT_CODE, EFFORT) %>%
  pivot_wider(names_from = EFFORT_UNIT_CODE, values_from = EFFORT, names_prefix = "effort_")

elasmo_species <- elasmo_species %>%
  left_join(effort_wide, by = keys, relationship = "many-to-one")

# **CHECK** All rows have effort
effort_cols <- grep("^effort_", names(elasmo_species), value = TRUE)
elasmo_species$any_effort <- rowSums(!is.na(elasmo_species[effort_cols])) > 0

elasmo_species %>% count(any_effort) # Should all be TRUE

# Are some reports recorded in both number and weight?
# A report is one stratum + species
rep_keys <- c(keys, "SPECIES_CODE")

# Same report and same unit more than once? Should be 0
elasmo_species %>%
  count(across(all_of(c(rep_keys, "CATCH_UNIT_CODE")))) %>%
  filter(n > 1) %>%
  nrow()

# How many reports have one unit vs both?
report_units <- elasmo_species %>%
  group_by(across(all_of(rep_keys))) %>%
  summarise(n_units = n_distinct(CATCH_UNIT_CODE), .groups = "drop")

table(report_units$n_units)

# Which fleets, gears and species do the double reports come from?
both <- elasmo_species %>%
  group_by(across(all_of(rep_keys))) %>%
  filter(n_distinct(CATCH_UNIT_CODE) == 2) %>%
  ungroup()

both %>% 
  count(GEAR_CODE, FLEET_CODE, SPECIES_CODE, sort = TRUE) %>% 
  print(n = 40)

# For double reports, what mean weight per fish do they imply?
both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  mutate(kg_per_fish = MT * 1000 / NO) %>%
  group_by(SPECIES_CODE) %>%
  summarise(n = n(),
            q25_kg = quantile(kg_per_fish, 0.25, na.rm = TRUE),
            median_kg = median(kg_per_fish, na.rm = TRUE),
            q75_kg = quantile(kg_per_fish, 0.75, na.rm = TRUE)) %>%
  arrange(desc(n))

# Beyond single reports: which fleet, gear and species series use each unit over time?
elasmo_species %>%
  group_by(FLEET_CODE, GEAR_CODE, SPECIES_CODE) %>%
  summarise(rows_MT = sum(CATCH_UNIT_CODE == "MT"), rows_NO = sum(CATCH_UNIT_CODE == "NO"),
            years_MT = n_distinct(YEAR[CATCH_UNIT_CODE == "MT"]),
            years_NO = n_distinct(YEAR[CATCH_UNIT_CODE == "NO"]),
            .groups = "drop") %>%
  filter(rows_MT > 0, rows_NO > 0) %>%
  arrange(desc(rows_MT + rows_NO)) %>% 
  print(n = 43)


# 1. Weight-only vs number-only vs both, across all reports
report_units2 <- elasmo_species %>%
  group_by(across(all_of(rep_keys))) %>%
  summarise(units = paste(sort(unique(CATCH_UNIT_CODE)), collapse = "+"), .groups = "drop")
table(report_units2$units)

# 2. By fleet and gear
report_units2 %>%
  count(FLEET_CODE, GEAR_CODE, units) %>%
  pivot_wider(names_from = units, values_from = n, values_fill = 0) %>%
  mutate(total = MT + NO + `MT+NO`, pct_both = round(100 * `MT+NO` / total, 1)) %>%
  arrange(desc(total)) %>% print(n = 92)

# 3. Is the weight measured, or computed from the number?
both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  mutate(kg = MT * 1000 / NO) %>%
  filter(FLEET_CODE == "TWN") %>%
  group_by(SPECIES_CODE, YEAR) %>%
  summarise(n = n(), distinct_kg = n_distinct(round(kg, 1)), sd_kg = sd(kg), .groups = "drop") %>%
  print(n = 59)


# Does kg per fish differ between fleets for the same species? (if similar, one factor is defensible)
both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  mutate(kg = MT * 1000 / NO) %>%
  filter(SPECIES_CODE %in% c("BSH", "MAK", "FAL", "SMA", "OCS")) %>%
  group_by(SPECIES_CODE, FLEET_CODE, GEAR_CODE) %>%
  summarise(n = n(), q25 = quantile(kg, 0.25, na.rm = TRUE),
            median_kg = median(kg, na.rm = TRUE),
            q75 = quantile(kg, 0.75, na.rm = TRUE), .groups = "drop") %>%
  filter(n >= 30) %>% arrange(SPECIES_CODE, desc(n)) %>% print(n = 40)

# Are some reports using fixed conversion factors? Most common kg/fish values, by fleet
both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  mutate(kg = round(MT * 1000 / NO, 1)) %>%
  count(FLEET_CODE, SPECIES_CODE, kg, sort = TRUE) %>% print(n = 20)



both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  mutate(kg = round(MT * 1000 / NO, 1),
         round5 = abs(kg / 5 - round(kg / 5)) < 1e-6) %>%
  group_by(FLEET_CODE, GEAR_CODE) %>%
  summarise(n = n(), pct_multiple_of_5 = round(100 * mean(round5, na.rm = TRUE), 1), .groups = "drop") %>%
  filter(n >= 30) %>% arrange(desc(n))

table(IOTC_catch$CATCH_UNIT_CODE, useNA = "ifany")

both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  filter(NO >= 50) %>%
  mutate(kg = round(MT * 1000 / NO, 1),
         round5 = abs(kg / 5 - round(kg / 5)) < 1e-6) %>%
  group_by(FLEET_CODE, GEAR_CODE) %>%
  summarise(n = n(), pct_multiple_of_5 = round(100 * mean(round5, na.rm = TRUE), 1), .groups = "drop") %>%
  filter(n >= 20) %>% arrange(desc(n))

elasmo_species %>% count(CATCH_UNIT_CODE, zero = CATCH == 0)

elasmo_species %>%
  filter(CATCH_UNIT_CODE == "MT") %>%
  summarise(min_t = min(CATCH),
            pct_at_0.01 = round(100 * mean(CATCH <= 0.01), 1),
            pct_under_0.05 = round(100 * mean(CATCH < 0.05), 1))


kg_table <- both %>%
  select(all_of(rep_keys), CATCH_UNIT_CODE, CATCH) %>%
  pivot_wider(names_from = CATCH_UNIT_CODE, values_from = CATCH) %>%
  filter(NO >= 50, MT > 0) %>%
  mutate(kg = MT * 1000 / NO)

tier1 <- kg_table %>%                       # species x fleet x gear
  group_by(SPECIES_CODE, FLEET_CODE, GEAR_CODE) %>%
  summarise(n1 = n(), kg1 = median(kg), .groups = "drop") %>% filter(n1 >= 30)

tier2 <- kg_table %>%                       # species x gear, fleets pooled
  group_by(SPECIES_CODE, GEAR_CODE) %>%
  summarise(n2 = n(), kg2 = median(kg), .groups = "drop") %>% filter(n2 >= 30)

tier1; tier2

no_only <- elasmo_species %>%
  group_by(across(all_of(rep_keys))) %>%
  filter(all(CATCH_UNIT_CODE == "NO")) %>%
  ungroup() %>%
  left_join(tier1, by = c("SPECIES_CODE", "FLEET_CODE", "GEAR_CODE")) %>%
  left_join(tier2, by = c("SPECIES_CODE", "GEAR_CODE")) %>%
  mutate(tier = case_when(!is.na(kg1) ~ 1, !is.na(kg2) ~ 2, TRUE ~ NA_real_))

no_only %>%
  count(FLEET_CODE, GEAR_CODE, SPECIES_CODE, tier) %>%
  arrange(desc(n)) %>% print(n = 30)

no_only %>% count(tier, useNA = "ifany")









# Adding stratum ID - a way to identify each row in each dataset
IOTC_catch_latlon <- IOTC_catch_latlon %>%
  unite(stratum_id, YEAR, QUARTER, MONTH_START, MONTH_END, FISHING_GROUND_CODE,
        FLEET_CODE, FISHERY_TYPE_CODE, FISHERY_GROUP_CODE, FISHERY_CODE, 
        GEAR_CODE, SCHOOL_TYPE_CODE, sep = "_", remove = FALSE)

IOTC_effort <- IOTC_effort %>%
  unite(stratum_id, YEAR, QUARTER, MONTH_START, MONTH_END, FISHING_GROUND_CODE,
        FLEET_CODE, FISHERY_TYPE_CODE, FISHERY_GROUP_CODE, FISHERY_CODE,
        GEAR_CODE, SCHOOL_TYPE_CODE, sep = "_", remove = FALSE) %>%
  mutate(grid_code = sprintf("%07d", as.integer(FISHING_GROUND_CODE))) %>% # Effort restricted to 1x1 and 5x5 cells
  filter(substr(grid_code, 1, 1) %in% c("5", "6"))

# **CHECK** Why does effort have several rows per stratum?
keys <- c("YEAR", "QUARTER", "MONTH_START", "MONTH_END", "FISHING_GROUND_CODE",
          "FLEET_CODE", "FISHERY_TYPE_CODE", "FISHERY_GROUP_CODE", "FISHERY_CODE",
          "GEAR_CODE", "SCHOOL_TYPE_CODE")

# A. Exact duplicates: rows with the same stratum AND same unit (expect 0)
IOTC_effort %>% 
  count(across(all_of(c(keys, "EFFORT_UNIT_CODE")))) %>% 
  filter(n > 1) %>% 
  nrow()

# B. Rows per stratum, and how many different units each stratum has
IOTC_effort %>%
  group_by(across(all_of(keys))) %>%
  summarise(n_rows = n(), n_units = n_distinct(EFFORT_UNIT_CODE), .groups = "drop") %>%
  count(n_rows, n_units)

# Result: Where there are multiple rows per stratum, it is because there are different units for the
# same stratum. The duplicated strata differ only in effort unit - there are no hidden duplicates

# C. Which units each gear uses
IOTC_effort %>% count(GEAR_CODE, EFFORT_UNIT_CODE) %>% arrange(GEAR_CODE, desc(n))

# D. Do school type codes agree between catch and effort?
IOTC_catch_latlon %>% count(GEAR_CODE, SCHOOL_TYPE_CODE)
IOTC_effort       %>% count(GEAR_CODE, SCHOOL_TYPE_CODE)


# Matching the effort to the catch
keys_noschool <- setdiff(keys, "SCHOOL_TYPE_CODE")

# One row per stratum and one column per unit in the effort dataframe
effort_wide <- IOTC_effort %>%
  select(all_of(keys), EFFORT_UNIT_CODE, EFFORT) %>%
  pivot_wider(names_from = EFFORT_UNIT_CODE, values_from = EFFORT, names_prefix = "effort_")

# In catch database, purse is never recorded as school type "UNCL" (unknown) but in effort the large
# majority are recorded as this. A join including school type would leave most purse seine without effort
# Therefore, first match includes school type and then the second match sends the leftover catch
# rows to the UNCL effort for the same stratum

effort_exact <- effort_wide %>% 
  mutate(effort_match = "exact")

effort_uncl  <- effort_wide %>% 
  filter(SCHOOL_TYPE_CODE == "UNCL") %>%
  select(-SCHOOL_TYPE_CODE) %>% 
  mutate(effort_match = "uncl_fallback")

step1 <- IOTC_catch_latlon %>%
  left_join(effort_exact, by = keys, relationship = "many-to-one")
done <- step1 %>%
  filter(!is.na(effort_match))
todo <- step1 %>%
  filter(is.na(effort_match)) %>%
  select(-starts_with("effort"))

# Leftover rows take the UNCL effort for the same stratum
step2 <- todo %>%
  left_join(effort_uncl, by = keys_noschool, relationship = "many-to-one")



# USE ONE UNIT PER SERIES - number where possible (because not rounded?)




















