# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data")

# Load packages
library(dplyr)
library(readr)

# Load timeseries data
fishglob_ts <- read.csv("Fishglob_timeseries.csv") %>%
  mutate(database = "FISHGLOB")
biotime_ts <- read.csv("Biotime_timeseries.csv") %>%
  mutate(database = "BioTIME")
IMOS_ts <- read.csv("IMOS_timeseries.csv") %>%
  mutate(database = "IMOS")
timefish_ts <- read.csv("TimeFISH_timeseries.csv") %>%
  mutate(database = "TimeFISH")
soviet_trawl_ts <- read.csv("SovietTrawl_timeseries.csv") %>%
  mutate(database = "SovietTrawl")
NSW_ts <- read.csv("NSW_timeseries.csv") %>%
  mutate(database = "NSW")
MEDITS_nw_ts <- read.csv("MEDITS_nw_timeseries.csv") %>%
  mutate(database = "MEDITS_nw")
MEDITS_spain_ts <- read.csv("MEDITS_spain_timeseries.csv") %>%
  mutate(database = "MEDITS_spain")
MEDITS_blacksea_ts <- read.csv("MEDITS_blacksea_timeseries.csv") %>%
  mutate(database = "MEDITS_blacksea")
IATTC_ts <- read.csv("IATTC_timeseries.csv") %>%
  mutate(database = "IATTC")
sharkipedia_ts <- read.csv("sharkipedia_timeseries.csv") %>%
  mutate(database = "sharkipedia")


# Bind rows
all_timeseries <- bind_rows(fishglob_ts, biotime_ts, IMOS_ts, timefish_ts, 
                            soviet_trawl_ts, NSW_ts, MEDITS_nw_ts, MEDITS_spain_ts,
                            MEDITS_blacksea_ts, IATTC_ts, sharkipedia_ts) %>%
  select(database, survey, year, lat_cell, lon_cell, accepted_name, order, class, mean_num_cpue,
         mean_num_cpua, mean_wgt_cpue, mean_wgt_cpua, other, source_observation)


# Add fisheries independent/fisheries information to source_observation column
all_timeseries <- all_timeseries %>%
  mutate(source_observation = case_when(
    database %in% c("FISHGLOB", "BioTIME", "IMOS", "TimeFISH", "SovietTrawl",
                    "NSW", "MEDITS_nw", "MEDITS_spain", "MEDITS_blacksea") & is.na(source_observation) ~ "independent",
    database %in% c("IATTC") & is.na(source_observation) ~ "fisheries",
    database == "sharkipedia" & source_observation == "fishery_independent" ~ "independent",
    database == "sharkipedia" & source_observation %in% c("stock_assessment", "fisheries") ~ "fisheries",
    database == "sharkipedia" & survey %in% c("sharkipedia_4012_everett2015",
                                              "sharkipedia_4096_macena2016",
                                              "sharkipedia_3551_dudley2006b",
                                              "sharkipedia_3552_dudley2006b",
                                              "sharkipedia_3555_dudley2006b",
                                              "sharkipedia_3575_dudley2006b",
                                              "sharkipedia_3584_dudley2006b",
                                              "sharkipedia_3613_dudley2006b",
                                              "sharkipedia_3626_dudley2006b",
                                              "sharkipedia_3644_dudley2006b",
                                              "sharkipedia_3648_dudley2006b",
                                              "sharkipedia_3766_dudley2006b",
                                              "sharkipedia_3851_dudley2006b",
                                              "sharkipedia_4145_dudley2006b",
                                              "sharkipedia_4159_dudley2006b",
                                              "sharkipedia_4173_dudley2006b",
                                              "sharkipedia_3642_Roff2018",
                                              "sharkipedia_3643_Roff2018",
                                              "sharkipedia_3646_Roff2018",
                                              "sharkipedia_3647_Roff2018",
                                              "sharkipedia_3764_Roff2018",
                                              "sharkipedia_3765_Roff2018",
                                              "sharkipedia_4153_Roff2018",
                                              "sharkipedia_4154_Roff2018",
                                              "sharkipedia_3645_reid2011",
                                              "sharkipedia_3649_reid2011",
                                              "sharkipedia_3767_reid2011",
                                              "sharkipedia_3825_reid2011",
                                              "sharkipedia_3972_reid2011",
                                              "sharkipedia_4164_reid2011",
                                              "sharkipedia_4174_reid2011",
                                              "sharkipedia_4240_reid2011",
                                              "sharkipedia_3758_holmes2012",
                                              "sharkipedia_3759_holmes2012",
                                              "sharkipedia_3760_holmes2012",
                                              "sharkipedia_3761_holmes2012",
                                              "sharkipedia_3774_simpfendorfer2011b",
                                              "sharkipedia_4155_simpfendorfer2011b",
                                              "sharkipedia_3961_carlson2017",
                                              "sharkipedia_4011_carlson2012b",
                                              "sharkipedia_4094_bradshaw2008", 
                                              "sharkipedia_4095_bradshaw2008",
                                              "sharkipedia_4144_noriega2011",
                                              "sharkipedia_4177_dicken2018") ~ "independent",
    database == "sharkipedia" & survey %in% c("sharkipedia_3453_baum2010",
                                              "sharkipedia_3592_baum2010",
                                              "sharkipedia_3757_baum2010",
                                              "sharkipedia_3850_baum2010",
                                              "sharkipedia_3976_baum2010",
                                              "sharkipedia_4156_baum2010",
                                              "sharkipedia_3652_curtis2014",
                                              "sharkipedia_3860_coelho2017a",
                                              "sharkipedia_3993_coelho2017a",
                                              "sharkipedia_3874_cortes2002b",
                                              "sharkipedia_3994_cortes2002b",
                                              "sharkipedia_4185_rago1998",
                                              "sharkipedia_4297_Latour2020",
                                              "sharkipedia_3496_wueringer2017",
                                              "sharkipedia_3495_wueringer2017",
                                              "sharkipedia_3617_walter2010",
                                              "sharkipedia_3632_walter2010",
                                              "sharkipedia_3653_curtis2014",
                                              "sharkipedia_4312_curtis2014",
                                              "sharkipedia_3868_skomal2005a",
                                              "sharkipedia_3869_skomal2005a",
                                              "sharkipedia_3870_skomal2005a",
                                              "sharkipedia_4003_skomal2005a",
                                              "sharkipedia_4004_skomal2005a",
                                              "sharkipedia_3984_campana2006b") ~ "fisheries"))


# Check that each row only has one metric column populated
value_cols <- c("mean_num_cpue", "mean_num_cpua", "mean_wgt_cpue", "mean_wgt_cpua", "other")

all_timeseries <- all_timeseries %>%
  mutate(n_filled = rowSums(!is.na(across(all_of(value_cols)))))

table(all_timeseries$n_filled)

all_timeseries <- all_timeseries %>%
  select(-n_filled)

# Collapse the value columns into one
all_timeseries <- all_timeseries %>%
  mutate(
    value  = coalesce(mean_num_cpue, mean_num_cpua, mean_wgt_cpue, mean_wgt_cpua, other),
    metric = case_when(
      !is.na(mean_num_cpue) ~ "mean_num_cpue",
      !is.na(mean_num_cpua) ~ "mean_num_cpua",
      !is.na(mean_wgt_cpue) ~ "mean_wgt_cpue",
      !is.na(mean_wgt_cpua) ~ "mean_wgt_cpua",
      !is.na(other)    ~ "other",
      TRUE             ~ NA_character_
    )
  )

# Check each survey uses a single metric
metric_check <- all_timeseries %>%
  filter(!is.na(metric)) %>%
  group_by(survey) %>%
  summarise(n_metrics = n_distinct(metric), .groups = "drop")

if (any(metric_check$n_metrics > 1)) {
  warning("Some surveys mix metrics: ",
          paste(metric_check$survey[metric_check$n_metrics > 1], collapse = ", "))
}

# Calculate z scores within each survey at species level 
# z = (x - mean) / sd, using the sample SD (n - 1), per survey.
all_timeseries_z <- all_timeseries %>%
  group_by(database, survey, lat_cell, lon_cell, accepted_name) %>%
  mutate(
    series_mean = mean(value, na.rm = TRUE),
    series_sd   = sd(value, na.rm = TRUE),
    z           = if_else(series_sd > 0,
                          (value - series_mean) / series_sd,
                          NA_real_)   # constant series -> NA, not NaN/Inf
  ) %>%
  ungroup() %>%
  select(-mean_num_cpue, -mean_num_cpua, -mean_wgt_cpue, -mean_wgt_cpua, -other)

# Determine the reasoning for NAs in Z score column
# Could be due to NA values, only series legnth of 1, or all the values are the same across the series
na_diagnosis <- all_timeseries_z %>%
  group_by(database, survey, lat_cell, lon_cell, accepted_name) %>%
  mutate(
    n_obs   = sum(!is.na(value)),
    n_years = n_distinct(year[!is.na(value)])
  ) %>%
  ungroup() %>%
  mutate(
    na_reason = case_when(
      !is.na(z)                  ~ "z calculated",
      is.na(value)               ~ "value is NA",
      n_obs == 1                 ~ "series length 1",
      series_sd == 0             ~ "constant series (sd = 0)",
      TRUE                       ~ "other"
    )
  )

# Overall breakdown
table(na_diagnosis$na_reason)

# Breakdown by database
na_diagnosis %>%
  filter(is.na(z)) %>%
  count(database, na_reason) %>%
  arrange(database, desc(n))

# Series length distribution of the NA rows
na_diagnosis %>%
  filter(is.na(z)) %>%
  count(n_obs) %>%
  arrange(n_obs)

# Remove the NA values if safe to do so
all_timeseries_z <- all_timeseries_z %>%
  filter(!is.na(z))

# **CHECK**
check <- all_timeseries_z %>%
  group_by(database, survey, lat_cell, lon_cell, accepted_name) %>%
  summarise(z_mean = mean(z, na.rm = TRUE),
            z_sd   = sd(z, na.rm = TRUE),
            .groups = "drop")

summary(check$z_mean)   # ~ 0 everywhere
summary(check$z_sd)     # 1 everywhere (NA for constant series)
sum(is.na(all_timeseries_z$z))   # rows with no z (constant series)

# Save CSV
write_csv(all_timeseries_z, "z_scores_all_timeseries.csv")





