# Set working directory
setwd("~/Library/CloudStorage/OneDrive-UniversityofBristol/PhD year 1/R code and data/Raw data")

# Load packages
library(dplyr)
library(tidyverse)
library(stringr)
library(rfishbase)
library(ggOceanMaps)
library(lme4) 
library(lmerTest) 
library(DHARMa) 
library(car)
library(ggeffects)
library(ggplot2)
library(sjPlot)
library(ggpubr)
library(patchwork)
library(sf)
library(ggnewscale)
library(rnaturalearth)
library(rnaturalearthdata)
library(broom)
library(glmmTMB)
library(emmeans)

#==========
# Load data
#==========

all_timeseries <- read.csv("z_scores_all_timeseries.csv")

# Adding superorder column
superorder <- all_timeseries %>%
  distinct(order) %>%
  mutate(
    superorder = case_when(
      order %in% c("Rajiformes",
                   "Myliobatiformes",
                   "Rhinopristiformes",
                   "Torpediniformes") ~ "Batoidea",
      order %in% c("Hexanchiformes",
                   "Squaliformes",
                   "Squatiniformes",
                   "Pristiophoriformes",
                   "Heterodontiformes",
                   "Orectolobiformes",
                   "Lamniformes",
                   "Carcharhiniformes",
                   "Echinorhiniformes") ~ "Selachimorpha",
      order %in% c("Chimaeriformes") ~ "Holocephalimorpha"))

all_timeseries <- left_join(all_timeseries, superorder, by = "order")

#=========================
# Add columns for analysis
#=========================

# Slope of trend per species per cell per survey added to dataframe
# Correlation coefficient between z score and year added per species per cell per survey
all_timeseries <- all_timeseries %>%
  group_by(survey, lat_cell, lon_cell, accepted_name) %>%
  mutate(slope = coef(lm(z ~ year)) [2],
         corr_coeff = cor(z, year)) %>%
  ungroup()

# Validate names
all_timeseries$accepted_name <- trimws(all_timeseries$accepted_name) # Remove leading/trailing spaces

species_data <- data.frame(accepted_name = unique(all_timeseries$accepted_name))

species_data$valid_species <- validate_names(species_data$accepted_name)

all_timeseries <- left_join(all_timeseries, species_data, by = "accepted_name")

# Adding length data from rfishbase
length_data <- fb_tbl("species") %>%
  select(Genus, Species, Length) %>%
  unite(valid_species, Genus, Species, sep = " ") %>%
  rename(max_total_length = Length)

all_timeseries <- left_join(all_timeseries, length_data, by = "valid_species")

# Species for which fishbase does not have length information
all_timeseries %>%
  filter(is.na(max_total_length) & !is.na(valid_species)) %>%
  distinct(accepted_name, valid_species)

# Adding distance to land data from dist2land in ggoceanmaps
all_timeseries %>%
  summarise(
    na_lat   = sum(is.na(lat_cell)),
    na_lon   = sum(is.na(lon_cell)),
    bad_lat  = sum(abs(lat_cell) > 90, na.rm = TRUE),
    bad_lon  = sum(abs(lon_cell) > 180, na.rm = TRUE)
  )

distance_data <- all_timeseries %>%
  distinct(lat_cell, lon_cell) %>%
  unite(latlon, lat_cell, lon_cell, sep = "_", remove = FALSE) %>%
  dist2land(lon = "lon_cell", lat = "lat_cell", binary = TRUE) %>%
  rename(on_sea = ldist) %>%
  dist2land(lon = "lon_cell", lat = "lat_cell") %>%
  select(latlon, ldist, on_sea)

all_timeseries <- all_timeseries %>%
  unite(latlon, lat_cell, lon_cell, sep = "_", remove = FALSE) %>%
  left_join(distance_data, by = "latlon") %>%
  rename(dist_to_land = ldist)

all_timeseries %>%
  summarise(count_true = sum(on_sea == TRUE),
            count_false = sum(on_sea == FALSE))

# Filtering for NAs in valid_species, max_total_length
# And FALSE for dist_to_land i.e. the point is on land
all_timeseries_filtered <- all_timeseries %>%
  filter(!is.na(valid_species)) %>%
  filter(!is.na(max_total_length)) %>%
  filter(on_sea == TRUE)


# Remove any timeseries which are shorter than 5 years
# **CHECK** is there more than one row per series per year?
all_timeseries_filtered %>%
  count(survey, lat_cell, lon_cell, accepted_name, year) %>%
  filter(n > 1)

all_timeseries_filtered <- all_timeseries_filtered %>%
  group_by(survey, lat_cell, lon_cell, accepted_name) %>%
  mutate(n_years = n_distinct(year[!is.na(value)])) %>%
  ungroup() %>%
  filter(n_years >= 5) %>%
  mutate(series_id = paste(database, survey, lat_cell, lon_cell, accepted_name, sep = "|"),
         group = ifelse(database == "sharkipedia", "sharkipedia", survey)) # Because all sharkipedia surveys are one series which causes problem for random effect



#==================================================
# Make slopes binomial (negative = 0, positive = 1)
#==================================================

all_timeseries_filtered <- all_timeseries_filtered %>%
  mutate(binomial_slope = as.integer(slope >0))

# The test subset groups the timeseries together into one row, so each series is represented by one row
# And has the information including survey, geography, species info and also the slope, correlation coefficient and binomial slope
# for those variables
test_subset <- all_timeseries_filtered %>%
  distinct(database, survey, source_observation, lat_cell, lon_cell, latlon, accepted_name, valid_species, 
           order, class, superorder, n_years, slope, binomial_slope, corr_coeff, max_total_length, dist_to_land) %>%
  mutate(group = ifelse(database == "sharkipedia", "sharkipedia", survey))
# %>%
#  filter(database != "IATTC")

# Subsetting data into superorders
shark_test_subset <- test_subset %>%
  filter(superorder == "Selachimorpha")
ray_test_subset <- test_subset %>%
  filter(superorder == "Batoidea")
chimaera_test_subset <- test_subset %>%
  filter(superorder == "Holocephalimorpha")

#===================================
# Overall abundance of biomass trend
#===================================

# One slope and uncertainty for each series
slopes <- all_timeseries_filtered %>%
  group_by(database, survey, lat_cell, lon_cell, accepted_name, order, superorder) %>% # Define a series
  group_modify(~ {
    fit <- lm(z ~ year, data = .x)                                  # Fit the slope for each series
    tibble(slope = coef(fit)["year"],                               # Extract slope and SE coefficients
                        se = summary(fit)$coefficients["year", 2])
  }) %>%
  ungroup() %>%
  filter(se > 0) %>%                                                # Remove series where the line fits perfectly otherwise get infinite weight
  mutate(series_id = row_number())                                  # Give each series their own ID number

# Intercept only model to show average overall trend (doesn't take SE into account)
averaged_slope_model <- lmer(slope ~ 1 + (1|survey), data = slopes)
summary(averaged_slope_model)

# PLOTTING
# Finding the slope for each series but also keeping the start and end year so that the slope
# is only predicted across the years that the survey actually took place
series_lines <- all_timeseries_filtered %>%
  group_by(database, survey, lat_cell, lon_cell, accepted_name, order, superorder) %>%
  filter(n() >= 5) %>%
  group_modify(~ {
    fit <- lm(z ~ year, data = .x)
    x0 <- min(.x$year)
    x1 <- max(.x$year)
    tibble(x0 = x0,
           x1 = x1,
           y0 = predict(fit, data.frame(year = x0)),
           y1 = predict(fit, data.frame(year = x1)))
  }) %>%
  ungroup()

# The average slope
estimate <- unname(fixef(averaged_slope_model)[1])
yr_mid <- mean(range(all_timeseries_filtered$year))

all_timeseries_trend_plot <- ggplot() +
  geom_segment(data = series_lines,
               aes(x = x0, xend = x1, y = y0, yend = y1),
               colour = "grey70", alpha = 0.08, linewidth = 0.2) +
  geom_abline(slope = estimate, intercept = -estimate*yr_mid, colour = "red", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_cartesian(ylim = c(-2, 2)) +
  labs(x = "Year", y = "Abundance (z-score)",
       title = "Trend of every time series, with the overall average",
       subtitle = "Grey = individual series, red = average trend") +
  theme_minimal()

all_timeseries_trend_plot

# Plotting the spatial and temporal trends
slope_summary <- slopes %>%
  group_by(lat_cell, lon_cell) %>%
  summarise(mean_slope = mean(slope), .groups = "drop")

# Creating general bathymetry map
bathymetry <- basemap(limits = c(-180, 180, -70, 80), bathymetry = TRUE, bathy.style = "rcb",
                      legends = FALSE)

# Global trends over full temporal scale
slope_map <- bathymetry +
  geom_point(data = slope_summary, aes(x = lon_cell, y = lat_cell, size = abs(mean_slope), 
             colour = ifelse(mean_slope < 0, "Negative", "Positive")), alpha = 0.5) +
  scale_size_continuous(range = c(2, 15), limits = c(0, 0.45)) +
  scale_colour_manual(
    values = c("Negative" = "red", "Positive" = "blue",
               name = "Trend direction")) +
  labs(title = "", x = "", y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, family = "Arial"),
        legend.position = "none")

slope_map


#==============
# Length models
#==============

# Plotting the correlation coefficient of z vs year against max_total_length
# Most data is coming from smaller species
ggplot() +
  geom_point(data = test_subset, aes(x = corr_coeff, y = max_total_length))

#================================
# glmer model with binomial slope
#================================

# Are larger species more or less likely to be increasing?
# Random intercept gives each species its own baseline probability
# (Several timeseries can come from the same species so they aren't independent)
length_model <- glmer(binomial_slope ~ max_total_length + (1|valid_species),
                      data = test_subset, family = binomial("logit"))
summary(length_model)
length_model_check <- simulateResiduals(length_model)
plot(length_model_check, quantreg = T) # quantreg used because dataset larger than 10,000

# Finding the 50% threshold
# The length at which the probability crosses 0.5 - the point at which species switch from
# mostly increasing to mostly declining (because estimate in summary is negative)
# The probability is 50% when the log-odds are zero - so the threshold is the negative intercept
# divided by the length coefficient
-fixef(length_model)[1] / fixef(length_model)["max_total_length"]

#Test for an interaction between length and superorder
length_super_model <- glmer(binomial_slope ~ max_total_length * superorder + (1|valid_species),
                            data=test_subset, family = binomial("logit"))
summary(length_super_model)
Anova(length_super_model)

# Sharks
shark_length_model <- glmer(binomial_slope ~ max_total_length + (1|valid_species),
                            data = shark_test_subset, family = binomial("logit"))
-fixef(shark_length_model)[1] / fixef(shark_length_model)["max_total_length"]

# Rays
ray_length_model <- glmer(binomial_slope ~ max_total_length + (1|valid_species),
                            data = ray_test_subset, family = binomial("logit"))
-fixef(ray_length_model)[1] / fixef(ray_length_model)["max_total_length"]

# Chimaera
chimaera_length_model <- glmer(binomial_slope ~ max_total_length + (1|valid_species),
                          data = chimaera_test_subset, family = binomial("logit"))
-fixef(chimaera_length_model)[1] / fixef(chimaera_length_model)["max_total_length"]


# PLOTTING
# SJ plot for binomial slope ~ length
length_plot <- plot_model(length_model, type = "pred", terms = c("max_total_length"), colors = "darkblue") +
  theme_classic() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(from = 0, to = 1800, by = 200), limits = c(0, 1800), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1), expand = c(0,0)) +
  labs(x = "Maximum total length", y = "Abundance trend", title = "All chondrichthyes") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))


length_plot

shark_length_plot <- plot_model(shark_length_model, type = "pred", terms = c("max_total_length"), colors = "darkblue") +
  theme_classic() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(from = 0, to = 1800, by = 200), limits = c(0, 1800), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1), expand = c(0,0)) +
  labs(x = "Maximum total length", y = "Abundance trend", title = "Sharks") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

shark_length_plot

ray_length_plot <- plot_model(ray_length_model, type = "pred", terms = c("max_total_length"), colors = "darkblue") +
  theme_classic() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(from = 0, to = 1800, by = 200), limits = c(0, 1800), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1), expand = c(0,0)) +
  labs(x = "Maximum total length", y = "Abundance trend", title = "Rays") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

ray_length_plot

chimaera_length_plot <- plot_model(chimaera_length_model, type = "pred", terms = c("max_total_length"), colors = "darkblue") +
  theme_classic() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(from = 0, to = 1800, by = 200), limits = c(0, 1800), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1), expand = c(0,0)) +
  labs(x = "Maximum total length", y = "Abundance trend", title = "Chimaera") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

chimaera_length_plot

length_plot + shark_length_plot + ray_length_plot + chimaera_length_plot


#=====================
# glmmTMB model with z
#=====================

year_length_model <- glmmTMB(z ~ year*max_total_length + (1|valid_species), data = all_timeseries_filtered)
summary(year_length_model)

year_length_model_lm <- lm(z ~ year*max_total_length, data = all_timeseries_filtered)
summary(year_length_model_lm)

# PLOTTING
# Estimate year and abundance slopes across maximum total body size
year_length_slopes <- emtrends(year_length_model, ~ max_total_length,
                           var = "year", at = list(max_total_length = seq(0, 1700, 10)))

# Extract estimates and standard errors
year_length_slopes_df <- as.data.frame(summary(year_length_slopes))

# Plot surface slopes with +/- 1 SE
year_length_slopes_plot <- ggplot(year_length_slopes_df,
       aes(x = max_total_length, y = year.trend)) +
  geom_ribbon(aes(ymin = year.trend - SE, ymax = year.trend + SE),
              fill = "darkblue", alpha = 0.2) +
  geom_line(colour = "darkblue", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 1700), breaks = seq(100, 1700, 200), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.055, 0.01), breaks = seq(-0.05, 0.01, 0.01)) +
  labs(x = "Maximum total length (cm)", y = "Abundance response to year (slope)",
       title = "Association between time (year) and abundance (z score) across maximum total length") +
  theme_classic()  +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        title = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

year_length_slopes_plot

#========================
# Distance to land models
#========================

# Plotting the correlation coefficient of z vs year against dist_to_land
# Most data is coming from smaller species
ggplot() +
  geom_point(data = test_subset, aes(x = corr_coeff, y = dist_to_land))

#================================
# glmer model with binomial slope
#================================

# Logistic mixed-effects model
# Are species further from land more or less likely to be increasing?
distance_model <- glmer(binomial_slope ~ dist_to_land + (1|valid_species),
                        data = test_subset, family = binomial("logit"))
summary(distance_model)
distance_model_check <- simulateResiduals(distance_model)
plot(distance_model_check, quantreg = T)

# Finding the 50% threshold where the species switch from mostly declining to mostly increasing
# Meaningless because doesn't cross 0.5 line within boundaries??
-fixef(distance_model)[1] / fixef(distance_model)["dist_to_land"]

# PLOTTING
# SJ plot for binomial slope ~ distance
distance_plot <- plot_model(distance_model, type = "pred", terms = c("dist_to_land"), colors = "darkblue") +
  theme_classic() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(from = 0, to = 2300, by = 200), limits = c(0, 2300), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1), expand = c(0,0)) +
  labs(x = "Distance to land (km)", y = "Abundance trend", title = "") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

distance_plot


#=====================
# glmmTMB model with z
#=====================

year_distance_model <- glmmTMB(z ~ year*dist_to_land + (1|valid_species), data = all_timeseries_filtered)
summary(year_distance_model)

# Standard errors, z values and p values all come out as NaN
# Can happen when the random effect variance is estimated to be near zero ( here is 9e-10)

#Linear model
year_distance_model_lm <- lm(z ~ year*dist_to_land, data = all_timeseries_filtered)

summary(year_distance_model_lm)

# PLOTTING
# Estimate year and abundance slopes across maximum total body size
year_distance_slopes <- emtrends(year_distance_model_lm, ~ dist_to_land,
                               var = "year", at = list(dist_to_land = seq(0, 2300, 10)))

# Extract estimates and standard errors
year_distance_slopes_df <- as.data.frame(summary(year_distance_slopes))

# Plot surface slopes with +/- 1 SE
year_distance_slopes_plot <- ggplot(year_distance_slopes_df,
                                  aes(x = dist_to_land, y = year.trend)) +
  geom_ribbon(aes(ymin = year.trend - SE, ymax = year.trend + SE),
              fill = "darkblue", alpha = 0.2) +
  geom_line(colour = "darkblue", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red", linewidth = 0.5) +
  scale_x_continuous(limits = c(0, 2300), breaks = seq(100, 2300, 200), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.055, 0.01), breaks = seq(-0.05, 0.01, 0.01)) +
  labs(x = "Distance to land (km)", y = "Abundance response to year (slope)",
       title = "Association between time (year) and abundance (z score) across distance to land") +
  theme_classic()  +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        title = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

year_distance_slopes_plot


#============================
# Length and distance togther
#============================

length_dist_model <- glmer(binomial_slope ~ max_total_length + dist_to_land + (1|valid_species) + (1|group),
                           data = test_subset, family = binomial)
summary(length_dist_model)
length_dist_model_check <- simulateResiduals(length_dist_model)
plot(length_dist_model_check, quantreg = T)


#========================================================
# Sliding window results - looking for an uplift in trend
#========================================================

series_id_df <- all_timeseries_filtered %>%
  mutate(series_id = paste(database, survey, lat_cell, lon_cell, accepted_name, sep = "|"),
         group = ifelse(database == "sharkipedia", "sharkipedia", survey)) # Because all sharkipedia surveys are one series which causes problem for random effect

window_length <- 20
min_years <- 5 # how many years a series needs inside a window to contribute to that window's slope

# Define the range of years in the dataset
start_year <- min(series_id_df$year, na.rm = TRUE)
end_year   <- max(series_id_df$year, na.rm = TRUE)

sliding_window_results <- data.frame()

# Loop through the dataset with a sliding window of 20 years, moving by 1 year
for (start in seq(start_year, end_year - window_length + 1)) {
  end <- start + window_length - 1
  
  data_window <- series_id_df %>%
    filter(year >= start & year <= end) %>%
    group_by(series_id, group, database) %>%
    filter(n() >= min_years) %>%
    summarise(slope = coef(lmer(z ~ year + (1|valid_species)))[2], .groups = "drop")
  
  # How many series a window needs to run - too low and the average isn't trustworthy
  # How many surveys a window needs (needs to be more than 1 for (1|group) to work) - 5 is a common minimum for a random effect
  if (nrow(data_window) < 10 || n_distinct(data_window$group) < 5) next
    
  
  model <- lmer(slope ~ 1 + (1|group), data = data_window)
  fe <- summary(model)$coefficients
  
  sliding_window_results <- rbind(sliding_window_results, data.frame(
    start_year = start,
    end_year = end,
    mid_year = start + (window_length - 1)/2,
    estimate = fe[1, "Estimate"],
    std_err = fe[1, "Std. Error"],
    n_series = nrow(data_window),
    n_groups = n_distinct(data_window$group)
  ))
}

print(sliding_window_results)


# PLOTTING
# Sliding window analysis results
sliding_window_plot <- ggplot(sliding_window_results, aes(x = mid_year, y = estimate)) +
  geom_ribbon(aes(ymin = estimate - 1.96 * std_err,
                  ymax = estimate + 1.96 * std_err),
              fill = "darkblue", alpha = 0.2) +
  geom_line(colour = "darkblue") +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.5) +
  theme_classic() +
  scale_x_continuous(breaks = seq(from = 1955, to = 2025, by = 10), expand = c(0,0)) +
  labs(x = "Window midpoint (year)", y = "Average trend (SD per year)") +
  theme(axis.text.x = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, family = "Arial"),
        axis.title.x = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

sliding_window_plot


#=======================
# Data description plots
#=======================

# Plot the density of year frequency
year_freq <- ggplot(all_timeseries_filtered, aes(x = year)) +
  geom_density(alpha = 0.7, fill = "lightblue") +
  geom_vline(aes(xintercept = median(year)), colour = "red", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  labs(x = "Year", y = "Proportion of data") +
  scale_x_continuous(breaks = seq(from = 1900, to = 2020, by = 20), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.035)) +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"),
        plot.margin = margin(3, 1, 1, 1, "cm")) 

year_freq

# Plot the density of timeseries length
ts_freq <- ggplot(all_timeseries_filtered %>%
                    distinct(survey, latlon, accepted_name, n_years),
                  aes(x = n_years)) +
  geom_density(alpha = 0.7, fill = "lightblue") +
  geom_vline(aes(xintercept = median(n_years)), colour = "red", linetype = "dashed", linewidth = 1) +
  theme_classic() +
  labs(x = "Length of timeseries (years)", y = "Proportion of timeseries") +
  scale_x_continuous(breaks = seq(from = 0, to = 120, by = 20), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.05)) +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
        axis.text.y = element_text(size = 16, family = "Arial"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"),
        plot.margin = margin(2, 1, 1, 1, "cm"))

ts_freq

# Map of data
world <- ne_countries(scale = "medium", returnclass = "sf") # World map

surveyed_cells <- all_timeseries_filtered %>%
  distinct(database, lat_cell, lon_cell)

surveyed_cells_sf <- surveyed_cells %>%
  mutate(
    xmin = lon_cell,        # Define the boundaries of each cell
    xmax = lon_cell + 1,
    ymin = lat_cell,
    ymax = lat_cell + 1
  ) %>%
  rowwise() %>%             
  mutate(geometry = list(st_polygon(list(matrix(c(
    xmin, ymin,             # Create the four corners of the square
    xmax, ymin,
    xmax, ymax,
    xmin, ymax, 
    xmin, ymin),
    ncol = 2, byrow = TRUE))))) %>%
  ungroup() %>%
  st_as_sf(crs = 4326) %>%
  select(database, lat_cell, lon_cell, geometry)

survey_coverage_map <- ggplot() +
  geom_sf(data = world, fill = "grey90", colour = "grey50", linewidth = 0.2) +
  geom_sf(data = surveyed_cells_sf, aes(fill = database), colour = "white", linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Surveyed locations (with Chondrichthyes timeseries data)",
       fill = "Database") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 20, family = "Arial"),
        axis.title.y = element_text(size = 20, family = "Arial"),
        plot.margin = margin(2, 1, 1, 1, "cm"))

survey_coverage_map


# Data distribution plot - number of time series in each cell
latlon_counts <- all_timeseries_filtered %>%
  group_by(lat_cell, lon_cell, latlon) %>%
  summarise(count = n_distinct(survey, accepted_name), .groups = "drop")

data_map <- basemap(limits = c(-180, 180, -70, 80), bathymetry = TRUE, bathy.style = "rcb",
                    legends = FALSE) +
  geom_point(data = latlon_counts, aes(x = lon_cell, y = lat_cell, size = count), 
             colour = "darkblue", alpha = 0.4) +
  scale_size_continuous(range = c(1, 10)) +
  labs(size = "Time series") +
  theme_classic() +
  theme(legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"))

data_map


# Bar plot of timseseries per order
order_timeseries_count <- all_timeseries_filtered %>%
  distinct(survey, accepted_name, order, superorder) %>%
  count(order, superorder) %>%
  filter(!is.na(superorder))

order_barplot <- ggplot(data = order_timeseries_count, aes(x = n, y = reorder(order, n), fill = superorder)) +
  geom_col(width = 0.6) +
  theme_classic() +
  scale_fill_manual(values = c("Selachimorpha" = "darkblue", "Batoidea" = "turquoise", "Holocephalimorpha" = "violet")) +
  labs(x = "Number of timeseries", y = "Taxonomic Order", fill = "Superorder") +
  scale_x_continuous(limits = c(0, 800), breaks = seq(from = 0, to = 800, by = 200), expand = c(0, 0)) +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        legend.text = element_text(size = 12, family = "Arial"),
        legend.title = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

order_barplot


























