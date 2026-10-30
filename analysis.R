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
all_timeseries <- all_timeseries %>%
  group_by(survey, lat_cell, lon_cell, accepted_name) %>%
  mutate(slope = coef(lm(z ~ year)) [2]) %>%
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
distance_data <- all_timeseries %>%
  distinct(lat_cell, lon_cell) %>%
  unite(latlon, lat_cell, lon_cell, sep = "_", remove = FALSE) %>%
  dist2land(lon = "lon_cell", lat = "lat_cell", binary = TRUE) %>%
  filter(ldist == TRUE) %>%                               # Removing any latlon on land
  select(-ldist) %>%
  dist2land(lon = "lon_cell", lat = "lat_cell") %>%
  select(latlon, ldist)
  

all_timeseries <- all_timeseries %>%
  unite(latlon, lat_cell, lon_cell, sep = "_", remove = FALSE) %>%
  left_join(distance_data, by = "latlon") %>%
  rename(dist_to_land = ldist)

# Filtering for NAs in valid_species, max_total_length, or dist_to_land
all_timeseries <- all_timeseries %>%
  filter(!is.na(valid_species)) %>%
  filter(!is.na(max_total_length)) %>%
  filter(!is.na(dist_to_land))


#==================================================
# Make slopes binomial (negative = 0, positive = 1)
#==================================================

all_timeseries <- all_timeseries %>%
  mutate(binomial_slope = as.integer(slope >0))

test_subset <- all_timeseries %>%
  distinct(database, survey, lat_cell, lon_cell, latlon, accepted_name, valid_species, 
           order, class, superorder, n_years, slope, binomial_slope, max_total_length, dist_to_land)
# %>%
#  filter(database != "IATTC")


#==============
# Length models
#==============

# Main model (logistic mixed-effects model)
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


#========================
# Distance to land models
#========================

hist(test_subset$slope)

distance_model <- glmer(binomial_slope ~ dist_to_land + (1|valid_species),
                        data = test_subset, family = binomial("logit"))
summary(distance_model)
distance_model_check <- simulateResiduals(distance_model)
plot(distance_model_check, quantreg = T)

# Finding the 50% threshold where the species switch from mostly declining to mostly increasing
-fixef(distance_model)[1] / fixef(distance_model)["dist_to_land"]

# Not binomial model?
distance_model_slope <- lmer(slope ~ dist_to_land + (1|valid_species),
                           data = test_subset)
summary(distance_model_slope)

#========================================================
# Sliding window results - looking for an uplift in trend
#========================================================

# Create a data frame to store results
sliding_window_results <- data.frame(
  start_year = integer(),
  end_year = integer(),
  slope = numeric(),
  p_value = numeric(),
  r_squared = numeric(),
  std_error = numeric(),
  stringsAsFactors = FALSE
)

# Define the range of years in the dataset
start_year <- min(all_timeseries$year, na.rm = TRUE)
end_year   <- max(all_timeseries$year, na.rm = TRUE)

# Loop through the dataset with a sliding window of 20 years, moving by 1 year
for (start in seq(start_year, end_year - 19)) {
  end <- start + 19
  
  data_window <- all_timeseries %>%
    filter(year >= start & year <= end)
  
  model <- lm(z ~ year, data = data_window)
  
  slope <- unname(coef(model)["year"])
  p_value <- summary(model)$coefficients[2,4]
  r_squared <- summary(model)$r.squared
  std_error <- summary(model)$coefficients[2,2]
  
  sliding_window_results <- rbind(sliding_window_results, data.frame(
    start_year = start,
    end_year = end,
    slope = slope,
    p_value = p_value,
    r_squared = r_squared,
    std_error = std_error
  ))
}

#========================================
# Abundance and abundance trend over time
#========================================

zscore_model <- lmer(z ~ year + (1|valid_species), data = all_timeseries)
summary(zscore_model)

slope_year_model <- glmer(binomial_slope ~ year + (1|valid_species),
                          data = all_timeseries, family = binomial("logit"))
summary(slope_year_model)


#==================
# Plots and figures
#==================

# Plot the density of year frequency
year_freq <- ggplot(all_timeseries, aes(x = year)) +
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
ts_freq <- ggplot(all_timeseries %>%
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

surveyed_cells <- all_timeseries %>%
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
latlon_counts <- all_timeseries %>%
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
order_timeseries_count <- all_timeseries %>%
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


# SJ plot for binomial slope ~ length
length_plot <- plot_model(length_model, type = "pred", terms = c("max_total_length"), colors = "darkblue") +
  theme_classic() +
  ylim(0, 1) +
  geom_hline(yintercept = 0.5, colour = "red", linewidth = 0.5) +
  scale_x_continuous(breaks = seq(from = 0, to = 1800, by = 200), limits = c(0, 1800), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(from = 0, to = 1, by = 0.25), limits = c(0, 1), expand = c(0,0)) +
  labs(x = "Maximum total length", y = "Abundance trend", title = "") +
  theme(axis.text.x = element_text(size = 12, family = "Arial"),
        axis.text.y = element_text(size = 12, family = "Arial"),
        axis.title.x = element_text(size = 14, family = "Arial"),
        axis.title.y = element_text(size = 14, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))
  

length_plot


#SJ plot for binomial slope ~ distance
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


# Sliding window analysis results
sliding_window_plot <- ggplot(sliding_window_results, aes(x = start_year, y = slope)) +
  geom_ribbon(aes(ymin = (slope - std_error), ymax = (slope + std_error)), fill = "darkblue", alpha = 0.2) +
  geom_line(colour = "darkblue") +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.5) +
  theme_classic() +
  labs(x = "Year",
       y = "Abundance trend")+
  theme(axis.text.x = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, family = "Arial"),
        axis.title.x = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))


sliding_window_plot

# Since 1950
sliding_plot_1950 <- ggplot(sliding_window_results %>% filter(start_year >= 1950),
                            aes(x = start_year, y = slope)) +
  geom_ribbon(aes(ymin = (slope - std_error), ymax = (slope + std_error)), fill = "darkblue", alpha = 0.2) +
  geom_line(colour = "darkblue") +
  geom_hline(yintercept = 0, colour = "red", linewidth = 0.5) +
  theme_classic() +
  labs(x = "Year",
       y = "Abundance trend")+
  theme(axis.text.x = element_text(size = 14, family = "Arial"),
        axis.text.y = element_text(size = 14, family = "Arial"),
        axis.title.x = element_text(size = 18, family = "Arial"),
        axis.title.y = element_text(size = 18, family = "Arial"),
        plot.margin = margin(1, 1, 1, 1, "cm"))

sliding_plot_1950






