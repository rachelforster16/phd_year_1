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

# Loading timeseries data which calculated z scores
# NO filtering has happened on this data e.g. length of timeseries or available species data etc.
maps_data <- read.csv("z_scores_all_timeseries.csv")

#============================
# Map with surveyed locations
#============================

# World map
world <- ne_countries(
  scale = "medium",
  returnclass = "sf")

# Creating a dataframe with each lat cell and lon cell combination which has been surveyed
surveyed_cells_all <- maps_data %>%
  distinct(database, lat_cell, lon_cell)

# Turning this into geographic polygons
surveyed_cells_sf_all <- surveyed_cells_all %>%
  mutate(
    xmin = lon_cell,        # Define the boundaries of each cell
    xmax = lon_cell + 1,
    ymin = lat_cell,
    ymax = lat_cell + 1
  ) %>%
  rowwise() %>%       # Do this for each cell
  mutate(geometry = list(st_polygon(list(matrix(c(
    xmin, ymin,       # Create the four corners of the square
    xmax, ymin,
    xmax, ymax,
    xmin, ymax, 
    xmin, ymin),
    ncol = 2, byrow = TRUE))))) %>%
  ungroup() %>%
  st_as_sf(crs = 4326) %>%
  select(database, lat_cell, lon_cell, geometry)


survey_coverage_map_all <- ggplot() +
  geom_sf(
    data = world,
    fill = "grey90",
    colour = "grey50",
    linewidth = 0.2
  ) +
  geom_sf(
    data = surveyed_cells_sf_all,
    aes(fill = database),
    colour = "white",
    linewidth = 0.1
  ) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(
    title = "Surveyed locations (with Chondrichthyes timeseries data)",
    fill = "Database"
  )

survey_coverage_map_all











