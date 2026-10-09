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

# Loading timeseries
fishglob_ts <- read.csv("Fishglob_timeseries.csv") %>%
  mutate(database = "FishGLOB")
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
MEDITS_ts <- read.csv("MEDITS_timeseries.csv") %>%
  mutate(database = "MEDITS")

#============================
# Map with surveyed locations
#============================

# World map
world <- ne_countries(
  scale = "medium",
  returnclass = "sf")

# Creating a dataframe with each lat cell and lon cell combination which has been surveyed
surveyed_cells <- bind_rows(fishglob_ts, biotime_ts, IMOS_ts, timefish_ts, soviet_trawl_ts, NSW_ts,
                            MEDITS_ts) %>%
  distinct(survey, lat_cell, lon_cell)

# Turning this into geographic polygons
surveyed_cells_sf <- surveyed_cells %>%
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
  select(survey, lat_cell, lon_cell, geometry)


survey_coverage_map <- ggplot() +
  geom_sf(
    data = world,
    fill = "grey90",
    colour = "grey50",
    linewidth = 0.2
  ) +
  geom_sf(
    data = surveyed_cells_sf,
    aes(fill = survey),
    colour = "white",
    linewidth = 0.1
  ) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Surveyed locations (with Chondrichthyes timeseries data)"
  )

survey_coverage_map

