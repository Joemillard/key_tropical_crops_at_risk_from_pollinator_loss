# bind the pesticide raster data to the predicts climate data
library(dplyr)
library(raster)
library(rworldmap)
library(rworldxtra)
library(ggplot2)

# source in additional functions
source("R/00_functions.R")

# read in the predict climate data
predicts_pollinators <- readRDS("PREDICT_pollinators_5_climate.rds")

# read in the pesticide data
pesticide_low <- raster("outputs/Pesticide_totalAPR_Low_cropped.tif")
pesticide_high <- raster("outputs/Pesticide_totalAPR_High_cropped.tif")

# convert the predicts coordinates to spatial points
predicts_points <- predicts_pollinators %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Latitude)) %>%
  SpatialPoints()

# use the extract function to extract the pesticide values from both the low and hihg estimate
pesticide_low_points <- extract(pesticide_low, predicts_points)
pesticide_high_points <- extract(pesticide_high, predicts_points)

# build data frame of predicts climate data with the fertiliser data
climate_pest_predicts <- data.frame(predicts_pollinators, "low_estimate" = pesticide_low_points, "high_estimate" = pesticide_high_points) %>%
  dplyr::select(-Latitude.x, -Latitude.y, -Longitude.x, -Longitude.y)

# bring in basemap for climate site plot 
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# plot map for pesticide application
climate_pest_predicts %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = x, y = y, colour = high_estimate), alpha = 0.5) + 
  facet_wrap(~Order) +
  #scale_colour_manual("Standardised climate anomaly", values = c("#000000", "darkred", "#D55E00", "#E69F00", "#F0E442", "#56B4E9")) +
  coord_map(projection = "mollweide") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "bottom")
