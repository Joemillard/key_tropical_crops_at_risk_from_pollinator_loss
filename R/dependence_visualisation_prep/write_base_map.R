# script to save basemap to be read in for shiny app

library(rworldmap)
library(rworldxtra)
library(ggplot2)

# download full basemap
base_map <- getMap(resolution = "high")

# convert to correction projection
proj4string(base_map) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# fortify the main map
map_fort <- fortify(base_map)

saveRDS(map_fort, "plot_base_map.rds")
