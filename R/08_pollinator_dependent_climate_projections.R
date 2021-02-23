# read in packages
library(dplyr)
library(raster)
library(rworldmap) 
library(rworldxtra)
library(ggplot2)
library(cowplot)

# source in additional functions
source("R/00_functions.R")

# read in the raster for the historical data to get the baseline
## standardised climate anomaly script
# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# take names of values for 1901 to 1931 - 30 year baseline
tmp1901_1931 <- tmp[[names(tmp)[1:361]]]

# calculate the mean and sd of the baseline values
tmp1901_1931mean <- calc(tmp1901_1931, mean)
tmp1901_1931sd <- calc(tmp1901_1931, stats::sd)

## read in the rasters for the future data, start with SSP585
# set up the starting directory
SSP_directory <- ("G:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies")

# set up historical anomaly to be added on
months.1979.2013 <- 937:1356

# calculate the average temperature for 1979-2013 onto which anomaly is added
hist.mean.temp.1979.2013 <- stack(stackApply(x = tmp[[months.1979.2013]],
                                             indices = (rep(1:35,each=12)),fun = mean))
hist.mean.temp.1979.2013 <- stackApply(x = hist.mean.temp.1979.2013,indices = rep(1,35),
                                       fun = mean)

# list the files in that directory
SSP_folders <- list.files(SSP_directory)

# select rcp26 for MIROC5
SSP_folders <- SSP_folders[grepl("MIROC5_rcp26", SSP_folders)]

# setp up empty vector for file paths
SSP_file_path <- c()

# iterate through each of the folders
for(i in 1:length(SSP_folders)){
  SSP_file_path[i] <- paste(SSP_directory, SSP_folders[i], sep = "/")
}

# future projection list
future_projection <- list()
future_projection_anomaly <- list()

# for each file path, read in the climate model
for(i in 1:length(SSP_file_path)){
  future_projection[[i]] <- raster(SSP_file_path[i]) / 10
}

# set up list for rolling average
rolling_average <- list()

# calculate an average for each 3 years
step <- 1
for(i in 1:12){
  rolling_average[[i]] <- stack(future_projection[step:(step+2)])
  rolling_average[[i]] <- calc(rolling_average[[i]], mean) + hist.mean.temp.1979.2013
  print(i)
  step <- step + 3
}

# set up climate anomaly lists
climate_anomaly <- list()
std_climate_anomaly <- list()

# calculate the climate anomaly
for(i in 1:length(rolling_average)){
  climate_anomaly[[i]] <- (rolling_average[[i]] - tmp1901_1931mean)
  std_climate_anomaly[[i]] <- climate_anomaly[[i]] / tmp1901_1931sd
}

map_anom <- function(std_climate_anom){
  
  # reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
  tmp2004_6std_climate_anomaly <- projectRaster(std_climate_anom, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  # convert the raster to a format amenable to ggplot
  # convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
  anom_spdf <- as(tmp2004_6std_climate_anomaly, "SpatialPixelsDataFrame")
  anom_df <- as.data.frame(anom_spdf)
  colnames(anom_df) <- c("value", "x", "y")
  
  # group the categories of climate anomaly into factors
  anom_df$value_group[anom_df$value >= 2] <- ">= 4"
  anom_df$value_group[anom_df$value > 1 & anom_df$value <= 2] <- "1 - 2"
  anom_df$value_group[anom_df$value > 0.5 & anom_df$value <= 1] <- "0.5 - 1"
  anom_df$value_group[anom_df$value > 0.25 & anom_df$value <= 0.5] <- "0.25 - 0.5"
  anom_df$value_group[anom_df$value >= 0 & anom_df$value <= 0.25] <- "0 - 0.25"
  anom_df$value_group[anom_df$value < 0] <- "< 0"
  
  # order the levels of those factors
  anom_df$value_group <- factor(anom_df$value_group, levels = c(">= 4", "1 - 2", "0.5 - 1", "0.25 - 0.5", "0 - 0.25", "< 0"))
  
  # bring in basemap for climate site plot 
  base_map <- get_basemap()
  
  base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  # fortify the main map
  map_fort <- fortify(base_map)
  
  # plot the ggplot map for climate anomaly
  anom <- anom_df %>%
    ggplot() +
    geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "grey", alpha = 0.3) +
    
    geom_tile(aes(x = x, y = y, fill = value_group)) +
    scale_fill_manual("Standardised climate anomaly", values = c("#000000", "grey", "grey", "grey", "grey", "grey")) +
    coord_equal() +
    theme(panel.background = element_blank(),
          panel.bord = element_rect(fill = NA),
          panel.grid = element_blank(), 
          axis.text = element_blank(),
          axis.ticks = element_blank(), 
          axis.title = element_blank())
  
  return(anom)
  
}

anom_maps <- lapply(std_climate_anomaly, map_anom)

plot_grid(anom_maps[[1]], anom_maps[[12]])

