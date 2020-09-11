#### Pesticide data processing
# Original author: Charlotte L. Outhwaite

# This script pulls together the PEST CHEM-GRID data into one raster of
# total pesticide application per grid cell.
# NOTE - use of raster in this script fills up the hard drive

# The data:
# 400 maps for 2015 for top 20 active ingredients in pesticides
# Resolution = 5 arc-minute (10km at the equator)

# files are per crop, per active ingredient, per year, for a high or low estimate
# load libraries
library(raster)
library(ggplot2)
library(viridis)
library(dplyr)

# set up the load directory
data_dir <- "data/ferman-v1-pest-chemgrids_geotiff/ApplicationRate/GEOTIFF"
out_dir <- "outputs"

# list the files
files <- list.files(data_dir) # 1200 files

# aggregate the high and low files for 2015 and 2025 with a for loop
pest_dates <- c("2015", "2025")

for(i in 1:length(pest_dates)){
  
  print(i)
  
  # just select the 2015 maps
  files <- list.files(data_dir, pattern = pest_dates[i]) # 400 files, high and low estimates
  
  print(files)
  
  # summarise the data across crops and active ingredients
  # separate out the low and high estimates
  files_H <- files[grep(files, pattern = "5_H")] # 200 files

  # stack the raster files and summarise the information
  pest_H <- stack(paste0(data_dir, "/", files_H))
  
  # convert all NA values to NA
  pest_H <- reclassify(pest_H, matrix(data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)
  
  print("pest_H")
  
  # get the sum of app rates across the crops/pesticides
  pest_H_total <- calc(x = pest_H, fun = sum, na.rm = TRUE)
  
  # remove large raster stack
  rm(pest_H)
  
  # save the raster of low and high estimate totals
  writeRaster(pest_H_total, filename = paste0(out_dir, "/", pest_dates[i], "Pesticide_totalAPR_High.tif"))
  
  # remove large raster stack
  rm(pest_H_total)
  
}

# plot of the pesticide datasets for 2015 and 2025, just for the high application
# set up empty list, one element for 2015 and one for 2025
pest_H_total_2015 <- raster(paste0(out_dir, "/", "2015Pesticide_totalAPR_High.tif"))
pest_H_total_2025 <- raster(paste0(out_dir, "/", "2025Pesticide_totalAPR_High.tif"))
  
# mask out the sea area that is set to 0
data("wrld_simpl", package = 'maptools')
pest_H_total_2015 <- mask(pest_H_total_2015$X2015Pesticide_totalAPR_High, wrld_simpl)
pest_H_total_2015 <- crop(pest_H_total_2015, wrld_simpl)
pest_H_total_2025 <- mask(pest_H_total_2025$X2025Pesticide_totalAPR_High, wrld_simpl)
pest_H_total_2025 <- crop(pest_H_total_2025, wrld_simpl)
  
# resave the 2015 and 2015 rasters for pesticide application before projecting into mollweide projection
writeRaster(pest_H_total_2015, filename = paste0(out_dir, "/", "2015_Pesticide_totalAPR_High_cropped.tif"))
writeRaster(pest_H_total_2025, filename = paste0(out_dir, "/", "2025_Pesticide_totalAPR_High_cropped.tif"))

# read the written rasters back in again
pest_H_total_2015 <- raster(paste0(out_dir, "/", "2015_Pesticide_totalAPR_High_cropped.tif"))
pest_H_total_2025 <- raster(paste0(out_dir, "/", "2025_Pesticide_totalAPR_High_cropped.tif"))

# raster change
pesticide_change <- pest_H_total_2025 - pest_H_total_2015

# reproject rasters as mollweide projection
pest_H_total_2015 <- projectRaster(pest_H_total_2015, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
pest_H_total_2025 <- projectRaster(pest_H_total_2025, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
# convert for use in ggplot
pest_H_total_2015 <- as.data.frame(pest_H_total_2015, xy = TRUE)
pest_H_total_2025 <- as.data.frame(pest_H_total_2025, xy = TRUE)
  
# remove NAs
pest_H_total_2015 <- pest_H_total_2015[!is.na(pest_H_total_2015$X2015_Pesticide_totalAPR_High_cropped), ]
pest_H_total_2025 <- pest_H_total_2025[!is.na(pest_H_total_2025$X2025_Pesticide_totalAPR_High_cropped), ]
  
# assign column for high or low estimate
pest_H_total_2015$est <- "2015"
pest_H_total_2025$est <- "2025"
  
# set APR column name
names(pest_H_total_2015)[3] <- "APR"
names(pest_H_total_2025)[3] <- "APR"
  
# bind together the high and low estimate pesticide application dataa, and then sort the factors as high and low for plot
pest_data <- rbind(pest_H_total_2015, pest_H_total_2025)
pest_data$est <- factor(pest_data$est, levels = c("2015", "2025"))

# set breaks for a facetted plot of high and low pesticide application
brk <- c(0, 10, 100, 300)

# plot of high and low pesticide application
pest_data %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = APR)) +
  facet_grid(est~., switch = "y") +
  scale_fill_gradientn(name = "Total application rate\n(kg/ha)", breaks = brk, trans = "log1p", colours = viridis(10), labels = brk) + 
  theme_bw() +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text = element_blank(),
        strip.text.y.left = element_text(size = 11, angle = 45)) 

# save the map of pesticides for 2015 and 2025
ggsave("pesticide_data_maps.png", dpi = 300, scale = 1)

# calculate change in pesticide application for 2015 and 2025