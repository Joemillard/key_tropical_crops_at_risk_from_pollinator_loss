#### Pesticide data processing
# Original author: Charlotte L. Outhwaite

# This script pulls together the PEST CHEM-GRID data into one raster of
# total pesticide application per grid cell.

# The data:
# 400 maps for 2015 for top 20 active ingredients in pesticides
# Resolution = 5 arc-minute (10km at the equator)

# files are per crop, per active ingredient, per year, for a high or low estimate
# load libraries
library(raster)
library(ggplot2)
#library(maps)
library(viridis)
library(dplyr)

# set up the load directory
data_dir <- "data/ferman-v1-pest-chemgrids_geotiff/ApplicationRate/GEOTIFF"
out_dir <- "outputs"

# list the files
files <- list.files(data_dir) # 1200 files

# just select the 2015 maps
files <- list.files(data_dir, pattern = "2015") # 400 files, high and low estimates

# summarise the data across crops and active ingredients
# separate out the low and high estimates
files_H <- files[grep(files, pattern = "5_H")] # 200 files
files_L <- files[grep(files, pattern = "5_L")] # 200 files

# stack the raster files and summarise the information
pest_H <- stack(paste0(data_dir, "/", files_H))

# convert all NA values to NA
pest_H <- reclassify(pest_H, matrix(data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# get the sum of app rates across the crops/pesticides
pest_H_total <- calc(x = pest_H, fun = sum, na.rm = TRUE)

# take a look at the total application rates
plot(pest_H_total$layer)

# same for the low estimates
pest_L <- stack(paste0(data_dir, "/", files_L))

#  convert all negative values to NA
pest_L <- reclassify(pest_L, matrix(data = c(-2, -1.5, -1, NA, NA, NA), byrow = F, ncol = 2), right = F)

# get the sum of app rates across the crops/pesticides
pest_L_total <- calc(x = pest_L, fun = sum, na.rm = TRUE)

# take a look at the total application rates
plot(pest_L_total$layer)

# save the raster of low and high estimate totals
writeRaster(pest_H_total, filename = paste0(outdir, "/Pesticide_totalAPR_High.tif"))
writeRaster(pest_L_total, filename = paste0(outdir, "/Pesticide_totalAPR_Low.tif"))

# plot of the pesticide datasets
pest_L_total <- raster("data/Pesticide_totalAPR_Low.tif")
pest_H_total <- raster("data/Pesticide_totalAPR_High.tif")

# mask out the sea area that is set to 0
data("wrld_simpl", package = 'maptools')
pest_L_total_crop <- mask(pest_L_total$Pesticide_totalAPR_Low, wrld_simpl)
pest_L_total_crop <- crop(pest_L_total_crop, wrld_simpl)
pest_H_total_crop <- mask(pest_H_total$Pesticide_totalAPR_High, wrld_simpl)
pest_H_total_crop <- crop(pest_H_total_crop, wrld_simpl)

# resave the high and low rasters for pesticide application before projecting into mollweide projection
writeRaster(pest_L_total_crop, filename = paste0(out_dir, "/Pesticide_totalAPR_Low_cropped.tif"))
writeRaster(pest_H_total_crop, filename = paste0(out_dir, "/Pesticide_totalAPR_High_cropped.tif"))

# reproject rasters as mollweide projection
pest_L_total_crop <- projectRaster(pest_L_total_crop, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
pest_H_total_crop <- projectRaster(pest_H_total_crop, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# convert for use in ggplot
pest_H_data <- as.data.frame(pest_H_total_crop, xy = TRUE)
pest_L_data <- as.data.frame(pest_L_total_crop, xy = TRUE)

# remove NAs
pest_H_data <- pest_H_data[!is.na(pest_H_data$Pesticide_totalAPR_High), ]
pest_L_data <- pest_L_data[!is.na(pest_L_data$Pesticide_totalAPR_Low), ]

# assign column for high or low estimate
pest_H_data$est <- "High"
pest_L_data$est <- "Low"

# set APR column name
names(pest_H_data)[3] <- "APR"
names(pest_L_data)[3] <- "APR"

# bind together the high and low estimate pesticide application dataa, and then sort the factors as high and low for plot
pest_data <- rbind(pest_H_data, pest_L_data)
pest_data$est <- factor(pest_data$est, levels = c("High", "Low"))

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

# save the map of pesticides
ggsave("pesticide_data_maps.png", dpi = 300, scale = 1)

