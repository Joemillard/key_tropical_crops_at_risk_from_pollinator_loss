# script for pollinator dependent production

library(raster)
library(ggplot2)
library(rasterVis)

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("G:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# loop through each directory and create a list of all files
all.ras <- NULL
crop.files <- list()

for(i in 1:length(cropdirs)){
  crop.files[[i]] <- list.files(cropdirs[i], pattern = "Production.tif$")
  crop.files[[i]] <- paste0(cropdirs[i], "/", crop.files[[i]])
}

# file paths for each of per hectare application
unlisted_crops <- unlist(crop.files)
rate_rasters <- list()

# subset the file paths for just those that are pollination dependent to some extent

# read in each of the rasters
for(i in 1:length(unlisted_crops)){
  rate_rasters[[i]] <- raster(unlisted_crops[i])
  print(i)
}

# multiple each raster by the percentage attributable to pollination

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
stacked_rasters <- stack(rate_rasters)
fert.total <- sum(stacked_rasters, na.rm = F)

plot(log1p(fert.total))

