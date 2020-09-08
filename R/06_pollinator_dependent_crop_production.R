# script for pollinator dependent production

library(raster)
library(ggplot2)
library(rasterVis)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)

# source in additional functions
source("R/00_functions.R")

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("G:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

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
# subset as strings to filter from klein_cleaned
pollinated_crops <- unique(grep(paste(klein_cleaned$MonfredaCrop, collapse = "|"), unlisted_crops, value = TRUE))
pollinat_crops_simp <- gsub("G:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff/", "", pollinated_crops)
pollinat_crops_simp <- gsub('([^/]+$)', "", pollinat_crops_simp)
pollinat_crops_simp <- gsub('/', "", pollinat_crops_simp)

# read in each of the rasters
for(i in 1:length(pollinated_crops)){
  rate_rasters[[i]] <- raster(pollinated_crops[i])
  print(i)
}

# multiple each raster by the percentage attributable to pollination
# pollination dependent ratios
# each crop in the monfreda data can have a different pollination dependence
# 0 = no increase
# 0.05 = little
# 0.25 = modest
# 0.65 = great
# 0.95 = essential

# dependence ratio assignment
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "no increase"] <- 0
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "little"] <- 0.05
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "modest"] <- 0.25
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "great"] <- 0.65
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "essential"] <- 0.95

# subset klein_cleaned for those with crop data
klein_cleaned_filt <- klein_cleaned %>%
  filter(MonfredaCrop %in% pollinat_crops_simp)

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
stacked_rasters <- stack(rate_rasters)
crop.total <- sum(stacked_rasters, na.rm = F)

# test plot
plot(crop.total)

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
crop.total <- projectRaster(crop.total, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# convert the raster to a dataframe to plot with ggplot
crop_spdf <- as(crop.total, "SpatialPixelsDataFrame")
crop_df <- as.data.frame(crop_spdf)

# assign 0 in layer as NA
crop_df$layer[crop_df$layer == 0] <- NA

# bring in basemap for climate site plot 
base_map <- get_basemap()

base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# fortify the main map
map_fort <- fortify(base_map)

crop_df %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_tile(aes(x = x, y = y, fill = log10(layer))) +
  scale_fill_viridis("Pollination dependent production \n (kg)", na.value = "transparent") +
  coord_equal() +
  guides(fill = guide_colourbar(ticks = FALSE)) +
  theme(panel.background = element_blank(),
        panel.bord = element_blank(),
        panel.grid = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())

ggsave("pollinated_crop_production.png", scale = 1, dpi = 350)
