library(dplyr)
library(ggplot2)
library(raster)
library(rworldmap)
library(rworldxtra)

# source in additional functions
source("R/00_functions.R")

# read in text file from Felipe without isolation term
trade_flow_new <- read.csv("data/trade_flow/Virtual_pollination_flow_no_isolation.csv", stringsAsFactors = FALSE)

# work out total exported per country
total_exported <- trade_flow_new %>%
  filter(year == 2001) %>%
  group_by(reporter_countries) %>%
  summarise(total_exported = sum(potential_poll_trade, na.rm = TRUE))

# bring in basemap for climate site plot 
base_map <- get_basemap()

# reproject basemap
base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# read in the total pollination production data

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

## calculate pollination dependence production
# select those with semi colon into a multiple rows
semi_colon_crop <- klein_cleaned$MonfredaCrop[grepl(";", klein_cleaned$MonfredaCrop)]

# subset those with semi colon and bind as duplicates onto the bottom of dataframe
semi_colon_crop <- klein_cleaned %>% filter(MonfredaCrop %in% semi_colon_crop) %>%
  mutate(MonfredaCrop = gsub("lemonlime; ", "", MonfredaCrop)) %>%
  mutate(MonfredaCrop = gsub("rasberry; ", "", MonfredaCrop)) %>%
  mutate(MonfredaCrop = gsub("cashew; ", "", MonfredaCrop))

# bind the rows for semi colon crops back onto the original klein data
klein_cleaned <- rbind(klein_cleaned, semi_colon_crop) %>% 
  mutate(MonfredaCrop = gsub("; citrusnes", "", MonfredaCrop)) %>%
  mutate(MonfredaCrop = gsub("; berrynes", "", MonfredaCrop)) %>%
  mutate(MonfredaCrop = gsub("; cashewapple", "", MonfredaCrop))

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
pollinated_crops <- grep(paste(unique(paste("/", klein_cleaned$MonfredaCrop, "_", sep = "")), collapse = "|"), unlisted_crops, value = TRUE)
pollinated_crops <- sort(pollinated_crops)
pollinat_crops_simp <- gsub("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff/", "", pollinated_crops)
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
# 0.45 = modest/great
# 0.65 = great
# 0.95 = essential

# dependence ratio assignment
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "no increase"] <- 0
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "little"] <- 0.05
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "modest"] <- 0.25
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "great"] <- 0.65
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "essential"] <- 0.95
klein_cleaned$dependence_ratio[klein_cleaned$Positive.impact.by.animal.pollination == "modest/great"] <- 0.45

# calculate average and standard deviation of pollination dependence for each Monfreda crop
av_dependence <- function(klein_cleaned){
  klein_cleaned_av <- klein_cleaned %>%
    group_by(MonfredaCrop) %>%
    mutate(av = mean(dependence_ratio, na.rm = TRUE)) %>%
    mutate(standard_dev = sd(dependence_ratio, na.rm = TRUE)) %>%
    ungroup() %>%
    select(MonfredaCrop, av, standard_dev) %>%
    unique()
  
  return(klein_cleaned_av)
}

# subset klein_cleaned for those with crop data
subset_klein <- function(klein_cleaned_av){
  klein_cleaned_filt <- klein_cleaned_av %>%
    filter(MonfredaCrop %in% pollinat_crops_simp) %>%
    arrange(MonfredaCrop)
  
  return(klein_cleaned_filt)
}

# run function for average pollination dependence and subset klein for those with crop data
klein_cleaned_filt <- subset_klein(av_dependence(klein_cleaned)) %>%
  arrange(MonfredaCrop)

# multiply each raster by its pollination dependence for that crop
rate_rasters_adj <- list()
for(i in 1:length(rate_rasters)){
  rate_rasters_adj[[i]] <- rate_rasters[[i]] * klein_cleaned_filt$av[i]
  print(i)
}

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
crop.total <- stack(rate_rasters_adj) %>% sum(na.rm = T)

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
crop.total <- projectRaster(crop.total, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

crop_points <- as(crop.total, "SpatialPixelsDataFrame")
crop_points <- as.data.frame(crop_points)

# convert spatial dataframe to coordinates
crop_points_spatial <- crop_points %>%
  dplyr::select(x, y) %>%
  unique() %>%
  SpatialPoints(proj4string = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# assign each set of coordinates to a country
country_coords <- over(crop_points_spatial, base_map, by = "ISO2", returnList = FALSE)

# assign the column for country, latitude and longitude back onto pollinator vulnerability data
crop_points <- cbind(crop_points, country_coords[c("SOVEREIGNT", "ISO3", "REGION", "LON", "LAT")])

# sum pollination dependent production per country
country_sums <- crop_points %>%
  group_by(ISO3) %>%
  summarise(total_country_production = sum(layer, na.rm = TRUE))

# merge the total pollination dependent production in with the import data
production_minus_export <- inner_join(country_sums, total_exported, by = c("SOVEREIGNT" = "reporter_countries")) %>%
  mutate(production_in_country = total_country_production - total_exported)

# subtract exported from produced to work out how much production stays in country


# sum the trade data for each year and convert to proportion of each countries exports
prop_flow <- trade_flow_new %>%
  filter(year == 2001) %>%
  group_by(reporter_countries, partner_countries) %>%
  summarise(total_flow = sum(potential_poll_trade, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(reporter_countries) %>%
  mutate(country_flow = sum(total_flow, na.rm = TRUE)) %>%
  mutate(percent_flow = (total_flow/country_flow) * 100)

# save proportional trade flow as an rds
saveRDS(prop_flow, "data/trade_flow/proportional_pollination_flow_non_isolation.rds")

  
