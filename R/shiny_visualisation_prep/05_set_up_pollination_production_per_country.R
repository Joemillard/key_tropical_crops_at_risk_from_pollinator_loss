# script for country level pollination dependent production


# read in packages
library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)
library(data.table)

# source in additional functions
source("R/00_functions.R")

# bring in basemap for climate site plot 
base_map <- get_basemap()

# reproject basemap
base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

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
klein_cleaned_filt <- subset_klein(av_dependence(klein_cleaned))

# multiply each raster by its pollination dependence for that crop
rate_rasters_adj <- list()
for(i in 1:length(rate_rasters)){
  rate_rasters_adj[[i]] <- rate_rasters[[i]] * klein_cleaned_filt$av[i]
  print(i)
}

# iteration vector for raster
adj_it_vec <- c()
adj_sum_vec <- c()
iteration <- c()

# select the top 5 crops for total pollination dependent production
for(i in 1:length(rate_rasters_adj)){
  adj_sum_vec[i] <- sum(rate_rasters_adj[[i]][])
  adj_it_vec[i] <- klein_cleaned_filt$MonfredaCrop[i]
  iteration[i] <- i
}

# sum of pollination dependent production production for each crop
pollination_production_sum <- data.frame(adj_it_vec, adj_sum_vec, iteration) %>%
  arrange(desc(adj_sum_vec)) %>%
  slice(1:22)

# subset original rate rasters for those in the top 20 pollination dependent production
rate_rasters_adj <- rate_rasters_adj[pollination_production_sum$iteration]

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
for(i in 1:length(rate_rasters_adj)){
  rate_rasters_adj[[i]] <- projectRaster(rate_rasters_adj[[i]], crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  print(i)
}

production_values <- list()

# convert raster of production to dataframe
for(i in 1:length(rate_rasters_adj)){
  production_values[[i]] <- as(rate_rasters_adj[[i]], "SpatialPixelsDataFrame")
  production_values[[i]] <- as.data.frame(production_values[[i]])
  print(i)
}

# set up list for spatial coordinates
spatial_production <- list()

# convert spatial dataframe to coordinates
spatial_production <- production_values[[1]] %>%
  dplyr::select(x, y) %>%
  unique() %>%
  SpatialPoints(proj4string = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# assign each set of coordinates to a country
country_coords <- over(spatial_production, base_map, by = "ISO2", returnList = FALSE)

production_values_coords <- list()

# assign the column for country, latitude and longitude back onto pollinator vulnerability data
for(i in 1:length(production_values)){
  production_values_coords[[i]] <- cbind(production_values[[i]], country_coords[c("SOVEREIGNT", "LON", "LAT")])
}

# bind the crop name onto each element of list
for(i in 1:length(production_values_coords)){
  production_values_coords[[i]]$crop <- pollination_production_sum$adj_it_vec[i]
  colnames(production_values_coords[[i]]) <- c("production", "x", "y", "SOVEREIGNT", "Lon", "Lat", "crop")
}

crop_sums <- list()

# too big to bind all together, so sum for each crop first
for(i in 1:length(production_values_coords)){
  crop_sums[[i]] <- production_values_coords[[i]] %>%
    group_by(crop, SOVEREIGNT) %>%
    summarise(total_production = sum(production))
}

# test for specific countries to check looks sensible
rbindlist(crop_sums) %>%
  filter(SOVEREIGNT == "Ivory Coast") %>%
  filter(total_production != 0) %>%
  mutate(crop = forcats::fct_reorder(crop, -total_production)) %>%
  ggplot() +
  geom_bar(aes(x = crop, y = total_production), stat = "identity")

# bind together and convert crop to characters and add in pollination dependence
bound_crop <- crop_sums %>%
  rbindlist() %>%
  inner_join(klein_cleaned_filt, by = c("crop" = "MonfredaCrop")) %>%
  dplyr::select(-standard_dev) %>%
  rename("pollination_dependence" = "av") %>%
  mutate(crop = factor(crop, levels = sort(unique(crop))))

# change the levels
levels(bound_crop$crop) <- c("Apple", "Bean", "Cocoa", "Coconut", "Coffee", "Cucumber",
                             "Eggplant", "Fruits (NE)", "Mango", "Melon", "Oilpalm", "Oilseed", 
                             "Peach", "Pear", "Plum", "Pumpkin","Rapeseed",  
                             "Soybean",  "Sunflower", "Tomato","Trop. fruits (NE)", "Watermelon")

saveRDS(bound_crop, "country_pollination_dependent_production.rds")


