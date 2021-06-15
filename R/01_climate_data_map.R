# script for pollinator dependent production

library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)
library(cowplot)

# source in additional functions
source("R/00_functions.R")

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

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
klein_cleaned_av <- klein_cleaned %>%
  group_by(MonfredaCrop) %>%
  mutate(av = mean(dependence_ratio, na.rm = TRUE)) %>%
  mutate(standard_dev = sd(dependence_ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  select(MonfredaCrop, av, standard_dev) %>%
  unique()

# subset klein_cleaned for those with crop data
klein_cleaned_filt <- klein_cleaned_av %>%
  filter(MonfredaCrop %in% pollinat_crops_simp) %>%
  arrange(MonfredaCrop)

# multiply each raster by its pollination dependence for that crop
rate_rasters_adj <- list()
for(i in 1:length(rate_rasters)){
  rate_rasters_adj[[i]] <- rate_rasters[[i]] * klein_cleaned_filt$av[i]
  print(i)
}

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
stacked_rasters <- stack(rate_rasters_adj)
crop.total <- sum(stacked_rasters, na.rm = T)

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

# split the pollination dependent production into a new variable of factors
crop_df$layer_group[crop_df$layer == 0] <- NA
crop_df$layer_group[crop_df$layer > 0 & crop_df$layer < 10^1] <- "0-10"
crop_df$layer_group[crop_df$layer >= 10^1 & crop_df$layer < 10^3] <- "10-1,000"
crop_df$layer_group[crop_df$layer >= 10^3 & crop_df$layer < 10^5] <- "1,000-100,000"

# create spatial data for crop values
# subset for unqiue sites
crop_df_locs <- crop_df %>%
  filter(!is.na(layer_group)) %>%
  dplyr::select(x, y) %>%
  unique() %>%
  SpatialPoints()

# assign new factor labels for 6 factors of production weight
crop_df$layer_group <- factor(crop_df$layer_group, levels = c("0-10", "10-1,000", "1,000-100,000"))

## standardised climate anomaly script
# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# take names of values for 1901 to 1931 - 30 year baseline
tmp1901_1905 <- tmp[[names(tmp)[1:360]]]

# calculate the mean and sd of the baseline values
tmp1901_1905mean <- calc(tmp1901_1905, mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)

# extract data for the years 2004-2006
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

### Calculate the standardised anomaly ###
# calc the mean for present time period
tmp2004_6mean <- calc(tmp[[names(tmp)[1237:1272]]], mean)

# calc mean for baseline
tmp2004_6_climate_anomaly <- (calc(tmp2004_6, mean) - tmp1901_1905mean)

# standardise the baseline
tmp2004_6std_climate_anomaly <- (calc(tmp2004_6, mean) - tmp1901_1905mean) / tmp1901_1905sd

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
tmp2004_6std_climate_anomaly <- projectRaster(tmp2004_6std_climate_anomaly, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

## combine the climate anomaly and the pollination dependence
# convert the raster to a format amenable to ggplot
# convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
climate_poll_values <- extract(tmp2004_6std_climate_anomaly, crop_df_locs, na.rm = FALSE)

# merge the climate standardised values onto the pollinator dependence data
climate_poll_data <- cbind((crop_df %>% filter(!is.na(layer_group))), climate_poll_values)

# create new dataframe for separate plot
climate_poll_data_map <- climate_poll_data

# group the categories of climate anomaly into factors
climate_poll_data_map$value_group[climate_poll_data_map$climate_poll_values > 2] <- "> 2"
climate_poll_data_map$value_group[climate_poll_data_map$climate_poll_values > 1 & climate_poll_data_map$climate_poll_values <= 2] <- "1 - 2"
climate_poll_data_map$value_group[climate_poll_data_map$climate_poll_values > 0.5 & climate_poll_data_map$climate_poll_values <= 1] <- "0.5 - 1"
climate_poll_data_map$value_group[climate_poll_data_map$climate_poll_values > 0.25 & climate_poll_data_map$climate_poll_values <= 0.5] <- "0.25 - 0.5"
climate_poll_data_map$value_group[climate_poll_data_map$climate_poll_values >= 0 & climate_poll_data_map$climate_poll_values <= 0.25] <- "0 - 0.25"
climate_poll_data_map$value_group[climate_poll_data_map$climate_poll_values < 0] <- "< 0"

# order the levels of those factors
climate_poll_data_map$value_group <- factor(climate_poll_data_map$value_group, levels = c("> 2", "1 - 2", "0.5 - 1", "0.25 - 0.5", "0 - 0.25", "< 0"))

# # plot the ggplot map for climate anomaly
climate_poll_data_map %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "grey", alpha = 0.3) +
  geom_tile(aes(x = x, y = y, fill = value_group)) +
  scale_fill_viridis_d("Standardised climate anomaly",
                    na.value = "transparent",
                    labels = c("> 2", "1 - 2", "0.5 - 1", "0.25 - 0.5", "0 - 0.25", "< 0", "No pollinator dependent prod.")) +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.bord = element_blank(),
        panel.grid = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())

# save the map for climate anomaly
ggsave("pollination_dependent_climate_anomaly_3.png", scale = 1.2, dpi = 350)
