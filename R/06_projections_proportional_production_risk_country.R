# script for country level proportion production risk projections at the cell level

# read in packages
library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)
library(cowplot)
library(data.table)
library(lme4)
library(yarg)
library(forcats)
library(ggforce)
library(ggrepel)

# source in additional functions
source("R/00_functions.R")

# read in local population sizes
population_size <- read.csv("data/population-past-future.csv") %>%
  filter(Year %in% c(2015, 2016, 2017, 2018, 2019)) %>%
  group_by(Code) %>%
  summarise(average_pop = mean(Population..historical.estimates.and.future.projections.))

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds")

# set up the starting directory for future climate data
SSP_directory <- ("D:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies")

# read in fao to monfreda conversion
fao_monfreda <- read.csv("data/trade_flow/FAO_Monfreda_conv.csv", stringsAsFactors = FALSE) %>%
  mutate(Cropname_FAO = gsub(",", "", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("corian.", "coriander", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Berries Nes", "Berries nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Stone fruit nes", "Fruit stone nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Chestnuts", "Chestnut", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Citrus fruit nes", "Fruit citrus nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Fruit Fresh Nes", "Fruit fresh nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Leguminous vegetables nes", "Vegetables leguminous nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Kolanuts", "Kola nuts", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Other melons", "Melons other", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Oilseeds Nes", "Oilseeds nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Pepper \\(Piper ", "Pepper \\(piper ", Cropname_FAO))

# calculate per country average total production value
fao_prod_value <- read.csv("data/trade_flow/FAOSTAT_data_3-28-2022_total_value.csv", stringsAsFactors = FALSE) %>%
  filter(Flag.Description == "Calculated data") %>%
  filter(!grepl("Meat", Item)) %>%
  filter(!grepl("Eggs", Item)) %>%
  filter(!grepl("China, mainland", Area)) %>%
  filter(!grepl("China, Hong Kong SAR", Area)) %>%
  filter(!grepl("Milk", Item)) %>%
  rename(economic_value = Value) %>%
  mutate(economic_value = economic_value * 1000)

# calculate per country average total production value
fao_prod <- read.csv("data/trade_flow/FAOSTAT_data_3-28-2022_total_production.csv", stringsAsFactors = FALSE) %>%
  filter(Flag.Description == "Official data") %>%
  filter(!grepl("Meat", Item)) %>%
  filter(!grepl("Eggs", Item)) %>%
  filter(!grepl("China, mainland", Area)) %>%
  filter(!grepl("China, Hong Kong SAR", Area)) %>%
  filter(!grepl("Milk", Item)) %>%
  rename(production = Value) %>%
  mutate(production = production * 1000)

# calc average price per crop
joined_prod_value <- inner_join(fao_prod, fao_prod_value, by = c("Year", "Area", "Item")) %>%
  mutate(price_per_kg = economic_value / production) %>% 
  group_by(Area, Item) %>%
  summarise(mean_price_kg = mean(price_per_kg, na.rm = TRUE)) %>%
  group_by(Item) %>%
  summarise(overall_price_kg = median(mean_price_kg, na.rm = TRUE)) %>% 
  mutate(Item = gsub(",", "", Item)) %>%
  inner_join(fao_monfreda, by = c("Item" = "Cropname_FAO"))

# PREDICTS data compilation
# filter for main pollinating taxa
PREDICTS_pollinators <- PREDICTS_pollinators_orig %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  dplyr::filter(Phylum %in% "Arthropoda") %>%
  droplevels()

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
order.sites.div <- SiteMetrics(diversity = PREDICTS_pollinators,
                               extra.cols = c("SSB", "SSBS", "Predominant_land_use", "UN_region"),
                               sites.are.unique = TRUE,
                               srEstimators = TRUE)

# set id column for merging back into correct place
order.sites.div$id_col <- 1:nrow(order.sites.div)

# PREDICTS sites with the month of the recording
PRED_sites <- order.sites.div %>% select(id_col, Latitude, Longitude, Sample_end_latest) %>%
  mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest)) %>%
  filter(!is.na(Latitude))

# calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# extract the points for each the predicts coordinates
PRED_sites_sp <- PRED_sites %>%
  select(Longitude, Latitude) %>%
  filter(!is.na(Latitude)) %>%
  SpatialPoints()

# calc baseline (mean and sd)
calc_baseline <- function(data_file, func, pred_points, pred_points_sp){
  
  # calcualte either the mean or standard error for baseline, then extract points for predicts sites
  data_fin <- calc(data_file, func) %>%
    extract(pred_points_sp)
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)
  
  return(data_fin)
  
}

# calculate the mean baseline, and convert to character for merging
climate_start_mean <- calc_baseline(tmp1901_1930, 
                                    func = mean, 
                                    pred_points = PRED_sites, 
                                    pred_points_sp = PRED_sites_sp) %>%
  mutate(Latitude = as.character(Latitude)) %>%
  mutate(Longitude = as.character(Longitude))

# calculate the sd baseline, and convert to character for merging
climate_start_sd <- calc_baseline(tmp1901_1930, 
                                  func = stats::sd, 
                                  pred_points = PRED_sites, 
                                  pred_points_sp = PRED_sites_sp) %>%
  mutate(Latitude = as.character(Latitude)) %>%
  mutate(Longitude = as.character(Longitude))

# calculate the mean temperatures for each predicts site, 11 months previously
# set up empty list for each dataframe
raster_means <- list()

# time the length of the loop over each unique date in predicts sites
system.time(
  # for each unique site month, select the pixels for that date, and 11 months previously, and then convert each raster set to a dataframe
  for(i in 1:length(unique(PRED_sites$Sample_end_latest))){
    
    # create unique list of the end dates for the predicts sites
    pred_dates <- unique(PRED_sites$Sample_end_latest)
    
    # select the raster pixels for the end date of the predicts site, and then the index of that name
    date_raster <- grepl(pred_dates[i], names(tmp))
    site_index <- which(date_raster == TRUE)
    
    # select the previous 11 indices - i.e. 11 months previous worth of pixels
    ind_raster <- tmp[[names(tmp)[(site_index - 11): site_index]]]
    
    # filter the raster for only coordinates we have PREDICTS sites for that date
    PRED_sites_filt <- PRED_sites %>%
      filter(Sample_end_latest == pred_dates[i]) %>%
      select(Longitude, Latitude) %>%
      SpatialPoints()
    
    # identify site ids for merging
    site_ids <- PRED_sites %>%
      filter(Sample_end_latest == pred_dates[i]) %>%
      select(id_col, Longitude, Latitude)
    
    # filter the raster for that date for the locations we have predicts sites
    PRED_coords <- cbind(site_ids, extract(ind_raster, PRED_sites_filt, na.rm = FALSE))
    
    # convert that set of dates to a dataframe
    ind_raster_frame <- as.data.frame(PRED_coords)
    
    # remove the extra coordinate columns for calculating the row means
    ind_raster_values <- ind_raster_frame %>% select(-id_col, -Longitude, -Latitude)
    
    # calculate the mean values for each coordinate and bind back onto the coordinates
    raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% select(id_col, Longitude, Latitude)), rowMeans(ind_raster_values))
    colnames(raster_means[[i]]) <- c("end_date", "id_col", "x", "y", "mean_value")
    
    # print the iteration number
    print(i)
  }
)

# checking the average temperate 11 months previous to each predicts site
rbindlist(raster_means) %>%
  ggplot() +
  geom_point(aes(x = x, y = y, colour = mean_value)) + 
  scale_colour_viridis()

## adjust the mean value for each site for the baseline at that site
# first, merge the baseline sd and mean by coordinate for each site
adjusted_climate <- rbindlist(raster_means) %>%
  select(-end_date) %>%
  unique() %>%
  inner_join(climate_start_mean, by = "id_col") %>%
  rename("mean_base" = "data_fin") %>%
  inner_join(climate_start_sd, by = "id_col") %>%
  rename("sd_base" = "data_fin") %>%
  mutate(anomaly = mean_value - mean_base) %>%
  mutate(standard_anom = anomaly / sd_base)

# bind the adjusted climate data back onto the predicts sites
predicts_climate <- inner_join(order.sites.div, adjusted_climate, by = "id_col")

# bring in basemap for climate site plot 
base_map <- get_basemap()

# reproject basemap
base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# add 1 for abundance and simpson diversity
predicts_climate$Total_abundance <- predicts_climate$Total_abundance + 1

# run model for total abundance for insect pollinators
model_2c_abundance <- lmerTest::lmer(log(Total_abundance) ~ standard_anom * Predominant_land_use + (1|SS) + (1|SSB), data = predicts_climate) 

# run predictions for the model of standard anomaly
abundance_model <- predict_continuous(model = model_2c_abundance,
                                      model_data = predicts_climate,
                                      response_variable = "Total_abundance",
                                      categorical_variable = c("Predominant_land_use"),
                                      continuous_variable = c("standard_anom"),
                                      continuous_transformation = "",
                                      random_variable = c("SS", "SSB", "SSBS"))

# plot for standardised anomaly and land-use for abundance
main_plot_abundance <- abundance_model %>% 
  ggplot() +
  geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = standard_anom, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
  scale_y_continuous("Total abundance", breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321), labels = c(5, 10, 20, 40, 80, 160, 320)) +
  scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  xlab("Standardised climate anomaly") +
  ylab("Total abundance") +
  theme_bw() +
  theme(panel.grid = element_blank())

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
rate_rasters_adj <- list()
rate_rasters_all <- list()

# subset the file paths for just those that are pollination dependent to some extent
# subset as strings to filter from klein_cleaned
pollinated_crops <- grep(paste(unique(paste("/", klein_cleaned$MonfredaCrop, "_", sep = "")), collapse = "|"), unlisted_crops, value = TRUE)
pollinat_crops_simp <- gsub("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff/", "", pollinated_crops)
pollinat_crops_simp <- gsub('([^/]+$)', "", pollinat_crops_simp)
pollinat_crops_simp <- gsub('/', "", pollinat_crops_simp)

# read in each of the rasters for pollinated crops
for(i in 1:length(pollinated_crops)){
  rate_rasters[[i]] <- terra::rast(pollinated_crops[i])
  print(i)
}

# read in each of the rasters for pollinated crops
for(i in 1:length(crop.files)){
  rate_rasters_all[[i]] <- terra::rast(crop.files[[i]])
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

rm(rate_rasters)

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
crop.total <- terra::rast(rate_rasters_adj) %>% terra::app(fun = "sum", na.rm = TRUE) %>% raster()
crop.total_all <- terra::rast(rate_rasters_all) %>% terra::app(fun = "sum", na.rm = TRUE) %>% raster()

# baseplot of proportion needing pollination
prop_poll <- crop.total / crop.total_all
plot(prop_poll)

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
crop.total_all <- projectRaster(crop.total_all, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
crop.total <- projectRaster(crop.total, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# calculate total pollination dependent production
total_production <- sum(crop.total[])

## standardised climate anomaly script
# calculate the mean and sd of the baseline values
tmp1901_1930mean <- calc(tmp1901_1930, mean)
tmp1901_1930sd <- calc(tmp1901_1930, stats::sd)

## read in the rasters for the future data, start with SSP585
# set up historical change to be added on
months.1979.2013 <- 937:1356

# calculate the average temperature for 1979-2013 onto which anomaly is added
hist.mean.temp.1979.2013 <- stack(stackApply(x = tmp[[months.1979.2013]],
                                             indices = (rep(1:35,each=12)),fun = mean))
hist.mean.temp.1979.2013 <- stackApply(x = hist.mean.temp.1979.2013,indices = rep(1,35),
                                       fun = mean)

# selection of years and empty year list
years <- 2048:2050
years_list <- list()

# set up list of years
for(i in 1:33){
  years <- years - 1
  years_list[[i]] <- years
}

# need to run for the average of climate models 
# need to run for each RCP scenario

# set up list for climate anomalies
tmp2069_71std_climate_anomaly <- list()

# average the set of climate models and calculate climate anomaly for the average
average_clim_models <- function(yr, RCP, clim_models){
  
  # print the set of years for that iteration
  print(yr)
  
  # subset for all files for rcp85, the set of years for that iteration, for the models of that iteration
  all.model.files <- all.files[grepl(RCP, all.files) & grepl(yr, all.files)]
  all.model.files <- all.model.files[grepl(clim_models, all.model.files)]
  
  meant.anom <- mean(stack(lapply(X = all.model.files, function(f){
    
    ras <- stack(f)$"X0.1"
    
  })), na.rm=TRUE)
  
  # add predicted anomaly onto the baseline average from historical temperatures
  meant <- hist.mean.temp.1979.2013 + (meant.anom/10)
  
  return(meant)
  
}

# set up vector of climate models
RCP_scenarios <- c("rcp60")
climate_model_combs_adj <- c("GFDL|HadGEM2|IPSL|MIROC5")

# iterate through each set of years as a rolling average
for(i in 1:length(years_list)){
      
  # file path for ISIMIP data
  all.files <- dir(path = SSP_directory,recursive = TRUE, full.names = TRUE)
      
  # using RCP 8.5 calculate average of separate models
  mean.temp.2069.2071 <- stack(lapply(X = years_list[[i]], FUN = average_clim_models, RCP = RCP_scenarios, clim_models = climate_model_combs_adj))
  
  mean.temp.2069.2071 <- stackApply(x = mean.temp.2069.2071,indices = rep(1,3), fun = mean)
      
  # calc the anomalies for the future years
  tmp2069_71_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1930mean)
  tmp2069_71std_climate_anomaly[[i]] <- (mean.temp.2069.2071-tmp1901_1930mean) / tmp1901_1930sd
     
  # reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
  tmp2069_71std_climate_anomaly[[i]] <- projectRaster(tmp2069_71std_climate_anomaly[[i]], crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
}
    
# set up empty list for standardised climate anomaly
std_high_abun_adj <- list()
std_anom_high <- list()
    
# predict abundance at 0 warming on cropland
zero_data <- data.frame("standard_anom" = 0, Predominant_land_use = "Cropland")
zero_warming_abundance <- predict(model_2c_abundance, zero_data, re.form = NA)
zero_warming_abundance <- exp(zero_warming_abundance)
    
# for each set of climate anomaly data, predict abundance reduction for all climate anomaly values in each cell
# and then sum abundance adjusted pollination dependence
for(i in 1:length(tmp2069_71std_climate_anomaly)){
      
  # convert the raster to a dataframe to plot with ggplot
  std_high_abun_adj[[i]] <- as(tmp2069_71std_climate_anomaly[[i]], "SpatialPixelsDataFrame")
  std_high_abun_adj[[i]] <- as.data.frame(std_high_abun_adj[[i]])
      
  # set up prediction data on basis of that set of years
  new_data_pred <- data.frame("standard_anom" = std_high_abun_adj[[i]]$layer, Predominant_land_use = "Cropland")
      
  # predict abundance for climate anomaly and join to data frame
  predicted_abundance <- predict(model_2c_abundance, new_data_pred, re.form = NA)
  std_high_abun_adj[[i]]$abundance <- exp(predicted_abundance)
      
  # for any location that's cooled abundance is that at no warming
  std_high_abun_adj[[i]]$abundance[std_high_abun_adj[[i]]$layer <= 0] <- zero_warming_abundance
      
  # calculate percentage change from place with 0 warming, and convert to vulnerability
  std_high_abun_adj[[i]]$abundance_change <- 1 - (std_high_abun_adj[[i]]$abundance / zero_warming_abundance)
      
  # convert spatial dataframe to coordinates
  std_anom_high[[i]] <- std_high_abun_adj[[i]] %>%
    dplyr::select(x, y) %>%
    unique() %>%
    SpatialPoints(proj4string = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
}
    
# for each set of coordinates, extract the pollination dependent values and sum
for(i in 1:length(std_anom_high)){
      
  # convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
  std_high_abun_adj[[i]]$production <- extract(crop.total, std_anom_high[[i]], na.rm = FALSE)
  std_high_abun_adj[[i]]$total_production <- extract(crop.total_all, std_anom_high[[i]], na.rm = FALSE)
  std_high_abun_adj[[i]]$pollinator_vulnerability <- unlist(std_high_abun_adj[[i]]$production * std_high_abun_adj[[i]]$abundance_change) / std_high_abun_adj[[i]]$total_production
}

# assign each set of coordinates to a country
country_coords <- over(std_anom_high[[1]], base_map, by = "ISO2", returnList = FALSE)

# assign the column for country, latitude and longitude back onto pollinator vulnerability data
for(i in 1:length(std_high_abun_adj)){
  std_high_abun_adj[[i]] <- cbind(std_high_abun_adj[[i]], country_coords[c("SOVEREIGNT", "GBD", "ISO3", "SRES", "continent", "LON", "LAT")])
}

# check coordinates and countries are vaguely correct
std_high_abun_adj[[1]] %>% 
  filter(!is.na(SOVEREIGNT)) %>% 
  ggplot() +
    geom_point(aes(x = x, y = y, colour = SOVEREIGNT)) + theme(legend.position = "none")

country_sums <- list()

# for each yearly set, calculate the pollinator vulnerability for each country
for(i in 1:length(std_high_abun_adj)){
  country_sums[[i]] <-  std_high_abun_adj[[i]] %>%
    group_by(SOVEREIGNT, ISO3, SRES, continent, GBD) %>%
    summarise(total = quantile(pollinator_vulnerability, probs = c(0.5), na.rm = TRUE), 
              upp_conf = quantile(pollinator_vulnerability, probs = c(0.975), na.rm = TRUE), 
              lower_conf = quantile(pollinator_vulnerability, probs = c(0.025), na.rm = TRUE),
              pollination_production = sum(production),
              all_production = sum(total_production),
              percent_pollinated = (pollination_production/all_production) * 100) %>%
    ungroup() %>%
    mutate(year = 2049 - i)
} 

# plot for trends in pollination vulnerability
change_obj <- rbindlist(country_sums) %>% 
  group_by(SOVEREIGNT, ISO3, SRES, continent, GBD) %>%
  arrange(year) %>%
  mutate(av_total = mean(total)) %>%
  mutate(change = max(total) - min(total)) %>%
  ungroup()

# create object just with summary values
plot_obj <- change_obj %>%
  select(SOVEREIGNT, ISO3, SRES, continent, GBD, av_total, change, pollination_production, percent_pollinated) %>%
  unique() %>% 
  filter(!is.na(change)) %>%
  filter(!is.na(continent))

# add separate regions
plot_obj$main_region[plot_obj$continent %in% c("Eurasia") & plot_obj$SRES %in% c("Central and Eastern Europe (EEU)", 
                                                                                 "Western Europe (WEU)")] <- "North America & Europe"
plot_obj$main_region[plot_obj$continent %in% c("North America")] <- "North America & Europe"
plot_obj$main_region[plot_obj$continent %in% c("Eurasia") & plot_obj$SRES %in% c("Centrally planned Asia and China (CPA)", 
                                                                                 "Middle East and North Africa (MEA)",
                                                                                 "Other Pacific Asia (PAS)", 
                                                                                 "Pacific OECD (PAO)", 
                                                                                 "South Asia (SAS)")] <- "Asia & Australia"
plot_obj$main_region[plot_obj$continent %in% c("Australia")] <- "Asia & Australia"
plot_obj$main_region[plot_obj$continent %in% c("Africa")] <- "Africa"
plot_obj$main_region[plot_obj$continent %in% c("South America and the Caribbean")] <- "South America & the Caribbean"

# correcting for former soviet union states
plot_obj$main_region[plot_obj$SRES %in% c("Newly Independent States of FSU (FSU)") & plot_obj$GBD== "Asia, Central"] <- "Asia & Australia"
plot_obj$main_region[plot_obj$SRES %in% c("Newly Independent States of FSU (FSU)") & plot_obj$GBD == "Europe, Eastern"] <- "North America & Europe"



# merge the price per kg value with the pollination dependent monfreda
joined_crop_val <- left_join(klein_cleaned_filt, joined_prod_value, by = c("MonfredaCrop" = "CROPNAME"))

# calculate pollination dependent production value per crop globally
# raster layers adjusted for value
# multiply each raster by its pollination dependence for that crop
rate_rasters_adj_val <- list()
for(i in 1:length(rate_rasters_adj)){
  rate_rasters_adj_val[[i]] <- rate_rasters_adj[[i]] * joined_crop_val$overall_price_kg[i]
  print(i)
}

# sum total cell level price of pollination dependent crops geographically
crop_total_price <- terra::rast(rate_rasters_adj_val) %>% terra::app(fun = "sum", na.rm = TRUE) %>% raster()

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
crop_total_price <- projectRaster(crop_total_price, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# remove rasters to free up space
rm(rate_rasters_adj)
rm(rate_rasters_adj_val)

# convert the raster to a dataframe to plot with ggplot
crop_total_price_frame <- as(crop_total_price, "SpatialPixelsDataFrame")
crop_total_price_frame <- as.data.frame(crop_total_price_frame) %>%
  filter(sum > 0) %>%
  rename(layer = sum)

# convert spatial dataframe to coordinates
crop_total_price_points <- crop_total_price_frame %>%
  dplyr::select(x, y) %>%
  unique() %>%
  SpatialPoints(proj4string = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# assign each set of coordinates to a country
country_value <- over(crop_total_price_points, base_map, by = "ISO2", returnList = FALSE)

# merge countries back onto values
country_value_bound <- cbind(crop_total_price_frame, country_value[c("SOVEREIGNT", "SRES", "continent", "LON", "LAT", "GDP_MD_EST", "NAME_FORMA")])

# calc country totals
country_totals <- country_value_bound %>%
  group_by(SOVEREIGNT, GDP_MD_EST, NAME_FORMA, SRES, continent) %>%
  summarise(total_value = sum(layer))
  
# check coordinates and countries are vaguely correct
country_value_bound %>% 
  filter(!is.na(SOVEREIGNT)) %>% 
  ggplot() +
  geom_point(aes(x = x, y = y, colour = SOVEREIGNT)) + theme(legend.position = "none")

# join plot data onto the country total value
plot_obj <- inner_join(plot_obj, country_totals, by = c("SOVEREIGNT", "SRES", "continent"))

plot_obj_pop <- plot_obj %>%
  inner_join(population_size, by = c("ISO3" = "Code")) %>%
  mutate(per_GDP_pollination = total_value/GDP_MD_EST) %>%
  filter(!(SOVEREIGNT == "United Kingdom" & is.na(NAME_FORMA))) %>%
  filter(NAME_FORMA != "Hong Kong Special Administrative Region") %>%
  filter(NAME_FORMA != "Gaza Strip")

# select top for adding text
plot_obj_top <- plot_obj_pop %>%
  filter(SOVEREIGNT  %in% c("Ivory Coast", "Somalia", "Guinea Bissau",
                            "Philippines", "Indonesia", "Papua New Guinea",
                            "Turkey", 
                            "Suriname", "Haiti") | NAME_FORMA  %in% c("Commonwealth of Puerto Rico")) %>%
  mutate(SOVEREIGNT = gsub("United States of America", "Puerto Rico", SOVEREIGNT))
  
# assign sub regions
plot_obj_pop$sub_region[plot_obj_pop$SRES %in% c("Central and Eastern Europe (EEU)", "Western Europe (WEU)")] <- "Europe"
plot_obj_pop$sub_region[plot_obj_pop$continent %in% c("North America")] <- "North America"
plot_obj_pop$sub_region[plot_obj_pop$continent %in% c("Australia")] <- "Australia"
plot_obj_pop$sub_region[plot_obj_pop$continent %in% c("Africa")] <- "Africa"
plot_obj_pop$sub_region[plot_obj_pop$continent %in% c("South America and the Caribbean")] <- "South America"
plot_obj_pop$sub_region[plot_obj_pop$main_region %in% c("North America & Europe") & is.na(plot_obj_pop$sub_region)] <- "Europe"
plot_obj_pop$sub_region[plot_obj_pop$GBD == "Caribbean"] <- "The Caribbean"
plot_obj_pop$sub_region[is.na(plot_obj_pop$sub_region)] <- "Asia"

# assign nudge parameter
plot_obj_top$y_nudge[plot_obj_top$main_region == "Africa"] <- 0.09
plot_obj_top$y_nudge[plot_obj_top$main_region == "Asia & Australia"] <- 0.3
plot_obj_top$y_nudge[plot_obj_top$main_region == "North America & Europe"] <- 0.3
plot_obj_top$y_nudge[plot_obj_top$main_region == "South America & the Caribbean"] <- 0.09

# build export risk plot
ggplot(plot_obj_pop) +
    geom_point(aes(x = change, y = av_total, size = per_GDP_pollination, colour = sub_region),  alpha = 0.7) +
    geom_label_repel(aes(x = change, y = av_total, label = SOVEREIGNT), data = plot_obj_top, alpha = 0.5,
                    nudge_x = .01,
                     nudge_y = c(0.012, 0.01, 0.01, 0.01, -0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                     segment.curvature = -1e-20,) +
    scale_x_continuous("Change in crop pollination risk", expand = c(0, 0), limits = c(0, 0.05), breaks = c(0, 0.02, 0.04), 
                      labels = c("0", "0.02", "0.04")) +
    scale_y_continuous("Overall crop pollination risk", expand = c(0, 0), limits = c(0, 0.14), breaks = c(0, 0.04, 0.08, 0.12),
                       labels = c("0", "0.04","0.08", "0.12")) +
    scale_colour_manual("Geographic region", values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999")) +
    scale_fill_manual("Geographic region", values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999")) +
    scale_size_continuous("Pollination dependent production value per GDP (US$/annum)", breaks = c(25, 50, 75, 100, 125), labels = c(25, 50, 75, 100, 125)) +
    theme_bw() +
    facet_wrap(~main_region, ncol = 4) +
    guides(size = guide_legend(order = 2, nrow = 2), 
             colour = guide_legend(order = 1, nrow = 2)) +
    theme(panel.grid = element_blank(), legend.position = "bottom", legend.box="vertical") 

ggsave("top_change_country_11.png", scale = 1.1, dpi = 350)

# write file to csv for Silvia
write.csv(plot_obj_pop, "country_level_risk.csv")
