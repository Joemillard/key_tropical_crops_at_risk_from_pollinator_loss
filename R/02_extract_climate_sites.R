# script for climate anomaly - dataframe rewrite to try and avoid repeated raster trim/crops
### Note: 2005 is the mean year for insect data
# try split up into three panels: overall (1901-2006), earlier (1901-1950), latter (1950-2006)

# load required libraries
library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(viridis)

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# read in the predicts pollinators
PREDICTS_pollinators <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_5.rds")

# PREDICTS sites with the month of the recording
PRED_sites <- PREDICTS_pollinators %>% select(Latitude, Longitude, Sample_end_latest) %>%
  filter(!is.na(Latitude)) %>%
  unique() %>%
  mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest))

#### calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]

# extract the points for each the predicts coordinates
PRED_sites_sp <- PRED_sites %>%
  select(Longitude, Latitude) %>%
  SpatialPoints()

# calc baseline (mean and sd)
calc_baseline <- function(data_file, func, pred_points, pred_points_sp){
  
  # calcualte either the mean or standard error for baseline, then extract points for predicts sites
  data_fin <- calc(data_file, func) %>%
    extract(pred_points_sp)
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- data.frame(pred_points[,1:2 ], data_fin)
  
  return(data_fin)
  
}

# calculate the mean baseline, and convert to character for merging
climate_start_mean <- calc_baseline(tmp1901_1905, 
                                    func = mean, 
                                    pred_points = PRED_sites, 
                                    pred_points_sp = PRED_sites_sp) %>%
  mutate(Latitude = as.character(Latitude)) %>%
  mutate(Longitude = as.character(Longitude))

# calculate the sd baseline, and convert to character for merging
climate_start_sd <- calc_baseline(tmp1901_1905, 
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
    
    # filter the raster for that date for the locations we have predicts sites
    PRED_coords <- cbind(PRED_sites_filt@coords, extract(ind_raster, PRED_sites_filt, na.rm = FALSE))
  
    # convert that set of dates to a dataframe
    ind_raster_frame <- as.data.frame(PRED_coords)
    
    # remove the extra coordinate columns for calculating the row means
    ind_raster_values <- ind_raster_frame %>% select(-Longitude, -Latitude)
  
    # calculate the mean values for each coordinate and bind back onto the coordinates
    raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% select(Longitude, Latitude)), rowMeans(ind_raster_values))
    colnames(raster_means[[i]]) <- c("end_date", "x", "y", "mean_value")

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
rbindlist(raster_means) %>%
  mutate(x = as.character(x)) %>%
  mutate(y = as.character(y)) %>%
  select(-end_date) %>%
  unique() %>%
  inner_join(climate_start_mean, by = c("x" = "Longitude", "y" = "Latitude")) %>%
  rename("mean_base" = "data_fin") %>%
  inner_join(climate_start_sd, by = c("x" = "Longitude", "y" = "Latitude")) %>%
  rename("sd_base" = "data_fin") %>%
  mutate(anomaly = mean_value - mean_base) %>%
  mutate(standard_anom = anomaly / sd_base)
