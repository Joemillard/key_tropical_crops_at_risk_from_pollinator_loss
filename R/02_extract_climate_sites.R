# script for climate anomaly - dataframe rewrite to try and avoid repeated raster trim/crops
### Note: 2005 is the mean year for insect data
# try split up into three panels: overall (1901-2006), earlier (1901-1950), latter (1950-2006)

# load required libraries
library(raster)
library(ggplot2)
library(dplyr)

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

# calculate the mean and sd of the baseline values as rasters
tmp1901_1905mean <- calc(tmp1901_1905, mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)

# convert the climate anomaly raster for the beginning of the series to a spatial pixels data frame, and then rename the columns
climate_start <- as(tmp1901_1905sd, "SpatialPixelsDataFrame")
climate_start_df <- as.data.frame(climate_start)
colnames(climate_start_df) <- c("value", "x", "y")

####

# set up empty list for each dataframe
raster_means <- list()

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

    # print 
    print(i)
  }
)

# join the mean temperatures back onto the PREDICTS sites by minimising distance between coordinates











# extract data for the years 2004-2006
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

### Calculate the standardised anomaly ###
# calc the mean for present time period
tmp2004_6mean <- calc(tmp[[names(tmp)[1237:1272]]], mean)

# calc mean for baseline
tmp2004_6_climate_anomaly <- (calc(tmp2004_6, mean) - tmp1901_1905mean)

# standardise the baseline
tmp2004_6std_climate_anomaly <- (calc(tmp2004_6, mean) - tmp1901_1905mean) / tmp1901_1905sd