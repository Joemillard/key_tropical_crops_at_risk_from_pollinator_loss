# script for pollinating/non-pollinating insect abundance models

# load required libraries
library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(yarg)
library(rworldmap) 
library(rworldxtra)
library(lme4)
library(cowplot)
library(viridis)

# read in the original predicts database 
PREDICTS <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/Data/PREDICTS/database.rds") %>%
  mutate(COL_ID = as.character(COL_ID))

# source in additional functions
source("R/00_functions.R")

# set active temperature threshold 
temp_threshold <- 10

# number of each month
month_number <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
month_number_filt <- paste("\\.", month_number, "\\.", sep = "")

# calc baseline (mean and sd) - now for only months in which insects are considered active
calc_base_months <- function(X, pred_points, pred_points_sp, temp_threshold){

  # bind the extracted values back onto the predicts coordinates
  data_fin <- X %>%
    extract(pred_points_sp)
  
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)  %>%
    mutate(Latitude = as.character(Latitude)) %>%
    mutate(Longitude = as.character(Longitude))  %>%
    filter(data_fin >= temp_threshold)
  return(data_fin)
  
}

# calc baseline (mean and sd), with filter for required months
calc_baseline <- function(data_file, func, pred_points, pred_points_sp, month_filter){
  
  # calcualte either the mean or standard error for baseline, then extract points for predicts sites
  data_fin <- raster::subset(data_file, grep(month_filter, names(data_file)), value = T) %>%
    calc(func) %>%
    extract(pred_points_sp)
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)  %>%
    mutate(Latitude = as.character(Latitude)) %>%
    mutate(Longitude = as.character(Longitude)) %>%
    mutate(months = month_filter)
  
  return(data_fin)
  
}

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# function to calculate avaerge temp for each month
calc_baseline_months <- function(X, months, func){
  monthly_average <- raster::subset(X, grep(months, names(X)), value = T) %>%
    calc(fun = func) 
  
  return(monthly_average)
}

# subset the baseline period for each month, and then calculate the average for each month
monthly_mean <- list()
for(i in 1:length(month_number)){
  monthly_mean[[i]] <- calc_baseline_months(tmp1901_1930, month_number_filt[i], func = mean)
}

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds") %>%
  dplyr::select(-clade_rank, -confidence)  %>%
  filter(Class == "Insecta") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()

# filter the pollinators from the overall predicts database
PREDICTS_non_pollinating <- PREDICTS %>%
  filter(Class == "Insecta") %>%
  filter(!COL_ID %in% as.character(PREDICTS_pollinators_orig$COL_ID)) %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()

# bind together the two dataframes
pollinat_bound <- list(PREDICTS_pollinators_orig, PREDICTS_non_pollinating)

# set up vector for filtering for vertebrates and invertebrates
predict_climate_list <- list()

# PREDICTS data compilation
# filter for main pollinating taxa
PREDICTS_pollinators <- pollinat_bound[[1]]

# correct for sampling effort
PREDICTS_pollinators <- CorrectSamplingEffort(PREDICTS_pollinators)

# calculate site metrics including all species (confirmed and not confirmed pollinator)
order.sites.div <- SiteMetrics(diversity = PREDICTS_pollinators,
                               extra.cols = c("SSB", "SSBS", "Predominant_land_use", "UN_region", "Pollinating"),
                               sites.are.unique = TRUE,
                               srEstimators = TRUE)

# set id column for merging back into correct place
order.sites.div$id_col <- 1:nrow(order.sites.div)

# PREDICTS sites with the month of the recording
PRED_sites <- order.sites.div %>% select(id_col, Latitude, Longitude, Sample_end_latest) %>%
  mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest)) %>%
  filter(!is.na(Latitude))

# extract the points for each the predicts coordinates
PRED_sites_sp <- PRED_sites %>%
  select(Longitude, Latitude) %>%
  filter(!is.na(Latitude)) %>%
  SpatialPoints()

# calculate the mean baseline, and convert to character for merging
climate_start_mean <- lapply(X = monthly_mean, FUN = calc_base_months, 
                             pred_points = PRED_sites, pred_points_sp = PRED_sites_sp, temp_threshold = temp_threshold)

# for each coordinate, identify months that should be included in the active months
assign_month <- function(X, month_number){
  data_fin <- X %>%
    mutate(month = month_number)
  
  return(data_fin)
}

# use mapply to ouput the coordinate over 10 for each month, and then add in the month
baseline_active_months <- mapply(X = climate_start_mean, FUN = assign_month, month_number = month_number, SIMPLIFY = FALSE) %>%
  rbindlist() %>%
  group_by(id_col) %>%
  mutate(month_group = paste(month, collapse = "|")) %>%
  select(id_col, month_group) %>%
  unique()

# derive the unique active months for combinations of baseline to calculate
unique_month_combinations <- unique(baseline_active_months$month_group)

climate_start_all_mean <- list()

# calculate the mean baseline for each potential combination of active months
for(i in 1:length(unique_month_combinations)){
  climate_start_all_mean[[i]] <- calc_baseline(tmp1901_1930, 
                                    func = mean, 
                                    pred_points = PRED_sites, 
                                    pred_points_sp = PRED_sites_sp,
                                    month_filter = unique_month_combinations[i])
}

# bind the active months onto the average for those months (at the baseline)
fin_baseline_mean <- rbindlist(climate_start_all_mean) %>%
  inner_join(baseline_active_months, by = c("months" = "month_group", "id_col"))

climate_start_sd <- list()

# calculate the mean baseline  for each potential combination of active months
for(i in 1:length(unique_month_combinations)){
  
  # calculate the sd baseline, and convert to character for merging
  climate_start_sd[[i]] <- calc_baseline(tmp1901_1930, 
                                  func = stats::sd, 
                                  pred_points = PRED_sites, 
                                  pred_points_sp = PRED_sites_sp,
                                  month_filter = unique_month_combinations[i])
}

# bind the active months onto the average for those months (at the baseline)
fin_baseline_sd <- rbindlist(climate_start_sd) %>%
  inner_join(baseline_active_months, by = c("months" = "month_group", "id_col"))

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
    
    # remove any columns for months that aren't active from the baseline for those coordinates
    frame_months <- inner_join(ind_raster_frame, baseline_active_months) %>%
      melt(id = c("id_col", "Longitude", "Latitude", "month_group"))
    
    # derive active month string for filtering months
    month_string <- paste("\\.", strsplit(frame_months$month_group[1], split = "|", fixed = TRUE)[[1]], "\\.", collapse = "|", sep = "")

    frame_months <- cbind(frame_months, month_string) %>%
      filter(grepl(month_string, variable )) %>% reshape2::dcast(id_col ~ variable)
      
    
    # remove the extra coordinate columns for calculating the row means
    ind_raster_values <- ind_raster_frame %>% select(-id_col)
    
    # calculate the mean values for each coordinate and bind back onto the coordinates
    raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% select(id_col, Longitude, Latitude)), rowMeans(ind_raster_values))
    colnames(raster_means[[i]]) <- c("end_date", "id_col", "x", "y", "mean_value")
    
    # print the iteration number
    print(i)
  }
)
  
  
  ## adjust the mean value for each site for the baseline at that site
  # first, merge the baseline sd and mean by coordinate for each site
  adjusted_climate <- rbindlist(raster_means) %>%
    select(-end_date) %>%
    unique() %>%
    inner_join(fin_baseline_mean, by = "id_col") %>%
    rename("mean_base" = "data_fin") %>%
    inner_join(fin_baseline_sd, by = "id_col") %>%
    rename("sd_base" = "data_fin") %>%
    mutate(anomaly = mean_value - mean_base) %>%
    mutate(standard_anom = anomaly / sd_base)
  
  # bind the adjusted climate data back onto the predicts sites
  predicts_climate <- inner_join(order.sites.div, adjusted_climate, by = "id_col")
  
  # add 1 for abundance and simpson diversity
  predicts_climate$Total_abundance <- predicts_climate$Total_abundance + 1
  predicts_climate$Simpson_diversity <- predicts_climate$Simpson_diversity + 1
  
  # assign to list of predicts_climate and insects and vertebrates
  predict_climate_list[[j]] <- predicts_climate
  

# set up new lists for output
model_2c_abundance <- list()
abundance_model <- list()
main_plot_abundance <- list()

# run models for both species richness and total abundance
for(m in 1:length(pollinat_bound)){
  
  # species richness, standard anom as a factor
  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ standard_anom * Predominant_land_use + (1|SS), data = predict_climate_list[[m]]) 
  print(AIC(model_2c_abundance[[m]]))
  
  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ standard_anom * Predominant_land_use + (1|SS) + (1|SSB), data = predict_climate_list[[m]]) 
  print(AIC(model_2c_abundance[[m]]))
  
  # run predictions for the model of standard anomaly
  abundance_model[[m]] <- predict_continuous(model = model_2c_abundance[[m]],
                                             model_data = predict_climate_list[[m]],
                                             response_variable = "Total_abundance",
                                             categorical_variable = c("Predominant_land_use"),
                                             continuous_variable = c("standard_anom"),
                                             continuous_transformation = "",
                                             random_variable = c("SS", "SSB", "SSBS"))
  
  # plot for standardised anomaly and land-use for abundance
  main_plot_abundance[[m]] <- ggplot(abundance_model[[m]]) +
    geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use), size = 1.5) +
    geom_ribbon(aes(x = standard_anom, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
    scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
    scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
    xlab("Standardised climate anomaly") +
    ylab("Total abundance") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
}

# plot for the pollinating insects and non-pollinating insects - climate anomaly of 4 corresponds to ~100% abundance loss
plot_grid(main_plot_abundance[[1]] +
            ggtitle("Pollinating insects") + 
            scale_y_continuous(limits = c(1, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
            theme(legend.position = "bottom"), main_plot_abundance[[2]] + 
            ggtitle("Non-pollinating insects") +
            scale_y_continuous(limits = c(1, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
            theme(legend.position = "bottom"), ncol = 2)

ggsave("pollinating_non-pollinating_6.png", scale = 1, dpi = 350)
