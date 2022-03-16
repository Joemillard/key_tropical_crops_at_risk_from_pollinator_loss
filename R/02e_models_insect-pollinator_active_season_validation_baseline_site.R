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
library(parallel)

# number of each month
month_number <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
month_number_filt <- paste("\\.", month_number, "\\.", sep = "")

# set up cores for parallel processing
cl <- makeCluster(detectCores())

# read in packages and data for each parallel session
clusterEvalQ(cl, {
  
  library(raster)
  library(dplyr)
  library(data.table)
  library(parallel)
  
})

# read in the original predicts database 
PREDICTS <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/Data/PREDICTS/database.rds") %>%
  mutate(COL_ID = as.character(COL_ID))

# source in additional functions
source("R/00_functions.R")

# function to calculate average temp for each month
calc_baseline_month_av <- function(X, data_file, func){
  monthly_average <- raster::subset(data_file, grep(X, names(data_file)), value = T) %>%
    calc(fun = func) 
  
  return(monthly_average)
}

# subset the baseline for just month temps that are above the threshold
subset_base_months <- function(X, pred_points, pred_points_sp, temp_threshold){
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- X %>%
    extract(pred_points_sp)
  
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)  %>%
    mutate(Latitude = as.character(Latitude)) %>%
    mutate(Longitude = as.character(Longitude))  %>%
    filter(data_fin >= temp_threshold)
  return(data_fin)
  
}

# calc baseline (mean and sd) - now for only months in which insects are considered active
calc_baseline <- function(X, func, pred_points, pred_points_sp, data_file){
  
  # calcualte either the mean or standard error for baseline, then extract points for predicts sites
  data_fin <- raster::subset(data_file, grep(X, names(data_file)), value = T) %>%
    calc(func) %>%
    extract(pred_points_sp)
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)  %>%
    mutate(Latitude = as.character(Latitude)) %>%
    mutate(Longitude = as.character(Longitude)) %>%
    mutate(months = X)
  
  return(data_fin)
  
}

# for each coordinate, identify months that should be included in the active months
assign_month <- function(X, month_number){
  data_fin <- X %>%
    mutate(month = month_number)
  
  return(data_fin)
}

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# subset the baseline period for each month, and then calculate the average for each month
monthly_mean <- lapply(X = month_number_filt, FUN = calc_baseline_month_av, data_file = tmp1901_1930, func = mean)

# read in the predicts pollinators
PREDICTS_pollinators <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds") %>%
  dplyr::select(-clade_rank, -confidence)  %>%
  filter(Class == "Insecta") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()



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

# set active temperature threshold, and empty list to assign to
temp_threshold <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)
predict_climate_list <- list()

# record initial time to calculate length of process time
st1 <- Sys.time()

# iterate through each temperature threshold value
for(j in 1:length(temp_threshold)){
  
  # subset baseline (where there are PREDICTS sites) for places where active months greater than threshold, and then add month comb column
  baseline_active_months <- lapply(X = monthly_mean, FUN = subset_base_months, 
                                   pred_points = PRED_sites, pred_points_sp = PRED_sites_sp, temp_threshold = temp_threshold[j]) %>%
    mapply(FUN = assign_month, month_number = month_number, SIMPLIFY = FALSE) %>%
    rbindlist() %>%
    group_by(id_col) %>%
    mutate(month_group = paste(month, collapse = "|")) %>%
    select(id_col, month_group) %>%
    unique()
  
  # derive the unique active months for combinations of baseline to calculate
  unique_month_combinations <- unique(baseline_active_months$month_group)
  
  # calculate the mean baseline for each potential combination of active months
  fin_baseline_mean <- parLapply(cl, X = unique_month_combinations, fun = calc_baseline, 
                                 func = mean, 
                                 pred_points = PRED_sites, 
                                 pred_points_sp = PRED_sites_sp,
                                 data_file = tmp1901_1930) %>%
    rbindlist() %>%
    inner_join(baseline_active_months, by = c("months" = "month_group", "id_col"))
  
  # calculate the sd baseline for each potential combination of active months
  fin_baseline_sd <- parLapply(cl, X = unique_month_combinations, fun = calc_baseline, 
                               func = stats::sd, 
                               pred_points = PRED_sites, 
                               pred_points_sp = PRED_sites_sp,
                               data_file = tmp1901_1930) %>%
    rbindlist() %>%
    inner_join(baseline_active_months, by = c("months" = "month_group", "id_col"))
  
  # calculate the mean temperatures for each predicts site, 11 months previously
  # set up empty list for each dataframe
  raster_means <- list()
  
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
    ind_raster_frame <- as.data.frame(PRED_coords) %>%
      inner_join(baseline_active_months, by = "id_col")
    
    # remove any columns for months that aren't active from the baseline for those coordinates
    frame_months <- ind_raster_frame %>%
      melt(id = c("id_col", "Longitude", "Latitude", "month_group"))
    
    # derive active month string for filtering months
    month_string <- paste("\\.", strsplit(frame_months$month_group[1], split = "|", fixed = TRUE)[[1]], "\\.", collapse = "|", sep = "")
    
    # some locations have no active months, if so add an NA row
    if(nrow(frame_months) == 0){
      
      # create dataframe for loations with no active months and rename columns
      raster_means[[i]] <- data.frame( "end_date" = names(tmp)[site_index], site_ids, "mean_value" = NA)
      colnames(raster_means[[i]]) <- c("end_date", "id_col", "x", "y", "mean_value")
      
      # skip onto next iteration
      next()
    }
    
    # recast dataframe after subsetting months for each PREDICTS site for the active months
    frame_months <- cbind(frame_months, month_string) %>%
      filter(grepl(month_string, variable )) %>% reshape2::dcast(id_col ~ variable)
    
    # remove the extra coordinate columns for calculating the row means
    ind_raster_values <- frame_months %>% select(-id_col)
    
    # calculate the mean values for each coordinate and bind back onto the coordinates
    raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% select(id_col, Longitude, Latitude)), rowMeans(ind_raster_values))
    colnames(raster_means[[i]]) <- c("end_date", "id_col", "x", "y", "mean_value")
    
  }
  
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
  
  print(j)
}
  
# stop the cluster
stopCluster(cl)

st2 <- Sys.time()
print(st2 - st1) # Time difference of 1.165101 hours on Charlie's laptop

# set up new lists for output
model_2c_abundance <- list()
abundance_model <- list()
main_plot_abundance <- list()

# run models for both species richness and total abundance
for(m in 1:length(predict_climate_list)){

  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ standard_anom * Predominant_land_use + (1|SS) + (1|SSB), data = predict_climate_list[[m]]) 
  
  # run predictions for the model of standard anomaly
  abundance_model[[m]] <- predict_continuous(model = model_2c_abundance[[m]],
                                             model_data = predict_climate_list[[m]],
                                             response_variable = "Total_abundance",
                                             categorical_variable = c("Predominant_land_use"),
                                             continuous_variable = c("standard_anom"),
                                             continuous_transformation = "",
                                             random_variable = c("SS", "SSB", "SSBS"))
  

  
}

# reformat data for plotting
reformat_sensitivity <- function(X, temp_threshold){
  
  data_fin <- X %>%
    dplyr::filter(standard_anom >= 0) %>%
    dplyr::filter(standard_anom <= 1) %>%
    mutate(temp_threshold = temp_threshold) %>%
    group_by(Predominant_land_use) %>%
    arrange(standard_anom) %>%
    mutate(difference_change = (exp(y_value[1]) - exp(tail(y_value, n = 1))) / exp(y_value[1])) %>%
    mutate(upper_conf_change = (exp(y_value_plus[1]) - exp(tail(y_value_minus, n = 1))) / exp(y_value_plus[1])) %>%
    mutate(lower_conf_change = (exp(y_value_minus[1]) - exp(tail(y_value_plus, n = 1))) / exp(y_value_minus[1])) %>%
    ungroup()
  
  return(data_fin)
}
 
# calculate difference in y value over 0-1 standard anomaly and plot
sensitive_object <- mapply(X = abundance_model, FUN = reformat_sensitivity, temp_threshold = temp_threshold, SIMPLIFY = FALSE) %>%
  rbindlist() %>%
  dplyr::select(Predominant_land_use, temp_threshold, difference_change, upper_conf_change, lower_conf_change) %>%
  unique()
  
# plot of change according to active season threshold
sensitive_object %>%
  filter(Predominant_land_use != "Primary vegetation") %>%
  mutate(Predominant_land_use = factor(Predominant_land_use, levels = c("Cropland"))) %>%
  ggplot() +
    geom_point(aes(x = temp_threshold, y = (difference_change * -1), colour = Predominant_land_use), position=position_dodge(width=1)) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    theme_bw() +
    scale_y_continuous(breaks = c(-0.5, -0.25, 0, 0.25, 0.5), labels = c("-50", "-25", "0", "25", "50")) +
    ylab("Mean abundance change (%)") +
    xlab("Active season temp. threshold (\u00B0C)") +
    scale_colour_manual("Land-use type", values = c( "#E69F00")) +
    theme(panel.grid = element_blank())

ggsave("pollinating-insects_active_season_sensitivity_3.png", scale = 0.9, dpi = 350)
