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

# calc baseline (mean and sd)
calc_baseline <- function(data_file, func, pred_points, pred_points_sp){
  
  # calcualte either the mean or standard error for baseline, then extract points for predicts sites
  data_fin <- calc(data_file, func) %>%
    extract(pred_points_sp)
  
  # bind the extracted values back onto the predicts coordinates
  data_fin <- data.frame(pred_points[,1:3 ], data_fin)
  
  return(data_fin)
  
}

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds") %>%
  dplyr::select(-clade_rank, -confidence)  %>%
  filter(Class == "Insecta") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()

# create table of number of unique species
PREDICTS_pollinators_orig %>%
  select(Class, Order, Best_guess_binomial) %>%
  filter(Best_guess_binomial != "") %>%
  unique() %>%
  group_by(Class, Order) %>%
  tally()

# number of unique sites per land use
PREDICTS_pollinators_orig %>%
  select(Predominant_land_use, SSBS) %>%
  unique() %>%
  group_by(Predominant_land_use) %>%
  tally()

# filter the pollinators from the overall predicts database
PREDICTS_non_pollinating <- PREDICTS %>%
  filter(Class == "Insecta") %>%
  filter(!COL_ID %in% as.character(PREDICTS_pollinators_orig$COL_ID)) %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()

# calculate variation in sampling effect
PREDICTS_non_pollinating %>%
  group_by(SS) %>%
  summarise(sampling_variation = sd(Sampling_effort, na.rm = TRUE)) %>%
  mutate(varies = ifelse(sampling_variation > 0, 1, 0)) %>%
  mutate(rows = 1) %>%
  summarise(total_vary = sum(varies, na.rm = TRUE), total = sum(rows)) %>%
  mutate(percentage = (total_vary / total) * 100)

# calculate variation in sampling effect
PREDICTS_pollinators_orig %>%
  group_by(SS) %>%
  summarise(sampling_variation = sd(Sampling_effort, na.rm = TRUE)) %>%
  mutate(varies = ifelse(sampling_variation > 0, 1, 0)) %>%
  mutate(rows = 1) %>%
  summarise(total_vary = sum(varies, na.rm = TRUE), total = sum(rows)) %>%
  mutate(percentage = (total_vary / total) * 100)

# create table of number of unique species
PREDICTS_non_pollinating %>%
  select(Class, Order, Best_guess_binomial) %>%
  filter(Best_guess_binomial != "") %>%
  filter(Order != "") %>%
  unique() %>%
  group_by(Class, Order) %>%
  tally() %>%
  arrange(desc(n))

# number of unique sites per land use
PREDICTS_non_pollinating %>%
  select(Predominant_land_use, SSBS) %>%
  unique() %>%
  group_by(Predominant_land_use) %>%
  tally()

# bind together the two dataframes
pollinat_bound <- list(PREDICTS_pollinators_orig, PREDICTS_non_pollinating)

# set up vector for filtering for vertebrates and invertebrates
predict_climate_list <- list()

# loop through each phylum
for(j in 1:length(pollinat_bound)){
  
  # PREDICTS data compilation
  # filter for main pollinating taxa
  PREDICTS_pollinators <- pollinat_bound[[j]]
  
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
  
  # calculate the means and standard deviation for the beginning of the series
  # take names of values for 1901 to 1930
  tmp1901_1930 <- tmp[[names(tmp)[1:360]]]
  
  # extract the points for each the predicts coordinates
  PRED_sites_sp <- PRED_sites %>%
    select(Longitude, Latitude) %>%
    filter(!is.na(Latitude)) %>%
    SpatialPoints()
  
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
  
  # add 1 for abundance and simpson diversity
  predicts_climate$Total_abundance <- predicts_climate$Total_abundance + 1
  predicts_climate$Simpson_diversity <- predicts_climate$Simpson_diversity + 1
  
  # assign to list of predicts_climate and insects and vertebrates
  predict_climate_list[[j]] <- predicts_climate
  
}

# set up new lists for output
model_2c_abundance <- list()
abundance_model <- list()
main_plot_abundance <- list()

# run models for both species richness and total abundance
for(m in 1:length(pollinat_bound)){
  
  # species richness, standard anom as a factor
  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ standard_anom * Predominant_land_use + (1|SS), data = predict_climate_list[[m]]) 
  #print(AIC(model_2c_abundance[[m]]))
  
  # retrive model R squared values
  print(StatisticalModels::R2GLMER(model_2c_abundance[[m]]))
  
  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ standard_anom * Predominant_land_use + (1|SS) + (1|SSB), data = predict_climate_list[[m]]) 
  #print(AIC(model_2c_abundance[[m]]))
  
  # retrive model R squared values
  print(StatisticalModels::R2GLMER(model_2c_abundance[[m]]))
  
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
   xlab("Standardised temperature anomaly") +
   ylab("Total abundance") +
   theme_bw() +
   theme(panel.grid = element_blank())
  
}

# predict abundance at 0 warming on cropland
zero_data <- data.frame("standard_anom" = 0, Predominant_land_use = "Cropland")
zero_warming_abundance <- predict(model_2c_abundance[[1]], zero_data, re.form = NA)
zero_warming_abundance <- exp(zero_warming_abundance)

# predict abundance at 1 STA warming on cropland
one_data <- data.frame("standard_anom" = 1, Predominant_land_use = "Cropland")
one_warming_abundance <- predict(model_2c_abundance[[1]], one_data, re.form = NA)
one_warming_abundance <- exp(one_warming_abundance)

# calculate percentage change 
percentage_change <- (zero_warming_abundance - one_warming_abundance) / zero_warming_abundance

# write percentage change to disk for active season
write.csv(percentage_change, "outputs/no_active_season_change.csv")

# plot for the pollinating insects and non-pollinating insects - climate anomaly of 4 corresponds to ~100% abundance loss
plot_grid(main_plot_abundance[[1]] +
            ggtitle("(A) Pollinating insects") + 
            scale_y_continuous(limits = c(1, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
            scale_x_continuous(limits = c(-0.6, 2.88)) +
            theme(legend.position = "bottom"), main_plot_abundance[[2]] + 
            ggtitle("(B) Non-pollinating insects") +
            scale_y_continuous(limits = c(1, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
            scale_x_continuous(limits = c(-0.6, 2.88)) +
            guides(color = guide_legend(override.aes = list(color = NA)), 
                   fill = guide_legend(override.aes = list(fill = NA))) +
            theme(legend.position = "bottom", legend.key = element_rect(fill = "white"), legend.text = element_text(color = "white"), legend.title = element_text(color = "white")), ncol = 2) 

ggsave("pollinating_non-pollinating_7.png", scale = 1.1, dpi = 350)

# subset data for just 4 rows, with 0 and 1 for cropland and primary vegetation
new_data <- predict_climate_list[[m]] %>%
  group_by(Predominant_land_use) %>%
  slice(0:2) %>%
  ungroup() %>%
  mutate("standard_anom" = c(0, 1, 0, 1))

# difference between cropland and primary, between standardised temperature anomaly of 0 and 1
standard_anom_1_prediction <- predict_continuous(model = model_2c_abundance[[1]],
                   model_data = new_data,
                   response_variable = "Total_abundance",
                   categorical_variable = c("Predominant_land_use"),
                   continuous_variable = c("standard_anom"),
                   continuous_transformation = "",
                   random_variable = c("SS", "SSB", "SSBS")) %>%
  mutate(y_value = exp(y_value)) %>%
  select(Predominant_land_use, standard_anom, y_value)

# calculate percentage change
(standard_anom_1_prediction$y_value[1] - standard_anom_1_prediction$y_value[4]) / standard_anom_1_prediction$y_value[1]
(standard_anom_1_prediction$y_value[4] / standard_anom_1_prediction$y_value[1])  * 100

