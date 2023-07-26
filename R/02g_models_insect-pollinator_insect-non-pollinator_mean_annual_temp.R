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
library(patchwork)

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

# PREDICTS data compilation
# filter for main pollinating taxa
PREDICTS_pollinators <- PREDICTS_pollinators_orig

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

# rescale abundance
RescaleAbundance <- function(sites){
  
  sites <- droplevels(sites)
  
  StudyMaxAbund <- suppressWarnings(tapply(
    X = sites$Total_abundance,INDEX = sites$SS,
    FUN = max,na.rm=TRUE))
  
  StudyMaxAbund[StudyMaxAbund == -Inf] <- NA
  
  AllMaxAbund <- StudyMaxAbund[match(sites$SS,names(StudyMaxAbund))]
  
  sites$Total_abundance_RS <- sites$Total_abundance/AllMaxAbund
  
  sites$LogAbund <- log(sites$Total_abundance_RS+0.01)
  
  return(sites)
}

predicts_climate <- RescaleAbundance(predicts_climate)

# plot rescaled abundance in relation to standardised climate anomaly on primary vegetation and cropland
predicts_climate %>%
  ggplot() +
    geom_point(aes(x = standard_anom, y = log(Total_abundance_RS), colour = Predominant_land_use)) +
    geom_smooth(aes(x = standard_anom, y = log(Total_abundance_RS), colour = Predominant_land_use), method = "lm") +
    facet_wrap(~Predominant_land_use)



# add 1 for abundance and simpson diversity
predicts_climate$Total_abundance <- predicts_climate$Total_abundance + 1
predicts_climate$Simpson_diversity <- predicts_climate$Simpson_diversity + 1

# correlation of std anomaly and mean temperature on cropland and primary vegetation
predicts_climate %>%
  ggplot() +
    geom_point(aes(x = mean_value, y = standard_anom, colour = Predominant_land_use)) +
    facet_wrap(~Predominant_land_use) +
    scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
    ylab("Standardised temperature anomaly") +
    xlab("Mean annual temperature") +
    theme_bw() +
    theme(panel.grid = element_blank())

ggsave("correlation_mean_temp_stf_anomaly.png", scale = 1, dpi = 350)

# correlation of std anomaly and mean temperature on cropland and primary vegetation
predicts_climate %>%
  ggplot() +
  geom_point(aes(x = mean_value, y = anomaly, colour = Predominant_land_use)) +
  facet_wrap(~Predominant_land_use) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  ylab("Temperature anomaly") +
  xlab("Mean annual temperature") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave("correlation_mean_temp_anomaly.png", scale = 1, dpi = 350)

# resample sites to get distribution of anomaly and mean temp within lowest random effect nesting
SS_vec <- unique(predicts_climate$SS)
land_use_vec <- unique(predicts_climate$Predominant_land_use)

climate_samples <- list()
land_use_samples <- list()

# iterate through each of the land use classes
for(j in 1:length(land_use_vec)){
  # for each of the SSB, calculate the mean for cropland and primary for mean temp and anomaly
  for(i in 1:length(SS_vec)){
    climate_samples[[i]] <- predicts_climate %>%
      filter(Predominant_land_use == land_use_vec[j]) %>%
      filter(SS %in% SS_vec[i]) %>%
      droplevels() %>%
      group_by(Predominant_land_use) %>%
      summarise(mean_sample = sd(mean_value, na.rm = TRUE))
    
  }
  
  land_use_samples[[j]] <- climate_samples %>%
    rbindlist()
  
}

rbindlist(land_use_samples) %>%
  ggplot() +
    geom_boxplot(aes(y = log(mean_sample), x = Predominant_land_use, fill = Predominant_land_use, colour = Predominant_land_use), alpha = 0.2) +
      scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
      scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
      scale_y_continuous("Std. dev. of SS mean annual temp.", breaks = c(-3.465736, -2.7725887, -2.0794415, -1.3862944, -0.6931472, 0), 
                         labels = c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1)) +
      xlab("") +
      theme_bw() +
      theme(panel.grid = element_blank())

# save plot for mean annual temp variation
ggsave("Mean_temp_std_dev.png", scale = 1, dpi = 350)

# bring in basemap for climate site plot 
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# plot the climate anom for each taxonomic order
predicts_climate %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = x, y = y, colour = mean_value), alpha = 0.5) + 
  facet_wrap(~Predominant_land_use) +
  scale_colour_viridis("Mean annual temperature", option = "B") +
  coord_map(projection = "mollweide") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "bottom")

# plot the climate anom for each taxonomic order
predicts_climate %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = x, y = y, colour = anomaly), alpha = 0.5) + 
  facet_wrap(~Predominant_land_use) +
  scale_colour_viridis("Mean annual temperature", option = "B") +
  coord_map(projection = "mollweide") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "bottom")

ggsave("mean_annual_temp_map.png", scale = 1, dpi = 350)

# what's the distribution of cropland and primary vegetation annual mean temperatures
mean_temp_plot <- predicts_climate %>%
  ggplot() + 
  geom_density(aes(x = mean_value, fill = Predominant_land_use, colour = Predominant_land_use), alpha = 0.2) +
    scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
    scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
    scale_y_continuous("Density", expand = c(0, 0), limits = c(0, 6.1)) +
    scale_x_continuous("Mean annual temperature",  expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")

std_anomaly_plot <- predicts_climate %>%
  ggplot() + 
  geom_density(aes(x = standard_anom, fill = Predominant_land_use, colour = Predominant_land_use), alpha = 0.2) +
  scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  scale_y_continuous("Density", expand = c(0, 0), limits = c(0, 6.1)) +
  scale_x_continuous("Standardised temperature anomaly", expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")

# save plots of dist of temp and anomaly
mean_temp_plot + std_anomaly_plot

ggsave("density_temp.png", scale = 1, dpi = 350)

# centre the predictors
predicts_climate$mean_value <- StdCenterPredictor(predicts_climate$mean_value)
predicts_climate$standard_anom <- StdCenterPredictor(predicts_climate$standard_anom)

model_abundance_mean_val <- lmerTest::lmer(log(Total_abundance_RS) ~ mean_value * Predominant_land_use + (1|SS) + (1|SSB), data = predicts_climate) 
model_abundance_anom <- lmerTest::lmer(log(Total_abundance_RS) ~ standard_anom * Predominant_land_use + (1|SS) + (1|SSB), data = predicts_climate) 
model_abundance_anom_non_stan <- lmerTest::lmer(log(Total_abundance_RS) ~ anomaly * Predominant_land_use + (1|SS) + (1|SSB), data = predicts_climate) 

AIC(model_abundance_mean_val, model_abundance_anom, model_abundance_anom_non_stan)

# run predictions for the model of standard anomaly
abundance_model <- predict_continuous(model = model_abundance_mean_val,
                                      model_data = predicts_climate,
                                      response_variable = "Total_abundance_RS",
                                      categorical_variable = c("Predominant_land_use"),
                                      continuous_variable = c("mean_value"),
                                      continuous_transformation = "",
                                      random_variable = c("SS", "SSB", "SSBS"))

# plot for standardised anomaly and land-use for abundance
mean_plot_abundance <- ggplot(abundance_model) +
  geom_line(aes(x = mean_value, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = mean_value, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
  scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
 # scale_y_continuous(limits = c(0.3, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  xlab("Mean temperature") +
  ylab("Total abundance") +
  theme_bw() +
  theme(panel.grid = element_blank())

# save plot for mean temperature
ggsave("mean_temp_plot_prediction.png", scale = 1, dpi = 350)

# run predictions for the model of standard anomaly
abundance_model <- predict_continuous(model = model_abundance_anom,
                                      model_data = predicts_climate,
                                      response_variable = "Total_abundance_RS",
                                      categorical_variable = c("Predominant_land_use"),
                                      continuous_variable = c("standard_anom"),
                                      continuous_transformation = "",
                                      random_variable = c("SS", "SSB", "SSBS"))

# plot for standardised anomaly and land-use for abundance
std_plot_abundance <- ggplot(abundance_model) +
  geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = standard_anom, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
  scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  #scale_y_continuous(limits = c(0.3, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  xlab("Std. temperature anomaly") +
  ylab("Total abundance") +
  theme_bw() +
  theme(panel.grid = element_blank())

# save plot for mean temperature
ggsave("std_anomaly_plot_prediction.png", scale = 1, dpi = 350)

# filter for just primary
predicts_climate_primary <- predicts_climate %>%
  filter(Predominant_land_use == "Primary vegetation")

# species richness, standard anom as a factor
pri_model_abundance_mean_anom <- lmerTest::lmer(log(Total_abundance_RS) ~ mean_value + standard_anom + (1|SS) + (1|SSB), data = predicts_climate_primary)
pri_model_abundance_mean <- lmerTest::lmer(log(Total_abundance_RS) ~ mean_value  + (1|SS) + (1|SSB), data = predicts_climate_primary) 
pri_model_abundance_anom <- lmerTest::lmer(log(Total_abundance_RS) ~ standard_anom + (1|SS) + (1|SSB), data = predicts_climate_primary) 
pri_model_abundance_anom_non_stan <- lmerTest::lmer(log(Total_abundance_RS) ~ anomaly + (1|SS) + (1|SSB), data = predicts_climate_primary) 

# extract model summaries
summary(pri_model_abundance_mean_anom)
summary(pri_model_abundance_mean)
summary(pri_model_abundance_anom)
summary(pri_model_abundance_anom_non_stan)

AIC(pri_model_abundance_mean_anom, pri_model_abundance_mean, pri_model_abundance_anom, pri_model_abundance_anom_non_stan)

# filter for just primary
predicts_climate_cropland <- predicts_climate %>%
  filter(Predominant_land_use == "Cropland")

# species richness, standard anom as a factor
cro_model_abundance_mean_anom <- lmerTest::lmer(log(Total_abundance_RS) ~ mean_value + standard_anom + (1|SS) + (1|SSB), data = predicts_climate_cropland)
cro_model_abundance_mean <- lmerTest::lmer(log(Total_abundance_RS) ~ mean_value  + (1|SS) + (1|SSB), data = predicts_climate_cropland) 
cro_model_abundance_anom <- lmerTest::lmer(log(Total_abundance_RS) ~ standard_anom + (1|SS) + (1|SSB), data = predicts_climate_cropland) 
cro_model_abundance_anom_non_stan <- lmerTest::lmer(log(Total_abundance_RS) ~ anomaly + (1|SS) + (1|SSB), data = predicts_climate_cropland) 

# extract model summaries
summary(cro_model_abundance_mean_anom)
summary(cro_model_abundance_mean)
summary(cro_model_abundance_anom)
summary(cro_model_abundance_anom_non_stan)

AIC(cro_model_abundance_mean_anom, cro_model_abundance_mean, cro_model_abundance_anom, cro_model_abundance_anom_non_stan)

# run predictions for the model of standard anomaly
abundance_model <- predict_continuous(model = model_abundance_mean_val,
                                           model_data = predicts_climate,
                                           response_variable = "Total_abundance",
                                           categorical_variable = c("Predominant_land_use"),
                                           continuous_variable = c("mean_value"),
                                           continuous_transformation = "",
                                           random_variable = c("SS", "SSB", "SSBS"))

# plot for standardised anomaly and land-use for abundance
main_plot_abundance <- ggplot(abundance_model) +
  geom_line(aes(x = mean_value, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = mean_value, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
  scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  scale_y_continuous(limits = c(0.3, 6.5), breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  xlab("Mean temperature") +
  ylab("Total abundance") +
  theme_bw() +
  theme(panel.grid = element_blank())

# code that randomsly samples sliding numbers from cropland and primary vegetation, and then checks how delta AIC changes
# hypothesis being that as the ratio moves towards primary vegetation, mean temperature better, explains, and then towards cropland standardised anomaly better explains
predicts_climate %>%
  group_by(Predominant_land_use) %>%
  tally()

# set up list for random filter values
filter_list <- list(c(100, 1100), 
     c(200, 1000), 
     c(300, 900), 
     c(400, 800), 
     c(500, 700),
     c(600, 600), 
     c(700, 500),
     c(800, 400),
     c(900, 300),
     c(1000, 200),
     c(1100, 100))

delta_AIC_list <- list()

# iterate through each iteration for the number randomly sampled
for(i in 1:length(filter_list)){
  
  # vector for delta AIC values
  delta_AIC <- c()
  
  for(j in 1:200){
    
    # filter for cropland and randomly sample rows for that simulation
    predicts_climate_cropland <- predicts_climate %>%
      filter(Predominant_land_use == "Cropland") %>%
      sample_n(filter_list[[i]][1]) %>%
      droplevels()
    
    # filter for primary vegetation and randomly sample rows for that simulation
    predicts_climate_prim <- predicts_climate %>%
      filter(Predominant_land_use == "Primary vegetation") %>%
      sample_n(filter_list[[i]][2]) %>%
      droplevels()
    
    # bind together the random samples for cropland and primary vegetation
    predicts_climate_sampled <- rbind(predicts_climate_cropland, predicts_climate_prim)
    
    # build models for both mean value and standardised temperature anomaly
    model_abundance_mean_val<- lmerTest::lmer(log(Total_abundance) ~ Predominant_land_use * mean_value + (1|SS) + (1|SSB), data = predicts_climate_sampled)
    model_abundance_stan_anom <- lmerTest::lmer(log(Total_abundance) ~ Predominant_land_use * standard_anom + (1|SS) + (1|SSB), data = predicts_climate_sampled)
    
    # calculate delta AIC per simulation
    delta_AIC[j] <- AIC(model_abundance_mean_val) - AIC(model_abundance_stan_anom)
    
  }
  
  delta_AIC_list[[i]] <- data.frame("delta_AIC" = delta_AIC, 
                                    "Cropland_sites" = filter_list[[i]][[1]], 
                                    "Primary_sites" = filter_list[[i]][[2]])
  
}

# bind together the simulations and calculate median delta AIC
rbindlist(delta_AIC_list) %>%
  group_by(Cropland_sites, Primary_sites) %>%
  summarise(lower = quantile(delta_AIC, 0.025), median_delta_AIC = median(delta_AIC), upper = quantile(delta_AIC, 0.975)) %>%
  ggplot() +
    geom_point(aes(x = Cropland_sites/(Cropland_sites+Primary_sites), y = median_delta_AIC)) +
    geom_errorbar(aes(x = Cropland_sites/(Cropland_sites+Primary_sites), ymin = lower, ymax = upper)) +
    geom_hline(yintercept = 0, alpha = 0.2) + 
    scale_x_continuous("Cropland site proportion", breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
    ylab("Delta AIC (Mean - Std. anomaly)") +
    theme_bw() +
    theme(panel.grid = element_blank())
 
ggsave("resampled_delta_AIC.png", scale = 1, dpi = 350)

# test for effect within one land use class and intensiy
predicts_cropland_light_use <- predicts_climate %>%
  filter(Predominant_land_use %in% "Cropland") %>%
  filter(Use_intensity == "Light use")

model_abundance_mean_val_light <- lmerTest::lmer(log(Total_abundance) ~ mean_value + (1|SS) + (1|SSB), data = predicts_cropland_light_use)
model_abundance_stan_anom_light <- lmerTest::lmer(log(Total_abundance) ~ standard_anom + (1|SS) + (1|SSB), data = predicts_cropland_light_use)

AIC(model_abundance_mean_val_light, model_abundance_stan_anom_light)

# test for effect in one sampling method
predicts_climate <- predicts_climate %>%
  filter()

model_abundance_mean_val_light <- lmerTest::lmer(log(Total_abundance) ~ mean_value + (1|SS) + (1|SSB), data = predicts_climate)
model_abundance_stan_anom_light <- lmerTest::lmer(log(Total_abundance) ~ standard_anom + (1|SS) + (1|SSB), data = predicts_climate)

