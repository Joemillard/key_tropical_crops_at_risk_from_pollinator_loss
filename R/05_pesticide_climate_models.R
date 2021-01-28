# script for climate anomaly - dataframe rewrite to try and avoid repeated raster trim/crops
### Note: 2005 is the mean year for insect data

# load required libraries
library(raster)
library(ggplot2)
library(dplyr)
library(data.table)
library(viridis)
library(yarg)
library(rworldmap) 
library(rworldxtra)
library(lme4)

# source in additional functions
source("R/00_functions.R")

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# read in the predicts pollinators
PREDICTS_pollinators <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds")

# PREDICTS data compilation
# filter for main pollinating taxa
PREDICTS_pollinators <- PREDICTS_pollinators %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  dplyr::filter(Class %in% c("Insecta")) %>%
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

# assign new variable for tropical/temperate, convert to factor, and filter out NA
#order.sites.div$zone <- ifelse(order.sites.div$Latitude >= -23.5 & order.sites.div$Latitude <= 23.5, "Tropics", "Temperate")
#order.sites.div$zone <- factor(order.sites.div$zone, levels = c("Temperate", "Tropics"))
#order.sites.div <- order.sites.div %>%
#  filter(zone == "Tropics") %>% 
#  droplevels()

# PREDICTS sites with the month of the recording
PRED_sites <- order.sites.div %>% select(id_col, Latitude, Longitude, Sample_end_latest) %>%
  mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest)) %>%
  filter(!is.na(Latitude))

# calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]

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

# plot the standardised climate anomaly
adjusted_climate %>%
  ggplot() +
  geom_point(aes(x = x, y = y, colour = standard_anom)) + 
  scale_colour_viridis()

# bind the adjusted climate data back onto the predicts sites
predicts_climate <- inner_join(order.sites.div, adjusted_climate, by = "id_col") #%>%
  #mutate(standard_anom  = standard_anom + 0.5)

# group the categories of climate anomaly into factors
predicts_climate$value_group[predicts_climate$standard_anom > 2] <- "> 2"
predicts_climate$value_group[predicts_climate$standard_anom > 1 & predicts_climate$standard_anom <= 2] <- "1 - 2"
predicts_climate$value_group[predicts_climate$standard_anom > 0.5 & predicts_climate$standard_anom <= 1] <- "0.5 - 1"
predicts_climate$value_group[predicts_climate$standard_anom > 0.25 & predicts_climate$standard_anom <= 0.5] <- "0.25 - 0.5"
predicts_climate$value_group[predicts_climate$standard_anom >= 0 & predicts_climate$standard_anom <= 0.25] <- "0 - 0.25"
predicts_climate$value_group[predicts_climate$standard_anom < 0] <- "< 0"

# order the levels of those factors
predicts_climate$value_group <- factor(predicts_climate$value_group, levels = c("> 2", "1 - 2", "0.5 - 1", "0.25 - 0.5", "0 - 0.25", "< 0"))

# bring in basemap for climate site plot 
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# plot the climate anom for each taxonomic order
predicts_climate %>%
  filter(!is.na(value_group)) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = x, y = y, colour = value_group), alpha = 0.5) + 
  scale_colour_manual("Standardised climate anomaly", values = c("#000000", "darkred", "#D55E00", "#E69F00", "#F0E442", "#56B4E9")) +
  coord_map(projection = "mollweide") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "bottom")

# save the plot for standardised climate anomaly by pollinating order
ggsave("stan-climate-anom_taxa.png", scale = 1.2, dpi = 350)

# read in the pesticide data
pesticide_low <- raster("outputs/Pesticide_totalAPR_Low_cropped.tif")
pesticide_high <- raster("outputs/Pesticide_totalAPR_High_cropped.tif")

# convert the predicts coordinates to spatial points
predicts_points <- predicts_climate %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Latitude)) %>%
  SpatialPoints()

# use the extract function to extract the pesticide values from both the low and hihg estimate
pesticide_low_points <- extract(pesticide_low, predicts_points)
pesticide_high_points <- extract(pesticide_high, predicts_points)

# build data frame of predicts climate data with the fertiliser data
climate_pest_predicts <- data.frame(predicts_climate, "low_estimate" = pesticide_low_points, "high_estimate" = pesticide_high_points) %>%
  dplyr::select(-Latitude.x, -Latitude.y, -Longitude.x, -Longitude.y)

# bring in basemap for climate site plot 
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# plot map for pesticide application
climate_pest_predicts %>%
  filter(!is.na(value_group)) %>%
  ggplot() +
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = x, y = y, colour = value_group, size = high_estimate), alpha = 0.3) + 
  scale_colour_manual("Standardised climate anomaly", values = c("#000000", "darkred", "#D55E00", "#E69F00", "#F0E442", "#56B4E9")) +
  scale_size_continuous("Pesticide application rate") +
  coord_map(projection = "mollweide") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        legend.position = "bottom")

# read in the fertiliser data 
fert_data <- raster("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/fertiliser_application_rate.tif")

# subset for unqiue sites
sites.sub_xy <- climate_pest_predicts %>%
  dplyr::select(Longitude, Latitude) %>%
  filter(!is.na(Longitude)) %>%
  filter(!is.na(Latitude)) %>%
  unique() %>%
  SpatialPoints()

# extract the fertislier values and join back onto the coordinates
sites.sub_xy$fert <- extract(fert_data, sites.sub_xy, na.rm = FALSE)

# turn the spatial points into a dataframe with the fertiliser data
fert_dat <- data.frame(coords = sites.sub_xy@coords, fert = sites.sub_xy@data$fert)

climate_pest_predicts <- inner_join(climate_pest_predicts, fert_dat, by = c("Longitude" = "coords.Longitude", "Latitude" = "coords.Latitude"))


## exploratory plots for relationship between climate, pesticide application, and taxa
climate_pest_predicts %>%
  ggplot() +
  geom_smooth(aes(x = high_estimate, y = log1p(Total_abundance), colour = value_group), method = "lm")

climate_pest_predicts %>%
  ggplot() +
  geom_smooth(aes(x = high_estimate, y = log1p(Total_abundance)), method = "lm")

climate_pest_predicts %>%
  ggplot() +
  geom_smooth(aes(x = standard_anom, y = log1p(Total_abundance)), method = "lm")

climate_pest_predicts %>%
  ggplot() +
  geom_smooth(aes(x = log10(fert), y = log1p(Total_abundance), colour = value_group), method = "lm")

climate_pest_predicts %>%
  filter(!is.na(value_group)) %>%
  ggplot() +
  geom_boxplot(aes(x = value_group, y = log10(Total_abundance)))

# add 1 for abundance and simpson diversity
climate_pest_predicts$Total_abundance <- climate_pest_predicts$Total_abundance + 1
climate_pest_predicts$Simpson_diversity <- climate_pest_predicts$Simpson_diversity + 1

# glmer plot for relationship between climate change and use intensity
# species richness, 2 continuous
model_1a <- glmer(Species_richness ~ log1p(standard_anom) * Use_intensity + (1|SS), data = climate_pest_predicts, family = poisson) 
model_1b <- glmer(Species_richness ~ log1p(standard_anom) * Use_intensity + (1|SS) + (1|SSB), data = climate_pest_predicts, family = poisson) 
model_1c <- glmer(Species_richness ~ log1p(standard_anom) * Use_intensity + (1|SS) + (1|SSB) + (1|SSBS), data = climate_pest_predicts, family = poisson) 

# check the AIC values
AIC(model_1a, model_1b, model_1c) # model_1c has the lowest AIC values

# species richness, standard anom as a factor
model_1c_1 <- glmer(Species_richness ~ log1p(standard_anom) * Use_intensity + (1|SS) + (1|SSB) + (1|SSBS), data = climate_pest_predicts, family = poisson) 
model_1c_2 <- glmer(Species_richness ~ log1p(standard_anom) + (1|SS) + (1|SSB) + (1|SSBS), data = climate_pest_predicts, family = poisson) 
model_1c_3 <- glmer(Species_richness ~ Use_intensity + (1|SS) + (1|SSB) + (1|SSBS), data = climate_pest_predicts, family = poisson) 
model_1c_4 <- glmer(Species_richness ~ 1 + (1|SS) + (1|SSB) + (1|SSBS), data = climate_pest_predicts, family = poisson) 


# check the AIC values
AIC(model_1c_1, model_1c_2, model_1c_3, model_1c_4) # model_1c_1 is the lowest again, but check unequal number of observations

# plot predicted values 
summary(model_1c_1)
anova(model_1c_1)

## total abundance, 2 continuous
model_2a <- lmer(log(Total_abundance) ~ log1p(standard_anom) * Use_intensity + (1|SS), data = climate_pest_predicts) 
model_2b <- lmer(log(Total_abundance) ~ log1p(standard_anom) * Use_intensity + (1|SS) + (1|SSB), data = climate_pest_predicts) 

# check the AIC values
AIC(model_2a, model_2b) # model_1c has the lowest AIC values

# species richness, standard anom as a factor
model_2c_1 <- lmerTest::lmer(log(Total_abundance) ~ log10(standard_anom + 1) * Predominant_land_use + (1|SS) + (1|SSB), data = climate_pest_predicts) 
model_2c_2 <- lmer(log(Total_abundance) ~ log1p(standard_anom) + (1|SS) + (1|SSB), data = climate_pest_predicts) 
model_2c_3 <- lmer(log(Total_abundance) ~ Use_intensity + (1|SS) + (1|SSB), data = climate_pest_predicts) 
model_2c_4 <- lmer(log(Total_abundance) ~ 1 + (1|SS) + (1|SSB), data = climate_pest_predicts) 

# check AIC values and summary
AIC(model_2c_1, model_2c_2, model_2c_3, model_2c_4)
summary(model_2c_1)
anova(model_2c_1)

# run predictions for the model of standard anomaly
abundance_model <- predict_continuous(model = model_2c_1,
                   model_data = climate_pest_predicts,
                   response_variable = "Total_abundance",
                   categorical_variable = c("Predominant_land_use"),
                   continuous_variable = c("standard_anom"),
                   continuous_transformation = log10,
                   random_variable = c("SS", "SSB", "SSBS"))

# plot for standardised anomaly and land-use for abundance
ggplot(abundance_model) +
  geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = standard_anom, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
  scale_fill_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  scale_colour_manual("Land-use type", values = c("#009E73", "#E69F00")) +
  xlab("Standardised climate anomaly") +
  ylab("Total abundance") +
  theme_bw() +
  theme(panel.grid = element_blank())

# save the insect pollinator anomaly plot
ggsave("insect_pollinator_anomaly.png", scale = 1, dpi = 350)
