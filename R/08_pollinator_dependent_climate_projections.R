# read in packages
library(raster)
library(ggplot2)
library(rasterVis)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)
library(cowplot)
library(data.table)
library(lme4)
library(yarg)
library(cowplot)

# source in additional functions
source("R/00_functions.R")

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("G:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds")

# set up the starting directory for future climate data
SSP_directory <- ("G:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies")

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
# take names of values for 1901 to 1931
tmp1901_1931 <- tmp[[names(tmp)[1:361]]]

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
climate_start_mean <- calc_baseline(tmp1901_1931, 
                                    func = mean, 
                                    pred_points = PRED_sites, 
                                    pred_points_sp = PRED_sites_sp) %>%
  mutate(Latitude = as.character(Latitude)) %>%
  mutate(Longitude = as.character(Longitude))

# calculate the sd baseline, and convert to character for merging
climate_start_sd <- calc_baseline(tmp1901_1931, 
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
main_plot_abundance <- ggplot(abundance_model) +
  geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = standard_anom, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
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

# subset the file paths for just those that are pollination dependent to some extent
# subset as strings to filter from klein_cleaned
pollinated_crops <- grep(paste(unique(paste("/", klein_cleaned$MonfredaCrop, "_", sep = "")), collapse = "|"), unlisted_crops, value = TRUE)
pollinat_crops_simp <- gsub("G:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff/", "", pollinated_crops)
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

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
crop.total <- stack(rate_rasters_adj) %>% sum(na.rm = T)

# calculate total pollination dependent production
total_production <- sum(crop.total[])

## standardised climate anomaly script
# take names of values for 1901 to 1931 - 30 year baseline
tmp1901_1931 <- tmp[[names(tmp)[1:361]]]

# calculate the mean and sd of the baseline values
tmp1901_1931mean <- calc(tmp1901_1931, mean)
tmp1901_1931sd <- calc(tmp1901_1931, stats::sd)

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
RCP_scenarios <- c("rcp85", "rcp60", "rcp26")
climate_model_combs <- c("GFDL","HadGEM2","IPSL","MIROC5")
climate_model_combs_adj <- c()
for(i in 1:length(climate_model_combs)){
  climate_model_combs_adj[i] <- paste(climate_model_combs[-i], collapse="|")
}

# add in the full model set for single vector of all climate models
climate_model_combs_adj <- append("GFDL|HadGEM2|IPSL|MIROC5", climate_model_combs_adj, after = 1)

# set up list for each RCP plot
RCP_plot <- list()

# iterate through each RCP scenario
for(k in 1:length(RCP_scenarios)){
  
  print(RCP_scenarios[k])
  
  # set up empty list for climate model combinations
  vulnerable_production_jack <- list()
  
  # iterate through each set of climate models to jack-knife by climate model  
  for(j in 1:length(climate_model_combs_adj)){   
    
    print(climate_model_combs_adj[j])
    
    # iterate through each set of years as a rolling average
    for(i in 1:length(years_list)){
      
      # file path for ISIMIP data
      all.files <- dir(path = SSP_directory,recursive = TRUE, full.names = TRUE)
      
      # using RCP 8.5 calculate average of separate models
      mean.temp.2069.2071 <- stack(lapply(X = years_list[[i]], FUN = average_clim_models, RCP = RCP_scenarios[k], clim_models = climate_model_combs_adj[j]))
      
      mean.temp.2069.2071 <- stackApply(x = mean.temp.2069.2071,indices = rep(1,3), fun = mean)
      
      # calc the anomalies for the future years
      tmp2069_71_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1931mean)
      tmp2069_71std_climate_anomaly[[i]] <- (mean.temp.2069.2071-tmp1901_1931mean) / tmp1901_1931sd
      
    }
    
    # set up empty list for standardised climate anomaly
    std_high_abun_adj <- list()
    std_anom_high <- list()
    
    # predict abundance at 0 warming on cropland
    zero_data <- data.frame("standard_anom" = 0, Predominant_land_use = "Cropland")
    zero_warming_abundance <- predict(model_2c_abundance, zero_data, re.form = NA)
    
    # for each set of climate anomaly data, predict abundance reduction for all climate anomaly values in each cell
    # and then sum abundance adjusted pollination dependence
    for(i in 1:length(tmp2069_71std_climate_anomaly)){
      
      # convert the raster to a dataframe to plot with ggplot
      std_high_abun_adj[[i]] <- as(tmp2069_71std_climate_anomaly[[i]], "SpatialPixelsDataFrame")
      std_high_abun_adj[[i]] <- as.data.frame(std_high_abun_adj[[i]])
      
      # set up prediction data on basis of that set of years
      new_data_pred <- data.frame("standard_anom" = std_high_abun_adj[[i]]$layer, Predominant_land_use = "Cropland")
      
      # predict abundance for climate anomaly and and to data frame
      predicted_abundance <- predict(model_2c_abundance, new_data_pred, re.form = NA)
      std_high_abun_adj[[i]]$abundance <- predicted_abundance
      
      # calculate percentage change from place with 0 warming, and convert to vulnerability
      std_high_abun_adj[[i]]$abundance_change <- 1 - (std_high_abun_adj[[i]]$abundance / zero_warming_abundance)
      
      # convert spatial dataframe to coordinates
      std_anom_high[[i]] <- std_high_abun_adj[[i]] %>%
        dplyr::select(x, y) %>%
        unique() %>%
        SpatialPoints()
      
    }
    
    # set up vector for total production
    vulnerable_production_list <- list()
    vulnerable_production <- c()
    
    # for each set of coordinates, extract the pollination dependent values and sum
    for(i in 1:length(std_anom_high)){
      
      # convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
      vulnerable_production_list[[i]] <- extract(crop.total, std_anom_high[[i]], na.rm = FALSE)
      vulnerable_production[i] <- unlist(vulnerable_production_list[[i]] * std_high_abun_adj[[i]]$abundance_change) %>% sum()
      vulnerable_production_jack[[j]] <- data.frame("vulnerability" = vulnerable_production, 
                                                    "model" = climate_model_combs_adj[j]) 
      
    }
  }
  
  # create dataframe for exposed production and build datafrmae
  RCP_plot[[k]] <- rbindlist(vulnerable_production_jack) %>%
    mutate(year = rep(c(seq(2048, 2016, -1)), 5)) %>%
    mutate(scenario = RCP_scenarios[k])
}

# bind together the outputs and plot as facetted plot for each scenario
rbindlist(RCP_plot) %>%
  mutate(model = factor(model, levels = c("GFDL|HadGEM2|IPSL|MIROC5", "HadGEM2|IPSL|MIROC5", "GFDL|IPSL|MIROC5", "GFDL|HadGEM2|MIROC5", "GFDL|HadGEM2|IPSL"),
                        labels = c("All 4 models", "Excluding GFDL", "Excluding HadGEM2", "Excluding IPSL", "Excluding MIROC5"))) %>%
  mutate(scenario = factor(scenario, levels = c("rcp85", "rcp60", "rcp26"),
                           labels = c("RCP 8.5", "RCP 6.0", "RCP 2.6"))) %>%
  ggplot() +
    geom_line(aes(x = year, y = vulnerability, colour = model, alpha = model)) +
    geom_point(aes(x = year, y = vulnerability, colour = model, alpha = model)) +
    facet_wrap(~scenario, ncol = 2) +
    scale_y_continuous(limits = c(700000, 2500000), expand = c(0, 0), breaks = c(1000000, 1500000, 2000000, 2500000), labels = c("1,000,000", "1,500,000", "2,000,000", "2,500,000")) +
    scale_x_continuous(limits = c(2015, 2050), expand = c(0, 0), breaks = c(2020, 2025, 2030, 2035, 2040, 2045, 2050)) +
    scale_colour_manual("Climate model", values = c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
    scale_alpha_manual("Climate model", values = c(1, 0.4, 0.4, 0.4, 0.4)) +
    ylab("Vulnerability weighted pollination dependent prod. (mt tonnes)") +
    xlab("") +
    theme_bw() +
    theme(panel.grid = element_blank())
  
# save facetted plot
ggsave("rcp_85_pollination_exposure_2.png", scale = 1, dpi = 350)
