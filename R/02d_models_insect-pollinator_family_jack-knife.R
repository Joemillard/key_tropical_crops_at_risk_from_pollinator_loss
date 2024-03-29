# script for pollinating/non-pollinating insects, either split into each insect order or jack-knifed (remove with replacement)

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
  dplyr::filter(Phylum %in% "Arthropoda") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels()

# set up vector for taxonomic orders
taxonomic_families <- c("", PREDICTS_pollinators_orig$Family %>% unique())
PREDICTS_pollinators_taxa <- list()

# loop through the taxonomic orders and remove each for jack-knife with replacement/select one for single order models
for(i in 1:length(taxonomic_families)){
  PREDICTS_pollinators_taxa[[i]] <- PREDICTS_pollinators_orig %>% 
    filter(Family != taxonomic_families[i]) %>% # amend to either jack-knife or select one order
    droplevels()
}

# set up vector for filtering for vertebrates and invertebrates
predict_climate_list <- list()

# loop through each phylum
for(j in 1:length(PREDICTS_pollinators_taxa)){
  
  # PREDICTS data compilation
  # filter for main pollinating taxa
  PREDICTS_pollinators <- PREDICTS_pollinators_taxa[[j]]
  
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
  
  print(j)
  
}

# set up new lists for output
model_2c_abundance <- list()
abundance_model <- list()
main_plot_abundance <- list()

# run models for both species richness and total abundance
for(m in 1:length(PREDICTS_pollinators_taxa)){
  
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
  
}

for(i in 1:length(abundance_model)){
  abundance_model[[i]]$taxa_jack <- taxonomic_families[i]
}

rbindlist(abundance_model) %>%
  filter(Predominant_land_use == "Cropland") %>%
  mutate(trend_type = ifelse(taxa_jack == "", "main", "jacked")) %>%
  ggplot() +
  geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use, group = taxa_jack, alpha = trend_type), size = 1) +
  scale_colour_manual("Land-use type", values = c("#E69F00")) +
  scale_alpha_discrete(range=c(0.15, 1)) +
  scale_linetype_discrete(c("solid", "dashed")) +
  coord_cartesian(xlim = c(0, 1.65), ylim = c(1.609438, 4.3), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.5, 1, 1.5), labels = c("0", "0.5", "1", "1.5")) +
  scale_y_continuous(breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321, 6.461468), labels = c(5, 10, 20, 40, 80, 160, 320, 640)) +
  xlab("Standardised temperature anomaly") +
  ylab("Total abundance") +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none")


ggsave("jack_knifed_crop_change.png", dpi = 350, scale  = 0.8)
