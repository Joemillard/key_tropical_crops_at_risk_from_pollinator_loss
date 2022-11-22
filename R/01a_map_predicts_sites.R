# map of site distribution

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
library(patchwork)
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

# build map
base_map <- get_basemap()

# fortify the main map
map_fort <- fortify(base_map)

# build map for distribution of sites
site_distribution_pollinator <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, fill = Predominant_land_use), data = predict_climate_list[[1]], shape = 21, colour = "black", alpha = 0.3) +
  scale_fill_manual("Predominant land-use", values = c("#009E73", "#E69F00"), guide = guide_legend(override.aes = list(fill = "white", colour = NA), ncol = 1)) +
  coord_map(projection = "mollweide") +
  labs(subtitle = "(A) Pollinators") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),
        legend.key=element_blank(), legend.position = "bottom", 
        legend.justification = "left")
  
# build map for distribution of sites
site_distribution_non_pollinator <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, fill = Predominant_land_use), data = predict_climate_list[[2]], shape = 21, colour = "black", alpha = 0.3) +
  scale_fill_manual("Predominant land-use", values = c("#009E73", "#E69F00"), guide = guide_legend(override.aes = list(alpha = 1, colour = NA), ncol = 1)) +
  coord_map(projection = "mollweide") +
  labs(subtitle = "(B) Non pollinators") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"), legend.position = "bottom", 
        legend.justification = "left")


# build map for distribution of sites
anomaly_pollinator <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, colour = standard_anom), data = predict_climate_list[[1]], alpha = 0.3) +
  scale_colour_viridis("Std. temp.\nanomaly", option = "B", guide = guide_legend(override.aes = list(fill = "white", colour = NA))) +
  coord_map(projection = "mollweide") +
  labs(subtitle = "(C) Pollinators") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"),
        legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),
        legend.key=element_blank(), legend.position = "bottom", 
        legend.justification = "left")

# build map for distribution of sites
anomaly_non_pollinator <- ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), data = map_fort, fill = "lightgrey") +
  geom_point(aes(x = Longitude, y = Latitude, colour = standard_anom), data = predict_climate_list[[2]], alpha = 0.3) +
  scale_colour_viridis("Std. temp.\nanomaly", option = "B") +
  coord_map(projection = "mollweide") +
  labs(subtitle = "(D) Non pollinators") +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        text = element_text(size = 13),
        panel.background = element_rect(fill = "white"), legend.position = "bottom", 
        legend.justification = "left")

# combine two maps for land use type of pollinators and non-pollinators
site_distribution_pollinator + site_distribution_non_pollinator + anomaly_pollinator + anomaly_non_pollinator + plot_layout(ncol = 2)



ggsave("pollinator_distribution_maps.png", scale = 1.3, dpi = 350)
