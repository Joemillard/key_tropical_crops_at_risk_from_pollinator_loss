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
#library(snow)


# read in the original predicts database 
PREDICTS <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/Data/PREDICTS/database.rds") %>%
  mutate(COL_ID = as.character(COL_ID))

# source in additional functions
source("R/00_functions.R")

# # calc baseline (mean and sd)
# calc_baseline <- function(data_file, func, pred_points, pred_points_sp){
#   
#   # calcualte either the mean or standard error for baseline, then extract points for predicts sites
#   data_fin <- calc(data_file, func) %>%
#     extract(pred_points_sp)
#   
#   # bind the extracted values back onto the predicts coordinates
#   data_fin <- data.frame(pred_points[,1:3 ], data_fin)
#   
#   return(data_fin)
#   
# }

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


# threshold temperature value for insect active months
thresh <- 10

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
  
  # CO: I have commented this out as the month in the original format is used in the 
  # code I have taken from the insects paper.
  # # PREDICTS sites with the month of the recording
  # PRED_sites <- order.sites.div %>% select(id_col, Latitude, Longitude, Sample_end_latest) %>%
  #   mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  #   mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest)) %>%
  #   filter(!is.na(Latitude))
  
  PRED_sites <- order.sites.div
  
  # calculate the means and standard deviation for the beginning of the series
  # take names of values for 1901 to 1931
  # CO: edited values taken to be those for 1901-1930 (30yr baseline)
  tmp1901_1930 <- tmp[[names(tmp)[1:360]]]
  
  # # extract the points for each the predicts coordinates
  # predicts_sp <- PRED_sites %>%
  #   select(Longitude, Latitude) %>%
  #   filter(!is.na(Latitude)) %>%
  #   SpatialPoints()
  

  #### Added by CO: code to determine anomaly values for each site ####
  wgs84 <- crs(tmp)
  
  predicts_sp <- SpatialPointsDataFrame(
    coords = cbind(PRED_sites$Longitude, PRED_sites$Latitude), 
    data = PRED_sites, proj4string = wgs84)
  
  ##Spatial points for rasterizing
  SP <- SpatialPoints(predicts_sp, proj4string=wgs84)
  
  ##Names of tmp layer, needed for subsettting
  names_tmp <- names(tmp)
  
  ##Create a list of a all names of tmp layers, that will be used for matching later on
  names_sub <- substr(names_tmp, 2, 8) 
  
  nCores <- parallel::detectCores()
  
  st1 <- Sys.time()
  
  cl <- snow::makeCluster(nCores-1)
  
  # export to clusters
  snow::clusterExport(
    cl = cl,
    list = c('predicts_sp','names_sub','names_tmp', 'values', 'names', 'length', 'mean', 'sd',
             'tmp', 'SP','rasterize','crop','trim', 'grep', 'sapply', 'strsplit',
             'cellStats', 'thresh', 'tmp1901_1930'),envir = environment())
  
  temperatureVars <- data.frame(t(parSapply(
    cl = cl,X = (1:nrow(predicts_sp)),FUN = function(i){
      
      # #for testing
      # temperatureVars <- NULL
      # for(i in 2306:2310){
      # print(i)
      
      #Get end sample date for sample in predicts
      sampDate <- predicts_sp$Sample_end_latest[i]
      
      #Reformat date for string matching
      sampDate <- substr(sampDate,1, 7)
      sampDate <- gsub("-", ".", sampDate, fixed = TRUE)
      
      #Match date in predicts with month in CRU climate data
      month_match <- which(names_sub==sampDate)
      
      # edit: use months from 5 year pre-sample, rather than 1 year
      surrounding_months <- names_tmp[(month_match-59):(month_match)]
      
      #Create a layer for average temperature in the year preceding end sample date
      temp <- tmp[[surrounding_months]]

      ## Mask to improve speed
      mask <- trim(rasterize(SP[i, ], temp[[1]]))
      mapCrop <- crop(temp, mask) 
      # plot(mapCrop) shows the one cell for each month
      
      
      # there are instances where there are no months above the threshold and
      # other instances where points do not line up with the tmp layers (in the sea?)
      # so this if statement is necessary to avoid errors in those instances.
      if(!length(names(mapCrop)[values(mapCrop) >= thresh]) == 0 & length(values(mapCrop)[!is.na(values(mapCrop))]) > 0 ){
        
        # Get the average temperature for each month across 5 years
        vals <- NULL
        
        # for each month, get the average temp over the 5 years
        for(k in 1:12){
          
          if(k < 10){ mon <- paste0(0, k) }else {mon <- k}
          
          monthmean <- values(mean(mapCrop[[grep(mon, sapply(strsplit(names(mapCrop), "[.]"), "[[", 2))  ]]))
          
          vals <- rbind(vals, c(mon, monthmean))
          
        }
        
        
        vals <- as.data.frame(vals)
        vals$V2 <- as.numeric(as.character(vals$V2))
        
        # which months are the 5 year average >= the threshold
        vals <- vals[vals$V2 >= thresh, ]
        
        # if the average monthly temp are below 10
        if(nrow(vals) == 0){
          
          avg_temp <- NA
          Anom <- NA
          StdAnom <- NA
          n_months = 0
          
          
          
          return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))    
          
          #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))
          
        }else{
        
        
        # which are the good months
        months <- vals$V1
        
        # how many months are at or above the threshold?
        n_months <- length(months)
        
        # calculate the "present day" mean 
        avg_temp <- mean(vals$V2)

        
        ### now work out the baseline mean and sd for the active months###
        
        # get the values for that grid cell across all years
        baseline <- crop(tmp1901_1930, mask)
        
        # subset the baseline to just the required months
        baseline <-  baseline[[names(baseline)[sapply(strsplit(names(baseline), "[.]"), "[[", 2) %in% months]]]
        
        # get the mean and sd
        mean_baseline <- mean(values(baseline))
        sd_mean_baseline <- sd(values(baseline))
        
        
        
        ### now calc the anomaly for that site, using the site specific baselines ###
        
        Anom <- avg_temp - mean_baseline
        StdAnom <-  Anom/sd_mean_baseline
        
        
        return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))
        
        #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))

      }}else{
        avg_temp <- NA
        Anom <- NA
        StdAnom <- NA
        n_months = 0

      
        
        return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))    
        
        #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))
      }
      
      
    }
    
  )))
  
  
  snow::stopCluster(cl)
  
  st2 <- Sys.time()
  
  print(st2 - st1) # Time difference of 1.165101 hours on Charlie's laptop
  
  
  # organise the anomaly info along with the predicts data
  temperatureVars <- as.data.frame(temperatureVars)
  
  
  ## add new values in temperatureVars into predicts dataset
  predicts_sp$avg_temp <- temperatureVars$avg_temp
  predicts_sp$TmeanAnomaly <- temperatureVars$Anom
  predicts_sp$StdTmeanAnomaly <- temperatureVars$StdAnom
  predicts_sp$n_months <- temperatureVars$n_months
  
  
  ## Joe's from here
  
  # # calculate the mean baseline, and convert to character for merging
  # climate_start_mean <- calc_baseline(tmp1901_1930, 
  #                                     func = mean, 
  #                                     pred_points = PRED_sites, 
  #                                     pred_points_sp = PRED_sites_sp) %>%
  #   mutate(Latitude = as.character(Latitude)) %>%
  #   mutate(Longitude = as.character(Longitude))
  # 
  # # calculate the sd baseline, and convert to character for merging
  # climate_start_sd <- calc_baseline(tmp1901_1930, 
  #                                   func = stats::sd, 
  #                                   pred_points = PRED_sites, 
  #                                   pred_points_sp = PRED_sites_sp) %>%
  #   mutate(Latitude = as.character(Latitude)) %>%
  #   mutate(Longitude = as.character(Longitude))
  # 
  # calculate the mean temperatures for each predicts site, 11 months previously
  # set up empty list for each dataframe
  
  # # unique to each site due to inclusion of insect active months only
  # raster_means <- list()
  # 
  # # time the length of the loop over each unique date in predicts sites
  # system.time(
  #   # for each unique site month, select the pixels for that date, and 11 months previously, and then convert each raster set to a dataframe
  #   for(i in 1:length(unique(PRED_sites$Sample_end_latest))){
  #     
  #     # create unique list of the end dates for the predicts sites
  #     pred_dates <- unique(PRED_sites$Sample_end_latest)
  #     
  #     # select the raster pixels for the end date of the predicts site, and then the index of that name
  #     date_raster <- grepl(pred_dates[i], names(tmp))
  #     site_index <- which(date_raster == TRUE)
  #     
  #     # select the previous 11 indices - i.e. 11 months previous worth of pixels
  #     ind_raster <- tmp[[names(tmp)[(site_index - 11): site_index]]]
  #     
  #     # filter the raster for only coordinates we have PREDICTS sites for that date
  #     PRED_sites_filt <- PRED_sites %>%
  #       filter(Sample_end_latest == pred_dates[i]) %>%
  #       select(Longitude, Latitude) %>%
  #       SpatialPoints()
  #     
  #     # identify site ids for merging
  #     site_ids <- PRED_sites %>%
  #       filter(Sample_end_latest == pred_dates[i]) %>%
  #       select(id_col, Longitude, Latitude)
  #     
  #     # filter the raster for that date for the locations we have predicts sites
  #     PRED_coords <- cbind(site_ids, extract(ind_raster, PRED_sites_filt, na.rm = FALSE))
  #     
  #     # convert that set of dates to a dataframe
  #     ind_raster_frame <- as.data.frame(PRED_coords)
  #     
  #     # remove the extra coordinate columns for calculating the row means
  #     ind_raster_values <- ind_raster_frame %>% select(-id_col, -Longitude, -Latitude)
  #     
  #     # calculate the mean values for each coordinate and bind back onto the coordinates
  #     raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% select(id_col, Longitude, Latitude)), rowMeans(ind_raster_values))
  #     colnames(raster_means[[i]]) <- c("end_date", "id_col", "x", "y", "mean_value")
  #     
  #     # print the iteration number
  #     print(i)
  #   }
  # )
  # 
  # ## adjust the mean value for each site for the baseline at that site
  # # first, merge the baseline sd and mean by coordinate for each site
  # adjusted_climate <- rbindlist(raster_means) %>%
  #   select(-end_date) %>%
  #   unique() %>%
  #   inner_join(climate_start_mean, by = "id_col") %>%
  #   rename("mean_base" = "data_fin") %>%
  #   inner_join(climate_start_sd, by = "id_col") %>%
  #   rename("sd_base" = "data_fin") %>%
  #   mutate(anomaly = mean_value - mean_base) %>%
  #   mutate(standard_anom = anomaly / sd_base)
  # 
  # bind the adjusted climate data back onto the predicts sites
  # predicts_climate <- inner_join(order.sites.div, adjusted_climate, by = "id_col")
  
  predicts_climate <- predicts_sp@data
  
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
  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ StdTmeanAnomaly * Predominant_land_use + (1|SS), data = predict_climate_list[[m]]) 
  print(AIC(model_2c_abundance[[m]]))
  
  model_2c_abundance[[m]] <- lmerTest::lmer(log(Total_abundance) ~ StdTmeanAnomaly * Predominant_land_use + (1|SS) + (1|SSB), data = predict_climate_list[[m]]) 
  print(AIC(model_2c_abundance[[m]]))
  
  # run predictions for the model of standard anomaly
  abundance_model[[m]] <- predict_continuous(model = model_2c_abundance[[m]],
                                             model_data = predict_climate_list[[m]],
                                             response_variable = "Total_abundance",
                                             categorical_variable = c("Predominant_land_use"),
                                             continuous_variable = c("StdTmeanAnomaly"),
                                             continuous_transformation = "",
                                             random_variable = c("SS", "SSB", "SSBS"))
  
  # plot for standardised anomaly and land-use for abundance
  main_plot_abundance[[m]] <- ggplot(abundance_model[[m]]) +
    geom_line(aes(x = StdTmeanAnomaly, y = y_value, colour = Predominant_land_use), size = 1.5) +
    geom_ribbon(aes(x = StdTmeanAnomaly, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
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
