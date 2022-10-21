# script for crop level proportion production risk projections at the cell level

# read in packages
#library(raster)
library(ggplot2)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)
library(cowplot)
library(data.table)
library(lme4)
library(yarg)
library(forcats)
library(ggrepel)
# library(StatisticalModels)

# source in additional functions
source("R/00_functions.R")

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds")

# set up the starting directory for future climate data
SSP_directory <- ("D:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies")

# read in fao to monfreda conversion
fao_monfreda <- read.csv("data/trade_flow/FAO_Monfreda_conv.csv", stringsAsFactors = FALSE) %>%
  mutate(Cropname_FAO = gsub(",", "", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("corian.", "coriander", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Stone fruit nes", "Other stone fruits", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Chestnuts", "Chestnuts in shell", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Fruit Fresh Nes", "Other fruits n.e.c.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Kolanuts", "Kola nuts", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Other melons \\(inc\\.cantaloupes\\)", "Cantaloupes and other melons", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Oilseeds Nes", "Other oil seeds n\\.e\\.c\\.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Pepper \\(Piper ", "Pepper \\(piper ", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Arecanuts", "Areca nuts", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Tung Nuts", "Tung nuts", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Roots and Tubers nes", "Roots and tubers nes", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("PyrethrumDried", "Pyrethrum dried", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Tangerines mandarins clem.", "Tangerines mandarins clementines satsumas", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Cinnamon \\(canella\\)", "Cinnamon \\(cannella\\)", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Hemp Tow Waste", "Hemp tow waste", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Natural rubber", "Rubber natural", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Plantains", "Plantains and others", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("MatŽ", "Mat?", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Onions \\(inc\\. shallots\\) green", "Onions shallots green", Cropname_FAO))%>%
  mutate(Cropname_FAO = gsub("Karite Nuts \\(Sheanuts\\)", "Karite nuts (sheanuts)", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Sour cherries", "Cherries sour", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Other Bastfibres", "Bastfibres other", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Agave Fibres Nes", "Agave fibres raw n.e.c.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Almonds with shell", "Almonds in shell", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Brazil nuts with shell", "Brazil nuts in shell", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Broad beans horse beans dry", "Broad beans and horse beans dry", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Cashew nuts with shell", "Cashew nuts in shell", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Anise badian fennel coriander", "Anise badian coriander cumin caraway fennel and juniper berries raw", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Chillies and peppers green", "Chillies and peppers green (Capsicum spp. and Pimenta spp.)", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Citrus fruit nes", "Other citrus fruit n.e.c.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Berries Nes", "Other berries and fruits of the genus vaccinium n.e.c.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Coconuts", "Coconuts in shell", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Mangoes mangosteens guavas", "Mangoes guavas and mangosteens", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Vanilla", "Vanilla raw", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Fruit tropical fresh nes", "Other tropical fruits n\\.e\\.c\\.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Spices nes", "Other stimulant spice and aromatic crops n\\.e\\.c\\.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Soybeans", "Soya beans", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Rapeseed", "Rape or colza seed", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Nutmeg mace and cardamoms", "Nutmeg mace cardamoms raw", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Pigeon peas", "Pigeon peas dry", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Nuts nes", "Other nuts \\(excluding wild edible nuts and groundnuts\\) in shell n\\.e\\.c\\.", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Leguminous vegetables nes", "Broad beans and horse beans green", Cropname_FAO)) %>%
  mutate(Cropname_FAO = gsub("Groundnuts with shell", "Groundnuts excluding shelled", Cropname_FAO))

# calculate per country average total production value
fao_prod_value <- read.csv("data/trade_flow/FAOSTAT_data_en_10-4-2022_producer_price.csv", stringsAsFactors = FALSE) %>%
  filter(!grepl("Meat", Item)) %>%
  filter(!grepl("Eggs", Item)) %>%
  filter(!grepl("China, mainland", Area)) %>%
  filter(!grepl("China, Hong Kong SAR", Area)) %>%
  filter(!grepl("Milk", Item)) %>%
  group_by(Area, Item) %>%
  summarise(mean_price_tonne = mean(Value, na.rm = TRUE)) %>%
  group_by(Item) %>%
  summarise(overall_price_tonne = median(mean_price_tonne, na.rm = TRUE)) %>% 
  mutate(Item = gsub(",", "", Item)) %>% 
  inner_join(fao_monfreda, by = c("Item" = "Cropname_FAO"))

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
PRED_sites <- order.sites.div %>% dplyr::select(id_col, Latitude, Longitude, Sample_end_latest) %>%
  mutate(Sample_end_latest = paste("X", substr(Sample_end_latest, start = 1, stop = 7), sep = "")) %>%
  mutate(Sample_end_latest = gsub("-", ".", Sample_end_latest)) %>%
  filter(!is.na(Latitude))

# calculate the means and standard deviation for the beginning of the series
# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# extract the points for each the predicts coordinates
PRED_sites_sp <- PRED_sites %>%
  dplyr::select(Longitude, Latitude) %>%
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
      dplyr::select(Longitude, Latitude) %>%
      SpatialPoints()
    
    # identify site ids for merging
    site_ids <- PRED_sites %>%
      filter(Sample_end_latest == pred_dates[i]) %>%
      dplyr::select(id_col, Longitude, Latitude)
    
    # filter the raster for that date for the locations we have predicts sites
    PRED_coords <- cbind(site_ids, extract(ind_raster, PRED_sites_filt, na.rm = FALSE))
    
    # convert that set of dates to a dataframe
    ind_raster_frame <- as.data.frame(PRED_coords)
    
    # remove the extra coordinate columns for calculating the row means
    ind_raster_values <- ind_raster_frame %>% dplyr::select(-id_col, -Longitude, -Latitude)
    
    # calculate the mean values for each coordinate and bind back onto the coordinates
    raster_means[[i]] <- cbind(names(tmp)[site_index], (ind_raster_frame %>% dplyr::select(id_col, Longitude, Latitude)), rowMeans(ind_raster_values))
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
  dplyr::select(-end_date) %>%
  unique() %>%
  inner_join(climate_start_mean, by = "id_col") %>%
  rename("mean_base" = "data_fin") %>%
  inner_join(climate_start_sd, by = "id_col") %>%
  rename("sd_base" = "data_fin") %>%
  mutate(anomaly = mean_value - mean_base) %>%
  mutate(standard_anom = anomaly / sd_base)

# bind the adjusted climate data back onto the predicts sites
predicts_climate <- inner_join(order.sites.div, adjusted_climate, by = "id_col")

# bring in basemap for climate site plot 
base_map <- get_basemap()

# reproject basemap
base_map <- spTransform(base_map, CRS = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

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
main_plot_abundance <- abundance_model %>% 
  ggplot() +
  geom_line(aes(x = standard_anom, y = y_value, colour = Predominant_land_use), size = 1.5) +
  geom_ribbon(aes(x = standard_anom, y = y_value, fill = Predominant_land_use, ymin = y_value_minus, ymax = y_value_plus), alpha = 0.4) +
  scale_y_continuous("Total abundance", breaks = c(1.609438, 2.302585, 2.995732, 3.6888795, 4.382027, 5.075174, 5.768321), labels = c(5, 10, 20, 40, 80, 160, 320)) +
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
pollinated_crops <- sort(pollinated_crops)
pollinat_crops_simp <- gsub("D:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff/", "", pollinated_crops)
pollinat_crops_simp <- gsub('([^/]+$)', "", pollinat_crops_simp)
pollinat_crops_simp <- gsub('/', "", pollinat_crops_simp)

# read in each of the rasters
for(i in 1:length(pollinated_crops)){
  rate_rasters[[i]] <- raster::raster(pollinated_crops[i])
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
    dplyr::select(MonfredaCrop, av, standard_dev) %>%
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
klein_cleaned_filt <- subset_klein(av_dependence(klein_cleaned)) %>%
  arrange(MonfredaCrop)

# multiply each raster by its pollination dependence for that crop and rebuffer
rate_rasters_adj <- list()
for(i in 1:length(rate_rasters)){
  rate_rasters_adj[[i]] <- rate_rasters[[i]] * klein_cleaned_filt$av[i]
  
  print(i)
}

# iteration vector for raster
adj_it_vec <- c()
adj_sum_vec <- c()
iteration <- c()

# select the top 5 crops for total pollination dependent production
for(i in 1:length(rate_rasters_adj)){
  adj_sum_vec[i] <- sum(rate_rasters_adj[[i]][])
  adj_it_vec[i] <- klein_cleaned_filt$MonfredaCrop[i]
  iteration[i] <- i
}

# sum of pollination dependent production production for each crop
pollination_production_sum <- data.frame(adj_it_vec, adj_sum_vec, iteration) %>%
  arrange(desc(adj_sum_vec)) %>%
  slice(1:22)

# pollination production barplot for top production
pollination_production_sum %>%
  filter(!adj_it_vec %in% c("fruitnes", "tropicalnes")) %>%
  mutate(adj_it_vec = factor(adj_it_vec, 
                       labels = c("Apple", "Bean", "Cocoa", "Coconut", "Coffee", "Cucumber",
                                  "Eggplant", "Mango", "Melon", "Oilpalm", "Oilseed", 
                                  "Peach", "Pear", "Plum", "Pumpkin","Rapeseed",  
                                  "Soybean",  "Sunflower", "Tomato", "Watermelon"))) %>%
  mutate(adj_it_vec = fct_reorder(adj_it_vec, -adj_sum_vec)) %>%
  ggplot() +
    geom_bar(aes(x = adj_it_vec, adj_sum_vec), stat =  "identity") +
    scale_y_continuous("Total pollination dependent production (million tonnes)", 
                       expand = c(0, 0),
                       limits = c(0, 73000000),
                       breaks = c(0, 20000000, 40000000, 60000000),
                       labels = c("0", "20", "40", "60")) +
    xlab("Crop") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("total_pollination_dependent_production_2.png", scale = 0.9, dpi = 350)

# subset original rate rasters for those in the top 20 pollination dependent production
rate_rasters_sub <- rate_rasters[pollination_production_sum$iteration]
rate_rasters_adj_sub <- rate_rasters_adj[pollination_production_sum$iteration]

# change the projection of each raster
for(i in 1:length(rate_rasters_sub)){
  
# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
  rate_rasters_sub[[i]] <- projectRaster(rate_rasters_sub[[i]], crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
  # reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
  rate_rasters_adj_sub[[i]] <- projectRaster(rate_rasters_adj_sub[[i]], crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
}

## standardised climate anomaly script
# calculate the mean and sd of the baseline values
tmp1901_1930mean <- calc(tmp1901_1930, mean)
tmp1901_1930sd <- calc(tmp1901_1930, stats::sd)

## read in the rasters for the future data, start with SSP585
# set up historical change to be added on
months.1979.2013 <- 937:1356

# calculate the average temperature for 1979-2013 onto which anomaly is added
hist.mean.temp.1979.2013 <- stack(stackApply(x = tmp[[months.1979.2013]],
                                             indices = (rep(1:35,each=12)),fun = mean))
hist.mean.temp.1979.2013 <- stackApply(x = hist.mean.temp.1979.2013,indices = rep(1,35),
                                       fun = mean)

# selection of years and empty year list
years <- 2049:2051
years_list <- list()

# set up list of years
for(i in 0:34){
  
  year <- years - i
  years_list[[i+1]] <- year
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
RCP_scenarios <- c("rcp60")
climate_model_combs_adj <- c("GFDL|HadGEM2|IPSL|MIROC5")

# iterate through each set of years as a rolling average
for(i in 1:length(years_list)){
  
  # file path for ISIMIP data
  all.files <- dir(path = SSP_directory,recursive = TRUE, full.names = TRUE)
  
  # using RCP 8.5 calculate average of separate models
  mean.temp.2069.2071 <- stack(lapply(X = years_list[[i]], FUN = average_clim_models, RCP = RCP_scenarios, clim_models = climate_model_combs_adj))
  
  mean.temp.2069.2071 <- stackApply(x = mean.temp.2069.2071,indices = rep(1,3), fun = mean)
  
  # calc the anomalies for the future years
  tmp2069_71_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1930mean)
  tmp2069_71std_climate_anomaly[[i]] <- (mean.temp.2069.2071-tmp1901_1930mean) / tmp1901_1930sd
  
  # reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
  tmp2069_71std_climate_anomaly[[i]] <- projectRaster(tmp2069_71std_climate_anomaly[[i]], crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
  
}

# set up empty list for standardised climate anomaly
std_high_abun_adj <- list()
std_anom_high <- list()

# predict abundance at 0 warming on cropland
zero_data <- data.frame("standard_anom" = 0, Predominant_land_use = "Cropland")
zero_warming_abundance <- predict(model_2c_abundance, zero_data, re.form = NA)
zero_warming_abundance <- exp(zero_warming_abundance)

# for each set of climate anomaly data, predict abundance reduction for all climate anomaly values in each cell
# and then sum abundance adjusted pollination dependence
for(i in 1:length(tmp2069_71std_climate_anomaly)){
  
  # convert the raster to a dataframe to plot with ggplot
  std_high_abun_adj[[i]] <- as(tmp2069_71std_climate_anomaly[[i]], "SpatialPixelsDataFrame")
  std_high_abun_adj[[i]] <- as.data.frame(std_high_abun_adj[[i]])
  
  # set up prediction data on basis of that set of years
  new_data_pred <- data.frame("standard_anom" = std_high_abun_adj[[i]]$layer, Predominant_land_use = "Cropland")
  
  # predict abundance for climate anomaly and join to data frame
  predicted_abundance <- predict(model_2c_abundance, new_data_pred, re.form = NA)
  std_high_abun_adj[[i]]$abundance <- exp(predicted_abundance)
  
  # for any location that's cooled abundance is that at no warming
  std_high_abun_adj[[i]]$abundance[std_high_abun_adj[[i]]$layer <= 0] <- zero_warming_abundance
  
  # calculate percentage change from place with 0 warming, and convert to vulnerability
  std_high_abun_adj[[i]]$abundance_change <- 1 - (std_high_abun_adj[[i]]$abundance / zero_warming_abundance)
  
  # convert spatial dataframe to coordinates
  std_anom_high[[i]] <- std_high_abun_adj[[i]] %>%
    dplyr::select(x, y) %>%
    unique() %>%
    SpatialPoints(proj4string = CRS("+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
}

# set up empty list for each of top pollination dependent crops
all_crop_list <- list()

for(j in 1:length(rate_rasters_sub)){
  
  # for each set of coordinates, extract the pollination dependent values and sum
  for(i in 1:length(std_anom_high)){
    # convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
    std_high_abun_adj[[i]]$production <- extract(rate_rasters_adj_sub[[j]], std_anom_high[[i]], na.rm = FALSE)
    std_high_abun_adj[[i]]$total_production <- extract(rate_rasters_sub[[j]], std_anom_high[[i]], na.rm = FALSE)
    std_high_abun_adj[[i]]$pollinator_vulnerability <- unlist(std_high_abun_adj[[i]]$production * std_high_abun_adj[[i]]$abundance_change) / std_high_abun_adj[[i]]$total_production
  }
  
  all_crop_list[[j]] <- std_high_abun_adj
}  

median_all_crop <- list()
all_crop_fin <- list()

# for each yearly set, calculate the pollinator vulnerability for each country
for(i in 1:length(all_crop_list)){
  for(j in 1:length(all_crop_list[[i]])){
    median_all_crop[[j]] <- all_crop_list[[i]][[j]] %>%
    summarise(total = quantile(pollinator_vulnerability, probs = c(0.5), na.rm = TRUE), 
              upp_conf = quantile(pollinator_vulnerability, probs = c(0.975), na.rm = TRUE), 
              lower_conf = quantile(pollinator_vulnerability, probs = c(0.025), na.rm = TRUE),
              vul_prod = sum(production * abundance_change, na.rm = TRUE),
              total_prod = sum(total_production, na.rm = TRUE),
              prop_prod = vul_prod/total_prod) %>%
    mutate(year = 2051 - j)
  }
  
  all_crop_fin[[i]] <- rbindlist(median_all_crop)
} 

for(i in 1:length(all_crop_fin)){
  all_crop_fin[[i]]$crop <- pollination_production_sum$adj_it_vec[i]
}

# build figure for variation with prop production at risk in brackets  
rbindlist(all_crop_fin) %>%
  group_by(crop) %>%
  mutate(change = max(total) - min(total)) %>%
  mutate(overall_val = mean(total)) %>%
  ungroup() %>%
  inner_join(fao_prod_value, by = c("crop" = "CROPNAME")) %>%
  dplyr::select(crop, overall_val, change, overall_price_tonne) %>%
  unique() %>%
  arrange(desc(overall_price_tonne)) %>%
  mutate(crop = ifelse(overall_val > 0.15, crop, "")) %>%
  mutate(crop = gsub("coffee", "Coffee", crop)) %>%
  mutate(crop = gsub("mango", "Mango", crop)) %>%
  mutate(crop = gsub("cocoa", "Cocoa", crop)) %>%
  mutate(crop = gsub("watermelon", "Watermelon", crop)) %>%
  mutate(crop = gsub("tropicalnes", "", crop)) %>%
  mutate(crop = gsub("melonetc", "", crop)) %>%
  mutate(crop = gsub("pumpkinetc", "", crop)) %>%
  ggplot() +
    geom_point(aes(x = change, y = overall_val, size = overall_price_tonne), pch=21, fill = "grey", colour = "black", alpha = 0.5) +
    geom_label_repel(aes(x = change, y = overall_val, label = crop), alpha = 0.7,
                     nudge_x = .04,
                     nudge_y = c(0.02),
                     segment.curvature = 0.1) +
    scale_size_continuous("Mean producer price \n(2015-2019; US$/tonne)") +
    scale_y_continuous(limits = c(0, 0.6), breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), labels = c("0", "0.1", "0.2", "0.3", "0.4", "0.5"), expand = c(0, 0)) +
    scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3), labels = c("0", "0.1", "0.2", "0.3"), limits = c(0, 0.3), expand = c(0, 0)) +
    xlab("Change in crop pollination risk") + 
    ylab("Overall crop pollination risk") +
    theme_bw() +
    theme(panel.grid = element_blank(), strip.text = element_text(size = 10.5), legend.position = "bottom")

ggsave("top_change_crop_8.png", scale = 0.85, dpi = 350)

# calculate crop with highest proportion of production at risk and write to csv
rbindlist(all_crop_fin) %>%
  group_by(crop) %>%
  mutate(change = max(total) - min(total)) %>%
  mutate(overall_val = mean(total)) %>%
  ungroup() %>%
  inner_join(joined_prod_value, by = c("crop" = "CROPNAME")) %>%
  dplyr::select(crop, overall_val, change, overall_price_kg) %>%
  unique() %>%
  arrange(desc(overall_price_kg)) %>%
  write.csv("crop_level_risk.csv")
