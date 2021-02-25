# read in packages
library(raster)
library(ggplot2)
library(rasterVis)
library(dplyr)
library(viridis)
library(rworldmap) 
library(rworldxtra)
library(cowplot)

# source in additional functions
source("R/00_functions.R")

# list the crop specific folders in the directory for external hard drive
cropdirs <- list.dirs("G:/Extra_data_files/HarvestedAreaYield175Crops_Geotiff/HarvestedAreaYield175Crops_Geotiff/Geotiff", recursive = FALSE)

# read in the Klein pollinator dependent crops
klein_cleaned <- read.csv(here::here("data/KleinPollinationDependentCrops.tar/KleinPollinationDependentCrops/data_cleaned.csv"))

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
klein_cleaned_av <- klein_cleaned %>%
  group_by(MonfredaCrop) %>%
  mutate(av = mean(dependence_ratio, na.rm = TRUE)) %>%
  mutate(standard_dev = sd(dependence_ratio, na.rm = TRUE)) %>%
  ungroup() %>%
  select(MonfredaCrop, av, standard_dev) %>%
  unique()

# subset klein_cleaned for those with crop data
klein_cleaned_filt <- klein_cleaned_av %>%
  filter(MonfredaCrop %in% pollinat_crops_simp) %>%
  arrange(MonfredaCrop)

# multiply each raster by its pollination dependence for that crop
rate_rasters_adj <- list()
for(i in 1:length(rate_rasters)){
  rate_rasters_adj[[i]] <- rate_rasters[[i]] * klein_cleaned_filt$av[i]
  print(i)
}

# sum the production for all the rasters
# organise all of the rasters into a stack and sum
stacked_rasters <- stack(rate_rasters_adj)
crop.total <- sum(stacked_rasters, na.rm = T)

# calculate total pollination dependent production
total_production <- sum(crop.total[])

# read in the raster for the historical data to get the baseline
## standardised climate anomaly script
# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# take names of values for 1901 to 1931 - 30 year baseline
tmp1901_1931 <- tmp[[names(tmp)[1:361]]]

# calculate the mean and sd of the baseline values
tmp1901_1931mean <- calc(tmp1901_1931, mean)
tmp1901_1931sd <- calc(tmp1901_1931, stats::sd)

## read in the rasters for the future data, start with SSP585
# set up the starting directory
SSP_directory <- ("G:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies")

# set up historical change to be added on
months.1979.2013 <- 937:1356

# calculate the average temperature for 1979-2013 onto which anomaly is added
hist.mean.temp.1979.2013 <- stack(stackApply(x = tmp[[months.1979.2013]],
                                             indices = (rep(1:35,each=12)),fun = mean))
hist.mean.temp.1979.2013 <- stackApply(x = hist.mean.temp.1979.2013,indices = rep(1,35),
                                       fun = mean)

# list the files in that directory
SSP_folders <- list.files(SSP_directory)

# select rcp26 for MIROC5
SSP_folders <- SSP_folders[grepl("MIROC5_rcp26", SSP_folders)]

# setp up empty vector for file paths
SSP_file_path <- c()

# iterate through each of the folders
for(i in 1:length(SSP_folders)){
  SSP_file_path[i] <- paste(SSP_directory, SSP_folders[i], sep = "/")
}

# future projection list
future_projection <- list()
future_projection_anomaly <- list()

# selection of years
years <- 2048:2050

years_list <- list()

# set up list of years
for(i in 1:33){
  years <- years - 1
  years_list[[i]] <- years
}

# set up list for climate anomalies
tmp2069_71std_climate_anomaly <- list()

for(i in 1:length(years_list)){

  # file path for ISIMIP data
  all.files <- dir(path = SSP_directory,recursive = TRUE,full.names = TRUE)
  
  # using RCP 8.5
  mean.temp.2069.2071 <- stack(lapply(X = years_list[[i]],FUN = function(yr){
    
    # print the set of years for that iteration
    print(yr)
    
    # subset for all files for rcp85 and the set of years for that iteration
    all.model.files <- all.files[grepl("rcp85",all.files) & grepl(yr,all.files)]
    
    # Check that there are the same files for each scenario-year combination
    stopifnot(all(sapply(
      X = gsub("G:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies/","",all.model.files),function(f) return(strsplit(x = f,split = "[-_]",fixed = FALSE)[[1]][1]))==
        c("GFDL","HadGEM2","IPSL","MIROC5")))
    
    meant.anom <- mean(stack(lapply(X = all.model.files,function(f){
      
      ras <- stack(f)$"X0.1"
      
    })),na.rm=TRUE)
    
    meant <- hist.mean.temp.1979.2013 + (meant.anom/10)
    
    return(meant)
    
  }))
  
  mean.temp.2069.2071 <- stackApply(x = mean.temp.2069.2071,indices = rep(1,3),fun = mean)
  
  # calc the anomalies for the future years
  tmp2069_71_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1931mean)
  tmp2069_71std_climate_anomaly[[i]] <- (mean.temp.2069.2071-tmp1901_1931mean)  / tmp1901_1931sd

}

std_anom_high <- list()

# for each set of climate anomaly data, subset for greater than equal to 4
# select all pollination dependence production at greater than 4 anomaly and sum
for(i in 1:length(tmp2069_71std_climate_anomaly)){
  
  # convert the raster to a dataframe to plot with ggplot
  std_anom_high[[i]] <- as(tmp2069_71std_climate_anomaly[[i]], "SpatialPixelsDataFrame")
  std_anom_high[[i]] <- as.data.frame(std_anom_high[[i]])
  
  # assign 0 in layer as NA
  std_anom_high[[i]]$layer[std_anom_high[[i]]$layer == 0] <- NA
  
  # convert spatial dataframe to coordinates
  std_anom_high[[i]] <- std_anom_high[[i]] %>%
    dplyr::filter(layer >= 4) %>%
    dplyr::select(x, y) %>%
    unique() %>%
    SpatialPoints()
  
  # print the number of coordinates
  print(length(std_anom_high[[i]]))
  
}

# set up vector for total production
vulnerable_production_list <- list()
vulnerable_production <- c()

# for each set of coordinates, extract the pollination dependent values and sum
for(i in 1:length(std_anom_high)){

  # convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
  vulnerable_production_list[[i]] <- extract(crop.total, std_anom_high[[i]], na.rm = FALSE)
  vulnerable_production[i] <- unlist(vulnerable_production_list[[i]]) %>% sum()
  
}

# create dataframe for exposed production and build datafrmae
data.frame("production" = vulnerable_production, "year" = c(seq(2048, 2016, -1))) %>%
  mutate("percentage" = (vulnerable_production / total_production) * 100) %>%
   ggplot() +
    geom_line(aes(x = year, y = percentage)) +
    geom_point(aes(x = year, y = percentage)) +
   # scale_y_continuous(limits = c(0, 265000), expand = c(0, 0)) +
    scale_x_continuous(limits = c(2015, 2050), expand = c(0, 0), breaks = c(2015, 2020, 2025, 2030, 2035, 2040, 2045, 2050)) +
    ylab("Un-viable pollination dependent prod. (mt tonnes)") +
    xlab("Year") +
    theme_bw() +
    theme(panel.grid = element_blank())

