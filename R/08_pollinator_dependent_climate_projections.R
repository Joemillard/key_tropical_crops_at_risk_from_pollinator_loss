# read in packages
library(dplyr)
library(raster)
library(rworldmap) 
library(rworldxtra)
library(ggplot2)
library(cowplot)

# source in additional functions
source("R/00_functions.R")

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

# set up historical anomaly to be added on
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
for(i in 1:11){
  years <- years - 3
  years_list[[i]] <- years
}

# set up list for climate anomalies
tmp2069_71std_climate_anomaly <- list()

for(i in 1:length(years_list)){

  # file path for ISIMIP data
  all.files <- dir(path = SSP_directory,recursive = TRUE,full.names = TRUE)
  
  # using RCP 8.5
  mean.temp.2069.2071 <- stack(lapply(X = years_list[[i]],FUN = function(yr){
    
    print(yr)
    
    all.model.files <- all.files[grepl("rcp85",all.files) & grepl(yr,all.files)]
    
    # Check that there are the same files for each scenario-year combination
    stopifnot(all(sapply(
      X = gsub("G:/Extra_data_files/climate_projections/ISIMIPAnomalies.tar/ISIMIPAnomalies/","",all.model.files),function(f) return(strsplit(x = f,split = "[-_]",fixed = FALSE)[[1]][1]))==
        c("GFDL","HadGEM2","IPSL","MIROC5")))
    
    # what are each of these files?
    
    # 
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

