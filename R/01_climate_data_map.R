# script for climate anomaly
### Note: 2005 is the mean year for insect data
# try split up into three panels: overall (1901-2006), earlier (1901-1950), latter (1950-2006)

# load required libraries
library(raster)
library(ggplot2)
library(dplyr)

# load in the mean temperature data from CRU
tmp <- raster::stack("data/cru_ts4.03.1901.2018.tmp.dat.nc", varname="tmp")

# take names of values for 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]

# calculate the mean and sd of the baseline values
tmp1901_1905mean <- calc(tmp1901_1905, mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)

# extract data for the years 2004-2006
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

### Calculate the standardised anomaly ###
# calc the mean for present time period
tmp2004_6mean <- calc(tmp[[names(tmp)[1237:1272]]], mean)

# calc mean for baseline
tmp2004_6_climate_anomaly <- (calc(tmp2004_6, mean) - tmp1901_1905mean)

# standardise the baseline
tmp2004_6std_climate_anomaly <- (calc(tmp2004_6, mean) - tmp1901_1905mean) / tmp1901_1905sd

# reproject on mollweide projection - note warning of missing points to check -- "55946 projected point(s) not finite"
tmp2004_6std_climate_anomaly <- projectRaster(tmp2004_6std_climate_anomaly, crs = "+proj=moll +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# convert the raster to a format amenable to ggplot
# convert the climate anomaly raster to a spatial pixels data frame, and then rename the columns
anom_spdf <- as(tmp2004_6std_climate_anomaly, "SpatialPixelsDataFrame")
anom_df <- as.data.frame(anom_spdf)
colnames(anom_df) <- c("value", "x", "y")

# group the categories of climate anomaly into factors
anom_df$value_group[anom_df$value > 2] <- "> 2"
anom_df$value_group[anom_df$value > 1 & anom_df$value <= 2] <- "1 - 2"
anom_df$value_group[anom_df$value > 0.5 & anom_df$value <= 1] <- "0.5 - 1"
anom_df$value_group[anom_df$value > 0.25 & anom_df$value <= 0.5] <- "0.25 - 0.5"
anom_df$value_group[anom_df$value >= 0 & anom_df$value <= 0.25] <- "0 - 0.25"
anom_df$value_group[anom_df$value < 0] <- "< 0"

# order the levels of those factors
anom_df$value_group <- factor(anom_df$value_group, levels = c("> 2", "1 - 2", "0.5 - 1", "0.25 - 0.5", "0 - 0.25", "< 0"))

# plot the ggplot map for climate anomaly
anom_df %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = value_group)) +
  scale_fill_manual("Standardised climate anomaly", values = c("#000000", "darkred", "#D55E00", "#E69F00", "#F0E442", "#56B4E9")) +
  coord_equal() +
  theme(panel.background = element_blank(),
        panel.bord = element_rect(fill = NA),
        panel.grid = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())

ggsave("global_climate_anomaly.png", scale = 1.2, dpi = 350)
