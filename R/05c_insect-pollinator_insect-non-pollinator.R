# script for pollinating/non-pollinating, need to think about species assigned as non-pollinators

# packages to read in
library(dplyr)
library(yarg)

# read in the original predicts database 
PREDICTS <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/Data/PREDICTS/database.rds") %>%
  filter(Class == "Insecta") %>%
  droplevels()

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds") %>%
  select(-clade_rank, -confidence)  %>%
  filter(Class == "Insecta") %>%
  droplevels() %>%
  mutate("Pollinating" = "Y")

# filter the pollinators from the overall predicts database
PREDICTS_non_pollinating <- PREDICTS %>%
  filter(!COL_ID %in% PREDICTS_pollinators_orig$COL_ID) %>%
  droplevels() %>%
  filter(Order != "Lepidoptera") %>%
  mutate("Pollinating" = "N")

# bind together the two dataframes
pollinat_bound <- rbind(PREDICTS_non_pollinating, PREDICTS_pollinators_orig)

# split diversity data into list of four for each order
diversityOrder <- split(x = pollinat_bound, f = pollinat_bound$Pollinating)

# drop unused levels from each list
diversityOrder <- lapply(diversityOrder, function(x) return(droplevels(x)))

# Remove empty rows
diversityOrder <- Filter(function(x) dim(x)[1] > 0, diversityOrder)

# calculate site metrics for each of the four order level subsets
order.sites.div <- do.call('rbind',lapply(X = diversityOrder, FUN = SiteMetrics,
                                          extra.cols = c("SSB", "SSBS","Biome", "Sampling_method",
                                                         "Study_common_taxon", "Sampling_effort",
                                                         "Sampling_effort_unit", "Realm",
                                                         "Use_intensity", "Order", "LUI", "zone"),
                                          sites.are.unique = TRUE, srEstimators = TRUE))

# add 1 to abundance and diversity 
order.sites.div$Total_abundance <- order.sites.div$Total_abundance + 1
order.sites.div$Simpson_diversity <- order.sites.div$Simpson_diversity + 1
