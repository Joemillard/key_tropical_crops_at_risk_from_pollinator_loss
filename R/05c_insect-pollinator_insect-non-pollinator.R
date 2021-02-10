# script for pollinating/non-pollinating, need to think about species assigned as non-pollinators

# packages to read in
library(dplyr)

# read in the original predicts database 
PREDICTS <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/Data/PREDICTS/database.rds") %>%
  filter(Class == "Insecta") %>%
  droplevels()

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds")  %>%
  filter(Class == "Insecta") %>%
  droplevels()

# filter the pollinators from the overall predicts database
PREDICTS_non_pollinating <- PREDICTS %>%
  filter(!COL_ID %in% PREDICTS_pollinators_orig$COL_ID) %>%
  droplevels() %>%
  filter(Order == "Hymenoptera") %>%
  droplevels()

table(PREDICTS_non_pollinating$Family)
