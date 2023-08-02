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

# read in the predicts pollinators
PREDICTS_pollinators_orig <- readRDS("C:/Users/joeym/Documents/PhD/Aims/Aim 2 - understand response to environmental change/outputs/PREDICTS_pollinators_8_exp.rds") %>%
  dplyr::select(-clade_rank, -confidence)  %>%
  filter(Class == "Insecta") %>%
  dplyr::filter(Predominant_land_use %in% c("Cropland", "Primary vegetation")) %>%
  droplevels() %>%
  mutate(Sampling_method = as.character(Sampling_method))

# calc longest gap between start and end date
data.frame("day" = as.numeric(as.Date(PREDICTS_pollinators_orig$Sample_end_latest) - as.Date(PREDICTS_pollinators_orig$Sample_start_earliest))) %>%
  arrange(day) %>%
  mutate(row_id = row_number()) %>%
  mutate(row_percentile = row_id / max(row_id)) %>%
  ggplot() + 
    geom_histogram(aes(day)) +
    theme_bw() +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 61000)) +
    ylab("Record count") +
    xlab("Sampling period (days)") +
    theme(panel.grid = element_blank(),
        text = element_text(size = 13),
        strip.text = element_text(size = 14))

ggsave("sampling_period_time.png", scale = 1, dpi = 350)

# grouping sampling methods
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "sweep nets"] <- "sweep net"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "sweep_nets"] <- "sweep net"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "sweep_netting"] <- "sweep net"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "various"] <- "multiple"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "soil cores along transects"] <- "soil core"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "pitfall trap transects"] <- "pit-fall traps"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "baited pit-fall traps"] <- "pit-fall traps"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "aerial flight-inception  trap"] <- "flight trap"
PREDICTS_pollinators_orig$Sampling_method[PREDICTS_pollinators_orig$Sampling_method == "aerial nets"] <- "flight trap"


# plot of sampling method frequency
PREDICTS_pollinators_orig %>% 
  group_by(Sampling_method) %>%
  tally() %>% 
  ungroup() %>%
  mutate(Sampling_method = stringr::str_to_sentence(Sampling_method)) %>%
  mutate(Sampling_method = forcats::fct_reorder(Sampling_method, -n)) %>%
  ggplot() +
    geom_bar(aes(x = Sampling_method, y = n), stat = "identity", fill = "black") +
    scale_y_continuous("Frequency", expand = c(0, 0), limits = c(0, 67000)) +
    scale_x_discrete("Sampling method") +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave("Sampling_method_distribution.png", scale = 1, dpi = 350)

# remove additional columns - number of unique species = 5100
pollinator_subset <- PREDICTS_pollinators_orig %>% 
  select(Kingdom, Phylum, Class, Order, Family, Best_guess_binomial) %>%
  filter(Best_guess_binomial != "") %>%
  filter(!grepl("\\d", Best_guess_binomial)) %>%
  unique()

# build table of counts
pollinator_subset_order <- pollinator_subset %>% 
  group_by(Order) %>%
  tally() %>%
  ungroup() %>%
  mutate(Order = forcats::fct_reorder(Order, -n))

# build table for number of orders
# plot of species count for PREDICTS pollinator subset by class and order
order_class <- pollinator_subset_order %>%
  ggplot() +
  geom_bar(aes(x = Order, y = n), stat = "identity", fill = "black") +
  theme_bw() +
  ylab("Species count") +
  scale_y_continuous(limits = c(0, 2500), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(), 
        text = element_text(size = 13), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

ggsave("Order_species_distribution.png", scale = 1, dpi = 350)

# build table of counts
pollinator_subset_family <- pollinator_subset %>% 
  group_by(Family, Order) %>%
  tally() %>%
  ungroup() %>%
  mutate(Family = forcats::fct_reorder(Family, -n)) %>%
  mutate(Order = factor(Order, levels = c("Lepidoptera", "Hymenoptera", "Coleoptera", "Diptera", "Thysanoptera"))) %>%
  ggplot() +
  geom_bar(aes(x = Family, y = n), stat = "identity", fill = "black") +
  theme_bw() +
  ylab("Species count") +
  facet_wrap(~Order, scales = "free_x") +
  scale_y_continuous(limits = c(0, 800), expand = c(0, 0)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(), 
        text = element_text(size = 12), 
        strip.text.x = element_text(size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(hjust = 1))
  
ggsave("Family_species_distribution.png", scale = 1.2, dpi = 350)
