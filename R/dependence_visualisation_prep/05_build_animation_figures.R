library(dplyr)
library(ggplot2)

# read in the total production data
pollination_production <- readRDS("pollinator_dependence_visualisation/data/global_change_production.rds")

# sort the years
sorted_years <- sort(unique(pollination_production$year))

# bind together the outputs and plot as facetted plot for each scenario
for(i in 1:length(sorted_years)){
  pollination_production %>%
    filter(year <= sorted_years[i]) %>%
      ggplot() +
      geom_line(aes(x = year, y = vulnerability, colour = model, alpha = model)) +
      geom_point(aes(x = year, y = vulnerability, colour = model, alpha = model)) +
      facet_wrap(~scenario, ncol = 2) +
      scale_y_continuous(limits = c(1700000, 4100000), expand = c(0, 0), breaks = c(2000000, 2500000, 3000000, 3500000, 4000000), labels = c("2,000,000", "2,500,000", "3,000,000", "3,500,000", "4,000,000")) +
      scale_x_continuous(limits = c(2015, 2050), expand = c(0, 0), breaks = c(2020, 2030, 2040, 2050)) +
      scale_colour_manual("Climate model", values = c("black", "#E69F00", "#56B4E9", "#009E73", "#F0E442")) +
      scale_alpha_manual("Climate model", values = c(1, 0.4, 0.4, 0.4, 0.4)) +
      ylab("Vulnerability-weighted pollination prod. (metric tonnes)") +
      xlab("") +
      theme_bw() +
      theme(panel.grid = element_blank())
  
  ggsave(paste("animated_production", sorted_years[i], ".png", sep = ""), scale = 0.9, dpi = 300)
}

