library(dplyr)
library(ggplot2)

# read in text file from Felipe without isolation term
trade_flow_new <- read.csv("data/trade_flow/Virtual_pollination_flow_no_isolation.csv", stringsAsFactors = FALSE)

# sum the trade data for each year and convert to proportion of each countries exports
prop_flow <- trade_flow_new %>%
  group_by(reporter_countries, partner_countries) %>%
  summarise(total_flow = sum(potential_poll_trade, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(reporter_countries) %>%
  mutate(country_flow = sum(total_flow, na.rm = TRUE)) %>%
  mutate(percent_flow = (total_flow/country_flow) * 100)

# save proportional trade flow as an rds
saveRDS(prop_flow, "data/trade_flow/proportional_pollination_flow_non_isolation.rds")

  
