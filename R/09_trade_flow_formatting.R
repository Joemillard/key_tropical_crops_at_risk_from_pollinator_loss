library(dplyr)
library(ggplot2)

# read in the trade flow data
trade_flow <- qs::qread(here::here("data/trade_flow/virtual-pollinators-flow.qs"))
str(trade_flow)  

# sum the trade data for each year and convert to proportion of each countries exports
prop_flow <- trade_flow %>%
  group_by(reporter_countries, partner_countries) %>%
  summarise(total_flow = sum(vp_flow)) %>%
  ungroup() %>%
  group_by(reporter_countries) %>%
  mutate(country_flow = sum(total_flow)) %>%
  mutate(percent_flow = (total_flow/country_flow) * 100)

# save proportional trade flow as an rds
saveRDS(prop_flow, "proportional_pollination_flow.rds")

  
