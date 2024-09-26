# Wolbachia frequency calculation
library(tidyverse)

wolb <- load_data("Frequency of Wolbachia and Release Number.xlsx")

wolb <- wolb %>%
  mutate(year_week = sub("Y", "", year_week)) %>%  # Remove the 'Y' prefix
  separate(year_week, into = c("Year", "Epid Week"), sep = "/") %>%
  mutate(Year = as.character(Year), Year = paste0("20", Year))

time_points <- dengue[1:2]

time_points <- time_points %>%
  mutate(n = row_number() + 999)

class(wolb$Year) = class(time_points$Year)
class(wolb$`Epid Week`) = class(time_points$`Epid Week`)

wolb <- left_join(wolb, time_points, by = c("Year" = "Year", "Epid Week" = "Epid Week"))
sites_to_exclude <- c("W04", "W10", "W26", "W05", "W27", "W15")

# overall 
wolb_overall <- wolb %>%
  filter(!(site %in% sites_to_exclude)) %>%
  summarise(
    total_infected = sum(infected, na.rm = TRUE),     # Sum of infected cases across all sites
    total_uninfected = sum(uninfected, na.rm = TRUE)  # Sum of uninfected cases across all sites
  ) %>%
  mutate(
    wolbachia_frequency = (total_infected * 100) / (total_infected + total_uninfected)  # Calculate Wolbachia frequency
  )

#calendar time
wolb_ct <- wolb %>%
  filter(!(site %in% sites_to_exclude)) %>%  # Exclude specified sites
  group_by(Year) %>%
  summarise(
    total_infected = sum(infected),
    total_uninfected = sum(uninfected)
  ) %>%
  mutate(
    wolbachia_frequency = (total_infected * 100) / (total_infected + total_uninfected)
  )

#event time 
# take each time - date$time and divide by 
wolb_et <- left_join(wolb, dates, by = c("site" = "Area"))

dengue_map <- time_points %>% 
  filter(`Epid Week` == 52) %>% 
  group_by(Year) %>% 
  summarise(rownum = n)

wolb_et <- left_join(wolb_et, dengue_map, by = "Year")

wolb_et <- wolb_et %>%
  mutate(n = if_else(is.na(n), as.integer(rownum), n)) %>%
  select(-rownum) %>%
  mutate(time_period = as.integer(((n - time - 1) / 26) + 1)) 

wolb_et$time_period[wolb_et$time_period > 4] <- 5

wolb_et <- wolb_et %>%
  filter(!(site %in% sites_to_exclude)) %>% 
  group_by(time_period) %>%
  summarise(
    total_infected = sum(infected),
    total_uninfected = sum(uninfected)
  ) %>%
  mutate(
    wolbachia_frequency = (total_infected * 100) / (total_infected + total_uninfected)
  )





