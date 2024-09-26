library(tidyverse)
library(stringr)
library(dplyr)
library(MSCMT)

# data cleaning for covariate data 
max <- read.csv("/Users/your-local-file/Malasia_variables_max_2010-2024.csv")
min <- read.csv("/Users/your-local-file/Malasia_variables_min_2010-2024.csv")
avg <- read.csv("/Users/your-local-file/Malasia_variables_mean_2010-2024.csv")

avg <- avg %>% # change here to min/max/avg
  mutate(year_epi = str_replace(year_epi, "_(\\d)$", "_0\\1")) %>% # Pad single-digit weeks with a leading zero
  arrange(year_epi) %>%
  filter(year_epi >= "2014_01" & year_epi <= "2022_26") %>% # trimming data to 2014_01 due to missing incidence data 
  mutate(sno = match(year_epi, unique(year_epi)) + 999) %>%
  mutate(Site_ID = case_when(Site_ID == "C76b" ~ "C76", Site_ID == "C78a" ~ "C78", Site_ID == "C79c" ~ "C79", TRUE ~ Site_ID))

# cleaning dengue incidence data 
dengue <- read_xlsx("/Users/your-local-file/data submission incidence.xlsx")
dengue <- dengue[-c(1:2, 18)]
dengue <- dengue[-c(1:52),]
dengue$sno <- 1000:(999 + nrow(dengue))

l.deng <- dengue %>%
  pivot_longer(cols = -sno, names_to = "site", values_to = "cases") %>%
  mutate(site.no = as.numeric(factor(site, levels = unique(site)))) %>%
  select(site.no, site, sno, cases)

# merging avg and l.deng 
m.deng <- as.data.frame(left_join(l.deng, avg, by = c("site" = "Site_ID", "sno" = "sno")))
m.deng <- select(m.deng, -year_epi)

# converting to list format 
Dengue <- listFromLong(m.deng, unit.variable="site.no", time.variable="sno", unit.names.variable="site")

# reset from here 
dweights_df <- data.frame()
sweights_df <- list()
processed_dates <- c()

type = "last" #mean, last 

for (i in dates$Area){
  start_date <- as.integer(dates$time[dates$Area == i]) - 3
  spill.size <- spillover.size(start_date)
  curr.ss <- ss[1:spill.size]
  
  #calculating direct effects
  direct <- covariate_scm(type, i, co2, start_date)
  w <- t(direct$cov_weights)
  colnames(w) = c("t", "t2m", "tp", "r", "cc", "cvl", "cvh", "lai_lv", "lai_hv")
  dweights_df <- rbind(dweights_df, data.frame(rmse = direct$rmse, w))
  
  if (start_date %in% processed_dates || !start_date %in% c(1198, 1288, 1315, 1354)) {
    next
  }
  
  sweights <- data.frame()
  processed_dates <- c(processed_dates, start_date)
  
  # spillover effects
  if (!is.null(s.date.map[[as.character(start_date)]])) {
    ind <- s.date.map[[as.character(start_date)]]
  } else {
    ind <- 1:spill.size
  }
  
  for (j in ind) {
    spillover <- covariate_scm(type, curr.ss[j], co2, start_date)
    
    sw <- t(spillover$cov_weights)
    colnames(sw) = c("t", "t2m", "tp", "r", "cc", "cvl", "cvh", "lai_lv", "lai_hv")
    sweights <- rbind(sweights, data.frame(rmse = spillover$rmse, sw))
  }
  
  sweights_df[[i]] <- as.data.frame(do.call(cbind, sweights))
}

sweights_df <- combined.sites(sweights_df)

write_xlsx(dweights_df, "mean_cov_direct.xlsx")
write_xlsx(sweights_df , "mean_cov_spilled.xlsx")


