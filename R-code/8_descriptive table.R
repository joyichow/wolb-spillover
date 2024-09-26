# average calculations
intervention_sites <- as.vector(intervention$Site)
intervention_sites <- setdiff(intervention_sites, c("W04", "W10", "W26", "W05", "W27", "W15"))
# spillover_sites <- c("C35", "C39", "C53", "C70", "C03") # for avg
spillover_sites <- c("C43", "C35", "C61", "C39", "C48", "C62", "C25", "C53", "C70", "C03", "C73") # per site
control_sites <- setdiff(co,ss)
start_times <- dates

#dengue incidence
# intervened sites 
dengue <- dengue
dengue_interventions <- dengue[, c("Year", "Epid Week", intersect(intervention_sites, names(dengue)))]
start_times$time <- start_times$time-999

dengue_long <- pivot_longer(dengue_interventions, cols = -c(Year, `Epid Week`), names_to = "Area", values_to = "Incidence")
merged_dengue <- merge(start_times, dengue_long, by = "Area")

avg_incidence <- merged_dengue %>%
  group_by(Area) %>%
  arrange(Area, Year, `Epid Week`) %>%
  mutate(Row = row_number()) %>%  # Row number within each Area
  summarise(
    Pre_intervention_mean = mean(Incidence[Row < time], na.rm = TRUE),
    Post_intervention_mean = mean(Incidence[Row >= time], na.rm = TRUE),
    Pre_intervention_sd = sd(Incidence[Row < time], na.rm = TRUE),
    Post_intervention_sd = sd(Incidence[Row >= time], na.rm = TRUE),
  )

incidence_int_site <- avg_incidence %>%
  mutate(
    Pre = sprintf("%.2f (%.2f)", Pre_intervention_mean, Pre_intervention_sd),
    Post = sprintf("%.2f (%.2f)", Post_intervention_mean, Post_intervention_sd)
  ) %>%
  select(Area, Pre, Post)

write.csv(incidence_int_site, file = "pre_post_intervention.csv", row.names = FALSE)

overall_avg_intervention <- avg_incidence %>%
  summarise(
    Overall_Pre_intervention_mean = mean(Pre_intervention_mean),
    Overall_Post_intervention_mean = mean(Post_intervention_mean),
    Overall_Pre_intervention_sd = sd(Pre_intervention_mean),
    Overall_Post_intervention_sd = sd(Post_intervention_mean),
  )

#spillover sites
spillover_times <- c(rep(289, 4), 316)
spillover_times <- c(199, rep(289, 8), 316, 355)

names(spillover_times) <- spillover_sites

intervention_times <- do.call(rbind, lapply(names(spillover_times), function(name) {
  data.frame(Area = name, time = spillover_times[[name]])
}))

dengue_long <- pivot_longer(dengue, cols = -c(Year, `Epid Week`), names_to = "Area", values_to = "Incidence")

df_incidence_long_spillover <- dengue_long %>%
  filter(Area %in% spillover_sites)

merged_data_spillover <- merge(df_incidence_long_spillover, intervention_times, by = "Area")

averages_per_spillover_site <- merged_data_spillover %>%
  group_by(Area) %>%
  arrange(Area, Year, `Epid Week`) %>%
  mutate(Row = row_number()) %>%
  summarise(
    Pre_intervention_mean = mean(Incidence[Row < time], na.rm = TRUE),
    Post_intervention_mean = mean(Incidence[Row >= time], na.rm = TRUE),
    Pre_intervention_sd = sd(Incidence[Row < time], na.rm = TRUE),
    Post_intervention_sd = sd(Incidence[Row >= time], na.rm = TRUE)
  )

incidence_spill_site <- averages_per_spillover_site %>%
  mutate(
    Pre = sprintf("%.2f (%.2f)", Pre_intervention_mean, Pre_intervention_sd),
    Post = sprintf("%.2f (%.2f)", Post_intervention_mean, Post_intervention_sd)
  ) %>%
  select(Area, Pre, Post)

write.csv(incidence_spill_site , file = "pre_post_spillover_actual.csv", row.names = FALSE)

overall_statistics_spillover <- averages_per_spillover_site %>%
  summarise(
    Overall_Pre_intervention_mean = mean(Pre_intervention_mean),
    Overall_Post_intervention_mean = mean(Post_intervention_mean),
    Overall_Pre_intervention_sd = sd(Pre_intervention_mean),
    Overall_Post_intervention_sd = sd(Post_intervention_mean)
  )

#control sites
df_incidence_long_control <- dengue_long  %>%
  filter(Area %in% control_sites)

control_intervention_time <- 176

df_incidence_long_control$time <- control_intervention_time 

control_statistics <- df_incidence_long_control %>%
  arrange(Area, Year, `Epid Week`) %>%
  summarise(
    Overall_Mean_Incidence = mean(Incidence, na.rm = TRUE),
    Overall_SD_Incidence = sd(Incidence, na.rm = TRUE)
  )

averages_per_control_site <- df_incidence_long_control %>%
  group_by(Area) %>%
  mutate(Row = row_number()) %>%
  summarise(
    Pre_intervention_mean = mean(Incidence[Row < time], na.rm = TRUE),
    Post_intervention_mean = mean(Incidence[Row >= time], na.rm = TRUE),
    Pre_intervention_sd = sd(Incidence[Row < time], na.rm = TRUE),
    Post_intervention_sd = sd(Incidence[Row >= time], na.rm = TRUE),
  )

overall_statistics_control <- averages_per_control_site %>%
  summarise(
    Overall_Pre_intervention_mean = mean(Pre_intervention_mean),
    Overall_Post_intervention_mean = mean(Post_intervention_mean),
    Overall_Pre_intervention_sd = sd(Pre_intervention_mean),
    Overall_Post_intervention_sd = sd(Post_intervention_mean)
  )

#covariates 
avg_cov <- read.csv("/Users/your-local-file-path/Malasia_variables_mean_2010-2024.csv")

intervention_data <- filter(avg_cov, Site_ID %in% intervention_sites)
spillover_data <- filter(avg_cov, Site_ID %in% spillover_sites)
control_data <- filter(avg_cov, Site_ID %in% control_sites)

intervention_stats <- cbind(Type = "Intervention", calculate_stats(intervention_data))
spillover_stats <- cbind(Type = "Spillover", calculate_stats(spillover_data))
control_stats <- cbind(Type = "Control", calculate_stats(control_data))

final_cov <- t(rbind(intervention_stats, spillover_stats, control_stats))

write_xlsx(as.data.frame(final_cov), "avg_cov.xlsx")

# synthetic control
variables <- c("t", "t2m", "tp", "r", "cc", "cvl", "cvh", "lai_lv", "lai_hv")
var_list <- list()

for (var in variables) {
  var_df <- avg %>%
    select(Site_ID, sno, !!sym(var)) %>%
    pivot_wider(names_from = Site_ID, values_from = !!sym(var))
  
  colnames(var_df) <- c("time", avg$Site_ID)
  
  var_list [[var]] <- var_df
}

colnames(sm_weights_df) <- spill.ie$site
sm_weights_df <- sm_weights_df[,-c(10)]
sm_weights_df <- sm_weights_df[,-c(2, 4, 5, 6, 10)]
sm_weights_df <- sm_weights_df[,-c(6)]

colnames(dm_weights_df) <- intervention$Site
dm_weights_df <- dm_weights_df[, !colnames(dm_weights_df) %in% c("W04", "W10", "W26", "W05", "W27", "W15")]

# Get the average 
result_list <- map(var_list, ~ calculate_weighted_avgs(.x, sm_weights_df)) # spillover

combined_df <- imap_dfr(result_list, function(df, i) {
  df %>%
    mutate(Variable = paste("Variable", i)) %>%  # Add a variable identifier
    pivot_longer(cols = -c(time, Variable), names_to = "Site", values_to = "Value")  # Reshape the dataframe
})

summary_df <- combined_df %>%
  group_by(Variable) %>%
  summarise(
    `Mean (SD)` = sprintf("%.2f (%.2f)", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE))
  )

# intervention covariates for SC 
intervention_list <- map(var_list, ~ calculate_weighted_avgs(.x, dm_weights_df)) # intervention

intervention_df <- imap_dfr(intervention_list, function(df, i) {
  df %>%
    mutate(Variable = paste("Variable", i)) %>%  # Add a variable identifier
    pivot_longer(cols = -c(time, Variable), names_to = "Site", values_to = "Value")  # Reshape the dataframe
})

summary_int_cov <- intervention_df  %>%
  group_by(Variable) %>%
  summarise(
    `Mean (SD)` = sprintf("%.2f (%.2f)", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE))
  )

#SC incidence
direct_sc <- as.data.frame(direct_sc)
direct_sc$time <- 1:nrow(direct_sc)
direct_sc_long <- pivot_longer(direct_sc, cols = -time, names_to = "Area", values_to = "Incidence")

# Merge with start_times
merged_direct_sc <- merge(start_times, direct_sc_long, by = "Area")

merged_direct_sc <- merged_direct_sc %>%
  rename(intervention_time = time.x, row_time = time.y) 
  # %>% filter(!Area %in% c("W04", "W10", "W26", "W05", "W27", "W15")) # for all 

# Calculate pre-intervention and post-intervention means and standard deviations
avg_incidence <- merged_direct_sc %>%
  group_by(Area) %>%
  summarise(
    Pre_intervention_mean = mean(Incidence[row_time < intervention_time], na.rm = TRUE),
    Post_intervention_mean = mean(Incidence[row_time >= intervention_time], na.rm = TRUE),
    Pre_intervention_sd = sd(Incidence[row_time < intervention_time], na.rm = TRUE),
    Post_intervention_sd = sd(Incidence[row_time >= intervention_time], na.rm = TRUE)
  )

incidence_intervened_site <- avg_incidence  %>%
  mutate(
    Pre = sprintf("%.2f (%.2f)", Pre_intervention_mean, Pre_intervention_sd),
    Post = sprintf("%.2f (%.2f)", Post_intervention_mean, Post_intervention_sd)
  ) %>%
  select(Area, Pre, Post)

write.csv(incidence_intervened_site , file = "pre_post_intervention_sc.csv", row.names = FALSE)

# Calculate overall statistics
overall_avg_intervention <- avg_incidence %>%
  summarise(
    Overall_Pre_intervention_mean = mean(Pre_intervention_mean, na.rm = TRUE),
    Overall_Post_intervention_mean = mean(Post_intervention_mean, na.rm = TRUE),
    Overall_Pre_intervention_sd = sd(Pre_intervention_mean, na.rm = TRUE),
    Overall_Post_intervention_sd = sd(Post_intervention_mean, na.rm = TRUE)
  )


# spillover sites
spillover_times_df <- data.frame(                        # for overall 
  Area = c("C35", "C39", "C53", "C70", "C03"),
  time = c(rep(289, 4), 316)
)

spillover_times_df <- data.frame(                        # per unit
  Area = c("C43", "C35", "C61", "C39", "C48", "C62", "C25", "C53", "C70", "C03", "C73"), 
  time = c(199, rep(289, 8), 316, 355)
)

spillover_sc$time <- 1:nrow(spillover_sc)
colnames(spillover_sc) <- c(unique(spill.ie$site), "time")
spillover_sc_long <- pivot_longer(spillover_sc, cols = -time, names_to = "Area", values_to = "Incidence")

merged_spillover_sc <- merge(spillover_times_df, spillover_sc_long, by = "Area")

merged_spillover_sc <- merged_spillover_sc %>%
  rename(intervention_time = time.x, row_time = time.y)

# Calculate pre-intervention and post-intervention means and standard deviations
avg_incidence <- merged_spillover_sc %>%
  group_by(Area) %>%
  summarise(
    Pre_intervention_mean = mean(Incidence[row_time < intervention_time], na.rm = TRUE),
    Post_intervention_mean = mean(Incidence[row_time >= intervention_time], na.rm = TRUE),
    Pre_intervention_sd = sd(Incidence[row_time < intervention_time], na.rm = TRUE),
    Post_intervention_sd = sd(Incidence[row_time >= intervention_time], na.rm = TRUE)
  )

incidence_spill_site <- avg_incidence  %>%
  mutate(
    Pre = sprintf("%.2f (%.2f)", Pre_intervention_mean, Pre_intervention_sd),
    Post = sprintf("%.2f (%.2f)", Post_intervention_mean, Post_intervention_sd)
  ) %>%
  select(Area, Pre, Post)

write.csv(incidence_spill_site , file = "pre_post_spillover_sc.csv", row.names = FALSE)

# Calculate overall statistics
overall_avg_intervention <- avg_incidence %>%
  summarise(
    Overall_Pre_intervention_mean = mean(Pre_intervention_mean, na.rm = TRUE),
    Overall_Post_intervention_mean = mean(Post_intervention_mean, na.rm = TRUE),
    Overall_Pre_intervention_sd = sd(Pre_intervention_mean, na.rm = TRUE),
    Overall_Post_intervention_sd = sd(Post_intervention_mean, na.rm = TRUE)
  )

# Print results
print(avg_incidence)
print(overall_avg_intervention)

# Average intervention time calculation for intervened sites 
calc_intervention_week <- function(val) {
  val <- 444-val+1
  return(val)
}

intervention_vals <- c(176, 199, 204, 271, rep(289, 7), rep(317, 2), rep(355,2))
count <- 0

for (i in intervention_vals){
  count <- count + calc_intervention_week(i)
  print(calc_intervention_week(i))
}

avg_intervention_time <- count/length(intervention_vals)

# Avg intervention time for spillover sites
spillover_vals <- c(rep(289, 4), 316)
count <- 0

for (j in spillover_vals){
  count <- count + calc_intervention_week(j)
}

avg_spillover_time <- count/length(spillover_vals)

# Average intervention time for controls 
control_vals <- c(rep(176, 65))
count <- 0

for (k in control_vals){
  count <- count + calc_intervention_week(k)
  print(count)
}

avg_control_time <- count/length(control_vals)

