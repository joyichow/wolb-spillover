# Load data 
load_data <- function(file_name) {
  file_path <- file.path(data_path, file_name)
  if (file.exists(file_path)) {
    return(read_xlsx(file_path))
  } else {
    stop("File does not exist: ", file_path)
  }
}

# Covariates
times.dep.func <- function(start, end){
  times.pred <- cbind("t" =  c(start,end),
                      "t2m" = c(start, end),
                      "tp" = c(start, end),
                      "r" = c(start, end),
                      "cc" = c(start, end),
                      "cvl" = c(start, end),
                      "cvh" = c(start, end),
                      "lai_lv" = c(start, end),
                      "lai_hv" = c(start, end)
  )
  return(times.pred)
}

# Finding spillover sites
prox.pairs <- function(distance, intervention, control) {
  results <- expand.grid(intervention = 1:nrow(intervention), control = 1:nrow(control))
  distances <- mapply(function(i, j) {
    geosphere::distHaversine(c(intervention$Longitude[i], intervention$Latitude[i]), 
                             c(control$Longitude[j], control$Latitude[j]))
  }, results$intervention, results$control)
  
  within_threshold <- distances <= distance
  data_frame(
    Study_PA = intervention[results$intervention[within_threshold], 1],
    Spillover = control[results$control[within_threshold], 1],
    distance = distances[within_threshold]
  )
}

# Getting pool size of spillover
spillover.size <- function(start_date) {
  
  spill.size <- ifelse(start_date >= 1198 & start_date < 1288, 1,
                       ifelse(start_date >=1288 & start_date <1315, 9,
                              ifelse(start_date >=1315 & start_date <1354, 10, 
                                     ifelse(start_date >= 1354, 12, 0))))
  
  return(spill.size)
}

# Analysis by event time
scm.et <- function(obs, sc, type){
  df <- data.frame(obs,sc)
  df$event_time <- (as.integer((seq_along(df$obs) - 1) / 26) + 1)
  df$event_time[df$event_time > 4] <- 5
  
  results <- df %>%
    group_by(event_time) %>%
    summarise(
      sum_obs = sum(obs, na.rm = TRUE),  # Sum of obs
      sum_sc = sum(sc, na.rm = TRUE),    # Sum of sc
      metric = case_when(
        type == "ie" ~ (sum_sc - sum_obs) * 100 / sum_sc,
        TRUE ~ sum(sc - obs)
      ),
      .groups = 'drop'
    ) %>%
    complete(event_time = 1:5, fill = list(obs = list(NULL), sc = list(NULL), sum_obs = 0, sum_sc = 0, metric = NA))
  
  return(results)
}

# Analysis by calendar time
scm.ct <- function(obs, sc, type, start){
  df <- data.frame(obs, sc)
  rownames(df) <- c(start:444)
  
  breaks <- c(175, 209, 261, 313, 366, 418, 444)
  labels <- c("2017", "2018", "2019", "2020", "2021", "2022")
  
  year_factor <- cut(as.numeric(rownames(df)), breaks=breaks, labels=labels, include.lowest=TRUE)
  df$year <- as.numeric(as.character(year_factor))
  
  results <- df %>%
    group_by(year) %>%
    summarise(
      sum_obs = sum(obs, na.rm = TRUE),  # Sum of obs
      sum_sc = sum(sc, na.rm = TRUE),    # Sum of sc
      metric = case_when(
        type == "ie" ~ (sum_sc - sum_obs) * 100 / sum_sc,
        TRUE ~ sum(sc - obs)
      ),
      .groups = 'drop'
    ) %>%
    complete(year = 2017:2022, fill = list(obs = list(NULL), sc = list(NULL), sum_obs = 0, sum_sc = 0, metric = NA))
  
  return(results)
}

# Main analysis 
partial_scm <- function(ind, data, co2, start_date, population) {
  times.dep  <- times.pred <- cbind("cases"  = c(1000, start_date))
  
  dat <- cbind(data[ind], data %>% select(co2))
  rownames(dat) <- c(1000:1443)
  
  l.dat <- list(cases = as.matrix(dat))
  
  # add to mscmt: agg.fns
  model <- mscmt(l.dat, ind, co2, times.dep, times.pred, seed = 1, verbose = FALSE)
  w <- model[["w"]]
  
  all_obs <- model[["data.synth"]][["cases"]]
  
  start <-  (as.integer(start_date) - 999)[1]
  
  obs <- (model$combined$cases[,1][(start):444])
  sc <- (model$combined$cases[,2][(start):444])
  dir <- (sum(sc)-sum(obs))*100/sum(sc)
  
  et <- scm.et(obs, sc, "ie")
  ct <- scm.ct(obs, sc, "ie", start)
  
  rmse <- sqrt((sum((model$combined$cases[,1][1:start] - model$combined$cases[,2][1:start])^2))/start)
  
  # absolute cases averted 
  pop <- (population$population[population$site == ind])/100000
  aca <- sum(round(sc*pop, digits = 0) - round(obs*pop, digits = 0))
  
  aca.et <- scm.et(round(obs*pop, digits = 0), round(sc*pop, digits = 0), "aca")
  aca.ct <- scm.ct(round(obs*pop, digits = 0), round(sc*pop, digits = 0), "aca", start)
  
  begin <- as.Date(paste0(2014, "-01-01"))
  date_seq <- seq.Date(begin, by = "week", length.out = nrow(model$combined$cases))
  
  plot_obj <- ggplot() +
    geom_line(aes(x = date_seq, y = model$combined$cases[,1], color = "Observed"), size = 0.5) +
    geom_line(aes(x = date_seq, y = model$combined$cases[,2], color = "Synthetic Control"), size = 0.5) +
    geom_vline(aes(xintercept = as.numeric(date_seq[start]), color = "Intervention Start"), size = 1, linetype = "dashed") +
    scale_x_date(breaks = seq.Date(from = begin, by = "year", length.out = length(unique(format(date_seq, "%Y")))),
                 labels = unique(format(date_seq, "%Y")),
                 expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = paste0(ind, " (RMSE = ", sprintf("%.2f", rmse), ")"),
         y = NULL, x = NULL) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(color = "black"),
          axis.ticks = element_line(color = "black"),       
          axis.ticks.length = unit(0.25, "cm"),               
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5),           
          legend.justification = c(0.05, 1),                    
          legend.position = c(0.05, 1),                    
          legend.title = element_blank(),
          legend.background = element_blank(),                
          legend.key = element_blank(),                     
          legend.text = element_text(size = 11),          
          legend.spacing.y = unit(0.5, "cm"),              
          legend.key.size = unit(1.1, "lines")) +            
    scale_color_manual(values = c("Observed" = "darkgrey", "Synthetic Control" = "red", "Intervention Start" = "black"))
  
  return(list(rmse = rmse, dir = dir, obs = sum(obs), sc = sum(sc), 
              aca = aca, et = et, aca.et = aca.et, ct = ct, aca.ct = aca.ct, w = w, all_obs = all_obs, plot = plot_obj))
}

# Covariate analysis 
covariate_scm <- function(type, ind, co2, start_date) {
  
  start_date <- as.integer(start_date) 
  
  times.dep  <- cbind("cases"  = c(1000, start_date))
  times.pred <- cbind(times.dep.func(1000, start_date))
  agg.fns <- rep(type, ncol(times.pred)) #mean, or id
  
  # add to mscmt: agg.fns
  model <- mscmt(Dengue, ind, co2, times.dep, times.pred, agg.fns, seed = 1, verbose = TRUE)
  cov_weights <- model[["v"]][,2]
  
  start <-  as.integer(start_date) - 999
  
  rmse <- sqrt((sum((model$combined$cases[,1][1:start] - model$combined$cases[,2][1:start])^2))/start)
  
  return(list(rmse = rmse, cov_weights = cov_weights))
}

# ie calculation 
calc.ie <- function(df){
  df$sc <- as.numeric(df$sc)
  df$obs <- as.numeric(df$obs)
  ie <- (sum(df$sc) - sum(df$obs)) * 100 / sum(df$sc)
  ie
}

# Summarising the spillover site dataframes 
combined.sites <- function(data_list) {
  combined_df <- lapply(names(data_list), function(name) {
    mutate(data_list[[name]], treated.site = name)
  }) %>%
    bind_rows()
  
  return(combined_df)
}

# Turning columns into numerical and getting the col sums for event time aggregations
con.sum <- function(df, n) {
  df[, 3:n] <- sapply(df[, 3:n], function(x) {
    if(is.factor(x) || is.character(x)) {
      as.numeric(as.character(x))
    } else {
      x
    }
  })
  
  col_sums <- colSums(df[, 3:n], na.rm = TRUE)
  return(col_sums)
}

# Filtering RMSE >100 in CI analysis 
dir.filter <- function(df) {
  filtered_data <- df[!df$site %in% c("W04", "W10", "W15" ,"W26"), ]
  filtered_data <- subset(filtered_data, as.numeric(rmse) < 100)
  
  return(filtered_data)
}
spill.filter <- function(df) {
  filtered_data <- df[!df$site %in% c("C43", "C48", "C25", "C61", "C73"), ]
  filtered_data <- subset(filtered_data, as.numeric(rmse) < 100)
  
  return(filtered_data)
}

# Getting CIs from each analysis 
calc.quantiles <- function(df) {
  quantiles <- apply(df, 2, function(column) quantile(column, probs = c(0.025, 0.975)))
  formatted_quantiles <- apply(quantiles, 2, function(q) {
    sprintf("(%.2f - %.2f)", q[1], q[2])
  })
  return(formatted_quantiles)
}

# Time placebo 
p.spillover.size <- function(start_date) {
  
  spill.size <- ifelse(start_date >= 1094 & start_date < 1184, 1,
                       ifelse(start_date >=1184 & start_date <1211, 9,
                              ifelse(start_date >=1211 & start_date <1250, 10, 
                                     ifelse(start_date >= 1250, 12, 0))))
  
  return(spill.size)
}

placebo.scm <- function(ind, data, co2, start_date, population) {
  times.dep  <- times.pred <- cbind("cases"  = c(1000, start_date))
  dat <- cbind(data[ind], data %>% select(co2))
  rownames(dat) <- c(1000:1443)
  
  l.dat <- list(cases = as.matrix(dat))
  # add to mscmt: agg.fns for covariates 
  model <- mscmt(l.dat, ind, co2, times.dep, times.pred, seed = 1, verbose = FALSE)
  
  start <- (as.integer(start_date) - 999)[1] #pre-intervention time-points 
  
  obs <- (model$combined$cases[,1][(start):(start+103)])
  sc <- (model$combined$cases[,2][(start):(start+103)])
  dir <- (sum(sc)-sum(obs))*100/sum(sc)
  
  rmse <- sqrt((sum((model$combined$cases[,1][1:start] - model$combined$cases[,2][1:start])^2))/start)
  p.error <- sum((model$combined$cases[,1][1:start] - model$combined$cases[,2][1:start])/model$combined$cases[,1][1:start])/start
  
  aca <- sum(round(sc*pop, digits = 0) - round(obs*pop, digits = 0))
  
  return(list(p.error = p.error, rmse = rmse, dir = dir, obs = sum(obs), sc = sum(sc), aca = sum(aca)))
}

# calculate averages and sd for covariates 
calculate_stats <- function(data) {
  data %>%
    summarise(
      t_mean = mean(t, na.rm = TRUE), t_sd = sd(t, na.rm = TRUE),
      t2m_mean = mean(t2m, na.rm = TRUE), t2m_sd = sd(t2m, na.rm = TRUE),
      tp_mean = mean(tp, na.rm = TRUE), tp_sd = sd(tp, na.rm = TRUE),
      r_mean = mean(r, na.rm = TRUE), r_sd = sd(r, na.rm = TRUE),
      cc_mean = mean(cc, na.rm = TRUE), cc_sd = sd(cc, na.rm = TRUE),
      cvl_mean = mean(cvl, na.rm = TRUE), cvl_sd = sd(cvl, na.rm = TRUE),
      cvh_mean = mean(cvh, na.rm = TRUE), cvh_sd = sd(cvh, na.rm = TRUE),
      lai_lv_mean = mean(lai_lv, na.rm = TRUE), lai_lv_sd = sd(lai_lv, na.rm = TRUE),
      lai_hv_mean = mean(lai_hv, na.rm = TRUE), lai_hv_sd = sd(lai_hv, na.rm = TRUE)
    ) %>%
    mutate(
      t = sprintf("%.2f (%.2f)", t_mean, t_sd),
      t2m = sprintf("%.2f (%.2f)", t2m_mean, t2m_sd),
      tp = sprintf("%.2f (%.2f)", tp_mean, tp_sd),
      r = sprintf("%.2f (%.2f)", r_mean, r_sd),
      cc = sprintf("%.2f (%.2f)", cc_mean, cc_sd),
      cvl = sprintf("%.2f (%.2f)", cvl_mean, cvl_sd),
      cvh = sprintf("%.2f (%.2f)", cvh_mean, cvh_sd),
      lai_lv = sprintf("%.2f (%.2f)", lai_lv_mean, lai_lv_sd),
      lai_hv = sprintf("%.2f (%.2f)", lai_hv_mean, lai_hv_sd)
    ) %>%
    select(t, t2m, tp, r, cc, cvl, cvh, lai_lv, lai_hv)
}

# calculate average for SC 
weighted_average <- function(clim_data, weights, site_name) {
  
  names(weights) <- rownames(sm_weights_df)
  
  weighted_avg <- clim_data %>% 
    mutate(across(starts_with("C"), ~ .x * weights[match(cur_column(), names(weights))])) %>%
    rowwise() %>%
    mutate(!!site_name := sum(c_across(starts_with("C")), na.rm = TRUE)) %>%
    select(all_of(site_name))
  
  return(weighted_avg)
}

