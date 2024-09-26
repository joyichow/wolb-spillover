# Method 2 for spillover without CIs
if (!require("pacman")) install.packages("pacman")
pacman::p_load("readxl", "geosphere", "MSCMT", "dplyr", "writexl", "tidyverse", 
               "ggplot2", "Synth", "glmnet", "osqp", "optimx", "meboot", "here", 
               "gridExtra", "cowplot")

base_path <- here::here("")
fxn_path <- file.path(base_path, "R")
data_path <- file.path(base_path, "Data")

# Load datasets
intervention <- load_data("intervention.xlsx")
control <- load_data("control.xlsx")
dengue <- as.data.frame(load_data("data submission incidence.xlsx")[-c(1:52), -c(18)]) # remove first 52 rows due to missing data, W20 missing 
dates <- load_data("dates.xlsx")
controls <- load_data("co.xlsx")
population <- as.data.frame(load_data("population.xlsx"))

# Finding intervention and control pairs within 1500m proximity
dist <- prox.pairs(1500, intervention, control) 

spillover <- as.vector(unlist(dist[,2])) # define spillover sites
co <- as.vector(unlist(controls[,1])) # define control pool 

# Define mapping for synthetic control estimation
mapT <- data.frame(sno = 1:444)
mapT$time <- sprintf("%dW%02d", dengue$Year, dengue$`Epid Week`)
mapT$ey.ew <- sprintf("%d.%02d", dengue$Year, dengue$`Epid Week`)
mapT$map <- seq(1000,1443)

# Match dates to serial numbers
dates$time <- sapply(dates$Date, function(x) mapT$map[mapT$time == x]) #scm

# Estimate synthetic control
ss <- c("C43", "C35", "C61", "C39", "C48", "C62", "C25", "C53", "C70", "C03", "C62", "C73") # sites with spillover ordered to start time 
co2 <- setdiff(co, ss)
#ss <- c() # without spillover 

s.date.map <- list(
  `1198` = c(1), 
  `1288` = c(2:9),       
  `1315` = c(10),      
  `1354` = c(11:12)       
)

#reset from here
processed_dates <- c()

# site-aggregated effects 
direct.ie <- data.frame()
spill.ie <- list()

# event time 
ae.sc <- data.frame()
ae.obs <- data.frame()
ae.aca <- data.frame()

spill.aca <- list()
spill.sc <- list()
spill.obs <- list()

# calendar time 
dct.sc <- data.frame()
dct.obs <- data.frame()
dct.aca <- data.frame()

sct.aca <- list()
sct.sc <- list()
sct.obs <- list()

#weights
dm_weights_df <- data.frame()
sm_weights_df <- list()

#plots 
plot_list <- list()
pdf("all-plots.pdf", width = 11.69, height = 8.27)

direct_sc <- data.frame(matrix(ncol = 0, nrow = 444))
spillover_sc <- data.frame(matrix(ncol = 0, nrow = 444))

for (i in dates$Area){
  start_date <- as.integer(dates$time[dates$Area == i])
  spill.size <- spillover.size(start_date)
  curr.ss <- ss[1:spill.size]
  
  #calculating direct effects
  direct <- partial_scm(i, dengue, co2, start_date, population)
  w <- t(direct$w)
  colnames(w) = co2
  dm_weights_df <- rbind(dm_weights_df, data.frame(w))

  direct.ie <- rbind(direct.ie, 
                     data.frame(site = i, dir = direct$dir, rmse = direct$rmse,
                                obs = direct$obs, sc = direct$sc, aca = direct$aca))
  
  direct_sc <- cbind(direct_sc, all_obs = direct$all_obs)
  
  plot_list <- c(plot_list, list(direct$plot))
  
  # event time 
  ae.sc <- rbind(ae.sc, data.frame(site = i, rmse = direct$rmse, t(direct$et$sum_sc)))
  ae.obs <- rbind(ae.obs, data.frame(site = i, rmse = direct$rmse, t(direct$et$sum_obs)))
  ae.aca <- rbind(ae.aca, data.frame(site = i, rmse = direct$rmse, t(direct$aca.et$metric)))
  
  dct.sc <- rbind(dct.sc, data.frame(site = i, rmse = direct$rmse, t(direct$ct$sum_sc)))
  dct.obs <- rbind(dct.obs, data.frame(site = i, rmse = direct$rmse, t(direct$ct$sum_obs)))
  dct.aca <- rbind(dct.aca, data.frame(site = i, rmse = direct$rmse, t(direct$aca.ct$metric)))
  
  if (start_date %in% processed_dates || !start_date %in% c(1198, 1288, 1315, 1354)) {
    next
  }
  
  spill.df <- data.frame(site = character(), dir = numeric(), rmse = numeric(),
                         obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)
  
  sp.sc.df <- data.frame()
  sp.obs.df <- data.frame()
  sp.aca.df <- data.frame()
  
  sct.sc.df <- data.frame()
  sct.obs.df <- data.frame()
  sct.aca.df <- data.frame()

  sm_weights <- data.frame()
  processed_dates <- c(processed_dates, start_date)
  
  # spillover effects
  if (!is.null(s.date.map[[as.character(start_date)]])) {
    ind <- s.date.map[[as.character(start_date)]]
  } else {
    ind <- 1:spill.size
  }
  
  for (j in ind) {
    spillover <- partial_scm(curr.ss[j], dengue, co2, start_date, population)
    spill.df <- rbind(spill.df,
                      data.frame(site = curr.ss[j], dir = spillover$dir, rmse = spillover$rmse, 
                                 obs = spillover$obs, sc = spillover$sc, aca = spillover$aca))
    
    sw <- t(spillover$w)
    colnames(sw) = co2
    sm_weights <- rbind(sm_weights, data.frame(sw))
    
    spillover_sc <- cbind(spillover_sc, all_obs = spillover$all_obs)
    
    sp.sc.df <- rbind(sp.sc.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$et$sum_sc)))
    sp.obs.df <- rbind(sp.obs.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$et$sum_obs)))
    sp.aca.df <- rbind(sp.aca.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$aca.et$metric)))
    
    sct.sc.df <- rbind(sct.sc.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$ct$sum_sc)))
    sct.obs.df <- rbind(sct.obs.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$ct$sum_obs)))
    sct.aca.df <- rbind(sct.aca.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$aca.ct$metric)))
    
    plot_list <- c(plot_list, list(spillover$plot))
  }
  
  spill.ie[[i]] <- as.data.frame(do.call(cbind, spill.df))
  
  spill.aca[[i]] <- as.data.frame(do.call(cbind, sp.aca.df))
  spill.obs[[i]] <- as.data.frame(do.call(cbind, sp.obs.df))
  spill.sc[[i]] <- as.data.frame(do.call(cbind, sp.sc.df))
  
  sct.aca[[i]] <- as.data.frame(do.call(cbind, sct.aca.df))
  sct.obs[[i]] <- as.data.frame(do.call(cbind, sct.obs.df))
  sct.sc[[i]] <- as.data.frame(do.call(cbind, sct.sc.df))
  
  sm_weights_df[[i]] <- as.data.frame(do.call(cbind, sm_weights))
}

spill_plots <- c(plot_list[2:9], plot_list[20], plot_list[25])
direct_plots <- plot_list[-c(2:9, 20, 24:25)]
all_plots <- c(direct_plots, spill_plots)

for (i in seq(1, length(all_plots), by = 6)) {
  do.call(grid.arrange, c(all_plots[i:min(i+5, length(all_plots))], ncol = 3, nrow = 2))
}
dev.off()
# getting sc
colnames(direct_sc) <- intervention$Site
colnames(spillover_sc) <- as.vector(spill.ie$site)
spillover_sc <- spillover_sc[,-c(10)]

# Combining results 
spill.ie <- combined.sites(spill.ie)
spill.aca <- combined.sites(spill.aca)
spill.sc <- combined.sites(spill.sc)
spill.obs <- combined.sites(spill.obs)
sm_weights_df <- combined.sites(sm_weights_df)
sm_weights_df <- t(sm_weights_df[1:65])

sct.aca <- combined.sites(sct.aca)
sct.sc <- combined.sites(sct.sc)
sct.obs <- combined.sites(sct.obs)

# Eliminating fits worse than RMSE 100 and those that failed placebo
dirdf <- dir.filter(direct.ie)
spilldf <- spill.filter(spill.ie)

# Aggregate effects
dir.ie <- calc.ie(dirdf)                  # direct IE
s.ie <- calc.ie(spilldf)                  # spillover IE

dir.aca <- sum(dirdf$aca)                 # direct ACA
spilled.aca <- sum(as.numeric(spilldf$aca)) # spillover ACA
tot.aca <- dir.aca + spilled.aca            # total ACA 

# Save files
write_xlsx(direct.ie, "directed.xlsx")
write_xlsx(spill.ie, "spillover.xlsx")
write_xlsx(as.data.frame(sm_weights_df), "weights_spillover.xlsx")
write_xlsx(as.data.frame(t(dm_weights_df)), "weights_direct.xlsx")

# event time aggregations
spill.aca <- spill.filter(spill.aca)
spill.obs <- spill.filter(spill.obs)
spill.sc <- spill.filter(spill.sc)

ae.sc <- dir.filter(ae.sc)
ae.obs <- dir.filter(ae.obs)
ae.aca <- dir.filter(ae.aca)

# direct absolute case aversions and aggregated ie 
con.sum(ae.aca, 7)
(con.sum(ae.sc, 7) - con.sum(ae.obs, 7))*100/con.sum(ae.sc, 7)

# spillover absolute case aversions and aggregated ie 
con.sum(spill.aca, 7)
(con.sum(spill.sc, 7) - con.sum(spill.obs, 7))*100/con.sum(spill.sc, 7)

# calendar time aggregations
sct.aca <- spill.filter(sct.aca)
sct.obs <- spill.filter(sct.obs)
sct.sc <- spill.filter(sct.sc)

dct.sc <- dir.filter(dct.sc)
dct.obs <- dir.filter(dct.obs)
dct.aca <- dir.filter(dct.aca)

# direct absolute case aversions and aggregated ie 
con.sum(dct.aca, 8)
(con.sum(dct.sc, 8) - con.sum(dct.obs, 8))*100/con.sum(dct.sc, 8)

# spillover absolute case aversions and aggregated ie 
con.sum(sct.aca, 8)
(con.sum(sct.sc, 8) - con.sum(sct.obs, 8))*100/con.sum(sct.sc, 8)
