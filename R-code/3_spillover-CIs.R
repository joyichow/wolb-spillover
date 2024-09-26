# Method 2 for spillover with CIs
library(meboot)

# Bootstrapping with meboot
s <- 1  # Seed number 
reps_chunk <- 2  # Number of reps to run in each chunk
total_reps <- 10000 # Total number of bootstrap samples
iterations <- total_reps / reps_chunk  

bsdf <- list()  

for (i in 1:iterations) {
  set.seed(s+1234)
  store <- meboot(as.ts(dengue[-c(1:2)]), trim = list(xmin = 0), reps = reps_chunk, elaps = TRUE)
  
  for (j in 1:reps_chunk) {
    time_series <- store$ensemble[, j]  
    data_frame <- cbind(dengue[c(1:2)], matrix(time_series, nrow = nrow(dengue), byrow = FALSE))
    
    colnames(data_frame) <- colnames(dengue)
    
    bsdf[[length(bsdf) + 1]] <- data.frame(data_frame)
  }
  s <- s + 1
}

# Define mapping for synthetic control estimation
mapT <- data.frame(sno = 1:444)
mapT$time <- sprintf("%dW%02d", dengue$Year, dengue$`Epid Week`)
mapT$ey.ew <- sprintf("%d.%02d", dengue$Year, dengue$`Epid Week`)
mapT$map <- seq(1000,1443)

# Match dates to serial numbers
dates$time <- sapply(dates$Date, function(x) mapT$map[mapT$time == x]) #scm

# Estimate synthetic control
ss <- c("C43", "C35", "C61", "C39", "C48", "C62", "C25", "C53", "C70", "C03", "C62", "C73") # sites with spillover ordered to start time
ss <- c() # without spillover
co2 <- setdiff(co, ss)


s.date.map <- list(
  `1198` = c(1), 
  `1288` = c(2:9),       
  `1315` = c(10),      
  `1354` = c(11:12)       
)

#reset from here each run
all.dir <- data.frame(matrix(ncol = length(unique(c(dates$Area, ss))), nrow = iterations))
all.rmse <- data.frame(matrix(ncol = length(unique(c(dates$Area, ss))), nrow = iterations))
names(all.dir) <- names(all.rmse) <- unique(c(dates$Area, ss))

det.ie <- det.aca <- set.ie <- set.aca <- matrix(ncol = 5, nrow = 0) #et 
dcal.ie <- dcal.aca <- scal.ie <- scal.aca <- matrix(ncol = 6, nrow = 0) #ct

for (iter in 1:length(bsdf)) {
  deng <- bsdf[[iter]]
  processed_dates <- c()
  
  direct.ie <- data.frame(site = character(), dir = numeric(), rmse = numeric(),
                          obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)
  
  spill.ie <- list()
  
  ## Event time ##
  # direct
  ae.sc <- data.frame()
  ae.obs <- data.frame()
  ae.aca <- data.frame()
  
  # spillover
  spill.aca <- list()
  spill.sc <- list()
  spill.obs <- list()
  
  ## Calendar time ##
  # direct
  dct.sc <- data.frame()
  dct.obs <- data.frame()
  dct.aca <- data.frame()
  
  # spillover 
  sct.aca <- list()
  sct.sc <- list()
  sct.obs <- list()
  
  for (i in dates$Area){
    start_date <- as.integer(dates$time[dates$Area == i])
    spill.size <- spillover.size(start_date)
    curr.ss <- ss[1:spill.size]
    curr.ss <- 0

    #calculating direct effects
    direct <- partial_scm(i, deng, co2, start_date, population)
    
    direct.ie <- rbind(direct.ie, 
                       data.frame(site = i, dir = direct$dir, rmse = direct$rmse,
                                  obs = direct$obs, sc = direct$sc, aca = direct$aca))
    
    ae.sc <- rbind(ae.sc, data.frame(site = i, rmse = direct$rmse, t(direct$et$sum_sc)))
    ae.obs <- rbind(ae.obs, data.frame(site = i, rmse = direct$rmse, t(direct$et$sum_obs)))
    ae.aca <- rbind(ae.aca, data.frame(site = i, rmse = direct$rmse, t(direct$aca.et$metric)))
    
    dct.sc <- rbind(dct.sc, data.frame(site = i, rmse = direct$rmse, t(direct$ct$sum_sc)))
    dct.obs <- rbind(dct.obs, data.frame(site = i, rmse = direct$rmse, t(direct$ct$sum_obs)))
    dct.aca <- rbind(dct.aca, data.frame(site = i, rmse = direct$rmse, t(direct$aca.ct$metric)))
    
    spill.df <- data.frame(site = character(), dir = numeric(), rmse = numeric(),
                           obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)

    sp.sc.df <- data.frame()
    sp.obs.df <- data.frame()
    sp.aca.df <- data.frame()

    sct.sc.df <- data.frame()
    sct.obs.df <- data.frame()
    sct.aca.df <- data.frame()

    if (start_date %in% processed_dates || !start_date %in% c(1198, 1288, 1315, 1354)) {
      next
    }

    processed_dates <- c(processed_dates, start_date)

    # spillover effects
    if (!is.null(s.date.map[[as.character(start_date)]])) {
      ind <- s.date.map[[as.character(start_date)]]
    } else {
      ind <- 1:spill.size
    }

    for (j in ind) {
      spillover <- partial_scm(curr.ss[j], deng, co2, start_date, population)
      spill.df <- rbind(spill.df,
                        data.frame(site = curr.ss[j], dir = spillover$dir, rmse = spillover$rmse,
                                   obs = spillover$obs, sc = spillover$sc, aca = spillover$aca))

      sp.sc.df <- rbind(sp.sc.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$et$sum_sc)))
      sp.obs.df <- rbind(sp.obs.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$et$sum_obs)))
      sp.aca.df <- rbind(sp.aca.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$aca.et$metric)))

      sct.sc.df <- rbind(sct.sc.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$ct$sum_sc)))
      sct.obs.df <- rbind(sct.obs.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$ct$sum_obs)))
      sct.aca.df <- rbind(sct.aca.df, data.frame(site = curr.ss[j], rmse = spillover$rmse, t(spillover$aca.ct$metric)))
    }

    spill.ie[[i]] <- as.data.frame(do.call(cbind, spill.df))

    spill.aca[[i]] <- as.data.frame(do.call(cbind, sp.aca.df))
    spill.obs[[i]] <- as.data.frame(do.call(cbind, sp.obs.df))
    spill.sc[[i]] <- as.data.frame(do.call(cbind, sp.sc.df))

    sct.aca[[i]] <- as.data.frame(do.call(cbind, sct.aca.df))
    sct.obs[[i]] <- as.data.frame(do.call(cbind, sct.obs.df))
    sct.sc[[i]] <- as.data.frame(do.call(cbind, sct.sc.df))
  }

  spill.ie <- combined.sites(spill.ie)
  spill.aca <- combined.sites(spill.aca)
  spill.sc <- combined.sites(spill.sc)
  spill.obs <- combined.sites(spill.obs)

  sct.aca <- combined.sites(sct.aca)
  sct.sc <- combined.sites(sct.sc)
  sct.obs <- combined.sites(sct.obs)
  
  ie effects - change to $aca for absolute case aversions
  direct.dirs <- setNames(direct.ie$dir, direct.ie$site)
  spill.dirs <- setNames(as.numeric(spill.ie$aca), spill.ie$site)

  c62_mean <- mean(spill.dirs["C62"], na.rm = TRUE)
  spill.dirs <- c(spill.dirs[names(spill.dirs) != "C62"], C62 = c62_mean)
  
  all.dir[iter, names(direct.dirs)] <- direct.dirs
  all.dir[iter, names(spill.dirs)] <- spill.dirs
  
  # rmse 
  direct.rmse <- setNames(direct.ie$rmse, direct.ie$site)
  spill.rmse <- setNames(as.numeric(spill.ie$rmse), spill.ie$site)

  c62 <- mean(spill.rmse["C62"], na.rm = TRUE)
  spill.rmse <- c(spill.rmse[names(spill.rmse) != "C62"], C62 = c62)

  all.rmse[iter, names(direct.rmse)] <- direct.rmse
  all.rmse[iter, names(spill.rmse)] <- spill.rmse

  # aggregate effects
  direct.ie <- dir.filter(direct.ie)
  spill.ie <- spill.filter(spill.ie)

  agg.ie.dir <- calc.ie(direct.ie)
  agg.ie.sc <- calc.ie(spill.ie)

  agg.dir.aca <- sum(direct.ie$aca)
  agg.spill.aca <- sum(as.numeric(spill.ie$aca))
  agg.total.aca <- agg.dir.aca+agg.spill.aca

  all.dir[iter, "AggIEDir"] <- agg.ie.dir
  all.dir[iter, "AggIESc"] <- agg.ie.sc
  all.dir[iter, "AggACADir"] <- agg.dir.aca
  all.dir[iter, "AggACASc"] <- agg.spill.aca
  all.dir[iter, "AggACATot"] <- agg.total.aca

  # aggregate by time
  dT <- map(list(dct.sc, dct.obs, dct.aca, ae.sc, ae.obs, ae.aca), ~filter(.x, !site %in% c("W05", "W27", "W04", "W10", "W15", "W26")))
  sT <- map(list(sct.sc, sct.obs, sct.aca, spill.sc, spill.obs, spill.aca), ~spill.filter(.x))

  # event time
  etd.aca <- con.sum(dT[[6]], 7)                                                        # direct ACA - et
  etd.ie <- (con.sum(dT[[4]], 7) - con.sum(dT[[5]], 7))*100/con.sum(dT[[4]], 7)        # direct IE - et

  ets.aca <- con.sum(sT[[6]], 7)                                                  # spillover ACA - et
  ets.ie <- (con.sum(sT[[4]], 7) - con.sum(sT[[5]], 7))*100/con.sum(sT[[4]], 7)   # spillover IE - et

  det.ie <- (rbind(det.ie, t(data.frame(etd.ie))))
  det.aca <- (rbind(det.aca, t(data.frame(etd.aca))))
  set.ie <- (rbind(set.ie, t(data.frame(ets.ie))))
  set.aca <- (rbind(set.aca, t(data.frame(ets.aca))))

  # calendar time
  ctd.aca <- con.sum(dT[[3]], 8)                                                    # direct ACA - ct
  ctd.ie <- (con.sum(dT[[1]], 8) - con.sum(dT[[2]], 8))*100/con.sum(dT[[1]], 8)     # direct IE - ct

  cts.aca <- con.sum(sT[[3]], 8)                                                       # spillover ACA - ct
  cts.ie <- (con.sum(sT[[1]], 8) - con.sum(sT[[2]], 8))*100/con.sum(sT[[1]], 8)       # spillover IE - ct

  dcal.ie <- (rbind(dcal.ie, t(data.frame(ctd.ie))))
  dcal.aca <- (rbind(dcal.aca, t(data.frame(ctd.aca))))
  scal.ie <- (rbind(scal.ie, t(data.frame(cts.ie))))
  scal.aca <- (rbind(scal.aca, t(data.frame(cts.aca))))

  save(direct.ie, spill.ie, file = paste0("iteration_", iter, "_results.RData"))
}

list.dfs <- list(all.dir, det.ie, det.aca, set.ie, set.aca, dcal.ie, dcal.aca, scal.ie[,3:6], scal.aca)

formatted.quantiles <- lapply(list.dfs, calc.quantiles)

# easy save
library(openxlsx)
list_of_columns <- list()

max_length <- max(sapply(formatted.quantiles, length))

for (i in seq_along(formatted.quantiles)) {
  col_data <- unlist(formatted.quantiles[[i]])
  
  col_data <- as.character(col_data)
  
  if (length(col_data) < max_length) {
    col_data <- c(col_data, rep(NA, max_length - length(col_data)))
  }

  col_name <- if (!is.null(names(formatted.quantiles))) names(formatted.quantiles)[i] else paste("Column", i)
  list_of_columns[[col_name]] <- col_data
}

combined_df <- data.frame(list_of_columns)
colnames(combined_df) <- c("all.dir", "det.ie", "det.aca", "set.ie", "set.aca", "dcal.ie", "dcal.aca", "scal.ie", "scal.aca")

write.xlsx(combined_df, "formatted_quantiles_10000_orig_ie.xlsx", rowNames = FALSE)
write.xlsx(all.dir, "test.xlsx")

