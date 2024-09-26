# placebo 
# time placebo - meaning intervention time gets pushed earlier 
dates$time <- sapply(dates$Date, function(x) mapT$map[mapT$time==x]) #scm
dates$time <- dates$time - 104

s.placebo.map <- list(
  `1094` = c(1), 
  `1184` = c(2:9),       
  `1211` = c(10),      
  `1250` = c(11:12)       
)

direct.ie <- data.frame(site = character(), dir = numeric(), rmse = numeric(), p.error = numeric(), 
                        obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)
spill.ie <- list()
processed_dates <- c()
co2 <- setdiff(co, ss)

for (i in dates$Area){
  start_date <- as.integer(dates$time[dates$Area == colnames(dengue[i])])
  spill.size <- p.spillover.size(start_date)
  
  # defining spillover (curr.ss) and control groups (co2)
  curr.ss <- ss[1:spill.size]
  
  #calculating direct effects
  direct <- placebo.scm(i, dengue, co2, start_date)
  
  direct.ie <- rbind(direct.ie, 
                     data.frame(site = i, dir = direct$dir, rmse = direct$rmse, p.error = direct$p.error,
                                obs = direct$obs, sc = direct$sc))
  
  
  if (start_date %in% processed_dates || !start_date %in% c(1094,1184,1211,1250)) {
    next
  }
  
  spill.df <- data.frame(site = character(), dir = numeric(), rmse = numeric(),
                         obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)
  
  processed_dates <- c(processed_dates, start_date)
  
  # spillover effects
  if (!is.null(s.placebo.map[[as.character(start_date)]])) {
    ind <- s.placebo.map[[as.character(start_date)]]
  } else {
    ind <- 1:spill.size
  }
  
  for (j in ind) {
    spillover <- placebo.scm(curr.ss[j], dengue, co2, start_date)
    spill.df <- rbind(spill.df,
                      data.frame(site = curr.ss[j], dir = spillover$dir, rmse = spillover$rmse, 
                                 obs = spillover$obs, sc = spillover$sc))
  }
  
  spill.ie[[i]] <- as.data.frame(do.call(cbind, spill.df))
}


spill.ie <- combined.sites(spill.ie)

dirdf <- subset(direct.ie, rmse <100)
spilldf <- subset(spill.ie, as.numeric(rmse) <100)

# Aggregate effects
dir.ie <- calc.ie(dirdf)                  # direct IE
s.ie <- calc.ie(spilldf)  

write_xlsx(direct.ie, "directed_tp.xlsx")
write_xlsx(spill.ie, "spillover_tp.xlsx")

# space placebo 
p.co <- subset(control, !control$Site %in% ss) # all sites that aren't treated or spillover

threshold_distance <- 1500

p.dist <- data.frame(intervention_site = integer(), 
                     control_site = integer(), 
                     stringsAsFactors = FALSE)

for (i in 1:(nrow(p.co) - 1)) {
  for (j in (i + 1):nrow(p.co)) {
    distance <- distHaversine(c(p.co$Longitude[i], p.co$Latitude[i]), 
                              c(p.co$Longitude[j], p.co$Latitude[j]))
    
    if (distance <= threshold_distance) {
      p.dist <- rbind(p.dist, 
                      data.frame(site1 = p.co$Site[i], 
                                 site2 = p.co$Site[j],
                                 distance = distance))
    }
  }
}

colnames(p.dist) <- colnames(dist[1:2]) 

direct.ie <- data.frame(site = character(), dir = numeric(), rmse = numeric(),
                        obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)

spill.df <- data.frame(site = character(), dir = numeric(), rmse = numeric(), p.error = numeric(), 
                       obs = numeric(), sc = numeric(), spillover.from = character(), stringsAsFactors = FALSE)

for (i in p.co$Site){
  start_date <- 1175-52
  times.dep  <- times.pred <- cbind("cases"  = c(1000, start_date))
  
  # defining spillover (curr.ss) and control groups (co2)
  curr.ss <- subset(p.dist, p.dist$Study_PA == i)$Within_1500
  co2 <- setdiff(c(p.co$Site), c(curr.ss, i))
  
  #calculating direct effects for placebo
  direct <- partial_scm(i, dengue, co2, start_date, population)
  
  direct.ie <- rbind(direct.ie, 
                     data.frame(site = i, dir = direct$dir, rmse = direct$rmse,
                                obs = direct$obs, sc = direct$sc))
  
  # spillover effects
  for (j in curr.ss) {
    spillover <- partial_scm(j, dengue, co2, start_date, population)
    spill.df <- rbind(spill.df,
                      data.frame(site = j, dir = spillover$dir, rmse = spillover$rmse, 
                                 obs = spillover$obs, sc = spillover$sc, spillover.from = i))
  }
}

# Eliminating fits worse than RMSE 100 and those that failed placebo
dirdf <- subset(direct.ie, rmse <100)
spilldf <- subset(spill.df, rmse <100)

# Aggregate effects
dir.ie <- calc.ie(dirdf)                  # direct IE
s.ie <- calc.ie(spilldf)  

write_xlsx(direct.ie, "directed_sp.xlsx")
write_xlsx(spill.df, "spillover_sp.xlsx")
write_xlsx(p.dist, "distance_sp.xlsx")
