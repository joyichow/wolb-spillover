# original SCM
ss <- c()
co2 <- setdiff(co, ss)
dates$time <- sapply(dates$Date, function(x) mapT$map[mapT$time == x]) #scm

# time placebo - meaning intervention time gets pushed earlier 
dates$time <- dates$time - 104

# site-aggregated effects 
direct.ie <- data.frame()

for (i in dates$Area){
  start_date <- as.integer(dates$time[dates$Area == i])
  spill.size <- spillover.size(start_date)
  curr.ss <- ss[1:spill.size]
  
  #calculating direct effects in main analysis 
  direct <- partial_scm(i, dengue, co2, start_date, population)

  direct.ie <- rbind(direct.ie,
                     data.frame(site = i, dir = direct$dir, rmse = direct$rmse,
                                obs = direct$obs, sc = direct$sc, aca = direct$aca))
  
  # time placebo
  direct <- placebo.scm(i, dengue, co2, start_date)

  direct.ie <- rbind(direct.ie,
                     data.frame(site = i, dir = direct$dir, rmse = direct$rmse, p.error = direct$p.error,
                                obs = direct$obs, sc = direct$sc))
  
}

# Eliminating fits worse than RMSE 100 and those that failed placebo
dirdf <- dir.filter(direct.ie)

# Aggregate effects
dir.ie <- calc.ie(dirdf)                  # direct IE
dir.aca <- sum(dirdf$aca)                 # direct ACA

write_xlsx(direct.ie, "directed.xlsx")
write_xlsx(direct.ie, "directed_tp.xlsx")

# space placebo 
p.co <- subset(control, !control$Site %in% ss) # all sites that aren't treated or spillover

direct.ie <- data.frame(site = character(), dir = numeric(), rmse = numeric(),
                        obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)

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
}

# Eliminating fits worse than RMSE 100 and those that failed placebo
dirdf <- dir.filter(direct.ie)

# Aggregate effects
dir.ie <- calc.ie(dirdf)                  # direct IE
dir.aca <- sum(dirdf$aca)                 # direct ACA

write_xlsx(direct.ie, "directed_sp.xlsx")


