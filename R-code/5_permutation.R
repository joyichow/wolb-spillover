packages <- c("readxl", "geosphere", "MSCMT", "dplyr", "writexl", "tidyverse", 
              "ggplot2", "Synth", "glmnet", "osqp", "optimx", "meboot")
lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
})

base_path <- "/Users/your-local-file-path/"
fxn_path <- file.path(base_path, "R")
data_path <- file.path(base_path, "Data")

# Load datasets
intervention <- read_xlsx(file.path(data_path, "intervention.xlsx"))
control <- read_xlsx(file.path(data_path, "control.xlsx"))
dengue <- read_xlsx(file.path(data_path, "data submission incidence.xlsx"))
dates <- read_xlsx(file.path(data_path, "dates.xlsx"))
controls <- read_xlsx(file.path(data_path, "co.xlsx"))

dist <- prox.pairs(1500, intervention, control) 

# Incidence data 
dengue <- dengue[-c(18)]
dengue <- dengue[-c(1:52),]

# Defining spillover sites 
spillover <- as.vector(unlist(dist[,2]))

mapT <- data.frame(sno = integer(444))
mapT$sno <- 1:444
mapT$time <- sprintf("%dW%02d", dengue$Year, dengue$`Epid Week`)
mapT$map <- 1000:1443
mapT$ey.ew <- sprintf("%d.%02d", dengue$Year, dengue$`Epid Week`)

dates <- read_xlsx("/Users/joyichow/Desktop/spillover_scm/Data/dates.xlsx")
dates$time <- sapply(dates$Date, function(x) mapT$map[mapT$time==x]) #scm

s.date.map <- list(
  `1198` = c(1), 
  `1288` = c(2:9),       
  `1315` = c(10),      
  `1354` = c(11:12)       
)

ss.l <- colnames(dengue[-c(1:2)]) # all sites 

direct.ie <- data.frame(itr = numeric(), site = character(), dir = numeric(), rmse = numeric(), 
                        obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)

spill.ie <- list()

# Estimate synthetic control
for (k in 1:1000) {
  ss <- sample(ss.l, 0, replace = FALSE)                   # pick 12 spillover sites from all sites
  treat <- sample(setdiff(ss.l, ss), 20, replace = FALSE)   # pick 20 intervention sites from all sites
  co <- setdiff(ss.l, c(ss, treat))                         # donour pool 
  
  processed_dates <- c()

  spill.df <- data.frame(itr = numeric(), site = character(), dir = numeric(), rmse = numeric(),
                         obs = numeric(), sc = numeric(), stringsAsFactors = FALSE)

  for (i in 1:length(treat)){
    start_date <- as.integer(dates$time[i])
    times.dep  <- times.pred <- cbind("cases"  = c(1000, start_date))
    spill.size <- spillover.size(start_date)
    
    # defining spillover (curr.ss) and control groups (co2)
    curr.ss <- ss[1:spill.size]
    
    #calculating direct effects
    direct <- partial_scm(treat[i], dengue, co, start_date, population)
    
    direct.ie <- rbind(direct.ie, 
                       data.frame(itr = k, site = treat[i], dir = direct$dir, rmse = direct$rmse,
                                  obs = direct$obs, sc = direct$sc))
    
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
      spillover <- partial_scm(curr.ss[j], dengue, co, start_date, population)
      spill.df <- rbind(spill.df,
                        data.frame(itr = k, site = curr.ss[j], dir = spillover$dir, rmse = spillover$rmse,
                                   obs = spillover$obs, sc = spillover$sc))
    }
  }
  spill.ie[[k]] <- as.data.frame(do.call(cbind, spill.df))
}


spill.comb <- do.call(rbind, spill.ie)

spill.ie <- lapply(names(spill.ie), function(name) {
  spill.ie[[name]] %>%
    mutate(itr = name)
})

dir.ie <- subset(direct.ie, as.numeric(direct.ie$rmse) <100)

spill.ie <- bind_rows(spill.ie)

agg.dir.ies <- dir.ie %>%
  mutate(obs = as.numeric(as.character(obs)), sc = as.numeric(as.character(sc))) %>% 
  group_by(itr) %>%
  summarise(sum.obs = sum(obs, na.rm = TRUE),
            sum.sc = sum(sc, na.rm = TRUE),
            ie = (sum.sc - sum.obs)*100/sum.sc)

x_range <- c(-65, 65)  # Define the desired range for the x-axis
breaks <- seq(from = -65, to = 65, by = 5)  # Define breaks at intervals of 5

pdf("ptests.pdf")

hist(agg.dir.ies$ie, main = "Permutation test for direct effects", 
     xlab = "Aggregated Direct Protective Efficacy (%)", xaxs = "i",yaxs="i", 
     xlim = x_range, ylim = c(0, 250), xaxt = "n" )
abline(v = 64.35, lty = "dashed")
text(64.35, 230, "Actual agg. direct = 64.35", adj = c(1.05, -0.5))
axis(1, at = seq(-65, 65, by = 5), labels = TRUE, las = 1)

dir.count <- agg.dir.ies %>%
  filter(ie < 64.19) %>%
  summarise(count = n()) %>%
  pull(count)

(1000 - dir.count)/1000

# spill analysis 
spilled.ie <- subset(spill.comb, as.numeric(spill.comb$rmse) <100)

agg.spill.ies <- spilled.ie %>%
  mutate(obs = as.numeric(as.character(obs)), sc = as.numeric(as.character(sc))) %>% 
  group_by(itr) %>%
  summarise(sum.obs = sum(obs, na.rm = TRUE),
            sum.sc = sum(sc, na.rm = TRUE),
            ie = (sum.sc - sum.obs)*100/sum.sc)

x_range <- c(-75, 70)  # Define the desired range for the x-axis
breaks <- seq(from = -75, to = 70, by = 5)  # Define breaks at intervals of 5

# Create the histogram with specified breaks and range
hist(agg.spill.ies$ie, breaks = breaks, main = "Permutation test for spillover effects",
     xlab = "Aggregated Spillover Protective Efficacy (%)", xaxs = "i",yaxs="i", xlim = x_range, ylim = c(0, 110), xaxt = "n")
abline(v = 37.69, lty = "dashed")
text(37.69, 100, "Actual agg. spillover = 37.69", adj = c(1.05, -0.5))
axis(1, at = seq(-75, 70, by = 5), labels = TRUE, las = 1)


spill.count <- agg.spill.ies %>%
  filter(ie < 37.68808) %>%
  summarise(count = n()) %>%
  pull(count)

(1000 - spill.count)/1000


