# plotting direct and spillover effects 
library(readxl)
library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(prettymapr)
library(ggrepel)
library(stringr)

base_path <- "/Users/your-local-file-path/"
fxn_path <- file.path(base_path, "R")
data_path <- file.path(base_path, "Data")

res <- read_xlsx(file.path(data_path, "ie-effects.xlsx"))
res$type <- as.factor(res$type)

res.map <- st_as_sf(res, coords = c("long", "lat"), crs = 4326, agr = "constant")

res.map <- res.map %>%
  mutate(type_ie_sign = interaction(type, ie >= 0, drop = TRUE),
         type_ie_sign = factor(type_ie_sign, labels = c("Spill-", "Direct-", "Spill+", "Direct+"))) %>%
  filter(!(unit %in% c("W04", "W05", "W10", "W26", "W27", 
                       "C62", "C25", "C43", "C61", "C48", "C73")))

res.map <- res.map %>%
  mutate(x = st_coordinates(geometry)[, 1],
         y = st_coordinates(geometry)[, 2])

units <- res.map %>% filter(type == 0)
bases <- res.map %>% filter(type == 1)

valid_units <- units %>% 
  filter(!is.na(base) & base != "") %>% 
  rowwise() %>%
  mutate(valid_base = any(str_split(base, ",\\s*")[[1]] %in% bases$unit)) %>%
  filter(valid_base) %>%
  ungroup()

lines <- list()

for (i in 1:nrow(valid_units)) {
  unit_row <- valid_units[i, ]
  base_ids <- str_split(unit_row$base, ",\\s*")[[1]]
  
  for (base_id in base_ids) {
    if (base_id %in% bases$unit) {
      base_row <- bases %>% filter(unit == base_id)
      # Assuming at least one valid base is found, generate the line
      if (nrow(base_row) > 0) {
        line_geom <- st_sfc(st_linestring(rbind(st_coordinates(unit_row), st_coordinates(base_row))), crs = 4326)
        lines <- c(lines, line_geom)
      }
    }
  }
}

lines_sfc <- st_sfc(lines, crs = 4326)

lines_sf <- st_sf(geometry = lines_sfc)

bases$ie <- as.numeric(bases$ie)
units$ie <- as.numeric(units$ie)

# version 2 
# with map
ggplot() +
  annotation_map_tile(type = "osm", zoom = 15) +  # Adding the map
  geom_sf(data = bases, aes(size = abs(ie)), color = "black", alpha = 0.8, show.legend = FALSE) +
  geom_sf(data = units, aes(size = abs(ie)), color = "black", alpha = 0.8, show.legend = FALSE) +

  geom_sf(data = bases, aes(size = abs(ie) * 0.8, color = type_ie_sign), alpha = 0.8) +
  geom_sf(data = units, aes(size = abs(ie) * 0.8, color = type_ie_sign), alpha = 0.8) +
  geom_sf(data = lines_sf, color = "gray50", size = 0.5, alpha = 1) +
  annotation_scale()

  geom_text_repel(data = res.map, aes(x = x, y = y, label = unit),
                  size = 2.5, box.padding = 0.2, point.padding = 0.2)+
  scale_size_continuous(name = "Abs IE", range = c(2, 9), limits = c(0, 100)) +
  scale_color_manual(values = c("Direct+" = "blue", "Direct-" = "purple", 
                                "Spill+" = "red", "Spill-" = "orange"),
                     name = "Effect Type",
                     labels = c("Direct-", "Direct+", "Spill-", "Spill+")) +
  labs(size = "IE Value") +
  ggtitle("Point Estimates, Effects, and Connections on Map") +
  theme_minimal() +
  coord_sf() +
  theme(legend.position = "bottom")

###### IE plots by site ######
library(forcats)
library(ggrepel)

df <- read_xlsx(file.path(data_path, "results_ie.xlsx"))

df <- df %>%
  mutate(
    IE_clean = gsub("[()]", "", IE),
    IE_numeric = as.numeric(sub(" .*", "", IE_clean)),
    CI = sub(".*\\((.*)\\).*", "\\1", IE),
    lower_IE = as.numeric(sub(" -.*", "", CI)),
    upper_IE = as.numeric(sub(".*- ", "", CI))
  ) %>%
  filter(!(Site %in% c("W04", "W05", "W10", "W15", "W26", "W27", 
                       "C62", "C25", "C43", "C61", "C48", "C73"))) %>%
  mutate(Type = factor(Type, levels = c("Direct", "Spill", "Agg. Direct", "Agg. Spillover"))) %>%
  arrange(desc(IE_numeric)) %>%
  mutate(Site = fct_inorder(Site)) %>%
  mutate(Site = fct_relevel(Site, "Agg. Spillover", "Agg. Direct", after = Inf))  # Move SE and DE to the last positions


p <- ggplot(df, aes(y = Site, x = IE_numeric, color = Type)) +
  geom_point(shape = 16, size = 2) +
  geom_errorbarh(aes(xmin = lower_IE, xmax = upper_IE), height = 0.2, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.2f", IE_numeric)), nudge_x = 0.5, hjust = 0.6, size = 3.5, vjust = 1.9, show.legend = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray", show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) +
  scale_color_manual(values = c("Direct" = "blue", "Spill" = "darkorange", 
                                "Agg. Direct" = "blue4", "Agg. Spillover" = "darkorange4")) +
  labs(y = "Site", x = "Protective Efficacy (%)", title = "Protective Efficacy (PE) by Site", color = "Effect type") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(angle = 0, hjust = 1, size = 12),
        axis.text.x = element_text(size = 12),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        axis.line = element_line(color = "black"),
        legend.position = c(0.17, 0.08),
        plot.margin = unit(c(0.6, 0.5, 0.5, 0.3), "cm"),
        legend.background = element_rect(colour = "black", size = 0.5, linetype = "solid")) +
  guides(color = guide_legend(keywidth = 0.5, keyheight = 0.8, title.position = "top", title.hjust = 0))


p

