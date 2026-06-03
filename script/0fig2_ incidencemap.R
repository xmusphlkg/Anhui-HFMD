rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggspatial)
  library(readxl)
  library(patchwork)
  library(lubridate)
})

if (dir.exists("../outcome")) {
  setwd("../outcome")
} else if (dir.exists("outcome")) {
  setwd("outcome")
} else {
  stop("Cannot find outcome/ directory from current working directory.")
}

data_dir <- if (dir.exists("../data_public")) "../data_public" else "../data"

load(file.path(data_dir, "datafile_city.RData"))

city_map <- st_read(file.path(data_dir, "安徽省_市.shp"), quiet = TRUE) %>%
  st_make_valid()

pop <- read_excel(file.path(data_dir, "city_pop.xlsx"))

city_counts <- datafiles_city %>%
  mutate(year = year(as.Date(date))) %>%
  filter(year >= 2008, year <= 2023) %>%
  count(name, year, name = "case_count")

pop_long <- pop %>%
  pivot_longer(cols = -name, names_to = "year", values_to = "population") %>%
  mutate(year = as.integer(year)) %>%
  filter(year >= 2008, year <= 2023)

incidence_df <- city_counts %>%
  left_join(pop_long, by = c("name", "year")) %>%
  mutate(incidence_rate = case_count / population * 10)

map_data <- city_map %>%
  left_join(incidence_df, by = "name")

breaks <- pretty(map_data$incidence_rate, n = 5)

city_labels <- city_map %>%
  filter(name %in% c("蚌埠市", "芜湖市", "黄山市")) %>%
  mutate(
    label = case_when(
      name == "蚌埠市" ~ "Bengbu",
      name == "芜湖市" ~ "Wuhu",
      name == "黄山市" ~ "Huangshan",
      TRUE ~ name
    )
  ) %>%
  st_point_on_surface()

plot_panel <- function(year_value, panel_label) {
  ggplot(map_data %>% filter(year == year_value)) +
    geom_sf(aes(fill = incidence_rate), color = "black", linewidth = 0.2) +
    geom_sf_text(
      data = city_labels,
      aes(label = label),
      size = 2.4,
      family = "Times New Roman",
      color = "black"
    ) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#C6DBEF", "#6BAED6", "#2171B5", "#08306B"),
      limits = range(breaks),
      breaks = breaks
    ) +
    labs(title = paste0(panel_label, ": ", year_value)) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering(),
      height = unit(0.8, "cm")
    ) +
    theme_classic(base_family = "Times New Roman") +
    theme(
      text = element_text(color = "black"),
      plot.title = element_text(size = 12.5, face = "bold"),
      legend.position = "bottom",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    ) +
    guides(
      fill = guide_colourbar(
        title = "Incidence rate (per 100,000)",
        title.position = "top",
        title.hjust = 0,
        barwidth = unit(5, "cm")
      )
    )
}

panels <- Map(
  f = plot_panel,
  year_value = 2008:2023,
  panel_label = LETTERS[1:16]
)

fig <- wrap_plots(panels, ncol = 4, guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "fig2.pdf",
  plot = fig,
  width = 14,
  height = 12,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  "fig2.png",
  plot = fig,
  width = 14,
  height = 12,
  dpi = 320
)
