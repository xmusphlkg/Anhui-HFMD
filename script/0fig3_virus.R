rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggspatial)
  library(readr)
  library(scales)
})

if (dir.exists("../outcome")) {
  setwd("../outcome")
} else if (dir.exists("outcome")) {
  setwd("outcome")
} else {
  stop("Cannot find outcome/ directory from current working directory.")
}

data_dir <- if (dir.exists("../data_public")) "../data_public" else "../data"

load(file.path(data_dir, "datafile_lo.RData"))

county_map <- st_read(file.path(data_dir, "安徽省_县.shp"), quiet = TRUE) %>%
  st_make_valid()

city_map <- st_read(file.path(data_dir, "安徽省_市.shp"), quiet = TRUE) %>%
  st_make_valid()

dir.create("supplementary_tables", showWarnings = FALSE, recursive = TRUE)

normalize_virus <- function(x) {
  dplyr::case_when(
    x %in% c("CV-A16", "Cox A16") ~ "CV-A16",
    x == "EV71" ~ "EV71",
    is.na(x) ~ NA_character_,
    TRUE ~ "Others"
  )
}

typed <- datafiles_lo %>%
  mutate(
    date = as.Date(date),
    year = lubridate::year(date),
    virus_std = normalize_virus(virus)
  ) %>%
  filter(year >= 2008, year <= 2023, !is.na(virus_std))

county_virus <- typed %>%
  count(name, year, virus_std, name = "typed_cases")

county_rank <- county_virus %>%
  group_by(name, year) %>%
  arrange(desc(typed_cases), virus_std, .by_group = TRUE) %>%
  mutate(rank_id = row_number())

county_conf <- county_rank %>%
  group_by(name, year) %>%
  summarise(
    typed_total = sum(typed_cases),
    top_virus = if_else(sum(typed_cases == max(typed_cases)) > 1, "Mixed", first(virus_std[typed_cases == max(typed_cases)])),
    top_cases = max(typed_cases),
    second_cases = if_else(any(rank_id == 2), typed_cases[rank_id == 2][1], 0L),
    top_share = top_cases / typed_total,
    dominance_gap = top_cases - second_cases,
    .groups = "drop"
  )

all_counties <- county_map %>%
  st_drop_geometry() %>%
  distinct(name)

county_conf <- all_counties %>%
  crossing(year = 2008:2023) %>%
  left_join(county_conf, by = c("name", "year")) %>%
  mutate(
    typed_total = coalesce(typed_total, 0L),
    top_virus = coalesce(top_virus, "No typed cases"),
    top_cases = coalesce(top_cases, 0L),
    second_cases = coalesce(second_cases, 0L),
    top_share = if_else(typed_total > 0, top_share, NA_real_),
    dominance_gap = coalesce(dominance_gap, 0L)
  )

# Preserve continuity across district boundary changes in the current shapefile.
county_conf <- county_conf %>%
  filter(!(name == "叶集区" & year %in% 2010:2015)) %>%
  filter(!(name == "博望区" & year %in% 2010:2011))

yeji_proxy <- county_conf %>%
  filter(name == "霍邱县", year %in% 2010:2015) %>%
  mutate(name = "叶集区")

bowang_proxy <- county_conf %>%
  filter(name == "当涂县", year %in% 2010:2011) %>%
  mutate(name = "博望区")

county_conf <- bind_rows(county_conf, yeji_proxy, bowang_proxy) %>%
  distinct(name, year, .keep_all = TRUE) %>%
  mutate(
    display_group = case_when(
      typed_total < 10 ~ "Typed n < 10",
      top_virus %in% c("EV71", "CV-A16", "Others", "Mixed") ~ top_virus,
      TRUE ~ "Typed n < 10"
    ),
    display_group = factor(
      display_group,
      levels = c("EV71", "CV-A16", "Others", "Mixed", "Typed n < 10")
    )
  ) %>%
  arrange(name, year)

write_csv(county_conf, "supplementary_tables/supplementary_data_D7_county_map_masking.csv")

color_map <- c(
  "EV71" = "#3CC8C0FF",
  "CV-A16" = "#F2EBBBFF",
  "Others" = "#04578CFF",
  "Mixed" = "#CB1724FF",
  "Typed n < 10" = "#B8B8B8"
)

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

plot_map_panel <- function(year_value, panel_label) {
  map_data <- county_map %>%
    left_join(
      county_conf %>% filter(year == year_value),
      by = "name"
    )

  ggplot(map_data) +
    geom_sf(aes(fill = display_group), color = "black", linewidth = 0.15) +
    geom_sf_text(
      data = city_labels,
      aes(label = label),
      size = 2.2,
      family = "Times New Roman",
      color = "black"
    ) +
    scale_fill_manual(values = color_map, drop = FALSE, na.value = "#B8B8B8", name = NULL) +
    labs(
      title = paste0(panel_label, ": ", year_value)
    ) +
    annotation_scale(location = "bl", width_hint = 0.2) +
    annotation_north_arrow(
      location = "tr",
      which_north = "true",
      style = north_arrow_fancy_orienteering(),
      height = unit(0.8, "cm")
    ) +
    theme_classic(base_family = "Times New Roman") +
    theme(
      plot.title = element_text(size = 12.5, face = "bold"),
      text = element_text(color = "black"),
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank()
    )
}

stacked_summary <- county_conf %>%
  count(year, display_group, name = "county_n") %>%
  complete(year = 2008:2023, display_group, fill = list(county_n = 0L))

fig_bar <- ggplot(stacked_summary, aes(x = year, y = county_n, fill = display_group)) +
  geom_col(position = "fill", color = "white", linewidth = 0.15) +
  scale_fill_manual(values = color_map, drop = FALSE) +
  scale_x_continuous(breaks = 2008:2023, expand = expansion(add = c(0.3, 0.3))) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = NULL,
    x = "Year",
    y = "Share of county-level units",
    fill = NULL
  ) +
  theme_bw(base_family = "Times New Roman") +
  theme(
    panel.grid = element_blank(),
    text = element_text(color = "black"),
    axis.text.x = element_text(size = 10, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.55, "cm"),
    plot.title.position = "plot"
  )

map_panels <- Map(
  f = plot_map_panel,
  year_value = 2008:2023,
  panel_label = LETTERS[2:17]
)

fig_maps <- wrap_plots(map_panels, ncol = 4)

fig <- fig_bar

ggsave(
  "fig3.pdf",
  plot = fig,
  width = 10,
  height = 5.6,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  "fig3.png",
  plot = fig,
  width = 10,
  height = 5.6,
  dpi = 320
)

ggsave(
  "figS8_county_typed_dominance_maps.pdf",
  plot = fig_maps,
  width = 14,
  height = 12,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  "figS8_county_typed_dominance_maps.png",
  plot = fig_maps,
  width = 14,
  height = 12,
  dpi = 320
)
