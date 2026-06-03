rm(list=ls())
#load package
library(spdep)
library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(readxl)
library(tidyr)
library(gridExtra)
library(patchwork)
library(paletteer)
library(lubridate)

#setwd
setwd('../outcome')
data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'
load(file.path(data_dir, 'datafile_city.RData'))

#导入各区县历年人口数据(/10000)
pop <- read_excel(file.path(data_dir, 'city_pop.xlsx'))

# 导入shapefile地图数据
anhui_map <- st_read(file.path(data_dir, '安徽省_市.shp'))
anhui_map <- st_make_valid(anhui_map)# 修复地图错误

# 计算每个地区每年病例数
datafiles_loo <-  datafiles_city%>%
  mutate(year = as.character(year(date))) %>%
  group_by(name, year) %>%
  summarise(case_count = n(), .groups = "drop")%>%
  filter(year>=2008&year<=2023)

# 将人口数据进行长格式转换，便于合并
pop_long <- pop %>%
  pivot_longer(cols = -name, names_to = "year", values_to = "population") %>%
  mutate(year = as.character(year))%>%
  filter(year>=2008&year<=2023)

# 合并病例数和人口数据
merged_data0 <- left_join(datafiles_loo, pop_long, by = c("name", "year"))

# 计算历年发病率（每10万人口的发病数）
merged_data <- merged_data0 %>%
  mutate(incidence_rate = (case_count / population) * 10)

# 1. 计算各城市2008-2023年平均发病率
avg_2008_2023 <- merged_data %>%
     group_by(name) %>%
     summarise(incidence_rate = mean(incidence_rate, na.rm = TRUE)) %>%
     mutate(period = "2008-2023")

# 2. 计算各城市2008-2016年平均发病率
avg_2008_2016 <- merged_data %>%
     filter(year %in% 2008:2016) %>%
     group_by(name) %>%
     summarise(incidence_rate = mean(incidence_rate, na.rm = TRUE)) %>%
     mutate(period = "2008-2016")

# 3. 计算各城市2017-2019年平均发病率
avg_2017_2019 <- merged_data %>%
     filter(year %in% 2017:2019) %>%
     group_by(name) %>%
     summarise(incidence_rate = mean(incidence_rate, na.rm = TRUE)) %>%
     mutate(period = "2017-2019")

# 4. 计算各城市2020-2023年平均发病率
avg_2020_2023 <- merged_data %>%
     filter(year %in% 2020:2023) %>%
     group_by(name) %>%
     summarise(incidence_rate = mean(incidence_rate, na.rm = TRUE)) %>%
     mutate(period = "2020-2023")

# 合并所有时期数据
all_periods <- bind_rows(
     avg_2008_2023,
     avg_2008_2016,
     avg_2017_2019,
     avg_2020_2023
) %>%
     mutate(period = factor(period, 
                            levels = c("2008-2023", "2008-2016", "2017-2019", "2020-2023"),
                            ordered = TRUE))

# 合并地图数据
map_data <- left_join(anhui_map, all_periods, by = "name")

# 创建统一的色标范围
breaks <- pretty(map_data$incidence_rate, n = 10)
color_palette <- paletteer::paletteer_d("dichromat::BluetoDarkOrange_12")

# 绘制四个时期的发病率地图
plot_period <- function(period_name) {
     map_data %>% 
          filter(period == period_name) %>% 
          ggplot() +
          geom_sf(aes(fill = incidence_rate), color = "black") +
          scale_fill_gradientn(
               colours = color_palette,
               limits = range(breaks),
               breaks = breaks,
               na.value = "grey90"  # 处理NA值
          ) +
          labs(title = period_name) +  # 使用时期作为标题
          theme_classic(base_family = "Times New Roman") +
          theme(
               text = element_text(color = "black"),
               plot.title = element_text(size = 14, hjust = 0, face = "bold"),
               legend.title = element_text(size = 12, face = "bold"),
               legend.text = element_text(size = 10),
               legend.position = "bottom",
               axis.text = element_blank(),
               axis.ticks = element_blank(),
               axis.line = element_blank()
          ) +
          guides(fill = guide_colourbar(
               title = "Incidence rate (per 100,000)",
               title.position = "top",
               barwidth = 15,
               label.position = "bottom"
          )) +
          annotation_scale(location = "bl", width_hint = 0.2) +
          annotation_north_arrow(location = "tr", 
                                 which_north = "true", 
                                 style = north_arrow_fancy_orienteering(),
                                 height = unit(1, "cm"))
}

# 生成四个时期的地图
p1 <- plot_period("2008-2023")
p2 <- plot_period("2008-2016")
p3 <- plot_period("2017-2019")
p4 <- plot_period("2020-2023")


# 组合图形
combined_plot <- (p1 + p2 + p3 + p4) + 
     plot_layout(ncol = 4, guides = "collect") &
     theme(legend.position = "bottom")

# 显示最终图形
ggsave("figS3.pdf",
       plot = combined_plot,
       width = 12, 
       height = 5, 
       device = cairo_pdf,
       family = "Times New Roman")

ggsave(
  "figS3.png",
  plot = combined_plot, 
  width = 12, 
  height = 5
)
