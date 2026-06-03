rm(list=ls())
##############空间自相关
library(spdep)
library(sf)
library(ggplot2)
library(dplyr)
library(ggspatial)
library(readxl)
library(tidyr)
library(gridExtra)
library(patchwork)
library(spatialreg)  
library(kableExtra)
library(writexl)

if (dir.exists('../outcome')) {
  setwd('../outcome')
} else if (dir.exists('outcome')) {
  setwd('outcome')
} else if (dir.exists('../AnhuiHFMD/outcome')) {
  setwd('../AnhuiHFMD/outcome')
} else {
  stop('Cannot find outcome/ directory from current working directory.')
}

data_dir <- if (dir.exists('../data_public')) '../data_public' else if (dir.exists('../data')) '../data' else '../AnhuiHFMD/data'

#加载历年发病率（含区县地图、区县级别发病率）
load(file.path(data_dir, 'datafile_lo.RData'))
load(file.path(data_dir, 'datafile_loo.RData'))
anhui_map <- st_read(file.path(data_dir, '安徽省_县.shp'))
anhui_map <- st_make_valid(anhui_map)# 修复地图错误
#导入各市级历年人口数据
pop <- read_excel(file.path(data_dir, 'city_pop.xlsx'))


# 求病例数每年和
datafiles_lo1 <- datafiles_city %>%
  filter(date < as.Date("2024-01-01")) %>%
  dplyr::select(date, name) %>%
  mutate(year = format(date, "%Y")) %>%
  group_by(year, name) %>%
  summarise(count = n(), .groups = "drop")

#将pop转化为长数据，得每年人口数
pop_long <- pop %>%
  pivot_longer(cols = -name, names_to = "year", values_to = "pop") %>%
  mutate(year = as.character(year))  # 确保 year 列是整数类型

# 合并病例数与人口数据，并计算发病率
datafiles_lo2 <- datafiles_lo1 %>%
  # 左连接人口数据（确保所有病例记录都保留）
  left_join(pop_long, by = c("name", "year")) %>%
  # 计算发病率（每10万人）
  mutate(
    incidence_rate = ifelse(pop > 0, (count / pop) * 10, NA_real_)
  ) %>%
  # 按年份和地区排序
  arrange(year, name)


# 确保使用正确的投影
anhui_neighbors <- poly2nb(anhui_map, queen = TRUE)

# 使用邻接矩阵创建空间权重矩阵，使用 W 加权
spatial_weights <- nb2listw(anhui_neighbors, style = "W")

# 检查邻接矩阵的基本信息
summary(spatial_weights)

# 获取区域名称
region_names <- anhui_map$name

# 为邻接矩阵添加区域名称
names(anhui_neighbors) <- region_names

###开始做市级的空间自相关
# 初始化一个空的结果数据框，包括 z 值和 p 值
global_moran <- data.frame(year = integer(), moran_i = numeric(), z_score = numeric(), p_value = numeric())


# 按年份逐个计算 Moran's I 和相关统计量
for (yr in unique(datafiles_lo2$year)) {
  # 获取当前年份的 incidence_rate 数据
  year_data <- datafiles_lo2 %>% filter(year == yr)
  
  # 检查当前年份的行数是否与邻接矩阵的区域数一致
  if (nrow(year_data) == length(anhui_neighbors)) {
    # 计算 Moran's I
    moran_result <- moran.test(year_data$incidence_rate, listw = spatial_weights)
    
    # 提取 Moran's I、z 值和 p 值
    moran_i <- moran_result$estimate["Moran I statistic"]
    z_score <- moran_result$statistic
    p_value <- moran_result$p.value
    
    # 将结果保存到数据框
    global_moran <- rbind(global_moran, data.frame(year = yr, moran_i = moran_i, z_score = z_score, p_value = p_value))
  } else {
    message(paste("Skipping year", yr, "due to mismatch in data length"))
  }
}

# 查看结果
print(global_moran)

# 导出
write_xlsx(global_moran, "tables2.xlsx")

