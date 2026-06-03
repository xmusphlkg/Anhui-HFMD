rm(list=ls())

# 加载必要的包
library(dplyr)
library(ggplot2)
library(lubridate)
library(viridis)

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
load(file.path(data_dir, 'datafile_lo.RData'))
data <- datafiles_lo%>%
  select(name,birth,date,address)%>%
  mutate(address = substr(address, 1, 4))


# 创建地址与名称对应关系
address_to_name <- data.frame(
  address = c("3401", "3402", "3403", "3404", "3405", "3406", "3407", "3408", 
              "3410", "3411", "3412", "3413", "3414", "3415", "3416", "3417", "3418"),
  new_name = c("Hefei", "Wuhu", "Bengbu", "Huainan", "Maanshan", "Huaibei", "Tongling", "Anqing", 
               "Huangshan", "Chuzhou", "Fuyang", "Suzhou", "Hefei", "Liuan", "Bozhou", "Chizhou", "Xuancheng")
)
# 数据处理
data0 <- data %>%
  left_join(address_to_name, by = "address") %>%
  mutate(name = new_name) %>%  
  select(-new_name)%>%
  select(date,name)%>%
  filter(date < as.Date("2024-01-01")&date >= as.Date("2010-01-01"))

# 数据处理
processed_data <- data0 %>%
  mutate(
    date = as.Date(date),
    # 处理闰日（保持与原始数据年一致）
    date = if_else(
      format(date, "%m-%d") == "02-29",
      as.Date(paste0(format(date, "%Y"), "-02-28")),  # 保持原始年份
      date
    ),
    # 计算年度周数（强制限制到52周）
    day_of_year = yday(date),
    week = pmin(ceiling(day_of_year / 7), 52)  # 关键修改：强制限制周数范围
  ) %>%
  # 合并城市计数
  group_by(name, week) %>%
  summarise(cases = n(), .groups = "drop") %>%
  # 计算百分比
  group_by(name) %>%
  mutate(
    total_cases = sum(cases),
    percent = cases / total_cases*100
  ) %>%
  ungroup()

#导入纬度
# 创建数据框（按纬度从高到低排序）
city_lat_data <- data.frame(
  name = c("Bozhou", "Suzhou", "Huaibei", "Fuyang", "Bengbu", 
           "Huainan", "Liuan", "Hefei", "Chuzhou", "Anqing", 
           "Maanshan", "Wuhu", "Tongling", "Xuancheng", "Chizhou", 
           "Huangshan"),
  lat = c(33.87, 33.63, 33.95, 32.90, 32.93, 
          32.62, 31.73, 31.87, 32.30, 30.52, 
          31.62, 31.19, 30.93, 30.95, 30.39, 
          29.71)
)
city_lat_data <- city_lat_data %>% arrange(desc(lat))


processed_data <- processed_data %>%
  left_join(city_lat_data, by = "name") %>%
  # 合并纬度数据
  mutate(
    name = factor(name, levels = rev(city_lat_data$name))  # 按纬度降序排列
  )


###绘图
# 生成精确的12等分刻度线位置（包含0和52）
month_breaks <- seq(0, 52, length.out = 13)

# 计算标签应显示的位置（两个刻度线之间的中点）
adjusted_breaks <- month_breaks[-1] - diff(month_breaks)/2

# 生成英文月份标签
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


figS1 <- ggplot(processed_data, aes(x = week - 0.5, y = name, fill = percent)) +  # 关键修改：x轴偏移
  geom_tile(color = NA, width = 1) +  # 确保宽度严格1个单位
  scale_fill_gradient(low = "white", high = "steelblue", name = "Percent (%)") +
  labs(x = "Month", y = "City", title = "") +
  theme_minimal(base_family = "Times New Roman" ) +
  theme(
    text = element_text(color="black",size = 13.5), #字体
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", size = 1, fill = NA), #边框
    axis.ticks = element_line(size = 1, color = "black"), #刻度线
    axis.ticks.length = unit(0.3, "cm"), #刻度线长度
    axis.title.x = element_text(margin = margin(t = 20),size= 14, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10),size= 14, face = "bold"),
    axis.text.y = element_text(
      color = "black",  # 强制纵坐标标签为纯黑
      margin = margin(r = 5)  # 右侧留白
    ),
    plot.title = element_text(face = "bold",family = "Times New Roman", size = 16, hjust = 0.5), #标题),
    legend.key.width = unit(1, "cm"),
    legend.key.height = unit(1, "cm"),
    legend.position = "right",
    plot.margin = margin(b = 25, t = 5, r = 5, l = 5, unit = "pt")
  ) +
  scale_x_continuous(
    breaks = month_breaks,  # 对齐新的坐标定位
    labels = NULL,
    limits = c(0, 52),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  annotate("text",
           x = adjusted_breaks,  # 对应新的坐标定位
           y = 0.2,
           label = month_labels,
           size = 4,
           family = "Times New Roman",
           color = "black",
           vjust = 1) +
  coord_cartesian(
    clip = "off",
    ylim = c(0.5, length(unique(processed_data$name)) + 0.5)
  )
figS1

##添加纬度排序

ggsave("figS1.pdf", plot=figS2,width = 15, height = 8, dpi = 600,device = cairo_pdf)
ggsave("figS1.png", plot=figS2,width = 15, height = 7, dpi = 600)
# tiff与eps
ggsave('figS1.tiff', plot = figS2, width = 15, height = 8, dpi = 700,compression = "lzw")
ggsave(
  'figS1.eps', 
  plot = figS2,  # 替换为你的图形对象
  width = 15, 
  height = 8, 
  dpi = 700,
  device = cairo_pdf,  # 关键：使用 Cairo 设备
  family = "Times New Roman"  # 强制指定字体
)

