rm(list=ls())
####小波分析
#load package
#install.packages("WaveletComp")
library(WaveletComp)
library(ggplot2)
library(dplyr)
library(lubridate)
library(tseries)
library(reshape2) 
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
datafiles_lo$date <- as.Date(datafiles_lo$date, format="%Y-%m-%d")

# 创建一个新列，将每个日期归类为所在的周
datafiles_looo <- datafiles_lo %>%
  mutate(week = floor_date(date, unit = "week"))  # 将每个日期归类为周一的日期

# 按周汇总病例数
data <- datafiles_looo %>%
  filter(date >= "2010-01-01"&date<="2024-01-01") %>%
  group_by(week) %>%
  summarise(cases = n())%>%
  mutate(time = seq_len(nrow(.)))%>%
  select(time,cases)

# 对病例数进行对数变换 + 标准化
data1 <- data %>%
  mutate(log_cases = log(cases))

data2 <- data1 %>%
  mutate(case = scale(log_cases))

# 填充零，确保数据长度为2的幂次
n <- length(data2$case)
next_power_of_2 <- 2^ceiling(log2(n))  # 计算最近的2的幂次
padding_zeros <- next_power_of_2 - n  # 需要填充的零的数量

# 填充零到数据的两端
padded_data <- c(rep(0, floor(padding_zeros / 2)), data2$case, rep(0, ceiling(padding_zeros / 2)))  # 两端填充零
padded_data <- as.numeric(padded_data) 

# 创建一个数据框，以便将数据传递给 analyze.wavelet 函数
padded_data_df <- data.frame(series = padded_data)
padded_data_df <- padded_data_df%>%
  mutate(time = 1:1024)

# 执行小波分析
wavelet_result <- analyze.wavelet(
  my.data = padded_data_df,   # 输入数据框
  my.series = "series",       # 指定要分析的列，假设列名为 "series"
  dt = 1,                     # 时间间隔：每周一个数据点
  dj = 1/128,                  # 小波变换分辨率
  lowerPeriod = 2,         # 最小周期
  upperPeriod = floor(nrow(padded_data_df)/3)*1, # 最大周期
  make.pval = TRUE,           # 计算 p 值
  method = "white.noise",     # 白噪声方法
  n.sim = 100,                # 模拟次数
  verbose = TRUE              # 显示过程信息
)


####全局小波功率
# 提取小波分析结果中的周期和功率
periods <- wavelet_result$Period    # 获取周期
power_matrix <- wavelet_result$Power # 获取每个时间点、每个周期的功率值

# 计算全局功率：对每个周期的所有时间点的功率求平均
global_power <- apply(power_matrix, 1, mean)  # 按行求均值

# 创建一个数据框用于绘图
global_power_df <- data.frame(period = periods, power = global_power)
global_power_p <- data.frame(period = periods,p=wavelet_result$Power.avg.pval)
#将两个表格合并，按照period列合并
global_power_gl <- merge(global_power_df,global_power_p,by="period")
global_power_gl <- global_power_gl%>%
  #添加一列pv，如果p<0.05,则pv=power,否则pv=NA
  mutate(pv=ifelse(p<0.05,power,NA))

# 可视化全局小波功率图
figS2 <- ggplot(data=global_power_gl) +
  geom_line( mapping=aes(x = period, y = global_power)) +   # 加粗线条，设置颜色
  #geom_point(mapping=aes(x = period, y = pv,color="p<0.05")) +       # 设置点的大小和颜色
  scale_x_continuous(trans = 'log10') +  # 将x轴转换为对数尺度
  labs(x = "Period (weeks)", y = "Global Power") +
  ggtitle("") +
  theme_minimal(base_family = "Times New Roman" ) +   # 设置为最小主题
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中，加粗，调整字体大小
    axis.title.x = element_text(size = 14, face = "bold"),  # X轴标题加粗
    axis.title.y = element_text(size = 14, face = "bold"),  # Y轴标题加粗
    axis.text.x = element_text(size = 12),   # X轴刻度字体大小
    axis.text.y = element_text(size = 12),   # Y轴刻度字体大小
    axis.ticks = element_line(color = "black", size = 0.5),  # 刻度线颜色和粗细
    axis.ticks.length = unit(0.2, "cm"),  # 刻度线长度
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  # 添加边框
    panel.grid.major = element_blank(),   # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次网格线
    plot.background = element_rect(fill = "white", color = NA)  # 设置白色背景，无框线
  )
figS2

#合并
ggsave("figS2.png",plot = figS3, width = 15, height = 8, dpi = 600)
ggsave("figS2.pdf",plot = figS3, width = 15, height = 8, dpi = 600,device = cairo_pdf, family = "Times New Roman")
# tiff与eps
ggsave('figS2.tiff', plot = figS3, width = 15, height = 8, dpi = 700,compression = "lzw")
ggsave(
  'figS2.eps', 
  plot = figS3,  # 替换为你的图形对象
  width = 15, 
  height = 8, 
  dpi = 700,
  device = cairo_pdf,  # 关键：使用 Cairo 设备
  family = "Times New Roman"  # 强制指定字体
)



