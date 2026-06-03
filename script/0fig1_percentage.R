#
rm(list=ls())
# 安装并运行 r 包 
library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggalluvial)
library(grid)
library(stringr) #首字母大写
library(patchwork)
library(lubridate)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
     dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
     getwd()
}
repo_root <- normalizePath(file.path(script_dir, ".."))
setwd(file.path(repo_root, 'outcome'))
data_dir <- if (dir.exists(file.path(repo_root, 'data_public'))) {
     file.path(repo_root, 'data_public')
} else {
     file.path(repo_root, 'data')
}
load(file.path(data_dir, 'clean2.RData'))

#查看列名
# 去除 2023 年 12 月 31 日之后的数据
datafile_all <- subset(datafiles_2, Date <= as.Date("2023-12-31"))
datafile_all$virus[datafile_all$virus == ''] <- NA
#将datafile_all 的 所有Cox A16 改为CV-A16
datafile_all$virus[datafile_all$virus == 'Cox A16'] <- 'CV-A16'
# 修改列名（首字母大写）并调整 Gender 变量值（为了画图大写）
datafile_all <- datafile_all %>%
     rename(
          Virus = virus,
          Gender = gender,
          Age = age
     ) %>%
     mutate(
          Gender = str_to_title(Gender)  # "male" → "Male", "female" → "Female"
     )

# 定义颜色
fill_color_setting <- paletteer::paletteer_d("futurevisions::atomic_clock")

color_palettes <- list(
     gender = fill_color_setting[c(2, 4)],
     age = fill_color_setting[c(1:2, 4:5)],
     virus = rev(fill_color_setting[1:3]),
     lab = fill_color_setting[c(4:5)]
)

# A-B ---------------------------------------------------------------------

theme_curve <- function(){
     theme_bw(base_family = 'Times New Roman')+
          theme(legend.position = 'none',
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                title = element_text(size = 16, vjust = .5, face = 'bold'),
                axis.title = element_text(size = 14, face = 'bold', color = 'black'),
                axis.text.y = element_text(size = 12, family = 'Times New Roman', color = 'black'), # Y 轴刻度标签字体
                axis.ticks = element_line(color = 'black', size = 0.5), # X 和 Y 轴刻度线的颜色和粗细
                axis.line = element_line(color = 'black', size = 0.8),   # 所有轴的轴线粗细
                legend.title = element_text(size = 14, vjust = .5, face = 'bold'),
                legend.text = element_text(size = 12, vjust = .5))
}

#筛选病毒列
epicurve_report_virus <- datafile_all %>% 
     filter(!is.na(Virus)) %>%
     mutate(virus = ifelse(Virus == "其他肠道病毒", "Other", Virus),
            month = month(Date),
            year = year(Date)) |> 
     group_by(year, month, virus) %>%
     summarise(count = n(), .groups = 'drop')

epicurve_report_all <- datafile_all %>%
     mutate(month = month(Date),
            year = year(Date),
            lab = if_else(!is.na(Virus), 'Typed', 'Not typed')) |> 
     group_by(year, month, lab) |> 
     summarise(count = n(), .groups = 'drop') |> 
     complete(year, month, lab, fill = list(count = 0))
     

#定义春秋
bg_df <- data.frame(
     start_date = seq.Date(from = as.Date('2008/2/1'), 
                           to = as.Date('2023/8/1'), by = '6 months'),
     end_date = seq.Date(from = as.Date('2008/2/1'), 
                         to = as.Date('2024/4/1'), by = '6 months')[-1],
     alpha = factor(c(rep(c("Spring", "Autumn"), times = 16)),
                    levels = c("Spring", "Autumn"))
)
bg_df$alpha_1 <- ifelse(bg_df$alpha == "Autumn", 1, 0)

## panel A --------------------------------------------------------------------

# 绘制流行曲线
figA1 <- ggplot()+
     geom_rect(data = bg_df, 
               aes(xmin = start_date, xmax = end_date, alpha = alpha), 
               ymin = -Inf, ymax = Inf, fill = "#E0E0E0",show.legend = T)+
     geom_col(data = epicurve_report_all, 
              aes(x = as.Date(paste(year, month, "01", sep = "-")), 
                  y = count/1000, fill = lab))+
     geom_vline(xintercept = as.Date(c('2016-06-01', '2019-12-31')))+
     coord_cartesian(xlim = c(as.Date("2008-01-01"), as.Date("2024-01-01")))+
     scale_fill_manual(values = color_palettes$lab)+
     scale_x_date(
          expand = c(0,0),
          date_breaks = "1 year",
          #date_minor_breaks = "month",
          date_labels = "%Y"
     )+
     scale_y_continuous(
          expand = c(0, 0),
          limits = c(0, 40)
     )+
     theme_curve()+
     theme(plot.margin = margin(r = 20),
          legend.position = c(0.01, 0.99),
           legend.justification = c(0, 1),
           legend.margin = margin(r = 5),
           legend.box = 'horizontal',
           legend.direction = 'vertical',
           plot.title.position = 'plot',
           axis.text.x = element_text(size = 14, angle = 0 , hjust = 0.5),  # 显示 x 轴刻度标签
           axis.text.y = element_text(vjust=0),
           axis.ticks.x = element_line(),  # 确保显示 x 轴刻度
           axis.title.x = element_blank()) +
     labs(x = 'year',
          # using scientific notation for y-axis
          y = expression(bold("Monthly cases") ~ ~"\u00D7" ~ 10^3),
          alpha = 'Semester',
          fill = 'Case type',
          title = 'A'
     )+
     guides(fill = guide_legend(title = "Case type", order = 1),
            alpha = guide_legend(title = "Semester", order = 2))

# pie plot A
df_percent_report <- epicurve_report_all |> 
     mutate(Date = as.Date(paste(year, month, "01", sep = "-")),
            lab = factor(lab, levels = c("Typed", "Not typed")),
            vaccine_period = case_when(Date < as.Date("2016-06-01") ~ 'Before 2016-06-01',
                                       Date >= as.Date("2016-06-01") & Date <= as.Date("2019-12-31") ~ '2016-06-01 to 2019-12-31',
                                       Date > as.Date("2019-12-31") ~ 'After 2019-12-31'),
            vaccine_period = factor(vaccine_period, 
                                         levels = c('Before 2016-06-01', 
                                                    '2016-06-01 to 2019-12-31', 
                                                    'After 2019-12-31'))) |>
     group_by(vaccine_period, lab)|> 
     summarise(n = sum(count), .groups = 'drop') |>
     group_by(vaccine_period) |> 
     mutate(percent = round(n/sum(n), 4),
            cumsum_prev = lag(cumsum(percent), default = 0),
            y_pos = cumsum_prev + percent / 2) |> 
     # adjust lab factor
     mutate(lab = as.character(lab)) 

df_text <- data.frame(label = c('Stage 1', 'Stage 2', 'Stage 3'),
                      vaccine_period = c('Before 2016-06-01', '2016-06-01 to 2019-12-31', 'After 2019-12-31'))

df_text$vaccine_period <- factor(df_text$vaccine_period, 
                                 levels = c('Before 2016-06-01', 
                                            '2016-06-01 to 2019-12-31', 
                                            'After 2019-12-31'))

# 绘制饼图
figA2 <- ggplot(df_percent_report,
                 mapping = aes(x = 2, y = percent, fill = lab)) +
     geom_col(color = "white",
              show.legend = FALSE) +
     # 百分比标签
     geom_text(aes(y = y_pos,
                   label = paste0(round(percent * 100, 1), "%")),
               family   = "Times New Roman",
               size     = 4,
               fontface = "bold",
               color    = "black") +
     # 中央文字
     geom_text(data = df_text,
               aes(x = 0.5, y = 0, label = label),
               inherit.aes = FALSE,
               family   = "Times New Roman",
               size     = 4,
               fontface = "bold",
               color    = "black") +
     facet_wrap(~ vaccine_period, nrow = 1) + # 使用指定顺序的 facet_wrap
     coord_polar(theta = "y", start = 0) +
     scale_fill_manual(values = color_palettes$lab) +
     scale_x_continuous(limits = c(0.5, 2.5)) +
     theme_void() +
     theme(plot.margin = margin(),
           strip.text = element_blank(),
           panel.spacing.x = unit(6, "lines")) # 调整饼图之间的水平间隔

figA <- figA1

## panel B --------------------------------------------------------------------

# 绘制病毒流行曲线
figB1 <- ggplot() +
     geom_rect(data = bg_df, 
               aes(xmin = start_date, xmax = end_date, alpha = alpha), 
               ymin = -Inf, ymax = Inf, fill = "#E0E0E0", show.legend = T) +
     geom_col(data = epicurve_report_virus, 
              aes(x = as.Date(paste(year, month, "01", sep = "-")), 
                  y = count/1000, fill = virus)) +
     geom_vline(xintercept = as.Date(c('2016-06-01', '2019-12-31'))) +
     coord_cartesian(xlim = c(as.Date("2008-01-01"), as.Date("2024-01-01"))) +
     scale_fill_manual(values = color_palettes$virus) +
     scale_x_date(
          expand = c(0, 0),
          date_breaks = "1 year",
          date_labels = "%Y"
     ) +
     scale_y_continuous(
          expand = c(0, 0),
          limits = c(0, 2.5)
     ) +
     theme_curve() +
     theme(plot.margin = margin(r = 20),
          legend.position = c(0.01, 0.99),
           legend.justification = c(0, 1),
           legend.margin = margin(r = 5),
           legend.box = 'horizontal',
           legend.direction = 'vertical',
           plot.title.position = 'plot',
           axis.text.x = element_text(size = 14, angle = 0 , hjust = 0.5),  # 显示 x 轴刻度标签
           axis.text.y = element_text(vjust=0),
           axis.ticks.x = element_line(),  # 确保显示 x 轴刻度
           axis.title.x = element_blank()) +
     labs(x = 'year',
          y = expression(bold("Monthly cases") ~ ~"\u00D7" ~ 10^3),
          alpha = 'Semester',
          fill = 'Virus type',
          title = 'B'
     )

# pie plot B
df_percent_report_virus <- epicurve_report_virus |> 
     mutate(Date = as.Date(paste(year, month, "01", sep = "-")),
            virus = factor(virus, levels = c("Other", "EV71", "CV-A16")),
            vaccine_period = case_when(Date < as.Date("2016-06-01") ~ 'Before 2016-06-01',
                                       Date >= as.Date("2016-06-01") & Date <= as.Date("2019-12-31") ~ '2016-06-01 to 2019-12-31',
                                       Date > as.Date("2019-12-31") ~ 'After 2019-12-31'),
            vaccine_period = factor(vaccine_period, 
                                         levels = c('Before 2016-06-01', 
                                                    '2016-06-01 to 2019-12-31', 
                                                    'After 2019-12-31'))) |>
     group_by(vaccine_period, virus)|> 
     summarise(n = sum(count), .groups = 'drop') |>
     group_by(vaccine_period) |> 
     mutate(percent = round(n/sum(n), 4),
            cumsum_prev = lag(cumsum(percent), default = 0),
            y_pos = cumsum_prev + percent / 2) |> 
     # adjust lab factor
     mutate(virus = as.character(virus),
            virus = factor(virus, levels = c("CV-A16", "EV71", "Other")))

df_text_virus <- data.frame(label = c('Stage 1', 'Stage 2', 'Stage 3'),
                            vaccine_period = c('Before 2016-06-01', '2016-06-01 to 2019-12-31', 'After 2019-12-31'))

df_text_virus$vaccine_period <- factor(df_text_virus$vaccine_period, 
                                       levels = c('Before 2016-06-01', 
                                                  '2016-06-01 to 2019-12-31', 
                                                  'After 2019-12-31'))

# 绘制饼图
figB2 <- ggplot(df_percent_report_virus,
                 mapping = aes(x = 2, y = percent, fill = virus)) +
     geom_col(color = "white", show.legend = F) +
     # 百分比标签
     geom_text(aes(y = y_pos,
                   label = paste0(round(percent * 100, 1), "%")),
               family   = "Times New Roman",
               size     = 4,
               fontface = "bold",
               color    = "white") +
     # 中央文字
     geom_text(data = df_text_virus,
               aes(x = 0.5, y = 0, label = label),
               inherit.aes = FALSE,
               family   = "Times New Roman",
               size     = 4,
               fontface = "bold",
               color    = "black") +
     facet_wrap(~ vaccine_period, nrow = 1) + # 使用指定顺序的 facet_wrap
     coord_polar(theta = "y", start = 0) +
     scale_fill_manual(values = color_palettes$virus) +
     scale_x_continuous(limits = c(0.5, 2.5)) +
     theme_void() +
     theme(plot.margin = margin(),
           strip.text = element_blank(),
           panel.spacing.x = unit(6, "lines"))+ # 调整饼图之间的水平间隔
     guides(fill = guide_legend(title = "Virus type", order = 1),
            alpha = guide_legend(title = "Semester", order = 2))

figB <- figB1

# C-H ---------------------------------------------------------------------

# 生成 data1
data1 <- datafile_all %>%
     select(Date, Age, Gender, Virus) %>%
     mutate(
          Age = case_when(
               Age == "＜1" ~ "≤1",       # 处理字符型"＜1"
               as.character(Age) == "1" ~ "≤1",  # 处理数值型或字符型的1
               Age %in% c(2, 3) ~ "2-3",
               Age %in% c(4, 5) ~ "4-5",
               TRUE ~ "6+"
          ),
          Age = factor(Age, levels = c("≤1", "2-3", "4-5", "6+"))  # 设置因子顺序
     )

# 生成 data2
data2 <- datafile_all %>%
     filter(severe == "Y") %>%
     select(Date, Age, Gender, Virus) %>%
     mutate(
          Age = case_when(
               Age == "＜1" ~ "≤1",       # 处理字符型"＜1"
               as.character(Age) == "1" ~ "≤1",  # 处理数值型或字符型的1
               Age %in% c(2, 3) ~ "2-3",
               Age %in% c(4, 5) ~ "4-5",
               TRUE ~ "6+"
          ),
          Age = factor(Age, levels = c("≤1", "2-3", "4-5", "6+"))  # 更新levels
     )

# 处理 data1 - age 列
data1_age_percent <- data1 %>%
  mutate(Year = year(Date)) %>%  # 提取年份
  group_by(Year, Age) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Year) %>%
  mutate(Total = sum(Count, na.rm = TRUE),  # 计算每年age的总数
         Percent = Count / Total * 100) %>%
  ungroup()

# 处理 data1 - gender 列
data1_gender_percent <- data1 %>%
  mutate(Year = year(Date)) %>%  # 提取年份
  group_by(Year, Gender) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Year) %>%
  mutate(Total = sum(Count, na.rm = TRUE),  # 计算每年gender的总数
         Percent = Count / Total * 100) %>%
  ungroup()

# 处理 data1 - virus 列
data1_virus_percent <- data1 %>%
  mutate(Year = year(Date)) %>%  # 提取年份
  filter(!is.na(Virus)) %>%  # 忽略 NA 值
  group_by(Year, Virus) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Year) %>%
  mutate(Total = sum(Count, na.rm = TRUE),  # 计算每年virus的总数
         Percent = Count / Total * 100) %>%
  #把virus列的"其他肠道病毒"改为"Other"
  mutate(Virus = case_when(
    Virus == "其他肠道病毒" ~ "Other",
    TRUE ~ Virus
  )) %>%
  ungroup()

# 对 data2 进行相同的处理
# 处理 data2 - age 列
data2_age_percent <- data2 %>%
     mutate(Year = as.character(format(as.Date(Date), "%Y"))) %>%  
  group_by(Year, Age) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Year) %>%
  mutate(Total = sum(Count, na.rm = TRUE),  # 计算每年age的总数
         Percent = Count / Total * 100) %>%
  ungroup()

# 处理 data2 - gender 列
data2_gender_percent <- data2 %>%
  mutate(Year = format(as.Date(Date), "%Y")) %>%  # 提取年份
  group_by(Year, Gender) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Year) %>%
  mutate(Total = sum(Count, na.rm = TRUE),  # 计算每年gender的总数
         Percent = Count / Total * 100) %>%
     ungroup()
          

# 处理 data2 - virus 列
data2_virus_percent <- data2 %>%
  mutate(Year = format(as.Date(Date), "%Y")) %>%  # 提取年份
  filter(!is.na(Virus)) %>%  # 忽略 NA 值
  group_by(Year, Virus) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Year) %>%
  mutate(Total = sum(Count, na.rm = TRUE),  # 计算每年virus的总数
         Percent = Count / Total * 100) %>%
  #把virus列的"其他肠道病毒"改为"Other"
  mutate(Virus = case_when(
    Virus == "其他肠道病毒" ~ "Other",
    TRUE ~ Virus
  )) %>%
  ungroup()

#' 绘图函数
plot_sankey <- function(data, fill_var, y_title = NULL, label_pos = "A", colors) {
  
  # 处理Age的特殊排序
  if(fill_var == "Age") {
    data[[fill_var]] <- factor(data[[fill_var]], 
                               levels = c("≤1", "2-3", "4-5", "6+"),
                               ordered = TRUE) %>% 
      forcats::fct_rev()
  }
  
  # 绘图主体
  p <- ggplot(data, aes_string(x = "Year", y = "Percent", 
                               fill = fill_var, 
                               stratum = fill_var, 
                               alluvium = fill_var)) +
    geom_stratum(width = 0.6, color = 'white') +
    geom_alluvium(alpha = 0.5, width = 0.6, 
                  color = 'white', linewidth = 0.5,
                  curve_type = "linear") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(mult=0, add=0),
                     breaks = seq(2008, 2024, by = 2))+
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
    theme_percentage() +  # 使用自定义主题
    theme(
      plot.margin = margin(20, 5, 5, 15),
      legend.position = "bottom",
      plot.title.position = "plot"
    )
  
  if (label_pos %in% c("C", "F")) {
    p <- p +
         labs(y = y_title, x = NULL, title = label_pos)
  } else {
    p <- p + 
         labs(y = NULL, x = NULL, title = label_pos)
  }
}

# 自定义主题
theme_percentage <- function(){
  theme_bw(base_family = 'Times New Roman') +
    theme(
      text = element_text(color = "black", size = 13.5,family = 'Times New Roman'),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12, color = "black",family = 'Times New Roman'),
      axis.text.x = element_text(size = 12, color = "black",family = 'Times New Roman'),
      axis.title.y = element_text(size = 14, face = "bold",family = 'Times New Roman'),
      axis.title.x = element_text(size = 14,family = 'Times New Roman'),
      plot.title = element_text(size = 16, face = "bold", family = 'Times New Roman'),
      strip.text = element_text(color = "black", size = 12,family = 'Times New Roman'),
      strip.background = element_rect(color = "black", fill = "grey90"),
      legend.text = element_text(),
      legend.title = element_text()
    )
}

# 绘制fig1的所有子图(A-F)
p_I_gender <- plot_sankey(data1_gender_percent, "Gender", "All reported cases (%)", "C", color_palettes$gender)
p_I_age <- plot_sankey(data1_age_percent, "Age", NULL, "D", color_palettes$age)
p_I_virus <- plot_sankey(data1_virus_percent, "Virus", NULL, "E", color_palettes$virus)
p_s_gender <- plot_sankey(data2_gender_percent, "Gender", "Severe cases (%)", "F", color_palettes$gender)
p_s_age <- plot_sankey(data2_age_percent, "Age", NULL, "G", color_palettes$age)
p_s_virus <- plot_sankey(data2_virus_percent, "Virus", NULL, "H", color_palettes$virus)


# 组合图形（2行3列）
combined_plot <- patchwork::wrap_plots(p_I_gender + p_s_gender + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = "bottom"),
                                       p_I_age + p_s_age + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = "bottom"),
                                       p_I_virus + p_s_virus + plot_layout(ncol = 1, guides = 'collect') & theme(legend.position = "bottom"),
                                       ncol = 3)

fig <- figA + 
     figB + 
     combined_plot + 
     plot_layout(ncol = 1)


# 保存comine_plot为 PDF
ggsave(
  'fig1.pdf',
  plot = fig,
  width = 16,
  height = 13,
  device = cairo_pdf,
  family = "Times New Roman"
)

ggsave(
  'fig1.png',
  plot = fig,
  width = 16,
  height = 13
)
