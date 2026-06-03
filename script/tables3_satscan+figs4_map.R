rm(list=ls())

library(sf)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(lubridate)
library(readxl) 
library(openxlsx) 
library(cowplot)
library(patchwork) 
library(gridExtra) 
library(ggspatial) 

setwd('../outcome')
data_dir <- if (dir.exists('../data_public')) '../data_public' else '../data'
anhui_map <- st_read(file.path(data_dir, '安徽省_市.shp'))
anhui_map <- st_make_valid(anhui_map)# 修复地图错误

# 读取 DBF 文件(来自satscan)
cluster_dbf2 <- st_read(file.path(data_dir, 'spt.col.dbf')) #sf文件
id <- read_excel(file.path(data_dir, 'ID.xls'))%>%
  mutate(id=as.character(id))%>%
  mutate(name=as.character(name)) # 把id和name都转为字符型

# 数据处理(画图)
cluster_dbf3 <- cluster_dbf2 %>% 
  filter(RADIUS != 0)%>%# 排除 RADIUS 为 0 的数据，数据量减少
  mutate(year = year(ymd(START_DATE)))%>% # 提取 START_DATE 的年份，并创建新的一列 'year'
  right_join(id, by = c("LOC_ID"="id")) %>% #把cluster_dbf2中的LOC_ID替换为name
  as.data.frame() %>%
  #给每一个name都填充2008-2023年的数据,并且填充REL_RISK、RADIUS为0
  complete(name, year = 2008:2023, fill = list(REL_RISK = 0, RADIUS = 0)) %>%
  filter(!is.na(year))%>%
  select(name, year, REL_RISK,RADIUS)

#plot
satscan_data <- left_join(anhui_map,cluster_dbf3,by = "name")%>%
     select(name, year, REL_RISK,RADIUS,geometry)


####figS3
satscan_map <- ggplot(data = satscan_data) +
     geom_sf(aes(fill = REL_RISK), color = "black") +
     scale_fill_viridis_c(option = "C", name = "space-time scan by Year", direction = -1) +  # 使用更优雅的渐变色 C
     labs(title = "") +
     facet_wrap(~ year, ncol = 7) +  # 调整每行显示7个图
     theme_classic(base_family = "Times New Roman") +
     theme(
          text = element_text(color = "black"),
          legend.position = "none",  # 去掉每年图层中的图例
          strip.text = element_text(size = 12),  # 调整标签大小
          axis.text = element_blank(),  # 去掉坐标轴标签
          axis.ticks = element_blank(),# 去掉坐标轴刻度
          axis.line = element_blank()    # 去掉坐标轴线
          #将左右边距调小,上下边距调大
          ,plot.margin = margin(1, 0, 1, 0, "cm"))+
     # 调整比例尺（更左下，更小）
     annotation_scale(
          location = "bl",
          pad_x = unit(0, "cm"),  # 向左移动
          pad_y = unit(0, "cm"),  # 向下移动
          width_hint = 0.2,         # 调整宽度
          text_cex = 0.8            # 调整文字大小
     ) +
     # 调整指北针（更右上，更小）
     annotation_north_arrow(
          location = "tr",
          pad_x = unit(0, "cm"),  # 向右移动
          pad_y = unit(0.2, "cm"),  # 向上移动
          which_north = "true",
          style = north_arrow_fancy_orienteering(),
          height = unit(1, "cm"),    # 调整大小
          width = unit(1, "cm"))
satscan_map

# 使用 patchwork 合并图层，确保图例位于底部
final_map <- satscan_map +
     scale_fill_viridis_c(option = "C", name = "Relative Risk", direction = -1) +
     theme(
          legend.position = "bottom",  # 设置图例位置
          legend.title = element_text(size = 10),  # 设置图例标题字体大小
          legend.text = element_text(size = 8),  # 设置图例文本字体大小
          legend.key.size = unit(0.5, "cm"),  # 调整图例键的大小
          legend.margin = margin(t = 1, r = 0, b = 0, l = 0),  # 调整图例的上边距
          plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")  # 上、右、下、左
)
final_map

#合并
ggsave("figS4.png",plot = final_map, width = 10, height = 7, dpi = 600)
ggsave("figS4.pdf",plot = final_map, width = 10, height = 7, dpi = 600,device = cairo_pdf, family = "Times New Roman")







# TableS4
cluster_dbf4 <- cluster_dbf2 %>% 
  filter(RADIUS != 0)%>%  # 排除 RADIUS 为 0 的数据，数据量减少
  mutate(year = year(ymd(START_DATE)))%>% # 提取 START_DATE 的年份，并创建新的一列 'year'
  left_join(id, by = c("LOC_ID"="id")) %>% #把cluster_dbf2中的LOC_ID替换为name
  as.data.frame()%>%
  select(RADIUS,START_DATE,END_DATE,P_VALUE,OBSERVED,EXPECTED,REL_RISK,year,name,LLR)%>%
  arrange(year)

write.xlsx(cluster_dbf4, "tableS3.xlsx", rownames = F)





