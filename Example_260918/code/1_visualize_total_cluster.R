####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
rm(list = ls())
gc()

####导入R包####
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(patchwork)
library(RColorBrewer)
library(ggsignif)
library(FSA)
library(ggforce)
library(viridis)
library(ggsci)

####构造函数####
# 设置CNS风格主题
theme_cns <- function() {
  theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", size = 0.5),
      axis.line = element_blank(),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11),
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      legend.position = "none"
    )
}

####保守基因簇统计####
#----读入数据----#
data <- read.csv("Example_260918/data/high_conserved_summary_annotated.csv")
colnames(data)[1:3] <- c("Name","Number_clusters","List")
head(data)

#----可视化----#
## 柱状图
data$Strain <- sapply(strsplit(data$Name, split="_"), function(x) paste0(x[1],"_",x[2]))
# 按属统计平均保守簇数量和总保守簇数量
genus_summary <- data %>%
  group_by(Strain) %>%
  summarise(
    Number_strains = n(),
    Total_clusters = sum(Number_clusters),
    Mean_clusters = mean(Number_clusters),
    Have_clusters = sum(Number_clusters > 0),
    Ratio_clusters = round(Have_clusters / Number_strains * 100, 1)
  ) %>%
  arrange(desc(Total_clusters))
# 打印统计表
print(genus_summary)
# 设置颜色（Nature/Science/Cell常用色板）
# 使用8种区分度高的颜色
cns_colors <- c(
  "#E64B35",      # 红色
  "#4DBBD5",          # 青色
  "#00A087",     # 绿色
  "#3C5488",  # 深蓝
  "#F39B7F",  # 橙色
  "#8491B4",     # 灰蓝
  "#91D1C2",  # 浅绿
  "#DC0000"          # 砖红
)

# 总保守簇数
p1 <- ggplot(data, aes(x = reorder(Strain, -Number_clusters),
                       y = Number_clusters,
                       fill = Number_clusters)) +      # 关键：fill 映射数值
  geom_bar(stat = "identity", width = 0.7) +
  # scale_fill_viridis_c(option = "D", direction = -1) +  # viridis 感知均匀
  # scale_fill_gradientn(colours = pal_npg("nrc")(10)) + # 或 NPG 渐变
  scale_fill_gradient(low = "#3C5488", high = "#E64B35") + # 或双色渐变
  theme_cns() +
  labs(x = "Strain", y = "Number",
       title = "Total conserved clusters",
       fill = "Clusters") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# 显示
pdf("Example_260918/result/保守簇的整体情况.pdf", width = 7, height = 4)
print(p1)
dev.off()


