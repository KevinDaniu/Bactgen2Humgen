# Author: Liu Zhehan
# 2026.9.13
####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Cluster_Visualization")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(ggplot2)
library(gggenes)

####可视化基因簇####
#----读入数据----#
gene_data <- read.csv("data/example/GCF_000007205.1_ASM720v1_genomic_Conserved_079_for_gggenes.csv")
colnames(gene_data)
gene_data <- gene_data[1:(length(unique(gene_data$homologous_id))*3),] # 抽取

#----可视化----#
pdf("result/Conserved_cluster_strand.pdf", width = 8, height = 5)
ggplot(gene_data, aes(xmin = start, xmax = end,
                      y = molecule,
                      fill = homologous_id,
                      label = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  geom_gene_label(align = "centre", size = 1.5) +
  theme_genes() +
  theme(
    strip.text.y = element_text(angle = 0, size = 6),
    legend.position = "right"
  ) +
  labs(title = "Conserved_150 Gene Cluster",
       x = "Position (bp)",
       y = "Genome")
dev.off()
