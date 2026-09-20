####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Example")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(ggplot2)
library(gggenes)
library(scales)
library(stringr)

####可视化基因簇####
#----读入数据----#
gene_data <- read.csv("data/GCF_000007205.1_ASM720v1_genomic_Conserved_079_for_gggenes.csv")
colnames(gene_data)
dim(gene_data)
head(gene_data, n = 18)

#----可视化----#
# 每个 cluster 选多少个 genome
n_genome_per_cluster <- 6

set.seed(123)  # 可复现

gene_sub <- gene_data %>%
  distinct(cluster, molecule) %>%
  group_by(cluster) %>%
  slice_sample(n = n_genome_per_cluster) %>%
  ungroup() %>%
  left_join(gene_data, by = c("cluster", "molecule"))

# 检查每个 cluster 抽到的 genome 数是否一致
gene_sub %>%
  distinct(cluster, molecule) %>%
  count(cluster)

# 固定 molecule 顺序，避免分面顺序乱跳
gene_sub <- gene_sub %>%
  mutate(molecule_label = str_wrap(gsub("_", " ", as.character(molecule)), width = 30))
gene_sub <- gene_sub %>%
  mutate(molecule_label = factor(molecule_label, levels = unique(molecule_label)))

# CNS / Nature 风格柔和配色
cns_colors <- c(
  "#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
  "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
  "#9C755F", "#BAB0AC"
)

pdf("result/Conserved_cluster_strand_subset.pdf", width = 8, height = 6)
ggplot(gene_sub, aes(xmin = start, xmax = end,
                     y = molecule_label,
                     fill = homologous_id,
                     label = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule_label, scales = "free", ncol = 1) +
  geom_gene_label(align = "centre", size = 1.5) +
  scale_fill_manual(values = cns_colors) +
  scale_x_continuous(
    labels = label_scientific(digits = 2),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  theme_genes() +
  theme(
    strip.text.y = element_text(angle = 0, size = 6),
    legend.position = "right",
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold")
  ) +
  labs(title = "Conserved_079 Gene Cluster (subset)",
       x = "Position (bp)",
       y = "Genome")

dev.off()

