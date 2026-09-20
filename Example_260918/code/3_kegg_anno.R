####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(KEGGREST)
library(networkD3)
library(tidyr)

####筛选基因簇####
#----读入数据----#
data <- read.csv("Example/data/final_homogroup_results.csv")
data <- data[!is.na(data$Human.Prot..by.Comprehensive.),]
conserved_id <- read.csv("Example/data/high_conserved_ge200.csv")
colnames(conserved_id) <- c("Cluster_ID","Number_genes","Number_strains","List")

#----筛选感兴趣的簇----#
colnames(data)
data$Conserved.ID %>% unique(.)
interest_id <- conserved_id$Cluster_ID
data_chosen <- subset(data, Conserved.ID %in% interest_id)

#----KEGG分析----#
if(!file.exists("Example/data/hsa_kegg.txt")){
  download.file("https://rest.kegg.jp/link/hsa/pathway", "Example/data/hsa_kegg.txt")
  download.file("https://rest.kegg.jp/list/pathway/hsa", "Example/data/hsa_pathway_names.txt")
}

# 读取KEGG注释文件
kegg_link <- read.table("Example/data/hsa_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(kegg_link) <- c("Pathway_ID", "ENTREZID")
kegg_link$ENTREZID <- gsub("hsa:", "", kegg_link$ENTREZID)
kegg_link$Pathway_ID <- gsub("path:", "", kegg_link$Pathway_ID)

# 读取通路名称
pathway_names <- read.table("Example/data/hsa_pathway_names.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(pathway_names) <- c("Pathway_ID", "Pathway_Name")
pathway_names$Pathway_Name <- gsub(" - Homo sapiens (human)", "", pathway_names$Pathway_Name, fixed = T)

# 合并
kegg_annotation <- merge(kegg_link, pathway_names, by = "Pathway_ID")

# 创建基因到KEGG的映射
gene_to_kegg <- split(kegg_annotation, kegg_annotation$ENTREZID)

# 快速注释函数
get_kegg_from_local <- function(gene_symbols){
  # 转换ENTREZID
  id <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  kegg_info <- data.frame()
  for(i in 1:nrow(id)){
    entrez <- as.character(id$ENTREZID[i])
    if(entrez %in% names(gene_to_kegg)){
      pathways <- gene_to_kegg[[entrez]]
      for(j in 1:nrow(pathways)){
        temp <- data.frame(
          SYMBOL = id$SYMBOL[i],
          ENTREZID = entrez,
          KEGG_ID = pathways$Pathway_ID[j],
          KEGG_Pathway = pathways$Pathway_Name[j],
          stringsAsFactors = FALSE
        )
        kegg_info <- rbind(kegg_info, temp)
      }
    }
  }
  return(kegg_info)
}

# 批量注释
all_genes_unique <- unique(data_chosen$Human.Prot..by.Comprehensive.)
all_genes_unique <- all_genes_unique[!is.na(all_genes_unique) & all_genes_unique != ""]

kegg_anno_local <- get_kegg_from_local(gene_symbols = all_genes_unique)

# 按基因簇分组
kegg_results_local <- list()
for(cluster in unique(data_chosen$Conserved.ID)){
  genes_in_cluster <- unique(subset(data_chosen, Conserved.ID == cluster)$Human.Prot..by.Comprehensive.)
  genes_in_cluster <- genes_in_cluster[!is.na(genes_in_cluster) & genes_in_cluster != ""]
  
  cluster_anno <- kegg_anno_local[kegg_anno_local$SYMBOL %in% genes_in_cluster, ]
  if(nrow(cluster_anno) > 0){
    kegg_results_local[[cluster]] <- cluster_anno
  }
}

# 合并所有KEGG注释
kegg_summary <- do.call(rbind, lapply(names(kegg_results_local), function(cluster){
  temp <- kegg_results_local[[cluster]]
  temp$Cluster_ID <- cluster
  return(temp)
}))

# ============ 计算每个簇中的总基因数 ============
cluster_total_genes <- data_chosen %>%
  group_by(Conserved.ID) %>%
  summarise(
    Cluster_Total_Genes = n_distinct(Human.Prot..by.Comprehensive.),
    .groups = "drop"
  )

colnames(cluster_total_genes)[1] <- "Cluster_ID"

# ============ 第一个图：每个簇中通路映射基因数 vs 该簇总基因数 ============
# 计算每个簇中每个通路的映射基因数
pathway_gene_counts <- kegg_summary %>%
  group_by(Cluster_ID, KEGG_Pathway) %>%
  summarise(
    Pathway_Genes = n_distinct(SYMBOL),
    Genes = paste(unique(SYMBOL), collapse = ";"),
    .groups = "drop"
  )

# 合并簇总基因数
first_plot_data <- pathway_gene_counts %>%
  left_join(cluster_total_genes, by = "Cluster_ID") %>%
  mutate(
    Proportion = Pathway_Genes / Cluster_Total_Genes,
    Proportion_Percent = Proportion * 100
  )

# 每个簇选择Pathway_Genes最大的通路作为代表（或者选择比例最高的）
# 这里选择比例最高的通路
top_pathway_first <- first_plot_data %>%
  group_by(Cluster_ID) %>%
  slice_max(Proportion, n = 1, with_ties = FALSE) %>%
  ungroup()

# 绘制第一个图 - 双轴图
# 为了双轴可视化，将总基因数缩放到与通路基因数相同的尺度
max_pathway_genes <- max(top_pathway_first$Pathway_Genes)
max_total_genes <- max(top_pathway_first$Cluster_Total_Genes)
scale_factor_first <- max_pathway_genes / max_total_genes

top_pathway_first$Total_Genes_Scaled <- top_pathway_first$Cluster_Total_Genes * scale_factor_first

p1 <- ggplot(top_pathway_first, aes(x = Cluster_ID)) +
  # 左侧柱子：通路映射基因数（红色）
  geom_bar(aes(y = Pathway_Genes, fill = "Pathway Genes"), 
           stat = "identity", width = 0.35, 
           position = position_nudge(x = -0.2)) +
  # 右侧柱子：簇总基因数（蓝色），缩放到左侧尺度
  geom_bar(aes(y = Total_Genes_Scaled, fill = "Total Cluster Genes"), 
           stat = "identity", width = 0.35, 
           position = position_nudge(x = 0.2), alpha = 0.8) +
  # 左侧柱子的标签（通路基因数）
  geom_text(aes(y = Pathway_Genes, label = Pathway_Genes), 
            position = position_nudge(x = -0.2),
            vjust = -0.3, size = 3.5, color = "#E74C3C") +
  # 右侧柱子的标签（簇总基因数，显示实际值）
  geom_text(aes(y = Total_Genes_Scaled, label = Cluster_Total_Genes), 
            position = position_nudge(x = 0.2),
            vjust = -0.3, size = 3.5, color = "#3498DB") +
  # 左侧Y轴：通路映射基因数
  scale_y_continuous(
    name = "Pathway Gene Count",
    sec.axis = sec_axis(
      trans = ~ . / scale_factor_first,
      name = "Total Cluster Genes"
    ),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  # 颜色设置
  scale_fill_manual(
    values = c("Pathway Genes" = "#E74C3C", 
               "Total Cluster Genes" = "#3498DB"),
    name = "Metric"
  ) +
  labs(
    title = "Top KEGG Pathway per Cluster",
    x = "Cluster ID"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.title.y.left = element_text(size = 12, color = "#E74C3C", margin = margin(r = 10)),
    axis.text.y.left = element_text(size = 10, color = "#E74C3C"),
    axis.title.y.right = element_text(size = 12, color = "#3498DB", margin = margin(l = 10)),
    axis.text.y.right = element_text(size = 10, color = "#3498DB"),
    legend.position = "top",
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank()
  )

print(p1)

ggsave("Example/result/Top_Pathway_per_Cluster.pdf", p1, width = 5, height = 4, dpi = 300)

# ============ 第二个图：覆盖率 = 通路映射基因数 / 该簇总基因数 ============
# 计算覆盖率（占簇总基因数的百分比）
second_plot_data <- first_plot_data %>%
  mutate(
    Coverage = (Pathway_Genes / Cluster_Total_Genes) * 100
  )

# 每个簇选择覆盖率最大的通路
top_pathway_second <- second_plot_data %>%
  group_by(Cluster_ID) %>%
  slice_max(Coverage, n = 1, with_ties = FALSE) %>%
  ungroup()

# 获取通路数量用于颜色
n_pathways <- length(unique(top_pathway_second$KEGG_Pathway))
colors <- brewer.pal(min(n_pathways, 12), "Set3")
if(n_pathways > 12) {
  colors <- rep(colors, length.out = n_pathways)
}

# 绘制第二个图 - 覆盖率柱状图
p2 <- ggplot(top_pathway_second, aes(x = Cluster_ID, y = Coverage, fill = KEGG_Pathway)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(round(Coverage, 1), "%\n(", Pathway_Genes, "/", Cluster_Total_Genes, ")")), 
            vjust = -0.3, size = 3.5, lineheight = 0.8) +
  labs(
    title = "KEGG Pathway Coverage per Cluster",
    x = "Cluster ID",
    y = "Coverage (%)",
    fill = "KEGG Pathway"
  ) +
  ylim(0, max(top_pathway_second$Coverage) * 1.15) +
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p2)
ggsave("Example/result/Pathway_Coverage_per_Cluster.pdf", p2, width = 8, height = 4, dpi = 300)

# ============ 导出结果表格 ============
# 第一个图的数据
write.csv(top_pathway_first, "Example/result/Top_Pathway_Data.csv", row.names = FALSE)

# 第二个图的数据
write.csv(top_pathway_second, "Example/result/Coverage_Data.csv", row.names = FALSE)

# 完整数据
write.csv(first_plot_data, "Example/result/All_Pathway_Data.csv", row.names = FALSE)

# ============ 计算KEGG数据库中每个通路的总基因数 ============
pathway_total_genes <- kegg_annotation %>%
  group_by(Pathway_ID, Pathway_Name) %>%
  summarise(Total_Genes = n_distinct(ENTREZID), .groups = "drop")

# ============ 继续原有分析 ============
# 创建基因到KEGG的映射
gene_to_kegg <- split(kegg_annotation, kegg_annotation$ENTREZID)

# 快速注释函数
get_kegg_from_local <- function(gene_symbols){
  # 转换ENTREZID
  id <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  
  kegg_info <- data.frame()
  for(i in 1:nrow(id)){
    entrez <- as.character(id$ENTREZID[i])
    if(entrez %in% names(gene_to_kegg)){
      pathways <- gene_to_kegg[[entrez]]
      for(j in 1:nrow(pathways)){
        temp <- data.frame(
          SYMBOL = id$SYMBOL[i],
          ENTREZID = entrez,
          KEGG_ID = pathways$Pathway_ID[j],
          KEGG_Pathway = pathways$Pathway_Name[j],
          stringsAsFactors = FALSE
        )
        kegg_info <- rbind(kegg_info, temp)
      }
    }
  }
  return(kegg_info)
}

# 批量注释
all_genes_unique <- unique(data_chosen$Human.Prot..by.Comprehensive.)
all_genes_unique <- all_genes_unique[!is.na(all_genes_unique) & all_genes_unique != ""]

kegg_anno_local <- get_kegg_from_local(gene_symbols = all_genes_unique)

# 按基因簇分组
kegg_results_local <- list()
for(cluster in unique(data_chosen$Conserved.ID)){
  genes_in_cluster <- unique(subset(data_chosen, Conserved.ID == cluster)$Human.Prot..by.Comprehensive.)
  genes_in_cluster <- genes_in_cluster[!is.na(genes_in_cluster) & genes_in_cluster != ""]
  
  cluster_anno <- kegg_anno_local[kegg_anno_local$SYMBOL %in% genes_in_cluster, ]
  if(nrow(cluster_anno) > 0){
    kegg_results_local[[cluster]] <- cluster_anno
  }
}

# 继续后续分析
kegg_summary <- do.call(rbind, lapply(names(kegg_results_local), function(cluster){
  temp <- kegg_results_local[[cluster]]
  temp$Cluster_ID <- cluster
  return(temp)
}))

# ============ 修改关键部分：计算覆盖率 ============
# 计算每个cluster中每个通路的映射基因数
result <- kegg_summary %>%
  group_by(Cluster_ID, KEGG_Pathway) %>%
  summarise(
    Genes = paste(unique(SYMBOL), collapse = ";"),
    Gene_Count = n_distinct(SYMBOL),
    .groups = "drop"
  )

# 添加KEGG数据库中的总基因数
result <- result %>%
  left_join(pathway_total_genes, by = c("KEGG_Pathway" = "Pathway_Name")) %>%
  # 计算覆盖率 = 映射基因数 / 通路总基因数
  mutate(
    Coverage = round(Gene_Count / Total_Genes * 100, 1),
    Proportion = Gene_Count / Total_Genes
  ) %>%
  dplyr::select(Cluster_ID, Pathway = KEGG_Pathway, Genes, Gene_Count, Total_Genes, Coverage, Proportion)

# ============ 提取每个cluster中覆盖率最大的通路 ============
top_pathway_df <- result %>%
  group_by(Cluster_ID) %>%
  slice_max(Coverage, n = 1, with_ties = FALSE) %>%
  ungroup()

# 查看结果
print(top_pathway_df)

# ============ 绘制覆盖率柱状图 ============
library(RColorBrewer)

# 获取通路数量
n_pathways <- length(unique(top_pathway_df$Pathway))

# 选择颜色方案
# 方案1: Set3 - 12种柔和的颜色
colors <- brewer.pal(min(n_pathways, 12), "Set3")
if(n_pathways > 12) {
  # 如果需要更多颜色，循环使用
  colors <- rep(colors, length.out = n_pathways)
}

# 或者使用Paired - 12种配对颜色
colors_paired <- brewer.pal(min(n_pathways, 12), "Paired")
if(n_pathways > 12) {
  colors_paired <- rep(colors_paired, length.out = n_pathways)
}

# 或者使用Spectral - 11种渐变色
colors_spectral <- brewer.pal(min(n_pathways, 11), "Spectral")
if(n_pathways > 11) {
  colors_spectral <- rep(colors_spectral, length.out = n_pathways)
}

# 绘图
p_coverage <- ggplot(top_pathway_df, aes(x = Cluster_ID, y = Coverage, fill = Pathway)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = paste0(Coverage, "%\n(", Gene_Count, "/", Total_Genes, ")")), 
            vjust = -0.3, size = 3.5, lineheight = 0.8) +
  labs(
    title = "Top KEGG Pathway Coverage per Cluster",
    x = "Cluster ID",
    y = "Coverage (%)",
    fill = "KEGG Pathway"
  ) +
  ylim(0, max(top_pathway_df$Coverage) * 1.15) +
  # 使用Set3调色板
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major.x = element_blank()
  )

print(p_coverage)
ggsave("Example/result/Pathway_Coverage_per_Pathway.pdf", p_coverage, width = 7, height = 4, dpi = 300)

