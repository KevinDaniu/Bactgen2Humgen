####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(KEGGREST)
library(networkD3)
library(tidyr)

####筛选基因簇####
#----读入数据----#
data <- read.csv("Example/data/final_homogroup_results.csv")
# data <- data[!is.na(data$Human.Gene..by.Comprehensive.) & !is.na(data$Human.Prot..by.Comprehensive.),]
data <- data[!is.na(data$Human.Prot..by.Comprehensive.),]
conserved_id <- read.csv("Example/data/high_conserved_ge200.csv")
colnames(conserved_id) <- c("Cluster_ID","Number_genes","Number_strains","List")

#----筛选感兴趣的簇----#
colnames(data)
data$Conserved.ID %>% unique(.)
# 可视化保守簇的物种共享情况
head(conserved_id[,1:3])

# 可视化
p1 <- ggplot(conserved_id, aes(x = Cluster_ID, y = Number_genes)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Number of genes in each cluster",
       x = "Cluster ID",
       y = "Number of genes")

# 保守物种数柱状图
p2 <- ggplot(conserved_id, aes(x = Cluster_ID, y = Number_strains)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Number of strains in each cluster",
       x = "Cluster ID",
       y = "Number of strains")

# 显示图形
print(p1)
print(p2)

# 方法2b：合并显示（分组柱状图）
# 先转换数据格式
# 按照保守物种数从高到低排序
# 按照保守物种数从高到低排序
conserved_id_sorted <- conserved_id[order(conserved_id$Number_strains, decreasing = TRUE), ]
conserved_id_sorted$Cluster_ID <- factor(conserved_id_sorted$Cluster_ID, 
                                    levels = conserved_id_sorted$Cluster_ID)

# 计算缩放比例，使两个指标在同一图表中显示
# 找到两个指标的最大值
max_gene <- max(conserved_id_sorted$Number_genes)
max_species <- max(conserved_id_sorted$Number_strains)
scale_factor <- max_species / max_gene

# 绘制双轴图
p3 <- ggplot(conserved_id_sorted, aes(x = Cluster_ID)) +
  # 柱状图：保守物种数（左轴）
  geom_bar(aes(y = Number_strains, fill = "Number_strains"), 
           stat = "identity", width = 0.8,
           position = position_dodge(0.9),
           alpha = 0.7) +
  # 柱状图或折线图：基因数量（右轴）
  geom_bar(aes(y = Number_genes * scale_factor, fill = "Number_genes"), 
           stat = "identity", width = 0.8,
           position = position_dodge(0.9),
           alpha = 0.7) +
  # 创建第二个纵轴
  scale_y_continuous(
    name = "Numbe of strains",
    sec.axis = sec_axis(~ . / scale_factor, name = "Number of genes")
  ) +
  scale_fill_manual(values = c("Number_strains" = "lightgreen", 
                               "Number_genes" = "skyblue")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.title = element_blank(),
        legend.position = "top") +
  labs(title = "Comparison of strains and genes in each cluster",
       x = "Cluster ID")

print(p3)

pdf("Example/result/高保守的簇整体情况.pdf", width = 5, height = 4)
print(p3)
dev.off()

data_chosen <- data[data$Conserved.ID %in% conserved_id$Cluster_ID,]

#----GO注释----#
annoGO <- function(x){
  x <- unique(x[!is.na(x) & x != ""])
  if(length(x) == 0) return(NA)
  
  id <- tryCatch({
    bitr(x, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  }, error = function(e) return(NA))
  
  if(length(id) == 1 && is.na(id)) return(NA)
  if(nrow(id) == 0) return(NA)
  
  ago_ALL <- groupGO(gene = id$ENTREZID, keyType = "ENTREZID", OrgDb = org.Hs.eg.db,
                     ont = "BP", level = 2, readable = T)
  return(as.data.frame(ago_ALL))
}

# 解析GO注释
anno_go <- lapply(unique(data_chosen$Conserved.ID), function(i) {
  proteins <- subset(data_chosen, Conserved.ID == i)$Human.Prot..by.Comprehensive.
  annoGO(proteins)
})
names(anno_go) <- unique(data_chosen$Conserved.ID)

anno_go <- anno_go[!is.na(anno_go)]
anno_sum <- data.frame(cluster_id = c(), GO_anno = c(), gene_count = c(), gene_list = c())
for(i in names(anno_go)){
  a_g <- anno_go[[i]]
  a_g <- subset(a_g, Count > 0)
  if(dim(a_g)[1] == 0){
    next
  }
  a_s <- data.frame(cluster_id = i, GO_anno = a_g$Description, gene_count = a_g$Count, gene_list = a_g$geneID)
  anno_sum <- rbind(anno_sum, a_s)
}

table(anno_sum$GO_anno)
plot_data <- as.data.frame(table(anno_sum$GO_anno))
colnames(plot_data) <- c("GO_Term", "Count")
plot_data <- plot_data[order(plot_data$Count), ]

plot_data$GO_Term_Wrapped <- str_wrap(plot_data$GO_Term, width = 30)
plot_data$ColorGroup <- cut(plot_data$Count, 
                            breaks = c(0, 50, 100, 150, 200),
                            labels = c("Low (0-50)", "Medium (50-100)", 
                                       "High (100-150)", "Very High (150+)"))
pdf("Example/result/基因GO注释柱状图.pdf", width = 5, height = 4)
ggplot(plot_data, aes(x = Count, y = reorder(GO_Term_Wrapped, Count), fill = Count)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3) +
  xlim(0, max(plot_data$Count) * 1.07) +
  scale_fill_gradient(low = "pink", high = "hotpink", name = "Gene Count") +
  labs(title = "GO Terms Distribution",
       x = "Gene Count",
       y = "") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14))
dev.off()

#----GO桑吉图----#
# 从anno_go重新整理桑基图数据
# 提取每个保守簇的GO注释
cluster_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#E5C494", "#B3B3B3", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#FDB462", "#B3DE69", "#A6D854",
  "#FFD92F", "#A65628", "#F781BF", "#999999", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#4DAF4A", "#984EA3"
)

# 获取所有Conserved簇
cluster_names <- grep("^Conserved_", names(anno_go), value = TRUE)
cat("Found", length(cluster_names), "clusters\n")

# 合并所有簇的GO数据
# 提取sankey
sankey_raw <- bind_rows(
  lapply(cluster_names, function(cluster) {
    df <- anno_go[[cluster]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      filter(Count > 0) %>%
      mutate(
        Cluster = cluster,
        GO_Term = Description,
        Genes = geneID
      ) %>%
      dplyr::select(Cluster, GO_Term, Genes)
  })
)
sankey_raw$GO_Term %>% unique(.)
# sankey_raw <- subset(sankey_raw, GO_Term %in% c("response to stimulus","immune system process"))

# 检查数据
if (is.null(sankey_raw) || nrow(sankey_raw) == 0) {
  stop("No data found!")
}

# 展开基因（GO项通常包含多个基因，用"/"分隔）
sankey_long <- sankey_raw %>%
  separate_rows(Genes, sep = "/") %>%
  filter(Genes != "")

cat("Total connections:", nrow(sankey_long), "\n")

sankey_long <- sankey_long %>%
  mutate(
    GO_Term_Short = case_when(
      GO_Term == "biological process involved in interspecies interaction between organisms" ~ "Interspecies interaction",
      GO_Term == "positive regulation of biological process" ~ "Positive regulation",
      GO_Term == "negative regulation of biological process" ~ "Negative regulation",
      GO_Term == "regulation of biological process" ~ "Regulation",
      GO_Term == "multicellular organismal process" ~ "Multicellular process",
      GO_Term == "response to stimulus" ~ "Response to stimulus",
      GO_Term == "cellular process" ~ "Cellular process",
      GO_Term == "biological regulation" ~ "Biological regulation",
      GO_Term == "immune system process" ~ "Immune response",
      TRUE ~ GO_Term
    )
  )


# 1 构建节点列表
all_clusters <- unique(sankey_long$Cluster)
all_genes <- unique(sankey_long$Genes)
all_go <- unique(sankey_long$GO_Term_Short)

nodes <- data.frame(
  name = c(all_clusters, all_genes, all_go),
  stringsAsFactors = FALSE
)

# 创建索引映射
idx <- setNames(0:(nrow(nodes) - 1), nodes$name)

# 2 构建链接：Cluster → Gene
links_1 <- sankey_long %>%
  distinct(Cluster, Genes) %>%
  mutate(
    source = idx[Cluster],
    target = idx[Genes],
    value = 1
  )

# 3 构建链接：Gene → GO_Term
links_2 <- sankey_long %>%
  distinct(Genes, GO_Term_Short) %>%
  mutate(
    source = idx[Genes],
    target = idx[GO_Term_Short],
    value = 1
  )

# 合并链接
links <- bind_rows(links_1, links_2) %>%
  group_by(source, target) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  filter(!is.na(source), !is.na(target))

cat("\n=== Data Summary ===\n")
cat("Nodes:", nrow(nodes), "\n")
cat("Links:", nrow(links), "\n")
cat("Clusters:", length(all_clusters), "\n")
cat("Genes:", length(all_genes), "\n")
cat("GO Terms:", length(all_go), "\n")
cat("Max source:", max(links$source), "(<", nrow(nodes), ")\n")
cat("Max target:", max(links$target), "(<", nrow(nodes), ")\n")

# 1 Cluster颜色
n_clusters <- length(all_clusters)
if (n_clusters <= length(cluster_colors)) {
  cluster_cols <- cluster_colors[1:n_clusters]
} else {
  cluster_cols <- colorRampPalette(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))(n_clusters)
}
names(cluster_cols) <- all_clusters

# 2 Gene颜色（统一为蓝色系）
gene_cols <- rep("#91BFDB", length(all_genes))
names(gene_cols) <- all_genes

# 3 GO Term颜色（统一为绿色系）
go_cols <- rep("#66C2A5", length(all_go))
names(go_cols) <- all_go

# 4 合并
all_colors <- c(cluster_cols, gene_cols, go_cols)
color_domain <- names(all_colors)
color_range <- unname(all_colors)

# 5 生成JS颜色代码
colour_scale <- JS(
  paste0(
    'd3.scaleOrdinal()',
    '.domain(["', paste(color_domain, collapse = '","'), '"])',
    '.range(["', paste(color_range, collapse = '","'), '"])'
  )
)

sankey <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "connections",
  fontSize = 15,
  nodeWidth = 30,
  nodePadding = 20,
  height = 600,
  width = 1000,
  colourScale = colour_scale
)

# 保存HTML
saveNetwork(sankey, file = paste0("Example/result/GO_Sankey_All_Clusters.html"))

# 显示
sankey

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

# 继续后续分析
kegg_summary <- do.call(rbind, lapply(names(kegg_results_local), function(cluster){
  temp <- kegg_results_local[[cluster]]
  temp$Cluster_ID <- cluster
  return(temp)
}))

kegg_summary

#----KEGG桑吉图----#
# 准备数据
# 1. 计算每个KEGG通路被多少个Cluster富集
colnames(kegg_summary)
pathway_cluster_count <- kegg_summary %>%
  distinct(Cluster_ID , KEGG_Pathway) %>%
  group_by(KEGG_Pathway) %>%
  summarise(n_clusters = n_distinct(Cluster_ID ), .groups = "drop") %>%
  arrange(desc(n_clusters))

cat("=== Pathway cluster distribution ===\n")
print(pathway_cluster_count)
pathway_cluster_count$n_clusters %>% table(.)

# 2. 筛选：只保留被 >= 2 个Cluster富集的通路
min_clusters <- 1  # 可调整为 2, 3, 等

selected_pathways <- pathway_cluster_count %>%
  filter(n_clusters >= min_clusters) %>%
  pull(KEGG_Pathway)

cat("\nSelected", length(selected_pathways), "pathways (shared by >= ", min_clusters, " clusters)\n")

# 3. 过滤数据
kegg_sankey_filtered <- kegg_summary %>%
  filter(KEGG_Pathway %in% selected_pathways) %>%
  mutate(
    KEGG_Pathway_Short = ifelse(nchar(KEGG_Pathway) > 35,
                                paste0(substr(KEGG_Pathway, 1, 32), "..."),
                                KEGG_Pathway)
  ) %>%
  distinct(Cluster_ID, SYMBOL, KEGG_Pathway_Short, KEGG_ID) %>%
  arrange(Cluster_ID, SYMBOL, KEGG_Pathway_Short)

cat("Filtered rows:", nrow(kegg_sankey_filtered), "\n")
cat("Pathways:", length(unique(kegg_sankey_filtered$KEGG_Pathway_Short)), "\n\n")

# 4. 构建桑基图（使用 filtered 数据）
nodes <- data.frame(
  name = unique(c(
    kegg_sankey_filtered$Cluster_ID,
    kegg_sankey_filtered$SYMBOL,
    kegg_sankey_filtered$KEGG_Pathway_Short
  ))
)

links_cluster_gene <- kegg_sankey_filtered %>%
  distinct(Cluster_ID, SYMBOL) %>%
  mutate(
    source = match(Cluster_ID, nodes$name) - 1,
    target = match(SYMBOL, nodes$name) - 1,
    value = 1
  )

links_gene_path <- kegg_sankey_filtered %>%
  distinct(SYMBOL, KEGG_Pathway_Short) %>%
  mutate(
    source = match(SYMBOL, nodes$name) - 1,
    target = match(KEGG_Pathway_Short, nodes$name) - 1,
    value = 1
  )

links <- bind_rows(links_cluster_gene, links_gene_path) %>%
  group_by(source, target) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  filter(!is.na(source), !is.na(target))

# 5. 颜色（自动）
all_clusters <- unique(kegg_sankey_filtered$Cluster_ID)
n_clusters <- length(all_clusters)

cluster_colors <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
  "#E5C494", "#B3B3B3", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#FDB462", "#B3DE69", "#A6D854"
)

if (n_clusters <= length(cluster_colors)) {
  cluster_cols <- cluster_colors[1:n_clusters]
} else {
  cluster_cols <- colorRampPalette(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))(n_clusters)
}
names(cluster_cols) <- all_clusters

all_genes <- unique(kegg_sankey_filtered$SYMBOL)
gene_cols <- rep("#91BFDB", length(all_genes))
names(gene_cols) <- all_genes

all_pathways <- unique(kegg_sankey_filtered$KEGG_Pathway_Short)
pathway_cols <- rep("#FDB462", length(all_pathways))
names(pathway_cols) <- all_pathways

all_colors <- c(cluster_cols, gene_cols, pathway_cols)
color_domain <- names(all_colors)
color_range <- unname(all_colors)

colour_scale <- JS(
  paste0(
    'd3.scaleOrdinal()',
    '.domain(["', paste(color_domain, collapse = '","'), '"])',
    '.range(["', paste(color_range, collapse = '","'), '"])'
  )
)

# 6. 绘制（放大尺寸）
sankey_kegg <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "name",
  units = "connections",
  fontSize = 15,
  nodeWidth = 30,
  nodePadding = 20,
  height = 600,
  width = 1000,
  colourScale = colour_scale
)

# 7. 保存
saveNetwork(sankey_kegg, file = paste0("Example/result/KEGG_Sankey_SharedPathways.html"))
sankey_kegg

# 8. 输出统计
cat("\n=== Summary ===\n")
cat("Clusters:", length(all_clusters), "\n")
cat("Genes:", length(all_genes), "\n")
cat("Pathways (shared):", length(all_pathways), "\n")
cat("Links:", nrow(links), "\n")

####保存桑基图####
library(webshot2)
Sys.setenv(CHROMOTE_CHROME =
             "C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe")

# 直接对 htmlwidget 截图存 PDF
webshot2::webshot(
  url    = "Example/result/KEGG_Sankey_SharedPathways.html",
  file   = "Example/result/KEGG_Sankey_SharedPathways.pdf",
  vwidth = 1000,   # 视口宽，和 sankeyNetwork 的 width 对应
  vheight = 600,   # 视口高
  zoom   = 2       # 放大倍数，提高分辨率
)
webshot2::webshot(
  url    = "Example/result/GO_Sankey_All_Clusters.html",
  file   = "Example/result/GO_Sankey_All_Clusters.pdf",
  vwidth = 1000,   # 视口宽，和 sankeyNetwork 的 width 对应
  vheight = 600,   # 视口高
  zoom   = 2       # 放大倍数，提高分辨率
)