####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Random/")
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
library(ggsci)   # CNS 风格配色；若未安装：install.packages("ggsci")

####筛选基因簇####
#----读入数据----#
data <- read.csv("data/final_homogroup_results.csv")
data <- data[!is.na(data$Human.Prot..by.Comprehensive.),]
conserved_id <- read.csv("data/high_conserved_ge200.csv")

#----筛选感兴趣的簇----#
colnames(data)
data$Conserved.ID %>% unique(.)
# 确认同时六个物种都包含
interest_id <- conserved_id$保守簇ID
data_chosen <- subset(data, Conserved.ID %in% interest_id)

#----KEGG分析----#
if(!file.exists("data/hsa_kegg.txt")){
  download.file("https://rest.kegg.jp/link/hsa/pathway", "data/hsa_kegg.txt")
  download.file("https://rest.kegg.jp/list/pathway/hsa", "data/hsa_pathway_names.txt")
}

# 读取KEGG注释文件
kegg_link <- read.table("data/hsa_kegg.txt", sep = "\t", stringsAsFactors = FALSE)
colnames(kegg_link) <- c("Pathway_ID", "ENTREZID")
kegg_link$ENTREZID <- gsub("hsa:", "", kegg_link$ENTREZID)
kegg_link$Pathway_ID <- gsub("path:", "", kegg_link$Pathway_ID)

# 读取通路名称
pathway_names <- read.table("data/hsa_pathway_names.txt", sep = "\t", stringsAsFactors = FALSE)
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

# 计算基因簇在通路的映射
cluster_gene_counts <- kegg_summary %>%
  distinct(Cluster_ID, SYMBOL) %>%
  group_by(Cluster_ID) %>%
  summarise(total_genes = n(), .groups = "drop")

# 生成通路-基因列表
result <- kegg_summary %>%
  group_by(Cluster_ID, KEGG_Pathway) %>%
  summarise(
    Genes = paste(unique(SYMBOL), collapse = ";"),
    Gene_Count = n_distinct(SYMBOL),
    .groups = "drop"
  ) %>%
  left_join(cluster_gene_counts, by = "Cluster_ID") %>%
  mutate(Proportion = Gene_Count / total_genes) %>%
  dplyr::select(Cluster_ID, Pathway = KEGG_Pathway, Genes, Gene_Count, Proportion)

# 分割为列表（每个Cluster一个元素）
result_list <- result %>%
  group_split(Cluster_ID) %>%
  setNames(unique(result$Cluster_ID)) %>%
  lapply(function(df) {
    df %>% arrange(desc(Proportion))
  })

# 从每个Cluster中提取Proportion最大的通路
top_pathway_df <- result_list %>%
  lapply(function(df) {
    df %>% slice_max(Proportion, n = 1, with_ties = FALSE)
  }) %>%
  bind_rows() %>%
  dplyr::select(Cluster_ID, Pathway, Gene_Count, Proportion, Genes)

# ============================================================
# 随机抽样：与实际簇按真实 Cluster_ID 对齐
# ============================================================
# 编码蛋白的基因
all_human_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
cat("Total human genes:", length(all_human_genes), "\n")
gene_info <- select(org.Hs.eg.db,
                    keys = all_human_genes,
                    columns = c("SYMBOL", "GENETYPE"),
                    keytype = "ENTREZID")
cat("\n=== Gene Type Distribution ===\n")
table(gene_info$GENETYPE)
protein_coding_genes <- gene_info %>%
  dplyr::filter(GENETYPE == "protein-coding") %>%
  pull(ENTREZID)
cat("\nProtein-coding genes:", length(protein_coding_genes), "\n")

# 用一个 data.frame 记录随机结果，直接带真实 Cluster_ID
random_records <- list()
not_interact_num <- 0

for(i in 1:nrow(top_pathway_df)){
  set.seed(i)
  cluster_id_i <- top_pathway_df$Cluster_ID[i]         # 真实簇 ID
  # 用 conserved_id 里对应簇的基因数量；按保守簇ID匹配，避免顺序错位
  n_gene_i <- conserved_id$基因数量[match(cluster_id_i, conserved_id$保守簇ID)]
  if(is.na(n_gene_i)) n_gene_i <- conserved_id$基因数量[i]  # 兜底
  
  random_protein_genes <- sample(protein_coding_genes, size = n_gene_i)
  random_protein_symbols <- bitr(random_protein_genes,
                                 fromType = "ENTREZID",
                                 toType = "SYMBOL",
                                 OrgDb = org.Hs.eg.db)
  kegg_anno_rand <- get_kegg_from_local(gene_symbols = random_protein_symbols$SYMBOL)
  
  if(nrow(kegg_anno_rand) > 0){
    table_path <- table(kegg_anno_rand$KEGG_Pathway)
    max_count  <- max(table_path)
    max_path   <- names(table_path)[which.max(table_path)]
  } else {
    max_count  <- 0
    max_path   <- "None"
    not_interact_num <- not_interact_num + 1
  }
  
  random_records[[i]] <- data.frame(
    Cluster_ID = cluster_id_i,
    Pathway    = max_path,
    Gene_Count = as.integer(max_count),
    Proportion = NA_real_,
    Type       = "Random",
    stringsAsFactors = FALSE
  )
}

random_df <- bind_rows(random_records)
cat("Number of random sets with no KEGG hit:", not_interact_num, "\n")

# ---- 整理实际簇结果 ----
observed_df <- top_pathway_df %>%
  dplyr::select(Cluster_ID, Pathway, Gene_Count, Proportion) %>%
  mutate(Type = "Observed")

# ---- 合并 ----
plot_df <- bind_rows(observed_df, random_df) %>%
  mutate(
    Type = factor(Type, levels = c("Observed", "Random")),
    Pathway_short = str_trunc(Pathway, 40)
  )

print(plot_df)

# ============================================================
# CNS 风格合并绘图
# ============================================================
pathway_levels <- unique(plot_df$Pathway_short)
n_path <- length(pathway_levels)

# NPG + D3 配色，保证通路多时也有足够颜色
base_cols <- c(pal_npg("nrc")(10), pal_d3("category20")(20))
path_cols <- setNames(base_cols[seq_len(n_path)], pathway_levels)

p_combined <- ggplot(plot_df,
                     aes(x = Cluster_ID,
                         y = Gene_Count,
                         fill = Pathway_short,
                         alpha = Type)) +
  geom_bar(stat = "identity",
           position = position_dodge(width = 0.8),
           width = 0.7,
           color = "black", linewidth = 0.25) +
  geom_text(aes(label = Gene_Count),
            position = position_dodge(width = 0.8),
            vjust = -0.35, size = 3, family = "sans") +
  scale_fill_manual(values = path_cols, name = "KEGG Pathway") +
  scale_alpha_manual(values = c("Observed" = 1, "Random" = 0.45),
                     name = "Gene set",
                     labels = c("Observed", "Random")) +
  labs(
    title = "Top shared KEGG pathway per conserved cluster",
    x = "Conserved cluster ID",
    y = "Number of genes"
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    plot.title      = element_text(face = "bold", size = 13, hjust = 0),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y     = element_text(size = 10, color = "black"),
    axis.title      = element_text(size = 11, color = "black"),
    axis.line       = element_line(linewidth = 0.4, color = "black"),
    axis.ticks      = element_line(linewidth = 0.4, color = "black"),
    legend.position = "right",
    legend.title    = element_text(size = 9, face = "bold"),
    legend.text     = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    plot.margin     = margin(8, 8, 8, 8)
  ) +
  guides(
    fill  = guide_legend(order = 1, ncol = 1),
    alpha = guide_legend(order = 2, override.aes = list(fill = "grey30"))
  )

print(p_combined)


ggsave("result/KEGG_observed_vs_random_CNS.pdf",
       p_combined, width = 8, height = 5, units = "in", device = cairo_pdf)
ggsave("result/KEGG_observed_vs_random_CNS.png",
       p_combined, width = 8, height = 5, units = "in", dpi = 300)
