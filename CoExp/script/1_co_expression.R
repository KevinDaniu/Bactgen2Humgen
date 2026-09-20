#!/usr/bin/env Rscript
# Author: Liu Zhehan
# Time: 2026.9.13
# ============================================================
# Co-Expression Analysis
# Input:
#       (1) Param 1: final_homogroup_results.csv
#       (2) Param 2: Conserved ID
# Data source: GTEx v8 (17,382 samples, 56,200 genes)
# Method: Spearman Correlation + Benjamini-Hochberg FDR-correction
# 
# dependency: data.table, curl, httr, jsonlite, psych, stats, ggplot2, pheatmap
# ============================================================

####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/CoExp")
rm(list = ls())
gc()

####加载R包####
library(data.table)
library(curl)
library(httr)
library(jsonlite)
library(psych)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(ggsci)
library(EnsDb.Hsapiens.v86)

####共表达分析####
#----参数设置----#
# 设置参数
INPUT_FILE <- "data/example/final_homogroup_results.csv"
data <- read.csv(INPUT_FILE)
CID <- "Conserved_079"
GENES <- subset(data, Conserved.ID == CID)$`Human.Prot..by.Comprehensive.`
GTEX_TPM_LOCAL <- "gtex_v8_tpm.gct"
TISSUE_LOCAL <- "sample_tissue_map.csv"
OUTPUT_DIR <- "result"

# 判断基因是否存在
if(sum(GENES != "") <= 1){
  stop("基因不足2个！")
}

# 创建输出目录
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 分析的组织列表
KEY_TISSUES <- c(
  "Liver", "Heart - Left Ventricle", "Heart - Atrial Appendage",
  "Muscle - Skeletal", "Kidney - Cortex", "Brain - Cortex",
  "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)",
  "Whole Blood", "Thyroid", "Lung", "Pancreas",
  "Adrenal Gland", "Testis", "Ovary"
)

# ============================================================
# 打印配置信息
# ============================================================
cat("=", rep("=", 60), "\n", sep="")
cat("Coexpression Analysis:", paste(GENES, collapse=", "), "\n")
cat("Data: GTEx v8 (17,382 samples)\n")
cat("Method: Spearman correlation + BH-FDR\n")
cat("=", rep("=", 60), "\n\n", sep="")


# ============================================================
# Step 2: 获取 GTEx 样本-组织映射
# ============================================================
cat("\nFetching GTEx sample-tissue mapping via API v2 ...\n")
if(!file.exists("data/GTEx_Data/sample_tissue_map.csv")){
  all_samples <- list()
  page <- 0
  while (TRUE) {
    url <- paste0("https://gtexportal.org/api/v2/dataset/sample?dataset=gtex_v8&pageSize=250&page=", page)
    response <- GET(url)
    data <- content(response, "parsed")$data
    if (length(data) == 0) break
    for (s in data) {
      all_samples <- c(all_samples, list(list(
        sampleId = s$tissueSampleId %||% "",
        tissue = s$tissueSiteDetail %||% "",
        subjectId = s$subjectId %||% ""
      )))
    }
    page <- page + 1
  }
  df_samples <- rbindlist(lapply(all_samples, as.data.table))
  cat("  Retrieved", nrow(df_samples), "samples across", uniqueN(df_samples$tissue), "tissues\n")
  write.csv(df_samples, file.path(OUTPUT_DIR, "sample_tissue_map.csv"))
} else {
  df_sample <- read.csv("data/GTEx_Data/sample_tissue_map.csv")
  df_sample <- df_sample[,-1]
  df_sample <- df_sample[!duplicated(df_sample),]
}


# ============================================================
# Step 3: 从 GCT 文件提取目标基因的 TPM
# ============================================================
cat("\nExtracting gene TPM from GCT file ...\n")

# 获取 Ensembl ID 映射
# gene_ensembl <- list()
# for (g in GENES) {
#   url <- paste0("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/", g, "?content-type=application/json")
#   response <- GET(url)
#   if (status_code(response) == 200) {
#     data <- content(response, "parsed")
#     for (item in data) {
#       if (!is.null(item$type) && item$type == "gene") {
#         gene_ensembl[[g]] <- item$id
#         break
#       }
#     }
#   }
# }
# cat("Gene -> Ensembl ID:", paste(names(gene_ensembl), unlist(gene_ensembl), sep="=", collapse=", "), "\n")
gene_ensembl <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = GENES,
  column = "GENEID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_ensembl <- as.list(gene_ensembl)
cat("Gene -> Ensembl ID:", 
    paste(names(gene_ensembl), unlist(gene_ensembl), sep = "=", collapse = ", "), 
    "\n")

# 读取GCT文件并提取目标基因
cat("\nExtracting gene TPM from GCT file...\n")

# 获取样本ID
con <- file(file.path("data/GTEx_Data", GTEX_TPM_LOCAL), "rt")
sample_ids <- strsplit(readLines(con, n = 3)[3], "\t")[[1]][-(1:2)]
close(con)

# 使用fread一次性读取全部，然后过滤（比逐行循环快10倍以上）
cat("Reading entire GCT file with fread...\n")
df_all <- fread(file.path("data/GTEx_Data", GTEX_TPM_LOCAL), 
                skip = 3, header = FALSE, 
                nThread = 4,  # 多线程
                showProgress = TRUE)

# 过滤目标基因（用基因符号或ensembl_id）
all_ensembl <- sapply(df_all$V1, function(x){
  chrs <- strsplit(x, split = ".", fixed = T)[[1]]
  return(chrs[1])
})
all_genes <- df_all$V2
df <- df_all[all_genes %in% GENES | all_ensembl %in% unlist(gene_ensembl)]

# 清理内存
rm(df_all); gc()

# 设置列名
setnames(df, c("ensembl_id", "gene_symbol", sample_ids))

# 如果有些基因没找到，用ensembl_id再匹配一次
found_genes <- unique(df$gene_symbol)
for(g in setdiff(GENES, found_genes)) {
  if(!is.null(gene_ensembl[[g]])) {
    idx <- which(grepl(paste0("^", gene_ensembl[[g]], "\\."), df$ensembl_id))
    if(length(idx) > 0) df[idx, gene_symbol := g]
  }
}

# 转换为样本x基因矩阵
df_expr <- melt(df, id.vars = c("ensembl_id", "gene_symbol"), 
                variable.name = "sample_id", value.name = "tpm")
df_expr <- dcast(df_expr, sample_id ~ gene_symbol, value.var = "tpm")

# 确保所有基因都存在
for(g in GENES) if(!g %in% names(df_expr)) df_expr[, (g) := NA_real_]

cat("Expression matrix:", dim(df_expr), "\n")
fwrite(df_expr, file.path(OUTPUT_DIR, "raw_tpm_matrix.csv"))

# ============================================================
# Step 4: 合并组织注释
# ============================================================
cat("\nMerging with tissue annotation ...\n")

# 提取样本ID前缀
df_expr[, tissue_sample_id := gsub("^([^-]+-[^-]+-[^-]+).*$", "\\1", sample_id)]
df_expr$sample_id %>% duplicated(.) %>% sum(.)

# 合并
df_merged <- merge(df_expr, df_sample, by.x = "tissue_sample_id", by.y = "sampleId", all.x = TRUE)
df_merged <- df_merged[!is.na(tissue)]
cat("Samples with tissue annotation:", nrow(df_merged), "\n")
fwrite(df_merged, file.path(OUTPUT_DIR, "tpm_with_tissue.csv"))

# ============================================================
# Step 5: 计算共表达
# ============================================================
cat("\n--- Pan-tissue coexpression ---\n")

# 泛组织相关分析
gene_pairs <- combn(GENES, 2, simplify = FALSE)

# 使用批量计算（更快且无警告）
expr_matrix <- as.matrix(df_merged[, ..GENES])
cor_matrix <- cor(expr_matrix, method = "spearman", use = "pairwise.complete.obs")
n <- nrow(df_merged)
t_stat <- cor_matrix * sqrt((n - 2) / (1 - cor_matrix^2))
p_matrix <- 2 * pt(-abs(t_stat), df = n - 2)

# 转换为长格式
pan_results <- rbindlist(lapply(1:(ncol(cor_matrix)-1), function(i) {
  rbindlist(lapply((i+1):ncol(cor_matrix), function(j) {
    data.table(gene1 = rownames(cor_matrix)[i],
               gene2 = rownames(cor_matrix)[j],
               spearman_rho = cor_matrix[i, j],
               pvalue = p_matrix[i, j],
               n_samples = n)
  }))
}))

# BH FDR校正 (注意参数是"BH"而不是"fdr_bh")
pan_results[, adj_pvalue := p.adjust(pvalue, method = "BH")]
fwrite(pan_results, file.path(OUTPUT_DIR, "coexpression_pantissue.csv"))

# 打印结果
for (i in 1:nrow(pan_results)) {
  r <- pan_results[i]
  sig <- ifelse(r$adj_pvalue < 0.001, "***", 
                ifelse(r$adj_pvalue < 0.01, "**", 
                       ifelse(r$adj_pvalue < 0.05, "*", "ns")))
  cat(sprintf("  %s vs %s: rho=%.4f, adj_p=%.2e %s\n", 
              r$gene1, r$gene2, r$spearman_rho, r$adj_pvalue, sig))
}

# 组织分层相关分析
cat("\n--- Tissue-stratified coexpression ---\n")

# 先查看df_merged中有哪些组织以及样本数
cat("Checking tissues in df_merged:\n")
tissue_counts <- df_merged[, .N, by = tissue]
print(tissue_counts[order(-N)], nrows = 30)

# 使用for循环，注意变量命名避免冲突
tissue_results <- data.table()  # 初始化为空data.table

for (tiss in KEY_TISSUES) {
  # 使用get()方法提取当前组织的数据
  df_t <- df_merged[get("tissue") == tiss, ]
  
  cat("  Processing", tiss, ":", nrow(df_t), "samples\n")
  
  if (nrow(df_t) < 20) {
    cat("    Skipped (too few samples)\n")
    next
  }
  
  # 批量计算
  expr_t <- as.matrix(df_t[, ..GENES])
  cor_t <- cor(expr_t, method = "spearman", use = "pairwise.complete.obs")
  n_t <- nrow(df_t)
  t_stat_t <- cor_t * sqrt((n_t - 2) / (1 - cor_t^2))
  p_t <- 2 * pt(-abs(t_stat_t), df = n_t - 2)
  
  # 转换为长格式，直接rbind到tissue_results
  for (i in 1:(length(GENES)-1)) {
    for (j in (i+1):length(GENES)) {
      tissue_results <- rbind(tissue_results, data.table(
        tissue = tiss,
        gene1 = GENES[i],
        gene2 = GENES[j],
        spearman_rho = cor_t[i, j],
        pvalue = p_t[i, j],
        n_samples = n_t
      ))
    }
  }
}

# 查看结果
cat("\nTotal results:", nrow(tissue_results), "\n")
print(head(tissue_results, 20))

# BH FDR校正 (按组织分层)
tissue_results[, adj_pvalue := p.adjust(pvalue, method = "BH"), by = tissue]
fwrite(tissue_results, file.path(OUTPUT_DIR, "coexpression_tissue_stratified.csv"))
cat("  Tissues analyzed:", uniqueN(tissue_results$tissue), "\n")

# 筛选存在数据的组织
tissues_present <- intersect(KEY_TISSUES, unique(tissue_results$tissue))

# ============================================================
# Step 6: 可视化
# ============================================================
cat("\n--- Generating figures ---\n")

# 6.1 共表达热图
# 准备泛组织矩阵
pan_matrix <- matrix(1, nrow = length(GENES), ncol = length(GENES))
rownames(pan_matrix) <- colnames(pan_matrix) <- GENES
for (i in 1:nrow(pan_results)) {
  r <- pan_results[i]
  idx1 <- which(GENES == r$gene1)
  idx2 <- which(GENES == r$gene2)
  pan_matrix[idx1, idx2] <- pan_matrix[idx2, idx1] <- r$spearman_rho
}

# 准备组织分层矩阵
pairs_names <- sapply(gene_pairs, function(p) paste(p[1], "vs", p[2]))
tissue_matrix <- matrix(NA, nrow = length(tissues_present), ncol = length(pairs_names))
rownames(tissue_matrix) <- tissues_present
colnames(tissue_matrix) <- pairs_names
tissue_sig <- matrix(FALSE, nrow = length(tissues_present), ncol = length(pairs_names))

for (ti in seq_along(tissues_present)) {
  for (pi in seq_along(gene_pairs)) {
    pair <- gene_pairs[[pi]]
    row <- tissue_results[tissue == tissues_present[ti] & 
                            gene1 == pair[1] & gene2 == pair[2]]
    if (nrow(row) > 0) {
      tissue_matrix[ti, pi] <- row$spearman_rho
      tissue_sig[ti, pi] <- row$adj_pvalue < 0.05
    }
  }
}

# 保存热图
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# 泛组织热图
png(file.path(OUTPUT_DIR, "fig_coexpression_heatmap_pan.png"), width = 4, height = 4, res = 300, units = "in")
pheatmap(pan_matrix, 
         display_numbers = TRUE,
         number_format = "%.2f",
         color = colors,
         main = "Pan-tissue coexpression",
         fontsize = 10,
         fontsize_number = 8)
dev.off()

# 组织分层热图（带显著性标记）
tissue_matrix_label <- matrix("", nrow = nrow(tissue_matrix), ncol = ncol(tissue_matrix))
for (i in 1:nrow(tissue_matrix)) {
  for (j in 1:ncol(tissue_matrix)) {
    if (!is.na(tissue_matrix[i, j])) {
      tissue_matrix_label[i, j] <- sprintf("%.2f", tissue_matrix[i, j])
      if (tissue_sig[i, j]) {
        tissue_matrix_label[i, j] <- paste0(tissue_matrix_label[i, j], "*")
      }
    }
  }
}

png(file.path(OUTPUT_DIR, "fig_coexpression_heatmap_tissue.png"), width = 4, height = 6, res = 300, units = "in")
pheatmap(tissue_matrix,
         display_numbers = tissue_matrix_label,
         number_format = "%s",
         color = colors,
         main = "Tissue-stratified coexpression (Spearman rho)\n(* FDR < 0.05)",
         fontsize = 9,
         fontsize_number = 7,
         labels_row = substr(rownames(tissue_matrix), 1, 15),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         number_color = "white")
dev.off()

cat("  Saved: fig_coexpression_heatmap_pan.png and fig_coexpression_heatmap_tissue.png\n")

# 6.2 组织表达条形图
cat("\nGenerating tissue expression barplots ...\n")

# 计算组织中位TPM
tissue_median <- df_merged[, lapply(.SD, median, na.rm = TRUE), 
                           by = tissue, .SDcols = GENES]
fwrite(tissue_median, file.path(OUTPUT_DIR, "gtex_median_tpm_by_tissue.csv"))

# 只保留目标组织
tissue_median_filtered <- tissue_median[tissue %in% tissues_present]
tissue_median_long <- melt(tissue_median_filtered, 
                           id.vars = "tissue", 
                           measure.vars = GENES,
                           variable.name = "gene", 
                           value.name = "median_tpm")


# 创建条形图
p <- ggplot(tissue_median_long, aes(x = reorder(tissue, -median_tpm), y = median_tpm, fill = gene)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), color = "black", linewidth = 0.2) +
  scale_fill_npg() +
  labs(title = "GTEx Tissue Median Expression (TPM)",
       x = "Tissue",
       y = "Median TPM",
       fill = "Gene") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank()) +
  facet_wrap(~ gene, scales = "free_y", ncol = 3)

ggsave(file.path(OUTPUT_DIR, "fig_tissue_expression.png"), p, width = 10, height = 5, dpi = 300)
cat("  Saved: fig_tissue_expression.png\n")

# ============================================================
# 完成
# ============================================================
cat("\n", "=", rep("=", 60), "\n", sep="")
cat("Analysis complete! All outputs saved to:", OUTPUT_DIR, "\n")
cat("=", rep("=", 60), "\n", sep="")

