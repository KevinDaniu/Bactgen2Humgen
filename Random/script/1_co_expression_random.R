####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Random/")
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
library(clusterProfiler)
library(org.Hs.eg.db)

####共表达分析####
#----参数设置----#
# 设置参数
INPUT_FILE <- "data/final_homogroup_results.csv"
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
gene_ensembl <- list()
for (g in GENES) {
  url <- paste0("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/", g, "?content-type=application/json")
  response <- GET(url)
  if (status_code(response) == 200) {
    data <- httr::content(response, "parsed")
    for (item in data) {
      if (!is.null(item$type) && item$type == "gene") {
        gene_ensembl[[g]] <- item$id
        break
      }
    }
  }
}
cat("Gene -> Ensembl ID:", paste(names(gene_ensembl), unlist(gene_ensembl), sep="=", collapse=", "), "\n")

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
# Step 5.5: 随机抽样对照分析
# ============================================================
cat("\n--- Random sampling control analysis ---\n")

set.seed(123)  # 可复现
N_RANDOM <- 1000          # 随机抽样次数（建议 1000，调试时可先用 100）
MIN_SAMPLES_PER_GENE <- 0.8  # 基因在样本中的最低表达比例（可选过滤）

# ---- 5.5.1 获取蛋白编码基因列表 ----
all_human_genes <- keys(org.Hs.eg.db, keytype = "ENTREZID")
gene_info <- select(org.Hs.eg.db,
                    keys = all_human_genes,
                    columns = c("SYMBOL", "GENETYPE"),
                    keytype = "ENTREZID")
protein_coding_genes <- gene_info %>%
  dplyr::filter(GENETYPE == "protein-coding") %>%
  pull(SYMBOL) %>%
  unique()
protein_coding_genes <- protein_coding_genes[!is.na(protein_coding_genes)]
cat("Protein-coding genes:", length(protein_coding_genes), "\n")

# ---- 5.5.2 从 GTEx 表达矩阵中筛选候选基因池 ----
# 只保留在 GTEx 中实际检测到的基因（避免抽到未表达基因导致 rho 全为 NA）
expr_pool <- fread(file.path("data/GTEx_Data", GTEX_TPM_LOCAL), 
                   skip = 3, header = FALSE, nThread = 4)

pool_symbols <- expr_pool$V2
pool_symbols <- pool_symbols[pool_symbols %in% protein_coding_genes]

# 去重
pool_symbols <- unique(pool_symbols)
cat("Genes available in GTEx & protein-coding:", length(pool_symbols), "\n")

# ---- 5.5.3 计算目标簇的共表达指标 ----
compute_pair_metrics <- function(expr_mat) {
  # expr_mat: samples x genes
  if (ncol(expr_mat) < 2) return(NULL)
  cor_m <- cor(expr_mat, method = "spearman", use = "pairwise.complete.obs")
  n <- nrow(expr_mat)
  t_stat <- cor_m * sqrt((n - 2) / (1 - cor_m^2))
  p_m <- 2 * pt(-abs(t_stat), df = n - 2)
  
  # 只取上三角
  idx <- upper.tri(cor_m)
  rho_vec <- cor_m[idx]
  p_vec <- p_m[idx]
  adj_p <- p.adjust(p_vec, method = "BH")
  
  list(
    mean_rho      = mean(rho_vec, na.rm = TRUE),
    mean_abs_rho  = mean(abs(rho_vec), na.rm = TRUE),
    max_abs_rho   = max(abs(rho_vec), na.rm = TRUE),
    prop_sig      = mean(adj_p < 0.05, na.rm = TRUE),
    n_pairs       = length(rho_vec),
    rho_vec       = rho_vec,
    adj_p         = adj_p
  )
}

# 目标簇指标
expr_target <- as.matrix(df_merged[, ..GENES])
target_metrics <- compute_pair_metrics(expr_target)
cat(sprintf("Target cluster (%d genes, %d pairs):\n", length(GENES), target_metrics$n_pairs))
cat(sprintf("  mean_rho=%.4f, mean_abs_rho=%.4f, max_abs_rho=%.4f, prop_sig=%.3f\n",
            target_metrics$mean_rho, target_metrics$mean_abs_rho,
            target_metrics$max_abs_rho, target_metrics$prop_sig))

# ---- 5.5.4 随机抽样 ----
k <- length(GENES)  # 每次抽样的基因数
random_metrics <- data.table(
  iter          = integer(),
  mean_rho      = numeric(),
  mean_abs_rho  = numeric(),
  max_abs_rho   = numeric(),
  prop_sig      = numeric()
)

# 预提取完整表达矩阵（一次性，避免循环中反复读取）
cat("Extracting full expression matrix for random sampling...\n")
# 注意：expr_pool 的行是基因，需要按 SYMBOL 取子集
expr_pool_mat <- expr_pool[expr_pool$V2 %in% pool_symbols, ]
setnames(expr_pool_mat, c("ensembl_id", "gene_symbol", sample_ids))
expr_pool_mat <- expr_pool_mat[!duplicated(gene_symbol)]
expr_mat_full <- as.matrix(expr_pool_mat[, ..sample_ids])
rownames(expr_mat_full) <- expr_pool_mat$gene_symbol

cat("Full expression matrix:", dim(expr_mat_full), "\n")

for (i in 1:N_RANDOM) {
  if (i %% 100 == 0) cat("  Random iteration", i, "/", N_RANDOM, "\n")
  
  random_genes <- sample(pool_symbols, size = k, replace = FALSE)
  # 确保都在矩阵中
  random_genes <- intersect(random_genes, rownames(expr_mat_full))
  if (length(random_genes) < k) next
  
  expr_rand <- t(expr_mat_full[random_genes, , drop = FALSE])
  m <- compute_pair_metrics(expr_rand)
  if (is.null(m)) next
  
  random_metrics <- rbind(random_metrics, data.table(
    iter         = i,
    mean_rho     = m$mean_rho,
    mean_abs_rho = m$mean_abs_rho,
    max_abs_rho  = m$max_abs_rho,
    prop_sig     = m$prop_sig
  ))
}

cat("Completed", nrow(random_metrics), "random iterations\n")
fwrite(random_metrics, file.path(OUTPUT_DIR, "random_sampling_metrics.csv"))

# ---- 5.5.5 经验 P 值 ----
empirical_p <- function(obs, null_vec, alternative = "greater") {
  if (alternative == "greater") {
    (sum(null_vec >= obs, na.rm = TRUE) + 1) / (length(null_vec) + 1)
  } else {
    (sum(null_vec <= obs, na.rm = TRUE) + 1) / (length(null_vec) + 1)
  }
}

emp_results <- data.table(
  metric        = c("mean_rho", "mean_abs_rho", "max_abs_rho", "prop_sig"),
  observed      = c(target_metrics$mean_rho, target_metrics$mean_abs_rho,
                    target_metrics$max_abs_rho, target_metrics$prop_sig),
  null_mean     = c(mean(random_metrics$mean_rho), mean(random_metrics$mean_abs_rho),
                    mean(random_metrics$max_abs_rho), mean(random_metrics$prop_sig)),
  null_sd       = c(sd(random_metrics$mean_rho), sd(random_metrics$mean_abs_rho),
                    sd(random_metrics$max_abs_rho), sd(random_metrics$prop_sig)),
  empirical_p   = c(empirical_p(target_metrics$mean_rho, random_metrics$mean_rho),
                    empirical_p(target_metrics$mean_abs_rho, random_metrics$mean_abs_rho),
                    empirical_p(target_metrics$max_abs_rho, random_metrics$max_abs_rho),
                    empirical_p(target_metrics$prop_sig, random_metrics$prop_sig))
)
emp_results[, z_score := (observed - null_mean) / null_sd]
emp_results[, fold_enrichment := observed / null_mean]

cat("\n=== Empirical P-value Results ===\n")
print(emp_results)
fwrite(emp_results, file.path(OUTPUT_DIR, "random_sampling_empirical_p.csv"))

# ---- 5.5.6 可视化：零分布 + 观测值 ----
plot_null_dist <- function(null_vec, obs, metric_name, out_file) {
  df_plot <- data.table(value = null_vec)
  p <- ggplot(df_plot, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                   fill = "grey70", color = "white") +
    geom_density(color = "steelblue", linewidth = 0.8) +
    geom_vline(xintercept = obs, color = "red", linewidth = 1, linetype = "dashed") +
    annotate("text", x = obs, y = Inf, label = sprintf("Observed = %.3f", obs),
             vjust = 1.5, hjust = -0.05, color = "red", size = 4) +
    labs(title = paste0("Random sampling null distribution: ", metric_name),
         subtitle = sprintf("N = %d random gene sets (k = %d genes)", 
                            length(null_vec), k),
         x = metric_name, y = "Density") +
    theme_bw(base_size = 12)
  ggsave(paste0(out_file,".png"), p, width = 6, height = 4, dpi = 300)
  ggsave(paste0(out_file,".pdf"), p, width = 6, height = 4, dpi = 300)
}

plot_null_dist(random_metrics$mean_rho, target_metrics$mean_rho,
               "mean_rho", file.path(OUTPUT_DIR, "fig_null_mean_rho"))
plot_null_dist(random_metrics$mean_abs_rho, target_metrics$mean_abs_rho,
               "mean_abs_rho", file.path(OUTPUT_DIR, "fig_null_mean_abs_rho"))
plot_null_dist(random_metrics$max_abs_rho, target_metrics$max_abs_rho,
               "max_abs_rho", file.path(OUTPUT_DIR, "fig_null_max_abs_rho"))
plot_null_dist(random_metrics$prop_sig, target_metrics$prop_sig,
               "prop_sig", file.path(OUTPUT_DIR, "fig_null_prop_sig"))

cat("  Saved 4 null distribution figures\n")

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
         number_color = "white",
         number_format = "%.2f",
         color = colors,
         main = "Pan-tissue coexpression",
         fontsize = 10,
         fontsize_number = 8)
dev.off()
pdf(file.path(OUTPUT_DIR, "fig_coexpression_heatmap_pan.pdf"), width = 4, height = 4)
pheatmap(pan_matrix, 
         display_numbers = TRUE,
         number_color = "white",
         number_format = "%.2f",
         color = colors,
         main = "Pan-tissue coexpression",
         fontsize = 10,
         fontsize_number = 8)
dev.off()

# 组织分层热图（带显著性标记）
# tissue_matrix <- t(tissue_matrix)
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
         main = "Tissue-stratified coexpression",
         fontsize = 9,
         fontsize_number = 7,
         labels_row = substr(rownames(tissue_matrix), 1, 15),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         number_color = "white")
dev.off()

pdf(file.path(OUTPUT_DIR, "fig_coexpression_heatmap_tissue.pdf"), width = 4, height = 6)
pheatmap(tissue_matrix,
         display_numbers = tissue_matrix_label,
         number_format = "%s",
         color = colors,
         main = "Tissue-stratified coexpression",
         fontsize = 9,
         fontsize_number = 7,
         labels_row = substr(rownames(tissue_matrix), 1, 15),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         number_color = "white")
dev.off()

cat("  Saved: fig_coexpression_heatmap_pan.png and fig_coexpression_heatmap_tissue.png\n")
