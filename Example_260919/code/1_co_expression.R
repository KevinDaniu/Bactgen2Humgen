####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Example_260919/")
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
library(EnsDb.Hsapiens.v86)

####共表达分析####
#----参数设置----#
# 设置参数
INPUT_FILE <- "data/final_homogroup_results.csv"
data <- read.csv(INPUT_FILE)
CID <- "Conserved_025"
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
    data <- httr::content(response, "parsed")$data
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
#     data <- httr::content(response, "parsed")
#     for (item in data) {
#       if (!is.null(item$type) && item$type == "gene") {
#         gene_ensembl[[g]] <- item$id
#         break
#       }
#     }
#   }
# }
# cat("Gene -> Ensembl ID:", paste(names(gene_ensembl), unlist(gene_ensembl), sep="=", collapse=", "), "\n")

# GENES 是你的基因符号向量
gene_ensembl <- mapIds(
  EnsDb.Hsapiens.v86,
  keys = GENES,
  column = "GENEID",
  keytype = "SYMBOL",
  multiVals = "first"
)

gene_ensembl <- as.list(gene_ensembl)
names(gene_ensembl)[1] <- "AC004832.3"
gene_ensembl[[1]] <- "ENSG00000249590"
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
