# Author: Liu Zhehan
# Time: 2026.9.13
# 功能：GO/KEGG 富集 + 基因共现网络图（点=基因，边=共享功能，颜色=功能）

#### 设置工作空间 ####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/GO_KEGG")
rm(list = ls())
gc()

#### 导入R包 ####
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(KEGGREST)
library(networkD3)
library(tidyr)
library(enrichplot)
library(igraph)

#### 参数设置 ####
INPUT_FILE <- "data/example/final_homogroup_results.csv"
CID        <- "Conserved_079"
OUTPUT_DIR <- "result"
TOP_N      <- 10        # 取前多少个功能（控制网络复杂度）
MIN_SHARED <- 1         # 两基因至少共享几个功能才连边

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#### 读入数据 ####
data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
GENES <- subset(data, Conserved.ID == CID)$`Human.Prot..by.Comprehensive.`
GENES <- GENES[GENES != "" & !is.na(GENES)]
GENES <- unique(GENES)

if (length(GENES) <= 1) stop("基因数小于2")

cat("输入基因数：", length(GENES), "\n")

#### 基因ID转换 ####
id <- tryCatch({
  bitr(GENES, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}, error = function(e) return(NA))

if (length(id) == 1 && is.na(id)) stop("基因转换失败")
if (nrow(id) <= 1) stop("基因转换后数量小于2")

cat("成功转换基因数：", nrow(id), "\n")

#### ============ GO 富集 ============ ####
ego <- enrichGO(gene          = id$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# 若富集为空，放宽阈值
if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
  warning("GO 严格阈值下无结果，放宽阈值")
  ego <- enrichGO(gene = id$ENTREZID, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "BP",
                  pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
}
cat("GO 富集条目数：", nrow(as.data.frame(ego)), "\n")

#### ============ KEGG 富集 ============ ####
KEGG_DIR <- "data/KEGG"
if (!dir.exists(KEGG_DIR)) dir.create(KEGG_DIR, recursive = TRUE)

kegg_link_file  <- file.path(KEGG_DIR, "hsa_kegg.txt")
kegg_name_file  <- file.path(KEGG_DIR, "hsa_pathway_names.txt")

if (!file.exists(kegg_link_file)) {
  download.file("https://rest.kegg.jp/link/hsa/pathway", kegg_link_file)
}
if (!file.exists(kegg_name_file)) {
  download.file("https://rest.kegg.jp/list/pathway/hsa", kegg_name_file)
}

# 读取 KEGG 注释
kegg_link <- read.table(kegg_link_file, sep = "\t", stringsAsFactors = FALSE)
colnames(kegg_link) <- c("Pathway_ID", "ENTREZID")
kegg_link$ENTREZID  <- gsub("hsa:", "", kegg_link$ENTREZID)
kegg_link$Pathway_ID <- gsub("path:", "", kegg_link$Pathway_ID)

pathway_names <- read.table(kegg_name_file, sep = "\t", stringsAsFactors = FALSE,
                            quote = "", fill = TRUE)
colnames(pathway_names) <- c("Pathway_ID", "Pathway_Name")
pathway_names$Pathway_Name <- gsub(" - Homo sapiens \\(human\\)", "",
                                   pathway_names$Pathway_Name)

kegg_annotation <- merge(kegg_link, pathway_names, by = "Pathway_ID")

# KEGG 富集
kk <- enricher(gene          = id$ENTREZID,
               TERM2GENE    = kegg_annotation[, c("Pathway_ID", "ENTREZID")],
               TERM2NAME    = kegg_annotation[, c("Pathway_ID", "Pathway_Name")],
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05)

if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
  warning("KEGG 严格阈值下无结果，放宽阈值")
  kk <- enricher(gene = id$ENTREZID,
                 TERM2GENE = kegg_annotation[, c("Pathway_ID", "ENTREZID")],
                 TERM2NAME = kegg_annotation[, c("Pathway_ID", "Pathway_Name")],
                 pvalueCutoff = 1, qvalueCutoff = 1)
}

kk_df <- as.data.frame(kk)

# 构建 ENTREZID -> SYMBOL 的映射表
id_map <- setNames(id$SYMBOL, id$ENTREZID)

# 把 geneID 列（"123/456/789" 格式）拆开逐个替换
replace_ids <- function(gene_str, map) {
  ids <- strsplit(gene_str, "/")[[1]]
  syms <- map[ids]
  syms[is.na(syms)] <- ids[is.na(syms)]   # 找不到的保留原ID
  paste(syms, collapse = "/")
}

kk_df$geneID <- sapply(kk_df$geneID, replace_ids, map = id_map)

# 写回 kk 对象（enricher 结果本质是 data.frame 子类）
kk@result <- kk_df

cat("KEGG 富集条目数：", nrow(as.data.frame(kk)), "\n")

#### ============ 通用：基因共现网络绘制函数 ============ ####
plot_gene_network <- function(enrich_res, title, out_file,
                              top_n = 10, min_shared = 1,
                              color_mode = c("main", "pie")) {
  color_mode <- match.arg(color_mode)
  
  df <- as.data.frame(enrich_res)
  if (is.null(df) || nrow(df) == 0) {
    warning("富集结果为空，跳过：", title)
    return(invisible(NULL))
  }
  df <- df[order(df$p.adjust), , drop = FALSE]
  df <- df[1:min(top_n, nrow(df)), , drop = FALSE]
  
  # 基因-功能长表
  gene_term <- do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    genes <- strsplit(df$geneID[i], "/")[[1]]
    genes <- genes[genes != ""]
    if (length(genes) == 0) return(NULL)
    data.frame(gene = genes,
               term = df$Description[i],
               rank = i,
               stringsAsFactors = FALSE)
  }))
  
  if (is.null(gene_term) || nrow(gene_term) == 0) {
    warning("无基因-功能映射，跳过：", title)
    return(invisible(NULL))
  }
  
  gene_list <- unique(gene_term$gene)
  terms     <- unique(gene_term$term)
  term_col  <- setNames(rainbow(length(terms)), terms)
  
  cat("  [", title, "] 基因数:", length(gene_list),
      " 功能数:", length(terms), "\n")
  
  # 构建边：两基因共享功能数 >= min_shared
  edges <- data.frame()
  if (length(gene_list) >= 2) {
    for (i in 1:(length(gene_list) - 1)) {
      for (j in (i + 1):length(gene_list)) {
        g1 <- gene_list[i]; g2 <- gene_list[j]
        t1 <- gene_term$term[gene_term$gene == g1]
        t2 <- gene_term$term[gene_term$gene == g2]
        shared <- intersect(t1, t2)
        if (length(shared) >= min_shared) {
          edges <- rbind(edges, data.frame(from = g1, to = g2,
                                           weight = length(shared),
                                           stringsAsFactors = FALSE))
        }
      }
    }
  }
  
  if (nrow(edges) == 0) {
    cat("  该网络无共享功能边，退化为散点节点图\n")
    edges <- data.frame(from = character(), to = character(),
                        weight = numeric())
  }
  
  g <- graph_from_data_frame(edges,
                             vertices = data.frame(name = gene_list,
                                                   stringsAsFactors = FALSE),
                             directed = FALSE)
  
  # 节点着色
  if (color_mode == "main") {
    # 每个基因取 rank 最小的功能作为主色
    main_term <- gene_term %>%
      group_by(gene) %>%
      slice_min(rank, n = 1, with_ties = FALSE) %>%
      ungroup()
    col_vec <- term_col[main_term$term[match(V(g)$name, main_term$gene)]]
    col_vec[is.na(col_vec)] <- "grey70"
    V(g)$color <- col_vec
    V(g)$shape <- "circle"
    V(g)$size  <- 10
  } else {
    # 饼图：每个基因各功能占比
    pie_list <- lapply(V(g)$name, function(gene) {
      tb <- table(gene_term$term[gene_term$gene == gene])
      as.numeric(tb)
    })
    pie_col_list <- lapply(V(g)$name, function(gene) {
      tb <- table(gene_term$term[gene_term$gene == gene])
      unname(term_col[names(tb)])
    })
    names(pie_list)     <- V(g)$name
    names(pie_col_list) <- V(g)$name
    V(g)$pie       <- pie_list
    V(g)$pie.color <- pie_col_list
    V(g)$shape <- "pie"
    V(g)$size  <- 15
  }
  
  # 绘图
  pdf(out_file, width = 12, height = 10)
  set.seed(123)
  layout_g <- layout_with_fr(g)
  
  if (color_mode == "main") {
    plot(g,
         layout             = layout_g,
         vertex.frame.color = "white",
         vertex.label.cex   = 0.7,
         vertex.label.dist  = 1.0,
         vertex.label.color = "black",
         edge.width         = E(g)$weight,
         edge.color         = "grey75",
         main               = title)
    legend("topleft",
           legend = terms,
           fill   = term_col,
           bty    = "n", cex = 0.7,
           title  = "Function (main)")
  } else {
    plot(g,
         layout             = layout_g,
         vertex.frame.color = "white",
         vertex.label.cex   = 0.7,
         vertex.label.dist  = 1.5,
         vertex.label.color = "black",
         edge.width         = E(g)$weight,
         edge.color         = "grey75",
         main               = title)
    legend("topleft",
           legend = terms,
           fill   = term_col,
           bty    = "n", cex = 0.7,
           title  = "Function (pie)")
  }
  dev.off()
  cat("  已保存：", out_file, "\n")
  
  # 同时输出边表和节点表，方便后续检查
  return(list(graph = g, edges = edges, gene_term = gene_term))
}

#### ============ 绘制 GO 网络图 ============ ####
cat("\n===== 绘制 GO 网络图 =====\n")
res_go_main <- plot_gene_network(
  ego,
  title      = paste0("GO-BP Gene Co-occurrence Network: ", CID),
  out_file   = file.path(OUTPUT_DIR, paste0(CID, "_GO_gene_net_main.pdf")),
  top_n      = TOP_N,
  min_shared = MIN_SHARED,
  color_mode = "main"
)

res_go_pie <- plot_gene_network(
  ego,
  title      = paste0("GO-BP Gene Co-occurrence Network (Pie): ", CID),
  out_file   = file.path(OUTPUT_DIR, paste0(CID, "_GO_gene_net_pie.pdf")),
  top_n      = TOP_N,
  min_shared = MIN_SHARED,
  color_mode = "pie"
)

#### ============ 绘制 KEGG 网络图 ============ ####
cat("\n===== 绘制 KEGG 网络图 =====\n")
res_kegg_main <- plot_gene_network(
  kk,
  title      = paste0("KEGG Gene Co-occurrence Network: ", CID),
  out_file   = file.path(OUTPUT_DIR, paste0(CID, "_KEGG_gene_net_main.pdf")),
  top_n      = TOP_N,
  min_shared = MIN_SHARED,
  color_mode = "main"
)

res_kegg_pie <- plot_gene_network(
  kk,
  title      = paste0("KEGG Gene Co-occurrence Network (Pie): ", CID),
  out_file   = file.path(OUTPUT_DIR, paste0(CID, "_KEGG_gene_net_pie.pdf")),
  top_n      = TOP_N,
  min_shared = MIN_SHARED,
  color_mode = "pie"
)

#### ============ 输出边表/节点表（可选） ============ ####
if (!is.null(res_go_main)) {
  write.csv(res_go_main$edges,
            file.path(OUTPUT_DIR, paste0(CID, "_GO_edges.csv")),
            row.names = FALSE)
}
if (!is.null(res_kegg_main)) {
  write.csv(res_kegg_main$edges,
            file.path(OUTPUT_DIR, paste0(CID, "_KEGG_edges.csv")),
            row.names = FALSE)
}

cat("\n全部完成！结果保存在：", normalizePath(OUTPUT_DIR), "\n")
