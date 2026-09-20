# Author: Liu Zhehan
# 2026.9.13
###############################################################################
# STRING PPI Network & Hub Gene Identification
# - Query STRING local database for protein-protein interactions of candidate genes
# - Build PPI network with igraph
# - Compute topological centrality: MCC, Degree, Betweenness, Closeness, Eigenvector
# - Identify hub genes
###############################################################################

##### 设置工作空间 ####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/PPI/")
rm(list = ls())
gc()

#### 导入R包 ####
library(STRINGdb)
library(igraph)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggraph)
library(patchwork)
library(data.table)
library(biomaRt)

#### Load candidate genes ####
INPUT_FILE <- "data/example/final_homogroup_results.csv"
CID        <- "Conserved_079"
OUTPUT_DIR <- "result"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#### 读入数据 ####
data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
GENES <- subset(data, Conserved.ID == CID)$`Human.Prot..by.Comprehensive.`
GENES <- GENES[GENES != "" & !is.na(GENES)]
GENES <- unique(GENES)

if (length(GENES) <= 1) stop("基因数小于2")

cat("输入基因数：", length(GENES), "\n")

#### Initialize STRINGdb ####
local_db_dir <- "data/stringdb/"

# 初始化 STRINGdb，指向本地目录
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = 400,
  input_directory = local_db_dir  # 关键参数：指定本地文件目录
)

##### Get interactions among candidate genes ####
# Map gene symbols to STRING IDs
genes_mapped <- string_db$map(data.frame(gene = GENES), 
                              "gene", 
                              removeUnmappedRows = TRUE)
candidate_string_ids <- genes_mapped$STRING_id

# Get interactions among candidate genes
interactions_candidate <- string_db$get_interactions(candidate_string_ids)

cat("STRING API returned", nrow(interactions_candidate), 
    "interactions among candidate genes\n")

# if(nrow(interactions_candidate) < 1) stop("没有互作关系")

#### Expand network with first-degree interactors ####
cat("\nExpanding network with first-degree interactors...\n")

all_interactors <- unique(unlist(lapply(candidate_string_ids, function(gene_id) {
  # 获取邻居
  neighbors <- string_db$get_neighbors(gene_id)
  if (length(neighbors) == 0) return(character(0))
  
  cat("  Gene", gene_id, "has", length(neighbors), "neighbors\n")
  
  # 因为get_interactions可能不支持大量查询
  top_neighbors <- neighbors[1:min(10, length(neighbors))]
  
  # 返回这些邻居的STRING ID
  return(top_neighbors)
})))

all_interactors <- setdiff(all_interactors, candidate_string_ids)
cat("Found", length(all_interactors), "unique interactors\n")

# 扩展网络
expanded_ids <- unique(c(candidate_string_ids, all_interactors))
cat("Total nodes after expansion:", length(expanded_ids), "\n")

# 获取扩展网络中所有基因之间的互作
cat("Fetching all interactions in expanded network...\n")
expanded_df <- data.frame(STRING_id = expanded_ids)
interactions_expanded <- string_db$get_interactions(expanded_df %>% unlist(.))

# 如果为空，用已有数据
if (is.null(interactions_expanded) || nrow(interactions_expanded) == 0) {
  cat("Using existing interactions_candidate data...\n")
  interactions_expanded <- interactions_candidate
}

if ("combined_score" %in% names(interactions_expanded)) {
  setnames(interactions_expanded, "combined_score", "score")
}

# 过滤score >= 400
if (nrow(interactions_expanded) > 0) {
  interactions_expanded <- interactions_expanded[interactions_expanded$score >= 400, ]
  cat("Total edges after expansion:", nrow(interactions_expanded), "\n")
} else {
  cat("Warning: No interactions found!\n")
}

#### Build igraph network ####
# Create edge list
edge_list <- data.frame(
  from = interactions_expanded$from,
  to = interactions_expanded$to,
  weight = interactions_expanded$score
)

# Create vertex data frame with labels
# Map back to gene symbols
all_ids <- unique(c(edge_list$from, edge_list$to))
ensp_ids <- sub("9606\\.", "", all_ids)
# 连接 Ensembl
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
# 获取基因符号
gene_symbols <- getBM(
  attributes = c("ensembl_peptide_id", "external_gene_name"),
  filters = "ensembl_peptide_id",
  values = ensp_ids,
  mart = ensembl
)
# 创建顶点数据框
vertices <- data.frame(
  name = all_ids,
  gene = gene_symbols$external_gene_name[match(ensp_ids, gene_symbols$ensembl_peptide_id)],
  stringsAsFactors = FALSE
)
# 如果某些没有匹配到，使用 ENSP 编号
vertices$gene[is.na(vertices$gene)] <- ensp_ids[is.na(vertices$gene)]

# 创建图
G <- graph_from_data_frame(edge_list, 
                           directed = FALSE, 
                           vertices = vertices)

cat("\nExpanded network:", vcount(G), "nodes,", ecount(G), "edges\n")


# Identify which nodes are candidate genes
V(G)$is_candidate <- V(G)$name %in% candidate_string_ids

# 直接用 gene 属性，不创建新属性
V(G)$gene_symbol <- V(G)$gene
V(G)$gene_symbol[!V(G)$is_candidate] <- NA

# Check candidate genes in network
candidate_in_net <- V(G)$name[V(G)$is_candidate]
cat("Candidate genes in network:", 
    paste(V(G)$gene_symbol[V(G)$is_candidate], collapse = ", "), "\n")

#### Compute topological metrics ####
# 1. Degree (度数)
degree_values <- igraph::degree(G)
degree_centrality <- igraph::centr_degree(G)$res

# 2. Betweenness centrality (介数中心性)
betweenness <- igraph::betweenness(G, weights = E(G)$weight)

# 3. Closeness centrality (紧密中心性)
closeness <- igraph::closeness(G, weights = E(G)$weight)

# 4. Eigenvector centrality (特征向量中心性)
eigenvector <- igraph::eigen_centrality(G, weights = E(G)$weight)$vector

# 5. Coreness (K-core 分解) - 比 MCC 更适合小网络
coreness <- igraph::coreness(G)

# 6. PageRank - 另一种衡量节点重要性的方法
pagerank <- igraph::page_rank(G, weights = E(G)$weight)$vector

# 7. MCC - 使用改进的版本，对小网络更友好
compute_mcc_improved <- function(graph) {
  # 获取所有团（包括大小为2的边）
  cliques <- igraph::cliques(graph, min = 2)
  
  if (length(cliques) == 0) {
    cat("没有找到任何团，使用 Degree * Closeness 作为 MCC 近似\n")
    return(igraph::degree(graph) * igraph::closeness(graph, weights = E(graph)$weight))
  }
  
  # 计算 MCC
  mcc <- rep(0, igraph::vcount(graph))
  names(mcc) <- igraph::V(graph)$name
  
  for (clique in cliques) {
    clique_size <- length(clique)
    if (clique_size >= 2) {
      # 贡献: (|clique| - 1)!
      contribution <- factorial(clique_size - 1)
      for (v in clique) {
        mcc[names(mcc) == v] <- mcc[names(mcc) == v] + contribution
      }
    }
  }
  
  # 如果没有找到团，使用近似
  if (sum(mcc) == 0) {
    cat("MCC 计算为 0，使用 Degree * Closeness 作为近似\n")
    mcc <- igraph::degree(graph) * igraph::closeness(graph, weights = E(graph)$weight)
  }
  
  return(mcc)
}

mcc_values <- compute_mcc_improved(G)

#### Compile results ####
topo_df <- data.frame(
  Gene = V(G)$gene_symbol,
  STRING_ID = V(G)$name,
  Degree = degree_values,
  DegreeCentrality = degree_centrality,
  MCC = mcc_values,
  Betweenness = betweenness,
  Closeness = closeness,
  Eigenvector = eigenvector,
  IsCandidate = V(G)$is_candidate,
  stringsAsFactors = FALSE
)

# Sort by MCC
topo_df <- topo_df[order(topo_df$MCC, decreasing = TRUE), ]

cat("\nTopological metrics (top 20):\n")
print(head(topo_df, 20))

write.csv(topo_df, file.path(OUTPUT_DIR, "PPI_topological_metrics.csv"), 
          row.names = FALSE)

#### Identify hub genes ####
# Hub genes: top-ranked by MCC (cytoHubba recommended method)
top_n <- min(10, nrow(topo_df))
hub_mcc_all <- topo_df[1:top_n, "Gene"]

# Hub genes among candidates: rank candidate genes by MCC
topo_candidate <- topo_df[topo_df$IsCandidate == TRUE, ]
topo_candidate <- topo_candidate[order(topo_candidate$MCC, decreasing = TRUE), ]

cat("\nCandidate gene rankings (by MCC):\n")
print(topo_candidate[, c("Gene", "Degree", "MCC", "Betweenness", "Closeness", "Eigenvector")])

# Hub genes = candidate genes that rank in top 10 of the full network by MCC
hub_genes <- intersect(topo_candidate$Gene[1:min(10, nrow(topo_candidate))], 
                       hub_mcc_all)

if (length(hub_genes) == 0) {
  # If no candidate is in top 10, take top 2 candidates by MCC
  hub_genes <- topo_candidate$Gene[1:min(2, nrow(topo_candidate))]
}

cat("\nHub genes (from candidates, by MCC ranking):", 
    paste(hub_genes, collapse = ", "), "\n")
cat("\nTop 10 network hubs (including interactors):", 
    paste(hub_mcc_all, collapse = ", "), "\n")

# Save hub genes
hub_df <- data.frame(Gene = hub_genes)
write.csv(hub_df, file.path(OUTPUT_DIR, "hub_genes.csv"), row.names = FALSE)

#### Save network for Cytoscape ####
# Save edge list
write.csv(edge_list, file.path(OUTPUT_DIR, "PPI_edges.csv"), row.names = FALSE)

# Save node attributes
node_attrs <- data.frame(
  id = V(G)$name,
  gene_symbol = V(G)$gene_symbol,
  degree = degree_values,
  mcc = mcc_values,
  is_candidate = V(G)$is_candidate
)
write.csv(node_attrs, file.path(OUTPUT_DIR, "PPI_nodes.csv"), row.names = FALSE)

# ---- Create PPI network visualization with ggplot2 ----
# Create layout
set.seed(42)
layout <- create_layout(G, layout = "fr", weights = E(G)$weight)

# Add node attributes to layout
layout$is_candidate <- V(G)$is_candidate
layout$gene_symbol <- V(G)$gene_symbol
layout$mcc <- mcc_values
layout$degree <- degree_values

# Network plot
p_network <- ggraph(layout) +
  geom_edge_link(aes(width = weight/1000), 
                 alpha = 0.3, 
                 color = "grey50") +
  geom_node_point(aes(size = mcc, 
                      color = is_candidate), 
                  alpha = 0.8) +
  geom_node_text(aes(label = ifelse(is_candidate | degree >= 3, 
                                    gene_symbol, 
                                    "")),
                 repel = TRUE,
                 size = 3,
                 fontface = "bold") +
  scale_color_manual(values = c("TRUE" = "#FF6B6B", 
                                "FALSE" = "#4ECDC4"),
                     labels = c("TRUE" = "Candidate genes", 
                                "FALSE" = "Interactors")) +
  scale_size_continuous(range = c(2, 10)) +
  scale_edge_width_continuous(range = c(0.2, 1.5)) +
  theme_void() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 14)) +
  guides(color = guide_legend(title = NULL),
         size = guide_legend(title = "MCC"),
         edge_width = guide_legend(title = "STRING score")) +
  labs(title = "PPI Network of Candidate Genes (STRING)\nNode size proportional to MCC")

ggsave(file.path(OUTPUT_DIR, "ppi_network.png"), 
       p_network, width = 6, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "ppi_network.pdf"), 
       p_network, width = 6, height = 5, dpi = 300)

cat("PPI network figure saved.\n")

# ---- Create centrality comparison bar plot (candidate genes only) ----
# Filter to candidate genes
topo_candidate_plot <- topo_df[topo_df$IsCandidate == TRUE, ]
topo_candidate_plot <- topo_candidate_plot[order(topo_candidate_plot$MCC, decreasing = TRUE), ]

# Melt data for faceted plotting
metrics <- c("Degree", "MCC", "Betweenness", "Closeness", "Eigenvector", "DegreeCentrality")

# Create individual bar plots
plot_list <- list()
for (i in seq_along(metrics)) {
  metric <- metrics[i]
  metric_df <- topo_candidate_plot
  metric_df$is_hub <- metric_df$Gene %in% hub_genes
  
  p <- ggplot(metric_df, aes(x = reorder(Gene, -!!sym(metric)), 
                             y = !!sym(metric))) +
    geom_bar(aes(fill = is_hub), stat = "identity") +
    scale_fill_manual(values = c("TRUE" = "#FF6B6B", "FALSE" = "#4ECDC4"),
                      labels = c("TRUE" = "Hub gene", "FALSE" = "Other")) +
    labs(x = "Gene", y = metric, title = paste(metric, "Centrality")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(hjust = 0.5, size = 10),
          legend.position = "none") +
    coord_flip()
  
  plot_list[[i]] <- p
}

# Combine plots
combined_plot <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(title = "Topological Centrality Metrics for Candidate Genes",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 14)))

ggsave(file.path(OUTPUT_DIR, "ppi_centrality_metrics.png"), 
       combined_plot, width = 9, height = 6, dpi = 300)

cat("Centrality metrics figure saved.\n")

# ---- Summary ----
cat("\n=== Stage 4 complete ===\n")
cat("PPI network:", vcount(G), "nodes,", ecount(G), "edges\n")
cat("Hub genes:", paste(hub_genes, collapse = ", "), "\n")

