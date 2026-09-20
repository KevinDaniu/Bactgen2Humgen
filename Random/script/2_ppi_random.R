####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Random/")
rm(list = ls())
gc()

####导入R包####
library(dplyr)
library(rbioapi)
library(igraph)
library(ggraph)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)

#### 单个基因集的 PPI 网络 ####
get_ppi <- function(genes, species = 9606) {
  genes <- unique(genes[!is.na(genes) & genes != ""])
  if (length(genes) < 2) return(NULL)
  
  ids <- tryCatch(rba_string_map_ids(genes, species = species),
                  error = function(e) NULL)
  if (is.null(ids) || nrow(ids) == 0) return(NULL)
  
  net <- tryCatch(
    rba_string_interactions_network(ids$stringId, species = species,
                                    add_nodes = 0, required_score = 400),
    error = function(e) NULL)
  if (is.null(net) || nrow(net) == 0) return(NULL)
  
  g <- graph_from_edgelist(
    as.matrix(net[, c("preferredName_A", "preferredName_B")]),
    directed = FALSE)
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
}

#### 多簇合并成一张图（节点带 cluster 属性；无互作也保留）####
build_multi <- function(genes_list) {
  edges <- list()
  nodes <- list()
  
  for (cn in names(genes_list)) {
    genes_cn <- unique(genes_list[[cn]][
      !is.na(genes_list[[cn]]) & genes_list[[cn]] != ""])
    if (length(genes_cn) == 0) next
    
    g <- get_ppi(genes_cn)
    Sys.sleep(1)
    
    if (!is.null(g)) {
      # 有互作：记录边 + 实际参与互作的节点
      el <- as_edgelist(g)
      edges[[cn]] <- data.frame(from = el[, 1], to = el[, 2], cluster = cn,
                                stringsAsFactors = FALSE)
      
      # 节点 = 网络里出现的节点 + 该簇输入但未进入网络的基因
      net_nodes  <- V(g)$name
      miss_nodes <- setdiff(genes_cn, net_nodes)
      nodes[[cn]] <- data.frame(
        name    = c(net_nodes, miss_nodes),
        cluster = cn,
        stringsAsFactors = FALSE
      )
    } else {
      # 无互作：全部作为孤立节点
      nodes[[cn]] <- data.frame(name = genes_cn, cluster = cn,
                                stringsAsFactors = FALSE)
    }
  }
  
  # 合并节点
  nodes <- bind_rows(nodes) %>% distinct(name, .keep_all = TRUE)
  
  # 合并边，并加空值保护
  edges <- bind_rows(edges)
  if (is.null(edges) || nrow(edges) == 0) {
    edges <- data.frame(from = character(0), to = character(0),
                        cluster = character(0), stringsAsFactors = FALSE)
  }
  edges <- edges[edges$from %in% nodes$name & edges$to %in% nodes$name, ]
  
  # 若完全没有边，用空边表 + 顶点表建图
  if (nrow(edges) == 0) {
    g <- make_empty_graph(n = nrow(nodes), directed = FALSE)
    V(g)$name    <- nodes$name
    V(g)$cluster <- nodes$cluster
    return(g)
  }
  
  graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
}

#### 绘图 ####
draw_net <- function(g, title, cols) {
  set.seed(42)
  lay <- create_layout(g, layout = "fr")
  ggraph(lay) +
    geom_edge_link(color = "grey70", alpha = 0.5, width = 0.5) +
    geom_node_point(aes(color = cluster), size = 3.2) +
    geom_node_text(aes(label = name), repel = TRUE, size = 2.2) +
    scale_color_manual(values = cols, name = "Cluster") +
    labs(title = sprintf("%s (%d nodes, %d edges)", title,
                         vcount(g), ecount(g))) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11))
}

#### 读入数据 ####
data <- read.csv("data/final_homogroup_results.csv", stringsAsFactors = FALSE)
data <- data[!is.na(data$Human.Prot..by.Comprehensive.),]
conserved_id <- read.csv("data/high_conserved_ge200.csv")
data_chosen <- subset(data, Conserved.ID %in% conserved_id$保守簇ID)

cluster_genes <- lapply(unique(data_chosen$Conserved.ID), function(cn) {
  g <- data_chosen$Human.Prot..by.Comprehensive.[data_chosen$Conserved.ID == cn]
  unique(g[!is.na(g) & g != ""])
})
names(cluster_genes) <- unique(data_chosen$Conserved.ID)

#### 每个簇随机抽一次（匹配基因数）####
gene_info <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = AnnotationDbi::keys(org.Hs.eg.db, "ENTREZID"),
  columns = c("SYMBOL", "GENETYPE"),
  keytype = "ENTREZID"
)
pcg <- gene_info$ENTREZID[gene_info$GENETYPE == "protein-coding"]

random_genes <- list()
for (i in seq_len(nrow(conserved_id))) {
  cn <- conserved_id$保守簇ID[i]
  if (!cn %in% names(cluster_genes)) next
  set.seed(i)
  random_genes[[cn]] <- bitr(
    sample(pcg, conserved_id$基因数量[i]),
    "ENTREZID", "SYMBOL", org.Hs.eg.db
  )$SYMBOL
}

#### 选取要画的簇：全部保留，最多 6 个 ####
# 如果想优先画有互作的，可以把 order 改成按边数降序
sel <- names(cluster_genes)
sel <- head(sel, 6)

#### 合并、绘图 ####
g_left  <- build_multi(cluster_genes[sel])
g_right <- build_multi(random_genes[sel])

cols <- setNames(
  colorRampPalette(c("#E64B35", "#4DBBD5", "#00A087",
                     "#3C5488", "#F39B7F", "#8491B4"))(length(sel)),
  sel)

p <- draw_net(g_left,  "Conserved clusters", cols) |
  draw_net(g_right, "Random (matched size)", cols)

dir.create("result", showWarnings = FALSE)
ggsave("result/multi_cluster_vs_random.pdf", p,
       width = 10, height = 5)
ggsave("result/multi_cluster_vs_random.png", p,
       width = 10, height = 5, dpi = 300)
print(p)