##### 设置工作空间 ####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Example_260919/")
rm(list = ls())
gc()

#### 导入R包 ####
library(data.table)

####数据分析####
RHO_CUT    <- 0.3    # |rho| threshold for C
PADJ_CUT   <- 0.05   # BH-adjusted p threshold for C
PPI_CUT    <- 400    # STRING combined_score threshold for P
FDR_CUT    <- 0.05   # enrichment FDR threshold for E

# ---- genes and display names -------------------------------------------------
GENES <- c("AC004832.3", "ETNPPL", "LONP1", "GCAT", "LIAS")
# GTEx v8 (GENCODE v26) legacy symbol for ENSG00000249590
GTEx_ALIAS <- c("AC004832.3" = "RP4-539M6.19")

norm_name <- function(x) {
  x <- as.character(x)
  x <- sub("\\.[0-9]+$", "", x)                 # strip Ensembl version
  x[x %in% c("ENSG00000249590", "RP4-539M6.19")] <- "AC004832.3"
  x
}

# ---- read inputs -------------------------------------------------------------
coexp <- fread("result/coexpression_pantissue.csv")
coexp[, gene1 := norm_name(gene1)]
coexp[, gene2 := norm_name(gene2)]

# PPI file: xlsx in disguise?
ppi <- fread("result/PPI_cluster.csv")
stopifnot(all(c("from_name", "to_name", "combined_score") %in% names(ppi)))
ppi[, from_name := norm_name(from_name)]
ppi[, to_name   := norm_name(to_name)]

enr <- fread("result/string_enrichment_shared_partners.csv")
enr[, gene_list := strsplit(genes, ",", fixed = TRUE)]

# ---- pairwise scoring --------------------------------------------------------
pairs <- combn(GENES, 2)
detail <- data.table(
  gene1   = pairs[1, ],
  gene2   = pairs[2, ],
  rho     = NA_real_, adj_p = NA_real_, C = NA_integer_,
  ppi_score = NA_real_, P = 0L,
  shared_terms = "", E = 0L,
  score = NA_integer_
)

for (i in seq_len(ncol(pairs))) {
  g1 <- pairs[1, i]; g2 <- pairs[2, i]

  # C: co-expression
  hit <- coexp[(gene1 == g1 & gene2 == g2) | (gene1 == g2 & gene2 == g1)]
  if (nrow(hit) > 0 && !is.na(hit$adj_pvalue[1]) && !is.na(hit$spearman_rho[1])) {
    detail$rho[i]   <- hit$spearman_rho[1]
    detail$adj_p[i] <- hit$adj_pvalue[1]
    detail$C[i] <- as.integer(hit$adj_pvalue[1] < PADJ_CUT &
                              abs(hit$spearman_rho[1]) >= RHO_CUT)
  }                                  # else stays NA -> "n.q."

  # P: direct STRING edge
  ph <- ppi[(from_name == g1 & to_name == g2) | (from_name == g2 & to_name == g1)]
  if (nrow(ph) > 0) {
    detail$ppi_score[i] <- ph$combined_score[1]
    detail$P[i] <- as.integer(ph$combined_score[1] >= PPI_CUT)
  }

  # E: shared enriched term
  terms <- enr[sapply(gene_list, function(gl) g1 %in% gl && g2 %in% gl) &
               fdr < FDR_CUT]
  if (nrow(terms) > 0) {
    detail$shared_terms[i] <- paste(terms$description, collapse = "; ")
    detail$E[i] <- 1L
  }
}

detail[, score := fifelse(is.na(C), P + E, C + P + E)]

# ---- matrices ----------------------------------------------------------------
mk <- function(col) {
  m <- matrix(NA_real_, length(GENES), length(GENES),
              dimnames = list(GENES, GENES))
  for (i in seq_len(nrow(detail))) {
    m[detail$gene1[i], detail$gene2[i]] <- detail[[col]][i]
    m[detail$gene2[i], detail$gene1[i]] <- detail[[col]][i]
  }
  diag(m) <- 0
  m
}
M_score <- mk("score")
M_C     <- mk("C")
M_P     <- mk("P")
M_E     <- mk("E")

# cell annotation: score, with "n.q." where C was not quantifiable
lab <- matrix("", length(GENES), length(GENES),
              dimnames = list(GENES, GENES))
for (g1 in GENES) for (g2 in GENES) {
  if (g1 == g2) next
  i <- which(detail$gene1 == g1 & detail$gene2 == g2 |
             detail$gene1 == g2 & detail$gene2 == g1)
  if (is.na(detail$C[i])) {
    lab[g1, g2] <- paste0(detail$score[i], "\n(n.q.)")
  } else {
    lab[g1, g2] <- as.character(detail$score[i])
  }
}

# ---- plot --------------------------------------------------------------------
# row/column labels: flag the not-quantified gene
labs <- GENES
names(labs) <- GENES
has_na <- GENES[sapply(GENES, function(g)
  any(is.na(detail$C[(detail$gene1 == g | detail$gene2 == g)])))]
names(labs)[names(labs) %in% has_na] <-
  paste0(has_na, "*")   # * = not quantified in GTEx v8 pipeline (see legend)

library(pheatmap)


pal <- c("#F5F5F0", "#FDF3D0", "#F9C46B", "#E8853D")  # 0,1,2,3
brk <- seq(-0.5, 3.5, by = 1)

pdf("result/Pairwise_evidence.pdf")
pheatmap(M_score,
  color            = pal,
  breaks           = brk,
  cluster_rows     = FALSE,
  cluster_cols     = FALSE,
  display_numbers  = lab,
  number_format    = "%.0f",
  legend_breaks    = 0:3,
  legend_labels    = c("0", "1", "2", "3"),
  border_color     = "white",
  cellwidth        = 46,
  cellheight       = 46,
  fontsize         = 12,
  fontsize_number  = 11,
  angle_col        = 45
)
dev.off()
