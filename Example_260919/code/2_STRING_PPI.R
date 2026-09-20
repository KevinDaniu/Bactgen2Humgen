##### 设置工作空间 ####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
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
library(EnsDb.Hsapiens.v86)

#### Load candidate genes ####
INPUT_FILE <- "Example_260919/data/final_homogroup_results.csv"
CID        <- "Conserved_025"
OUTPUT_DIR <- "Example_260919/result"

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#### 读入数据 ####
data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
GENES <- subset(data, Conserved.ID == CID)$`Human.Prot..by.Comprehensive.`
GENES <- GENES[GENES != "" & !is.na(GENES)]
GENES <- unique(GENES)

if (length(GENES) <= 1) stop("基因数小于2")

cat("输入基因数：", length(GENES), "\n")

#### Initialize STRINGdb ####
local_db_dir <- "Example_260919/data/stringdb/"

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
write.csv(interactions_candidate, file = "Example_260919/result/PPI_cluster.csv")
