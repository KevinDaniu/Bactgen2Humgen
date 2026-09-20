##### 设置工作空间 ####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/")
rm(list = ls())
gc()

#### 导入R包 ####
library(httr)
library(jsonlite)
library(dplyr)

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

## ---- 2. STRING enrichment API（POST）----
resp <- POST("https://version-12.string-db.org/api/json/enrichment",
             body = list(identifiers = paste(GENES, collapse = "\r"),
                         species = 9606),
             encode = "form", timeout(120))
stop_for_status(resp)
enr <- fromJSON(content(resp, as = "text"))
enr_out <- enr[order(enr$fdr), c("category", "term", "description", "fdr", "number_of_genes")]
enr_out$genes <- sapply(enr$inputGenes, function(x){paste0(unlist(x),collapse = ",")})
head(enr_out[, c("category", "term", "description", "fdr")], 15)
write.csv(enr_out, "Example_260919/result/string_enrichment_shared_partners.csv", row.names = FALSE)
cat("\n富集条目数:", nrow(enr_out), "；前 10 项:\n")
print(head(enr_out[, c("category", "description", "fdr")], 10))
