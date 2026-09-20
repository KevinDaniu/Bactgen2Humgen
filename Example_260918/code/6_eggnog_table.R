####设置工作空间####
setwd("C:/Users/yuxip/Desktop/LZH/Tool/Example")
rm(list = ls())
gc()

####导入R包####
library(dplyr)

####处理eggNOG注释####
#----读入数据----#
file_path <- dir("data/eggnog_table/", full.names = T)
for(i in 1:length(file_path)){
  if(i == 1){
    data <- read.table(file_path[i], fill = T, sep = "\t", header = T)
    data <- data[,c("Cluster_ID","Gene","Preferred_name","Description")]
  } else {
    d <- read.table(file_path[i], fill = T, sep = "\t", header = T)
    d <- d[,c("Cluster_ID","Gene","Preferred_name","Description")]
    data <- rbind(data, d)
  }

}

#----保存为表格----#
writexl::write_xlsx(data, path = "result/eggnog_all_result.xlsx")
