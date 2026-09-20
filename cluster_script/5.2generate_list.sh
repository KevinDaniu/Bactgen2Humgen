#!/bin/bash

FNA_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_fna_files"
FILE_FORMAT=".gz"
OUTPUT_LIST="genome_list.txt"

# 更高效的版本：使用 sed 一次性处理
find -L "$FNA_DIR" -maxdepth 1 -name "*.fna${FILE_FORMAT}" -type f | \
    grep -v "\.fna\.fai" | \
    sed "s|^$FNA_DIR/||" | \
    sed "s|\.fna$FILE_FORMAT$||" > "$OUTPUT_LIST"

echo "生成基因组列表: $(wc -l < $OUTPUT_LIST) 个"
