#!/bin/bash

SELECT_FILE="../select_bacteria_result/representative_list_0.25.txt"
BASE_DIR="../select_bacteria_result/linked_genomes"
FNA_DIR="${BASE_DIR}/selected_fna"
GFF_DIR="${BASE_DIR}/selected_gff"

mkdir -p "$FNA_DIR"
mkdir -p "$GFF_DIR"

# 从 representative_list.txt 读取每个GCF号
while read -r line; do
    # 跳过空行
    [[ -z "$line" ]] && continue
    
    # 提取纯 GCF 号（去掉 .fna.gz 或 .fna 后缀）
    gcf="${line%.fna.gz}"
    gcf="${gcf%.fna}"  # 如果没有 .gz 的情况

    echo "处理: $gcf"

    # 在 fna 目录中查找匹配的文件
    fna_file=$(find "${BASE_DIR}/fna" -name "${gcf}*.fna.gz" | head -n 1)
    if [[ -n "$fna_file" ]]; then
        cp -P "$fna_file" "$FNA_DIR/"
        echo "  ✅ 复制 fna: $(basename "$fna_file")"
    else
        echo "  ❌ 未找到 fna: $gcf"
    fi

    # 在 gff 目录中查找匹配的文件
    gff_file=$(find "${BASE_DIR}/gff" -name "${gcf}*.gff.gz" | head -n 1)
    if [[ -n "$gff_file" ]]; then
        cp -P "$gff_file" "$GFF_DIR/"
        echo "  ✅ 复制 gff: $(basename "$gff_file")"
    else
        echo "  ❌ 未找到 gff: $gcf"
    fi

done < "$SELECT_FILE"

echo "完成！"
