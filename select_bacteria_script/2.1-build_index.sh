#!/bin/bash
# 1-build_index.sh - 极简版

BASE_DIR="../../bacteria_genomes/20260201"
INFO_FILE="${BASE_DIR}/assembly_accessions.info"
INDEX_FILE="../select_bacteria_result/ncbi_db_list.txt"

echo "开始建立索引: $(date)"

# 直接从 assembly_accessions.info 提取
# fna 和 gff 只是后缀不同
# 使用 tail -n +3 跳过前两行注释
awk -F'\t' 'NR>1 && $2 ~ /fna.gz$/ {
    acc = $1
    base = acc
    gsub(/\.[0-9]+$/, "", base)
    # fna
    print base "\tfna\t" $2
    # gff（把 fna.gz 替换为 gff.gz）
    gff_file = $2
    gsub(/fna\.gz$/, "gff.gz", gff_file)
    print base "\tgff\t" gff_file
}' "$INFO_FILE" > "$INDEX_FILE"

total=$(wc -l < "$INDEX_FILE")
echo "索引完成: $(date)"
echo "索引文件: $INDEX_FILE"
echo "共索引 $total 个文件"
echo ""
echo "fna 数量: $(grep -c $'\tfna\t' "$INDEX_FILE")"
echo "gff 数量: $(grep -c $'\tgff\t' "$INDEX_FILE")"
echo ""
echo "索引示例:"
head -10 "$INDEX_FILE"
