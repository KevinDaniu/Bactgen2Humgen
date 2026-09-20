#!/bin/bash
# 2-link_genomes.sh - 使用绝对路径

RESULT_DIR="../select_bacteria_result"
ASSEMBLY_FILE="${RESULT_DIR}/matched_gcf_list.txt"
INDEX_FILE="${RESULT_DIR}/ncbi_db_list.txt"
LINK_DIR="${RESULT_DIR}/linked_genomes"
# 使用绝对路径
BASE_DIR="/home/liuzhh/bacteria_cluster/bacteria_genomes/20260201"
SOURCE_DIR="${BASE_DIR}/RefSeq/bacteria/GCF"

mkdir -p "$LINK_DIR/fna" "$LINK_DIR/gff"

echo "开始创建链接: $(date)"

# 先用脚本创建链接
awk -v acc_file="$ASSEMBLY_FILE" -v link_dir="$LINK_DIR" -v source_dir="$SOURCE_DIR" '
BEGIN {
    while ((getline acc < acc_file) > 0) {
        gsub(/\.[0-9]+$/, "", acc)
        keep[acc] = 1
    }
    close(acc_file)
}
$1 in keep {
    if ($2 == "fna") {
        print "ln -sf " source_dir "/" $3 " " link_dir "/fna/"
    } else if ($2 == "gff") {
        print "ln -sf " source_dir "/" $3 " " link_dir "/gff/"
    }
}
' "$INDEX_FILE" | sh

echo "完成: $(date)"
echo "FNA: $(ls -1 $LINK_DIR/fna/*.fna.gz 2>/dev/null | wc -l)"
echo "GFF: $(ls -1 $LINK_DIR/gff/*.gff.gz 2>/dev/null | wc -l)"
