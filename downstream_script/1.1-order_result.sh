#!/bin/bash
# 1.1-order_result.sh
# 对单个菌株的保守簇排序并删除第二列

SAVE_DIR=$1
OUTPUT_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result"
FILE_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene/${SAVE_DIR}"

INPUT="$FILE_DIR/gene_cluster/cluster_genes/gene_cooccurrence_conserve.csv"
OUTFILE="$OUTPUT_DIR/$SAVE_DIR/gene_cooccurrence_conserve_order.csv"

# 输入不存在 -> 跳过
if [[ ! -f "$INPUT" ]]; then
    echo "[SKIP] 输入不存在: $INPUT"
    exit 0
fi

# 已有输出且非空 -> 跳过
if [[ -s "$OUTFILE" ]]; then
    echo "[SKIP] 已完成: $OUTFILE"
    exit 0
fi

mkdir -p "$OUTPUT_DIR/$SAVE_DIR"

(head -n 1 "$INPUT" | cut -d',' -f1,3,4,5 && \
 tail -n +2 "$INPUT" | cut -d',' -f1,3,4,5 | sort -t',' -k3 -rn) \
 > "$OUTFILE"

echo "[OK] $SAVE_DIR -> $OUTFILE"
