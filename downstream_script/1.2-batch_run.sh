#!/bin/bash
# 1.2-batch_run.sh
# 批量并行对所有菌株运行排序脚本

RESULT_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ORDER_SCRIPT="$SCRIPT_DIR/1.1-order_result.sh"

# 并行核数，可按机器改
JOBS=16

# 收集所有菌株名（只要目录名含 GCF 的）
strains=$(for d in "$RESULT_DIR"/*/; do
    name=$(basename "$d")
    [[ "$name" =~ GCF ]] && echo "$name"
done)

total=$(echo "$strains" | grep -c . || true)
echo "共发现 $total 个菌株目录，并行核数: $JOBS"

echo "$strains" | xargs -P "$JOBS" -I {} bash -c '
    strain="{}"
    input="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene/${strain}/gene_cluster/cluster_genes/gene_cooccurrence_conserve.csv"
    out="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result/${strain}/gene_cooccurrence_conserve_order.csv"

    if [[ ! -f "$input" ]]; then
        echo "[SKIP-MISSING] $strain"
        exit 0
    fi
    if [[ -s "$out" ]]; then
        echo "[SKIP-DONE] $strain"
        exit 0
    fi

    bash "'"$ORDER_SCRIPT"'" "$strain"
'

echo "批量排序完成！"
