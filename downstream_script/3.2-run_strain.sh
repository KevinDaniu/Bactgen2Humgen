#!/bin/bash
# ============ 配置 ============
RESULT_DIR=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result
SCRIPT_DIR=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_script
PY_SCRIPT=$SCRIPT_DIR/3.1-extract_gene_withGID.py
LOGDIR=$SCRIPT_DIR/logs

mkdir -p "$LOGDIR"

# ============ 用法 ============
if [[ $# -lt 1 ]]; then
    echo "用法: bash run_one_strain.sh <strain_name>"
    echo "示例: bash run_one_strain.sh GCF_000007205.1_ASM720v1_genomic"
    exit 1
fi

strain="$1"

# ============ 找簇文件 ============
CSV=$RESULT_DIR/$strain/gene_cooccurrence_conserve_order.csv
if [[ ! -f "$CSV" ]]; then
    echo "[ERROR] 找不到簇文件: $CSV"
    exit 1
fi

echo "菌株: $strain"
echo "簇文件: $CSV"
echo "python 脚本: $PY_SCRIPT"

# ============ 循环跑 ============
# 跳过表头；逐个保守簇 ID
tail -n +2 "$CSV" | cut -d',' -f1 | while read -r cid; do
    [[ -z "$cid" ]] && continue

    done_flag="$LOGDIR/${strain}__${cid}.done"
    log_file="$LOGDIR/${strain}__${cid}.log"

    if [[ -f "$done_flag" ]]; then
        echo "[SKIP-DONE] $cid"
        continue
    fi

    echo "[RUN] $cid"
    if python "$PY_SCRIPT" \
        --base_dir "$strain" \
        --id "$cid"
    then
        touch "$done_flag"
        echo "[OK] $cid"
    else
        echo "[FAIL] $cid  见 $log_file"
    fi
done

echo "完成: $strain"
