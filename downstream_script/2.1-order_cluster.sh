#!/bin/bash
# 2.1-order_cluster.sh
# 为每个菌株提取保守物种数>=CUT_OFF 的结果，并生成总表

RESULT_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result"
SUMMARY_FILE="$RESULT_DIR/high_conserved_summary.csv"
CUT_OFF=200
JOBS=16

# 重新生成汇总表表头（最后再重建，这里先清空）
TMP_SUMMARY="$SUMMARY_FILE.tmp"
echo "菌株名称,保守簇数量,保守簇ID列表" > "$TMP_SUMMARY"

# 收集菌株名
strains=$(for d in "$RESULT_DIR"/*/; do
    [[ -d "$d" ]] && basename "$d"
done)

echo "$strains" | xargs -P "$JOBS" -I {} bash -c '
    strain="{}"
    RESULT_DIR="'"$RESULT_DIR"'"
    CUT_OFF="'"$CUT_OFF"'"
    TMP_SUMMARY="'"$TMP_SUMMARY"'"

    input_file="$RESULT_DIR/$strain/gene_cooccurrence_conserve_order.csv"
    output_file="$RESULT_DIR/$strain/high_conserved_ge${CUT_OFF}.csv"

    # 输入不存在 -> 跳过（只记录日志，不写入汇总表）
    if [[ ! -f "$input_file" ]]; then
        echo "[SKIP-MISSING] $strain"
        exit 0
    fi

    # 输出已有且非空 -> 跳过提取，但仍要写回汇总表
    if [[ -s "$output_file" ]]; then
        echo "[SKIP-DONE] $strain"
    else
        awk -F"," -v cutoff="$CUT_OFF" "NR==1 || \$3>=cutoff" "$input_file" > "$output_file"
    fi

    count=$(tail -n +2 "$output_file" | wc -l)
    cluster_ids=$(tail -n +2 "$output_file" | cut -d"," -f1 | tr "\n" ";" | sed "s/;$//")

    # 用 flock 防止多进程同时写同一文件
    (
        flock 9
        echo "$strain,$count,\"$cluster_ids\"" >> "$TMP_SUMMARY"
    ) 9>>"$TMP_SUMMARY.lock"

    echo "[OK] $strain: $count 个保守簇"
'

# 原子替换汇总表
mv "$TMP_SUMMARY" "$SUMMARY_FILE"
rm -f "$TMP_SUMMARY.lock"

echo ""
echo "完成！"
echo "汇总表已保存到: $SUMMARY_FILE"
echo ""
echo "汇总表预览:"
head -10 "$SUMMARY_FILE"
