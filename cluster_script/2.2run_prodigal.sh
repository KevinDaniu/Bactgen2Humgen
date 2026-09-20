#!/bin/bash

# 定义目录路径
FNA_DIR=~/bacteria_cluster/260811_allgenome_ani/result/all_fna_files
OUTPUT_GFF_DIR=~/bacteria_cluster/260811_allgenome_ani/result/predicted_gff
GCF_LIST=~/bacteria_cluster/260811_allgenome_ani/result/error_gff.txt

# 创建输出目录
mkdir -p "$OUTPUT_GFF_DIR"
mkdir -p "$OUTPUT_GFF_DIR/temp"

# 清空日志
> "${OUTPUT_GFF_DIR}/prodigal.log"

echo "开始批量基因预测..."
START_TIME=$(date +%s)

total=0
success=0
failed=0

while read -r GCF; do
    if [[ -z "$GCF" ]]; then
        continue
    fi
    
    total=$((total + 1))
    
    FNA_GZ="${FNA_DIR}/${GCF}.fna.gz"
    
    if [[ ! -f "$FNA_GZ" ]]; then
        echo "  ✗ 文件不存在: $GCF.fna.gz"
        failed=$((failed + 1))
        continue
    fi
    
    echo "正在处理 ($total): $GCF"
    
    # 解压到临时文件
    TEMP_FNA="${OUTPUT_GFF_DIR}/temp/${GCF}.fna"
    gunzip -c "$FNA_GZ" > "$TEMP_FNA" 2>/dev/null
    
    if [[ ! -s "$TEMP_FNA" ]]; then
        echo "  ✗ 解压失败: $GCF"
        failed=$((failed + 1))
        rm -f "$TEMP_FNA"
        continue
    fi
    
    # 运行 Prodigal（使用临时文件）
    prodigal -i "$TEMP_FNA" \
             -o "${OUTPUT_GFF_DIR}/${GCF}.gff" \
             -f gff \
             -p single \
             2>> "${OUTPUT_GFF_DIR}/prodigal.log"
    
    if [[ $? -eq 0 ]] && [[ -s "${OUTPUT_GFF_DIR}/${GCF}.gff" ]]; then
        echo "  ✓ 成功: ${GCF}.gff"
        success=$((success + 1))
    else
        echo "  ✗ 失败: $GCF"
        failed=$((failed + 1))
    fi
    
    # 删除临时文件（节省空间）
    rm -f "$TEMP_FNA"
    
done < "$GCF_LIST"

# 清理临时目录
rm -rf "${OUTPUT_GFF_DIR}/temp"

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

echo ""
echo "====================================="
echo "批量预测完成！"
echo "  总计: $total 个基因组"
echo "  成功: $success 个"
echo "  失败: $failed 个"
echo "  用时: ${DURATION}秒"
echo "结果保存在: $OUTPUT_GFF_DIR"
echo "====================================="
