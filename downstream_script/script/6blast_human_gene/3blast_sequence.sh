#!/bin/bash

# 设置参数
INPUT_DIR="$1"
OUTPUT_DIR="$2"
DB="$3"
THREADS=48

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/summary"

# 记录开始时间
start_time=$(date +%s)
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

# 获取所有fa文件
fa_files=($(ls "$INPUT_DIR"/Conserved_*.fa 2>/dev/null))
total_files=${#fa_files[@]}

if [ $total_files -eq 0 ]; then
    echo "错误: 在 $INPUT_DIR 中没有找到 Conserved_*.fa 文件"
    exit 1
fi

echo "找到 $total_files 个文件需要处理"
echo "========================================\n"

# 创建汇总文件
SUMMARY_FILE="$OUTPUT_DIR/summary/all_results.txt"
echo -e "Cluster\tTotal_Queries\tMatches\tMatched_Queries\tTop_Hit\tTop_Identity\tTop_Evalue\tTop_Description" > "$SUMMARY_FILE"

# 计数器
processed=0
total_matches=0

# 遍历每个文件
for fa_file in "${fa_files[@]}"; do
    # 获取文件名（不含路径和扩展名）
    basename=$(basename "$fa_file" .fa)
    output_file="${OUTPUT_DIR}/${basename}_blast.txt"
    
    # 计算进度
    processed=$((processed + 1))
    echo "[$processed/$total_files] 处理文件: $basename"
    echo "开始时间: $(date '+%H:%M:%S')"
    
    # 单个文件的开始时间
    file_start=$(date +%s)
    
    # 统计查询序列数
    query_count=$(grep -c "^>" "$fa_file")
    echo "  查询序列数: $query_count"
    
    # 运行 BLAST
    blastn -query "$fa_file" \
        -db "$DB" \
        -evalue 10 \
        -word_size 7 \
        -perc_identity 70 \
        -num_threads $THREADS \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
        -out "$output_file" \
        2> "${OUTPUT_DIR}/logs/${basename}.log"
    
    # 检查结果
    if [ -f "$output_file" ]; then
        match_count=$(wc -l < "$output_file")
        unique_queries=$(cut -f1 "$output_file" | sort -u | wc -l)
        total_matches=$((total_matches + match_count))
        
        echo "  总匹配数: $match_count"
        echo "  匹配到的唯一序列: $unique_queries/$query_count"
        
        # 获取最佳匹配
        if [ $match_count -gt 0 ]; then
            best_hit=$(sort -k12,12g "$output_file" | head -1)
            best_qseqid=$(echo "$best_hit" | cut -f1)
            best_sseqid=$(echo "$best_hit" | cut -f2)
            best_pident=$(echo "$best_hit" | cut -f3)
            best_evalue=$(echo "$best_hit" | cut -f11)
            best_stitle=$(echo "$best_hit" | cut -f13)
            
            echo "  最佳匹配: $best_sseqid (identity: $best_pident%, evalue: $best_evalue)"
            
            # 写入汇总
            echo -e "$basename\t$query_count\t$match_count\t$unique_queries\t$best_sseqid\t$best_pident\t$best_evalue\t$best_stitle" >> "$SUMMARY_FILE"
        else
            echo -e "$basename\t$query_count\t0\t0\t-\t-\t-\t-" >> "$SUMMARY_FILE"
        fi
    else
        echo "  警告: BLAST未生成结果文件"
        echo -e "$basename\t$query_count\t0\t0\t-\t-\t-\t-" >> "$SUMMARY_FILE"
    fi
    
    # 计算单个文件耗时
    file_end=$(date +%s)
    file_duration=$((file_end - file_start))
    echo "完成时间: $(date '+%H:%M:%S'), 耗时: $((file_duration / 60))分 $((file_duration % 60))秒"
    echo "----------------------------------------"
done

# 统计结果
echo "\n========================================"
echo "所有文件处理完成！"
echo "结果目录: $OUTPUT_DIR"
echo "总匹配数: $total_matches"
echo "汇总文件: $SUMMARY_FILE"

# 显示前10个最佳匹配的簇
echo "\n最佳匹配的簇（按identity排序）:"
sort -k6,6nr "$SUMMARY_FILE" | head -11 | column -t

# 总运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "\n总运行时间: $((duration / 60)) 分 $((duration % 60)) 秒"
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
