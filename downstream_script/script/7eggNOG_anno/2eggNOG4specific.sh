#!/bin/bash

# 设置参数
#INPUT_DIR="../../result/gene2human/sequence_file"
#INPUT_FILE="$1"
#OUTPUT_DIR="../../result/gene2human/gene4eggnog"
#EGGNOG_DATA_DIR="../../eggnog-mapper-data"
#THREADS=48

INPUT_DIR="$1"
INPUT_FILE="$2"
OUTPUT_DIR="$3"
EGGNOG_DATA_DIR="$4"
THREADS=48

# 创建输出目录
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"
mkdir -p "$OUTPUT_DIR/annos"

# 记录开始时间
start_time=$(date +%s)
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================"

# 获取所有fa文件
fa_file=($(ls "$INPUT_DIR"/"$INPUT_FILE" 2>/dev/null))
total_files=${#fa_file[@]}

if [ $total_files -eq 0 ]; then
    echo "错误: 在 $INPUT_DIR 中没有找到 $INPUT_FILE 文件"
    exit 1
fi

echo "找到 $INPUT_FILE 需要处理"
echo "========================================\n"

# 创建注释文件
filename=$(basename "$INPUT_FILE")
basename="${filename%.*}"
SUMMARY_FILE="$OUTPUT_DIR/annos/${basename}_eggnog_result.txt"
echo -e "Cluster\tTotal_Queries\tMatches\tMatched_Queries\tTop_Hit\tTop_Evalue\tTop_Score\tTop_Pident\tTop_Description" > "$SUMMARY_FILE"

# 计数器
processed=0
total_matches=0

# 开始分析
# 获取文件名（不含路径和扩展名）
output_prefix="${OUTPUT_DIR}/annos/${basename}"
log_file="${OUTPUT_DIR}/logs/${basename}.log"

# 计算进度
processed=$((processed + 1))
echo "[$processed/$total_files] 处理文件: $basename"
echo "开始时间: $(date '+%H:%M:%S')"

# 单个文件的开始时间
file_start=$(date +%s)

# 统计查询序列数
query_count=$(grep -c "^>" "$fa_file")
echo "  查询序列数: $query_count"

# 运行 eggNOG-mapper
emapper.py -i "$fa_file" \
	--itype CDS \
	-o "$output_prefix" \
        -m diamond \
        --cpu $THREADS \
        --data_dir "$EGGNOG_DATA_DIR" \
        --override \
        --sensmode fast \
        --evalue 0.1 \
        --score 30 \
        --query_cover 30 \
        --subject_cover 30 \
        --go_evidence all \
        --pfam_realign realign \
        --usemem \
        2> "$log_file"

# 检查结果
seed_file="${output_prefix}.emapper.seed_orthologs"
annotations_file="${output_prefix}.emapper.annotations"
    
if [ -f "$seed_file" ]; then
	# 统计匹配数（排除注释行）
	match_count=$(grep -v "^#" "$seed_file" | wc -l)
        unique_queries=$(grep -v "^#" "$seed_file" | cut -f1 | sort -u | wc -l)
        total_matches=$((total_matches + match_count))

        echo "  总匹配数: $match_count"
        echo "  匹配到的唯一序列: $unique_queries/$query_count"

        # 获取最佳匹配（基于e-value）
        if [ $match_count -gt 0 ]; then
            # 提取最佳匹配（按e-value排序）
            best_hit=$(grep -v "^#" "$seed_file" | sort -k4,4g | head -1)
            best_query=$(echo "$best_hit" | cut -f1)
            best_hit_name=$(echo "$best_hit" | cut -f2)
            best_evalue=$(echo "$best_hit" | cut -f4)
            best_score=$(echo "$best_hit" | cut -f5)
            
            # 从annotations文件中获取更多信息
            if [ -f "$annotations_file" ]; then
                best_line=$(grep -m1 "^$best_query" "$annotations_file")
                best_pident=$(echo "$best_line" | cut -f4)
                best_description=$(echo "$best_line" | cut -f10)
            else
                best_pident="-"
                best_description="-"
            fi

            echo "  最佳匹配: $best_hit_name (e-value: $best_evalue, score: $best_score)"

            # 写入汇总
            echo -e "$basename\t$query_count\t$match_count\t$unique_queries\t$best_hit_name\t$best_evalue\t$best_score\t$best_pident\t$best_description" >> "$SUMMARY_FILE"
        else
            echo -e "$basename\t$query_count\t0\t0\t-\t-\t-\t-\t-" >> "$SUMMARY_FILE"
        fi
        
        # 创建简化的同源基因结果
        simplified_file="${OUTPUT_DIR}/annos/${basename}_homologs.txt"
        echo -e "Query\tHit\tEvalue\tScore\tPident" > "$simplified_file"
        if [ -f "$annotations_file" ]; then
            paste <(grep -v "^#" "$seed_file" | cut -f1,2,4,5) \
                  <(grep -v "^#" "$annotations_file" | cut -f4) \
                  | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' >> "$simplified_file"
        else
            grep -v "^#" "$seed_file" | cut -f1,2,4,5 | awk '{print $1"\t"$2"\t"$3"\t"$4"\t-"}' >> "$simplified_file"
        fi
        
    else
        echo "  警告: eggNOG未生成结果文件"
        echo -e "$basename\t$query_count\t0\t0\t-\t-\t-\t-\t-" >> "$SUMMARY_FILE"
    fi

    # 计算单个文件耗时
    file_end=$(date +%s)
    file_duration=$((file_end - file_start))
    echo "完成时间: $(date '+%H:%M:%S'), 耗时: $((file_duration / 60))分 $((file_duration % 60))秒"
    echo "----------------------------------------"

# 统计结果
echo "\n========================================"
echo "所有文件处理完成！"
echo "结果目录: $OUTPUT_DIR"
echo "总匹配数: $total_matches"
echo "汇总文件: $SUMMARY_FILE"

# 显示前10个最佳匹配的簇（按e-value排序）
echo "\n最佳匹配的簇（按e-value排序）:"
if [ -f "$SUMMARY_FILE" ]; then
    # 排除标题行，按e-value排序，显示前10个
    tail -n +2 "$SUMMARY_FILE" | sort -t$'\t' -k6,6g | head -10 | while IFS=$'\t' read -r cluster queries matches unique best_hit evalue score pident desc; do
        printf "  %-15s %-30s e-value: %-10s score: %-5s\n" "$cluster" "$best_hit" "$evalue" "$score"
    done
fi

# 创建最终的同源基因汇总表
FINAL_SUMMARY="$OUTPUT_DIR/annos/${basename}_homology_summary.csv"
echo "Conserved ID,Cluster Gene,Human Gene,E-value,Score,Percent Identity" > "$FINAL_SUMMARY"

for fa_file in "${fa_files[@]}"; do
    basename=$(basename "$fa_file" .fa)
    seed_file="${OUTPUT_DIR}/${basename}.emapper.seed_orthologs"
    annotations_file="${OUTPUT_DIR}/${basename}.emapper.annotations"
    
    if [ -f "$seed_file" ]; then
        # 获取每个查询的最佳匹配
        grep -v "^#" "$seed_file" | sort -k4,4g | awk '!seen[$1]++' | while IFS=$'\t' read -r query hit evalue score rest; do
            # 从annotations文件中获取percent identity
            if [ -f "$annotations_file" ]; then
                pident=$(grep "^$query" "$annotations_file" | head -1 | cut -f4)
            else
                pident="-"
            fi
            # 提取cluster gene（从query中）
            cluster_gene=$(echo "$query" | cut -d'|' -f2)
            echo "$basename,$cluster_gene,$hit,$evalue,$score,$pident" >> "$FINAL_SUMMARY"
        done
    fi
done

echo "\n同源基因汇总表已保存到: $FINAL_SUMMARY"

# annotations文件的处理
annotations_file="${OUTPUT_DIR}/annos/${basename}.emapper.annotations"
clean_file="${OUTPUT_DIR}/annos/${basename}.clean.annotations"

# 去除以 ## 开头的行，保留以 # 开头的列名行和数据行
grep -v '^##' "$annotations_file" | \
awk 'BEGIN{FS=OFS="\t"}
     NR==1 {
         # 列名行：从第2列开始输出
         printf "Cluster_ID\tGene"
         for(i=2; i<=NF; i++) printf "\t%s", $i
         printf "\n"
         next
     }
     {
         # 数据行：拆分第一列，从第2列开始输出
         split($1, a, "|")
         printf "%s\t%s", a[1], a[2]
         for(i=2; i<=NF; i++) printf "\t%s", $i
         printf "\n"
     }' > "$clean_file"

# 总运行时间
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "\n总运行时间: $((duration / 60)) 分 $((duration % 60)) 秒"
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"

# 快速查看结果统计
echo "\n========================================"
echo "结果统计:"
if [ -f "$FINAL_SUMMARY" ]; then
    total_hits=$(tail -n +2 "$FINAL_SUMMARY" | wc -l)
    unique_genes=$(tail -n +2 "$FINAL_SUMMARY" | cut -d',' -f3 | sort -u | wc -l)
    echo "  总匹配数: $total_hits"
    echo "  唯一基因: $unique_genes"
    
    echo "\n最常见的基因:"
    tail -n +2 "$FINAL_SUMMARY" | cut -d',' -f3 | sort | uniq -c | sort -nr | head -10 | while read count gene; do
        printf "  %-5s %s\n" "$count" "$gene"
    done
fi

rm *.annotations
rm *.seed_orthologs
rm *.hits
rm -rf emappertmp_*
echo "已移除当前路径重复的结果文件与临时文件"
