#!/bin/bash

# ===== 配置路径 =====
base_dir="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result"
nextflow_script="bactgen2humgen.nf"
output_base="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene"

script_dir="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/script/script"
rawdata_fna="${base_dir}/all_fna_files"
rawdata_gff="${base_dir}/all_gff_files"
file_format=".gz"

run_blast=$1
run_eggnog=$2
run_mmseqs=$3
dir_name=$4
output_name=${5:-"spacedust"}
chosen_long=${6:-"false"}
chosen_genome=${7:-"none"}

result_dir="${output_base}/${output_name}/gene_cluster"

# 创建输出目录
mkdir -p ${output_base}
mkdir -p "${output_base}/${output_name}"

# 日志文件
log_file="${output_base}/${output_name}/run_${output_name}.log"
> "$log_file"

# ===== 记录开始时间 =====
start_time=$(date +%s)
echo "开始时间: $(date)" | tee -a "$log_file"
        
# 缓存文件路径
CACHE_FILE="${output_base}/${output_name}/query_genome.cache"

# 获取最长基因组（使用缓存）
if [ "$chosen_long" == "true" ]; then
    if [ -f "$CACHE_FILE" ] && [ -s "$CACHE_FILE" ]; then  # -s 检查文件非空
        query_genome=$(cat "$CACHE_FILE")
        echo "使用缓存的查询基因组: $query_genome"
    else
        echo "正在查找最大的基因组文件..."

        # 直接使用 ls -S 按大小排序（最简单可靠）
        largest_file=$(ls -S "$rawdata_fna"/*.fna.gz 2>/dev/null | head -1)

        if [ -n "$largest_file" ] && [ -f "$largest_file" ]; then
            # 提取文件名（去掉路径和 .fna.gz 后缀）
            query_genome=$(basename "$largest_file" .fna.gz)
            file_size=$(ls -lh "$largest_file" | awk '{print $5}')
            echo "最大的文件: $(basename "$largest_file") (大小: $file_size)"
            echo "提取的基因组名: $query_genome"
        else
            echo "错误: 在 $rawdata_fna 中没有找到 .fna.gz 文件"
            exit 1
        fi

        # 缓存结果
        echo "$query_genome" > "$CACHE_FILE"
        echo "已缓存查询基因组: $query_genome"
    fi
else
    query_genome="$chosen_genome"
    echo "使用指定基因组: $query_genome"
fi

# 运行 Nextflow 流程
nextflow run "$nextflow_script" \
	--input_dir "$base_dir" \
	--fna_dir "$rawdata_fna" \
	--gff_dir "$rawdata_gff" \
	--file_format "$file_format" \
        --output_dir "$result_dir" \
    	--script_dir "$script_dir" \
        --query_genome "$query_genome" \
	--run_blast "$run_blast" \
	--run_eggnog "$run_eggnog" \
	--run_mmseqs "$run_mmseqs" \
	--homo_pattern "strict" \
	--max_parallel 48 \
	-work-dir $dir_name \
	-resume \
	2>&1 | tee -a "$log_file"

# ===== 记录结束时间和运行时长 =====
end_time=$(date +%s)
echo "结束时间: $(date)" | tee -a "$log_file"
duration=$((end_time - start_time))
echo "运行时长: $((duration / 3600))小时 $(((duration % 3600) / 60))分钟 $((duration % 60))秒" | tee -a "$log_file"

echo "批量处理完成！" | tee -a "$log_file"
echo "结果保存在: $output_base" | tee -a "$log_file"
echo "========================================" | tee -a "$log_file"
