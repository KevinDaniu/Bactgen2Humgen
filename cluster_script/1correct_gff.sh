#!/bin/bash

# 设置路径
input_dir="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_gff_files"
output_dir="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_cgff_files"
file_format=".gff.gz"

# 创建输出目录
mkdir -p "$output_dir"

# 定义处理单个文件的函数
process_file() {
    local gz_file=$1
    local output_dir=$2
    local file_format=$3
    local base_name=$(basename "$gz_file" "$file_format")
    
    echo "[$(date +%H:%M:%S)] 处理文件: $base_name"
    
    # 使用 awk 修正 GFF 文件
    if [[ "$gz_file" == *.gff.gz ]]; then
        zcat "$gz_file" | awk 'BEGIN {FS="\t"; OFS="\t"}
          ##sequence-region 行记下新 seqid,重置"第一条"标志
          $1 ~ /^##sequence-region/ {seqid=$2; first[seqid]=1; print; next}
          # 注释行原样输出
          $1 ~ /^#/ {print; next}
          # 非 CDS 行原样输出
          $3!="CDS" {print; next}
          # 规则 1:第一条 CDS
          first[seqid] && $7=="-" { $7="+"; first[seqid]=0 }
          # 规则 2:start==1 的后续 CDS
          $4==1 && $7=="-" { $7="+" }
          # 输出
          {print}
        ' > "$output_dir/${base_name}_corrected.gff"
    else
        cat "$gz_file" | awk 'BEGIN {FS="\t"; OFS="\t"}
          ##sequence-region 行记下新 seqid,重置"第一条"标志
          $1 ~ /^##sequence-region/ {seqid=$2; first[seqid]=1; print; next}
          # 注释行原样输出
          $1 ~ /^#/ {print; next}
          # 非 CDS 行原样输出
          $3!="CDS" {print; next}
          # 规则 1:第一条 CDS
          first[seqid] && $7=="-" { $7="+"; first[seqid]=0 }
          # 规则 2:start==1 的后续 CDS
          $4==1 && $7=="-" { $7="+" }
          # 输出
          {print}
        ' > "$output_dir/${base_name}_corrected.gff"
    fi
    
    if [ $? -eq 0 ]; then
        echo "✓ 完成: $base_name"
        return 0
    else
        echo "✗ 失败: $base_name"
        return 1
    fi
}

# 导出函数
export -f process_file

# 获取CPU核心数（默认使用所有核心）
if command -v nproc &> /dev/null; then
    CPU_CORES=$(nproc)
else
    CPU_CORES=$(grep -c ^processor /proc/cpuinfo)
fi

echo "使用 $CPU_CORES 个并行核心"
echo "开始并行处理..."
echo "处理文件格式: *${file_format}"

# 使用parallel并行处理
find -L "$input_dir" -maxdepth 1 -name "*${file_format}" -type f | \
    parallel --jobs "$CPU_CORES" --bar \
        "process_file {} $output_dir $file_format"

# 统计结果
total=$(find -L "$input_dir" -maxdepth 1 -name "*${file_format}" -type f | wc -l)
success=$(find -L "$output_dir" -name "*_corrected.gff" -type f | wc -l)
failed=$((total - success))

echo "========================================"
echo "所有文件处理完成！"
echo "总计: $total, 成功: $success, 失败: $failed"
echo "结果保存在: $output_dir"
