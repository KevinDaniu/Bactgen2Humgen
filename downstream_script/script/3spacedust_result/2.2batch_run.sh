#!/bin/bash

INPUT_BASE=$1
BASE_DIR=$2
MAX_JOBS=${3:-96}  # 默认64个并行，可通过第3个参数指定

# 检查输入目录是否存在
if [ ! -d "$INPUT_BASE" ]; then
    echo "错误: 输入目录 $INPUT_BASE 不存在"
    exit 1
fi

# 获取所有文件夹列表
folders=()
for dir in "$INPUT_BASE"/*/; do
    if [ -d "$dir" ]; then
        folder=$(basename "$dir")
        folders+=("$folder")
    fi
done

total=${#folders[@]}
if [ $total -eq 0 ]; then
    echo "警告: 在 $INPUT_BASE 中没有找到子文件夹"
    exit 0
fi

echo "========================================"
echo "开始批量排序（并行版本）"
echo "输入目录: $INPUT_BASE"
echo "脚本目录: $BASE_DIR"
echo "总文件夹数: $total"
echo "并行数: $MAX_JOBS"
echo "========================================"

# 定义处理函数
process_folder() {
    local folder=$1
    local input_base=$2
    local base_dir=$3
    
    echo "[$(date '+%H:%M:%S')] 开始处理: $folder"
    
    # 运行排序脚本
    if python "${base_dir}/2.1order_result.py" --name "$folder" --input_dir "$input_base" 2>&1; then
        echo "[$(date '+%H:%M:%S')] ✓ 完成: $folder"
        return 0
    else
        echo "[$(date '+%H:%M:%S')] ✗ 失败: $folder"
        return 1
    fi
}

export -f process_folder

# 使用parallel并行处理
printf '%s\n' "${folders[@]}" | \
    parallel --jobs "$MAX_JOBS" --bar \
        "process_folder {} $INPUT_BASE $BASE_DIR"

# 统计结果
echo "========================================"
echo "全部完成！"
echo "========================================"
