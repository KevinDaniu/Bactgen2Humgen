#!/bin/bash

GENOME_LIST="genome_list.txt"
SPLIT_DIR="task_splits"
mkdir -p "$SPLIT_DIR"

total=$(wc -l < "$GENOME_LIST")
num_nodes=${1:-2000}
lines_per_node=$(( (total + num_nodes - 1) / num_nodes ))

# 获取CPU核心数（兼容Linux和macOS）
if command -v nproc &> /dev/null; then
    cores=$(nproc)
elif command -v sysctl &> /dev/null; then
    cores=$(sysctl -n hw.ncpu)
else
    cores=4  # 默认值
fi

echo "总基因组数: $total"
echo "节点数: $num_nodes"
echo "每个节点任务数: $lines_per_node"
echo "CPU核心数: $cores"

# 并行拆分：将文件分成 num_nodes 块，并行写入
if command -v parallel &> /dev/null; then
    # 使用 parallel，指定并行任务数
    cat "$GENOME_LIST" | parallel -j "$cores" --pipe -N "$lines_per_node" \
        "cat > ${SPLIT_DIR}/task_{#}.txt"
else
    # 降级方案：使用 split
    split -d -a 4 -l "$lines_per_node" "$GENOME_LIST" "${SPLIT_DIR}/task_"
    # 并行重命名（使用获取到的核心数）
    find "$SPLIT_DIR" -name "task_*" ! -name "*.txt" | \
        xargs -P "$cores" -I {} mv {} {}.txt
fi

echo "任务已拆分到 $SPLIT_DIR"
ls -lh "$SPLIT_DIR"/
