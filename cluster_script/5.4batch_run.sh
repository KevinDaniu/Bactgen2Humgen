#!/bin/bash
# 用法: bash 5.2batch_run_parallel.sh <任务文件> [节点编号]

TASK_FILE=${1:-"task_splits/task_1.txt"}
NODE_ID=${2:-"node_$(hostname)"}

# 设置路径
FNA_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_fna_files"
FILE_FORMAT=".gz"
OUTPUT_BASE="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene"

# 创建日志目录
mkdir -p logs

echo "========================================"
echo "节点: $NODE_ID"
echo "任务文件: $TASK_FILE"
echo "开始时间: $(date)"
echo "========================================"

# 读取任务列表
mapfile -t genomes < "$TASK_FILE"
total=${#genomes[@]}
current=0

for genome_name in "${genomes[@]}"; do
    # 跳过空行
    [ -z "$genome_name" ] && continue
    
    current=$((current + 1))
    
    # ===== 检查是否已完成 =====
    result_dir="${OUTPUT_BASE}/${genome_name}/gene_cluster/cluster_genes"
    done_file="${result_dir}/gene_cooccurrence_conserve.csv"
    
    if [ -f "$done_file" ] && [ -s "$done_file" ]; then
        echo "[$current/$total] ✓ 跳过已完成: $genome_name (找到 $done_file)"
        continue
    fi
    
    echo "========================================"
    echo "[$current/$total] 处理基因组: $genome_name"
    echo "开始时间: $(date)"
    echo "========================================"
    
    # 创建临时工作目录（使用节点ID避免冲突）
    work_dir="work_${NODE_ID}_${genome_name}"
    mkdir -p "$work_dir"
    
    # 运行处理
    bash 5.1run_bact2hum.sh false false true "$work_dir" "$genome_name" false "$genome_name"
    
    # 检查执行结果
    if [ $? -eq 0 ]; then
        echo "✓ 完成: $genome_name"
    else
        echo "✗ 失败: $genome_name"
        # 记录失败
        echo "$genome_name" >> "logs/failed_${NODE_ID}.log"
    fi
    
    # 清理工作目录（可选）
    rm -rf "$work_dir"
    
    echo ""
done

echo "========================================="
echo "节点 $NODE_ID 任务完成！结束时间: $(date)"
echo "========================================="
