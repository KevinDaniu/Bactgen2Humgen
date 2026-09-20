#!/bin/bash

# 设置路径
#HUMAN_SEQ_DIR="../../result/gene2human/mmseqs/human_sequences"
#CLUSTER_DIR="../../result/gene2human/mmseqs/human_clusters"

HUMAN_SEQ_DIR="$1"
CLUSTER_DIR="$2"

mkdir -p "$CLUSTER_DIR"

# 遍历每个 human_hits.fa 文件
for fasta_file in "$HUMAN_SEQ_DIR"/*_human_hits.fa; do
    # 获取 basename（去掉 _human_hits.fa 后缀）
    basename=$(basename "$fasta_file" _human_hits.fa)
    
    # 检查文件是否为空
    if [ ! -s "$fasta_file" ]; then
        echo "Skipping $basename (empty file)"
        continue
    fi
    
    echo "========================================="
    echo "Clustering: $basename"
    echo "========================================="
    
    # 统计输入的序列数
    seq_count=$(grep -c "^>" "$fasta_file")
    echo "Input sequences: $seq_count"
    
    # 使用 MMseqs2 聚类
    mmseqs easy-cluster "$fasta_file" \
        "$CLUSTER_DIR/${basename}_cluster" \
        "$CLUSTER_DIR/tmp_${basename}" \
        --min-seq-id 0.3 \
        -c 0.3 \
        --threads 48
    
    # 检查结果
    if [ -f "$CLUSTER_DIR/${basename}_cluster_cluster.tsv" ]; then
        cluster_count=$(wc -l < "$CLUSTER_DIR/${basename}_cluster_cluster.tsv")
        echo "Output clusters: $cluster_count"
        
        # 显示前几个簇的信息
        echo "Top clusters (representative and member count):"
        awk -F'\t' '{count[$1]++} END {for(rep in count) print rep, count[rep]}' \
            "$CLUSTER_DIR/${basename}_cluster_cluster.tsv" | \
            sort -k2,2nr | head -5
    else
        echo "ERROR: Clustering failed for $basename"
    fi
    
    # 只处理第一个文件后停止
    #break
    echo ""
done

# 清理临时目录
rm -rf "$CLUSTER_DIR"/tmp_*

echo "All clustering completed!"
echo "Results in: $CLUSTER_DIR"
