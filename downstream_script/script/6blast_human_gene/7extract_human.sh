#!/bin/bash
#BLAST_DIR="../../result/gene2human/blastp_result"
#OUTPUT_DIR="../../result/gene2human/mmseqs"
#FILE_DIR="human_sequences"
#HUMAN_DB="../../rawdata/gencode.v49.pc_translations.fa"  # 人源蛋白数据库路径

BLAST_DIR="$1"
OUTPUT_DIR="$2"
FILE_DIR="human_sequences"
HUMAN_DB="$3"  # 人源蛋白数据库路径

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/$FILE_DIR"

OUTPUT_DIR="$OUTPUT_DIR/$FILE_DIR"

# 遍历每个 BLAST 结果文件
for blast_file in "$BLAST_DIR"/Conserved_*_blast.txt; do
    basename=$(basename "$blast_file" _blast.txt)
    # 检查 BLAST 结果文件是否为空
    if [ ! -s "$blast_file" ]; then
        echo "Skipping $basename (BLAST result is empty)"
        continue
    fi
    echo "Processing $basename..."

    # 提取该文件中的所有 human gene ID（BLAST 结果第二列）
    cut -f2 "$blast_file" | sort -u > "$OUTPUT_DIR/${basename}_human_genes.txt"

    # 从人源蛋白数据库中提取这些基因的序列
    seqkit grep -f "$OUTPUT_DIR/${basename}_human_genes.txt" \
        "$HUMAN_DB" > "$OUTPUT_DIR/${basename}_human_hits.fa" 2>/dev/null

    # 统计
    gene_count=$(wc -l < "$OUTPUT_DIR/${basename}_human_genes.txt")
    seq_count=$(grep -c "^>" "$OUTPUT_DIR/${basename}_human_hits.fa" 2>/dev/null || echo 0)
    echo "  Found $gene_count unique human genes, extracted $seq_count sequences"
done
