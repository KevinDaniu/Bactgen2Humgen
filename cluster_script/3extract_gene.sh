#!/bin/bash

FNA_DIR="../result/all_fna_files"
GFF_DIR="../result/all_gff_files"
OUTPUT_DIR="../result/gff_meta"
FILE_FORMAT=".gz"
THREADS=$1  # 并行核心数

mkdir -p $OUTPUT_DIR

# 设置日志文件
LOG_FILE="${OUTPUT_DIR}/extract_cds.log"

# 清空旧日志
> "$LOG_FILE"

# 定义处理单个样本的函数
process_sample() {
    local fna_file="$1"
    local FNA_DIR="$2"
    local GFF_DIR="$3"
    local OUTPUT_DIR="$4"
    local FILE_FORMAT="$5"
    local LOG_FILE="$6"
    
    # 提取样本名称
    prefix=$(basename "$fna_file" ".fna${FILE_FORMAT}")
    sample="$prefix"
    
    # 为每个样本创建临时目录（使用$$和随机数避免冲突）
    TEMP_DIR="/tmp/temp_${sample}_$$_${RANDOM}"
    mkdir -p "$TEMP_DIR"
    
    # 重定向日志到临时文件，避免并行写入冲突
    local log_tmp="${TEMP_DIR}/log.tmp"
    
    {
        echo "======================================"
        echo "Processing sample: $sample"
        echo "======================================"
        
        # 解压文件
        if [[ ${FILE_FORMAT} == ".gz" ]]; then
            echo "Decompressing files..."
            gunzip -c "${GFF_DIR}/${sample}.gff${FILE_FORMAT}" > "${TEMP_DIR}/${sample}.gff" 2>/dev/null
            gunzip -c "${FNA_DIR}/${sample}.fna${FILE_FORMAT}" > "${TEMP_DIR}/${sample}.fna" 2>/dev/null
        fi
        
        # 检查解压是否成功
        if [[ ! -f "${TEMP_DIR}/${sample}.gff" ]] || [[ ! -f "${TEMP_DIR}/${sample}.fna" ]]; then
            echo "ERROR: Failed to decompress files for $sample"
            rm -rf "$TEMP_DIR"
            return 1
        fi
        
        # 从GFF提取CDS区域并转换为BED格式
        echo "Converting GFF to BED..."
        awk -F'\t' '$3 == "CDS"' "${TEMP_DIR}/${sample}.gff" | \
        awk -F'\t' -v OFS='\t' '{
            match($0, /ID=([^;]+)/, id_arr);
            if (id_arr[1] != "") {
                id = id_arr[1];
            } else {
                id = $1 "_" $4 "_" $5;
            }
            print $1, $4-1, $5, id, "0", $7
        }' > "${TEMP_DIR}/${sample}_cds.bed"
        
        # 使用bedtools提取序列
        echo "Extracting sequences with bedtools..."
        bedtools getfasta \
            -fi "${TEMP_DIR}/${sample}.fna" \
            -bed "${TEMP_DIR}/${sample}_cds.bed" \
            -fo "${OUTPUT_DIR}/${sample}_cds.fa" \
            -name \
            -s 2>/dev/null
        
        if [[ -f "${OUTPUT_DIR}/${sample}_cds.fa" ]]; then
            echo "Done! Output: ${OUTPUT_DIR}/${sample}_cds.fa"
        else
            echo "ERROR: bedtools failed for $sample"
        fi
        
        echo ""
        
    } >> "$LOG_FILE" 2>&1
    
    # 清理临时文件
    rm -rf "$TEMP_DIR"
}

export -f process_sample

echo "Starting parallel processing with $THREADS cores..." | tee -a "$LOG_FILE"
echo "Start time: $(date)" | tee -a "$LOG_FILE"

# 使用find + parallel处理所有fna文件
find -L "${FNA_DIR}" -maxdepth 1 -name "*.fna${FILE_FORMAT}" -type f | \
    parallel -j "$THREADS" \
        process_sample {} "$FNA_DIR" "$GFF_DIR" "$OUTPUT_DIR" "$FILE_FORMAT" "$LOG_FILE"

echo "======================================" | tee -a "$LOG_FILE"
echo "All samples processed!" | tee -a "$LOG_FILE"
echo "End time: $(date)" | tee -a "$LOG_FILE"
