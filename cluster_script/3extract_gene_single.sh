#!/bin/bash

FNA_DIR="../result/all_fna_files"
GFF_DIR="../result/all_gff_files"
OUTPUT_DIR="../result/gff_meta"
FILE_FORMAT=".gz"
mkdir -p $OUTPUT_DIR

# 设置日志文件
LOG_FILE="${OUTPUT_DIR}/extract_cds.log"

# 获取所有fna文件的basename
for fna_file in ${FNA_DIR}/*.fna${FILE_FORMAT}; do
    # 提取样本名称，例如从"GCF_000009705.1_ASM970v1.fna.gz"提取"GCF_000009705.1_ASM970v1"
    prefix=$(basename "$fna_file" ".fna${FILE_FORMAT}")
    #sample=$(echo "$prefix" | cut -d'_' -f1-2)
    sample=$(echo "$prefix")
    
    # 为每个样本创建临时目录
    TEMP_DIR="/tmp/temp_${sample}"
    mkdir -p $TEMP_DIR

    # 解压文件
    if [[ ${FILE_FORMAT} == ".gz" ]]; then
	    echo "Decompressing files..."
	    gunzip -c "${GFF_DIR}/${sample}.gff${FILE_FORMAT}" > "${TEMP_DIR}/${sample}.gff"
	    gunzip -c "${FNA_DIR}/${sample}.fna${FILE_FORMAT}" > "${TEMP_DIR}/${sample}.fna"
    fi
    echo "======================================" | tee -a "$LOG_FILE"
    echo "Processing sample: $sample" | tee -a "$LOG_FILE"
    echo "======================================" | tee -a "$LOG_FILE"

    # 从GFF提取CDS区域并转换为BED格式
    echo "Converting GFF to BED..." | tee -a "$LOG_FILE"
    awk -F'\t' '$3 == "CDS"' "${TEMP_DIR}/${sample}.gff" | \
	    awk -F'\t' -v OFS='\t' '{
                match($0, /ID=([^;]+)/, id_arr);
		if (id_arr[1] != "") {
		    id = id_arr[1];
	 	} else {
		    id = $1 "_" $4 "_" $5;
	        }
                # BED格式: chrom, start, end, name, score, strand
                # GFF的start是1-based，BED需要0-based
                print $1, $4-1, $5, id, "0", $7
	}' > "${GFF_DIR}/${sample}_cds.bed"

    # 使用bedtools提取序列
    echo "Extracting sequences with bedtools..." | tee -a "$LOG_FILE"
    bedtools getfasta \
        -fi "${TEMP_DIR}/${prefix}.fna" \
        -bed "${GFF_DIR}/${sample}_cds.bed" \
        -fo "${OUTPUT_DIR}/${sample}_cds.fa" \
        -name \
        -s

    # 清理临时文件
    rm -rf $TEMP_DIR

    echo "Done! Output: ${OUTPUT_DIR}/${sample}_cds.fa" | tee -a "$LOG_FILE"

    echo "" | tee -a "$LOG_FILE"
done

echo "All samples processed!" | tee -a "$LOG_FILE"
