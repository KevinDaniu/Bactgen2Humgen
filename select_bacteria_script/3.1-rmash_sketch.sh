#!/bin/bash
# 使用多核生成sketch

RESULT_DIR="../select_bacteria_result"
FNA_LINK_DIR="${RESULT_DIR}/linked_genomes/fna"
MASH_WORK_DIR="${RESULT_DIR}/mash_work"

mkdir -p "$MASH_WORK_DIR"

echo "步骤1: 生成Mash sketches (使用32线程)..."

# 使用 -P 32 和 -p 参数
./mash_avx512/mash sketch \
    -p 48 \
    -o "${MASH_WORK_DIR}/all_genomes" \
    -k 21 \
    -s 10000 \
    -a \
    "${FNA_LINK_DIR}"/*.fna.gz

echo "完成！"
