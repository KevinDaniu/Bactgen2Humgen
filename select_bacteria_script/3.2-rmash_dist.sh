#!/bin/bash
# 使用多核计算距离矩阵

RESULT_DIR="../select_bacteria_result/"
MASH_WORK_DIR="${RESULT_DIR}/mash_work"

echo "步骤2: 计算距离矩阵 (使用多线程)..."

# 使用 -p 启用多核
./mash_avx512/mash dist \
    -p 48 \
    "${MASH_WORK_DIR}/all_genomes.msh" \
    "${MASH_WORK_DIR}/all_genomes.msh" \
    > "${MASH_WORK_DIR}/aai_distances.txt"

echo "完成！"
