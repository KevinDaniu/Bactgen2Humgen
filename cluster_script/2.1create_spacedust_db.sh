#!/bin/bash
# 智能 spacedust数据库创建脚本
# 用法: ./2create_spacedust_db.sh [可选: 并行任务数]

# 设置路径
MAX_JOBS="${1:-96}"  # 并行任务数
OUTPUT_BASE="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/spacedust_db"
FILE_FORMAT=".fna.gz"

# 基础目录设置 - 使用绝对路径
FNA_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_fna_files"
GFF_DIR="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_cgff_files"

# 日志文件 - 放在输出目录
MAIN_LOG="${OUTPUT_BASE}/spacedustDB_creation_$(date +%Y%m%d_%H%M%S).log"

# 创建输出目录（如果不存在）
mkdir -p "${OUTPUT_BASE}"

{
    echo "=== 智能 spacedust target数据库创建脚本 ==="
    echo "开始时间: $(date)"
    echo "最大并行任务数: ${MAX_JOBS}"
    echo "FNA目录: ${FNA_DIR}"
    echo "GFF目录: ${GFF_DIR}"
    echo "输出目录: ${OUTPUT_BASE}"
    echo "----------------------------------------"
} > "${MAIN_LOG}"

# 获取所有fna文件列表
echo "正在扫描FNA目录: ${FNA_DIR}" | tee -a "${MAIN_LOG}"

# 获取所有基因组基础名称
all_genomes=()
# 使用绝对路径查找文件
for fna_file in "${FNA_DIR}"/*${FILE_FORMAT}; do
    # 检查文件是否存在（处理没有匹配文件的情况）
    if [ -f "${fna_file}" ]; then
        # 提取基础名称
        filename=$(basename "${fna_file}")
        base_name=${filename%${FILE_FORMAT}}
        all_genomes+=("${base_name}")
        echo "找到: ${filename} -> ${base_name}" | tee -a "${MAIN_LOG}"
    fi
done

# 如果没有找到文件，报错退出
if [ ${#all_genomes[@]} -eq 0 ]; then
    echo "错误: 在 ${FNA_DIR} 中没有找到 *${FILE_FORMAT} 文件" | tee -a "${MAIN_LOG}"
    exit 1
fi

echo "找到以下基因组文件:" | tee -a "${MAIN_LOG}"
printf '%s\n' "${all_genomes[@]}" | tee -a "${MAIN_LOG}"
echo "总共找到 ${#all_genomes[@]} 个基因组文件" | tee -a "${MAIN_LOG}"
echo "----------------------------------------" | tee -a "${MAIN_LOG}"

# 过滤掉query文件，得到target列表
target_genomes=()
for genome in "${all_genomes[@]}"; do
    target_genomes+=("${genome}")
done

echo "将创建以下spacedust数据库:" | tee -a "${MAIN_LOG}"
printf '%s\n' "${target_genomes[@]}" | tee -a "${MAIN_LOG}"
echo "总共 ${#target_genomes[@]} 个spacedust数据库需要创建" | tee -a "${MAIN_LOG}"
echo "----------------------------------------" | tee -a "${MAIN_LOG}"

# 函数：为单个target创建数据库
create_target_db() {
    local target=$1
    local main_log=$2
    local job_id=$3
    local file_format=$4
    
    # 设置文件路径
    local raw_file="${FNA_DIR}/${target}${file_format}"
    #PREFIX=$(echo ${target} | cut -d'_' -f1-2)
    PREFIX=$(echo ${target})
    local gff_file="${GFF_DIR}/${PREFIX}_corrected.gff"
    local output_dir="${OUTPUT_BASE}/${target}"
    # 临时目录也放在output_base下，避免路径混乱
    local temp_dir="${OUTPUT_BASE}/temp_${target}_${job_id}"
    local log_file="${output_dir}/spacedust_setup_$(date +%Y%m%d_%H%M%S).log"
    
    # 检查必要文件是否存在
    if [[ ! -f "${raw_file}" ]]; then
        echo "[${target}] 错误: fna文件不存在: ${raw_file}" | tee -a "${main_log}"
        return 1
    fi
    
    if [[ ! -f "${gff_file}" ]]; then
        echo "[${target}] 错误: gff文件不存在: ${gff_file}" | tee -a "${main_log}"
        return 1
    fi
    
    # 创建必要的目录
    mkdir -p "${temp_dir}"
    mkdir -p "${output_dir}"
    
    # 开始日志记录
    {
        echo "=== 为 ${target} 创建spacedust数据库 ==="
        echo "开始时间: $(date)"
        echo "raw_file: ${raw_file}"
        echo "gff_file: ${gff_file}"
        echo "output_dir: ${output_dir}"
        echo "temp_dir: ${temp_dir}"
        echo "----------------------------------------"
        
        # 1. 创建 GFF 文件列表
        echo "步骤1: 创建 GFF 文件列表..."
        gff_list="${temp_dir}/gff_files.txt"
        echo "${gff_file}" > "${gff_list}"
        echo "GFF 列表已创建: ${gff_list}"
        
        # 2. 创建 spacedust 数据库
        echo "步骤2: 创建 spacedust 数据库..."
        setdb_name="${output_dir}/setDB"
        tmp_folder="${output_dir}/tmp_$(date +%Y%m%d_%H%M%S)_${job_id}"
        
        echo "数据库名称: ${setdb_name}"
        echo "临时文件夹: ${tmp_folder}"
        
        # 执行 spacedust 命令 - 使用绝对路径！
        echo "开始执行 spacedust createsetdb..."
        spacedust createsetdb "${raw_file}" "${setdb_name}" "${tmp_folder}" \
            --gff-dir "${gff_list}" \
            --gff-type CDS
        
        # 检查执行结果
        if [ $? -eq 0 ]; then
            echo "----------------------------------------"
            echo "spacedust 数据库创建成功!"
            echo "数据库位置: ${setdb_name}"
            
            # 清理临时文件
            rm -rf "${tmp_folder}"
            rm -rf "${temp_dir}"
            echo "已清理临时文件夹"
            echo "----------------------------------------"
            echo "脚本执行完成: $(date)"
            return 0
        else
            echo "----------------------------------------"
            echo "错误: spacedust 数据库创建失败!"
            echo "临时文件夹保留在: ${tmp_folder}"
            echo "----------------------------------------"
            echo "脚本执行失败: $(date)"
            return 1
        fi
    } > "${log_file}" 2>&1
    
    # 获取函数返回值
    local result=$?
    
    # 输出简要信息到主日志
    if [ ${result} -eq 0 ]; then
        echo "[$(date +%H:%M:%S)] ${target}: 数据库创建成功" | tee -a "${main_log}"
    else
        echo "[$(date +%H:%M:%S)] ${target}: 数据库创建失败" | tee -a "${main_log}"
    fi
    
    # 显示最后几行日志
    echo "[${target}] 最后10行日志:" | tee -a "${main_log}"
    tail -10 "${log_file}" | sed "s/^/[${target}] /" | tee -a "${main_log}"
    echo "----------------------------------------" | tee -a "${main_log}"
    
    return ${result}
}

# 导出函数和变量，以便在并行环境中使用
export -f create_target_db
export FNA_DIR GFF_DIR OUTPUT_BASE

# 创建成功和失败的计数文件（放在输出目录）
SUCCESS_FILE="${OUTPUT_BASE}/success_$$.tmp"
FAIL_FILE="${OUTPUT_BASE}/fail_$$.tmp"
> "${SUCCESS_FILE}"
> "${FAIL_FILE}"

# 使用GNU parallel并行执行（如果可用）
if command -v parallel &> /dev/null; then
    echo "使用GNU parallel并行执行，最大任务数: ${MAX_JOBS}" | tee -a "${MAIN_LOG}"
    
    # 准备参数列表 - 使用绝对路径
    printf '%s\n' "${target_genomes[@]}" | \
    parallel -j "${MAX_JOBS}" \
        "create_target_db {} ${MAIN_LOG} {#} ${FILE_FORMAT} && echo {} >> ${SUCCESS_FILE} || echo {} >> ${FAIL_FILE}"
else
    # 如果没有parallel，使用串行执行
    echo "GNU parallel未安装，使用串行执行" | tee -a "${MAIN_LOG}"
    
    job_id=0
    for target in "${target_genomes[@]}"; do
        job_id=$((job_id + 1))
        echo "开始处理 [${job_id}/${#target_genomes[@]}] ${target}" | tee -a "${MAIN_LOG}"
        
        if create_target_db "${target}" "${MAIN_LOG}" "${job_id}" "${FILE_FORMAT}"; then
            echo "${target}" >> "${SUCCESS_FILE}"
        else
            echo "${target}" >> "${FAIL_FILE}"
        fi
        
        echo "----------------------------------------" | tee -a "${MAIN_LOG}"
    done
fi

# 统计结果
SUCCESS_COUNT=$(wc -l < "${SUCCESS_FILE}" 2>/dev/null || echo 0)
FAIL_COUNT=$(wc -l < "${FAIL_FILE}" 2>/dev/null || echo 0)

echo "========================================" | tee -a "${MAIN_LOG}"
echo "数据库创建完成统计:" | tee -a "${MAIN_LOG}"
echo "成功: ${SUCCESS_COUNT}" | tee -a "${MAIN_LOG}"
echo "失败: ${FAIL_COUNT}" | tee -a "${MAIN_LOG}"

if [ ${FAIL_COUNT} -gt 0 ]; then
    echo "失败的基因组:" | tee -a "${MAIN_LOG}"
    cat "${FAIL_FILE}" | tee -a "${MAIN_LOG}"
fi

echo "完成时间: $(date)" | tee -a "${MAIN_LOG}"
echo "主日志文件: ${MAIN_LOG}" | tee -a "${MAIN_LOG}"

# 清理临时文件
rm -f "${SUCCESS_FILE}" "${FAIL_FILE}"

# 根据执行结果返回状态
if [ ${FAIL_COUNT} -eq 0 ]; then
    echo "所有target数据库创建成功！"
    exit 0
else
    echo "有 ${FAIL_COUNT} 个数据库创建失败，请检查日志"
    exit 1
fi
