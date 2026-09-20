#!/bin/bash

QUERYDB=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
THREADS=$4

mkdir -p $OUTPUT_DIR

# 检查targetDB目录是否存在
if [ ! -d "${INPUT_DIR}" ]; then
    echo "错误：spacedustDB目录不存在"
    exit 1
fi

# 检查queryDB是否存在
if [ ! -d "${INPUT_DIR}/${QUERYDB}" ]; then
    echo "错误：queryDB '${QUERYDB}' 不存在"
    echo "请检查路径：${INPUT_DIR}/${QUERYDB}"
    exit 1
fi

# 创建日志文件
LOG_FILE="${OUTPUT_DIR}/batch_run_${QUERYDB}_$(date +%Y%m%d_%H%M%S).log"
echo "批量运行日志 - $(date)" > "$LOG_FILE"
echo "QueryDB: ${QUERYDB}" >> "$LOG_FILE"
echo "spacedustDB目录: ${INPUT_DIR}" >> "$LOG_FILE"
echo "使用并行核心数: 64" >> "$LOG_FILE"
echo "----------------------------------------" >> "$LOG_FILE"

# 获取所有targetDB列表
targets=()
for target in "${INPUT_DIR}"/*; do
    if [ -d "$target" ]; then
        TARGETDB=$(basename "$target")
        if [[ "$TARGETDB" != .* && "$TARGETDB" != "$QUERYDB" ]]; then
            targets+=("$TARGETDB")
        fi
    fi
done

total=${#targets[@]}
echo "发现 ${total} 个targetDB目录"
echo "总计需要处理: ${total} 个任务" >> "$LOG_FILE"

# 创建临时目录存放状态文件
TEMP_STATUS_DIR="${OUTPUT_DIR}/temp_status_${QUERYDB}_$$"
mkdir -p "$TEMP_STATUS_DIR"

# 定义单个任务函数
run_single_task() {
    local TARGETDB=$1
    local QUERYDB=$2
    local INPUT_DIR=$3
    local OUTPUT_DIR=$4
    local TEMP_STATUS_DIR=$5
    
    local FILENAME="${QUERYDB}_with_${TARGETDB}"
    local OUTPUT_FILE="${OUTPUT_DIR}/${FILENAME}/${FILENAME}_result.tsv"
    local STATUS_FILE="${TEMP_STATUS_DIR}/${TARGETDB}.status"
    
    # 创建输出目录
    mkdir -p "${OUTPUT_DIR}/${FILENAME}"
    
    # 创建临时目录（每个任务独立）
    local TMP_DIR="tmpFolder_${FILENAME}_$$_${RANDOM}"
    
    # 运行spacedust
    if spacedust clustersearch \
        "${INPUT_DIR}/${TARGETDB}/setDB" \
        "${INPUT_DIR}/${QUERYDB}/setDB" \
        "${OUTPUT_FILE}" \
        "${TMP_DIR}" 2> /dev/null; then
        
        if [ -f "${OUTPUT_FILE}" ]; then
            cluster_count=$(grep -c "^#" "${OUTPUT_FILE}" 2>/dev/null || echo "0")
            echo "SUCCESS:${TARGETDB}:${cluster_count}" > "$STATUS_FILE"
            rm -rf "$TMP_DIR"
            return 0
        else
            echo "FAILED:${TARGETDB}:no_output" > "$STATUS_FILE"
            rm -rf "$TMP_DIR"
            return 1
        fi
    else
        echo "FAILED:${TARGETDB}:execution_error" > "$STATUS_FILE"
        rm -rf "$TMP_DIR"
        return 1
    fi
}

# 导出函数和变量，以便parallel使用
export -f run_single_task
export QUERYDB INPUT_DIR OUTPUT_DIR TEMP_STATUS_DIR

# 使用parallel并行运行（64个核心）
echo "开始并行处理... (使用64个核心)"
echo "开始并行处理 - $(date)" >> "$LOG_FILE"

# 运行parallel，每个任务单独记录
printf '%s\n' "${targets[@]}" | \
    parallel --jobs $THREADS --bar \
        "run_single_task {} $QUERYDB $INPUT_DIR $OUTPUT_DIR $TEMP_STATUS_DIR"

# 收集结果
success=0
failed=0

echo "" >> "$LOG_FILE"
echo "任务完成情况:" >> "$LOG_FILE"

for status_file in "${TEMP_STATUS_DIR}"/*.status; do
    if [ -f "$status_file" ]; then
        status=$(cat "$status_file")
        if [[ $status == SUCCESS:* ]]; then
            IFS=':' read -r _ target cluster_count <<< "$status"
            echo "  ✓ ${target} - 簇数量: ${cluster_count}" >> "$LOG_FILE"
            ((success++))
        else
            IFS=':' read -r _ target error <<< "$status"
            echo "  ✗ ${target} - 错误: ${error}" >> "$LOG_FILE"
            ((failed++))
        fi
    fi
done

# 清理临时目录
rm -rf "$TEMP_STATUS_DIR"

# 输出统计信息
echo "========================================"
echo "批量运行完成！"
echo "QueryDB: ${QUERYDB}"
echo "总计处理: ${total} 个targetDB"
echo "成功: ${success}"
echo "失败: ${failed}"
echo "日志文件: ${LOG_FILE}"

echo "" >> "$LOG_FILE"
echo "运行统计 - $(date)" >> "$LOG_FILE"
echo "总计: ${total}, 成功: ${success}, 失败: ${failed}" >> "$LOG_FILE"
