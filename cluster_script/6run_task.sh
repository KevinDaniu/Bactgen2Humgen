#!/bin/bash
# submit_in_batches_smart.sh
# 智能分批提交：检测任务是否已完成，动态调整等待时间
# 支持断电恢复：启动时快速扫描，跳过已完成的批次

TASK_SPLITS_DIR="task_splits"
LOG_DIR="submission_logs"
mkdir -p "$LOG_DIR"

# 配置参数
BATCH_SIZE=10
MAX_INTERVAL=7200   # 最大等待2小时:7200
MIN_INTERVAL=300    # 最小等待5分钟
CHECK_INTERVAL=60   # 每60秒检查一次
# 如果未设置 OUTPUT_BASE，尝试从环境获取或设置默认值
OUTPUT_BASE="/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene"

# 颜色输出（可选）
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "========================================="
echo -e "${BLUE}智能批量任务提交脚本 (断电恢复版)${NC}"
echo "总任务数: $TOTAL_TASKS"
echo "每批提交: $BATCH_SIZE 个任务"
echo "最大等待间隔: 2小时"
echo "最小等待间隔: 5分钟"
echo "开始时间: $(date)"
echo "========================================="

# ============================================
# 核心函数：检查任务文件是否全部完成
# ============================================
check_task_completion() {
    local task_file=$1
    local completed_count=0
    local total_count=0
    
    # 读取任务列表
    mapfile -t genomes < "$task_file"
    
    for genome_name in "${genomes[@]}"; do
        [ -z "$genome_name" ] && continue
        total_count=$((total_count + 1))
        
        # 检查是否已完成
        result_dir="${OUTPUT_BASE}/${genome_name}/gene_cluster/cluster_genes"
        done_file="${result_dir}/gene_cooccurrence_conserve.csv"
        
        if [ -f "$done_file" ] && [ -s "$done_file" ]; then
            completed_count=$((completed_count + 1))
        fi
    done
    
    # 如果所有任务都已完成，返回0（成功）
    if [ $total_count -eq 0 ]; then
        return 1
    elif [ $completed_count -eq $total_count ]; then
        return 0
    else
        return 1
    fi
}

# ============================================
# 获取任务的完成统计（用于智能等待）
# ============================================
get_completion_stats() {
    local task_file=$1
    local completed_count=0
    local total_count=0
    
    mapfile -t genomes < "$task_file"
    
    for genome_name in "${genomes[@]}"; do
        [ -z "$genome_name" ] && continue
        total_count=$((total_count + 1))
        
        result_dir="${OUTPUT_BASE}/${genome_name}/gene_cluster/cluster_genes"
        done_file="${result_dir}/gene_cooccurrence_conserve.csv"
        
        if [ -f "$done_file" ] && [ -s "$done_file" ]; then
            completed_count=$((completed_count + 1))
        fi
    done
    
    echo "${completed_count}:${total_count}"
}

# ============================================
# 计算等待时间（基于完成率）
# ============================================
calculate_wait_time() {
    local completion_rate=$1
    
    if (( $(echo "$completion_rate >= 0.9" | bc -l) )); then
        echo $MIN_INTERVAL
    elif (( $(echo "$completion_rate >= 0.7" | bc -l) )); then
        echo 1800
    elif (( $(echo "$completion_rate >= 0.5" | bc -l) )); then
        echo 3600
    else
        echo $MAX_INTERVAL
    fi
}

# ============================================
# 启动扫描：快速判断哪些批次已完成
# ============================================
echo ""
echo -e "${BLUE}[启动扫描] 正在检查已有任务完成情况...${NC}"

# 获取所有任务文件并排序
TASK_FILES=($(ls -1 ${TASK_SPLITS_DIR}/task_*.txt 2>/dev/null | sort -V))
TOTAL_TASKS=${#TASK_FILES[@]}

if [ $TOTAL_TASKS -eq 0 ]; then
    echo -e "${RED}错误: 未找到任务文件 ${TASK_SPLITS_DIR}/task_*.txt${NC}"
    exit 1
fi

# 计算总批次
TOTAL_BATCHES=$(( (TOTAL_TASKS + BATCH_SIZE - 1) / BATCH_SIZE ))

# 存储每个批次的完成状态
declare -A BATCH_COMPLETED
declare -A BATCH_COMPLETION_RATE

for ((i=0; i<TOTAL_TASKS; i+=BATCH_SIZE)); do
    BATCH_NUM=$((i/BATCH_SIZE + 1))
    END=$((i + BATCH_SIZE))
    [ $END -gt $TOTAL_TASKS ] && END=$TOTAL_TASKS
    
    total_in_batch=0
    completed_in_batch=0
    
    for ((j=i; j<END; j++)); do
        TASK_FILE="${TASK_FILES[$j]}"
        stats=$(get_completion_stats "$TASK_FILE")
        comp=$(echo $stats | cut -d: -f1)
        total=$(echo $stats | cut -d: -f2)
        completed_in_batch=$((completed_in_batch + comp))
        total_in_batch=$((total_in_batch + total))
    done
    
    if [ $total_in_batch -eq 0 ]; then
        rate=0
    else
        rate=$(echo "scale=2; $completed_in_batch / $total_in_batch" | bc)
    fi
    
    BATCH_COMPLETION_RATE[$BATCH_NUM]=$rate
    
    if [ "$completed_in_batch" -eq "$total_in_batch" ] && [ $total_in_batch -gt 0 ]; then
        BATCH_COMPLETED[$BATCH_NUM]=true
	percentage=$(echo "$rate * 100" | bc)
        echo -e "  ✅ 批次 $BATCH_NUM: 已完成 (完成率 ${percentage}%)"
    else
        BATCH_COMPLETED[$BATCH_NUM]=false
	percentage=$(echo "$rate * 100" | bc)
        echo -e "  ⏳ 批次 $BATCH_NUM: 已完成 (完成率 ${percentage}%)"
    fi
done

# 找到第一个未完成的批次
FIRST_INCOMPLETE=1
for ((i=1; i<=TOTAL_BATCHES; i++)); do
    if [ "${BATCH_COMPLETED[$i]}" = "false" ]; then
        FIRST_INCOMPLETE=$i
        break
    fi
    FIRST_INCOMPLETE=$((i + 1))
done

if [ $FIRST_INCOMPLETE -gt $TOTAL_BATCHES ]; then
    echo ""
    echo -e "${GREEN}🎉 所有批次均已完成！无需提交任何任务。${NC}"
    echo "结束时间: $(date)"
    exit 0
fi

echo ""
echo -e "${YELLOW}📌 将从批次 $FIRST_INCOMPLETE 开始继续提交${NC}"
echo "========================================="

# ============================================
# 主循环：从第一个未完成的批次开始
# ============================================
for ((i=(FIRST_INCOMPLETE-1)*BATCH_SIZE; i<TOTAL_TASKS; i+=BATCH_SIZE)); do
    BATCH_NUM=$((i/BATCH_SIZE + 1))
    END=$((i + BATCH_SIZE))
    [ $END -gt $TOTAL_TASKS ] && END=$TOTAL_TASKS
    
    echo ""
    echo "========================================="
    echo -e "${BLUE}[$(date)] 处理批次 $BATCH_NUM/$TOTAL_BATCHES (任务 $((i+1)) 到 $END)${NC}"
    echo "========================================="
    
    # 检查整批是否已完成
    if [ "${BATCH_COMPLETED[$BATCH_NUM]}" = "true" ]; then
        echo -e "${GREEN}✅ 批次 $BATCH_NUM 已完成，跳过${NC}"
        continue
    fi
    
    # 提交当前批次中未完成的任务
    submitted_count=0
    skipped_count=0
    
    for ((j=i; j<END; j++)); do
        TASK_FILE="${TASK_FILES[$j]}"
        TASK_NAME=$(basename "$TASK_FILE" .txt)
        
        # 检查单个任务是否已完成
        if check_task_completion "$TASK_FILE"; then
            echo -e "  ${GREEN}⏭ 已完成: $TASK_NAME${NC}"
            skipped_count=$((skipped_count + 1))
            continue
        fi
        
        # 提交任务
        sbatch --job-name="${TASK_NAME}" \
               --output="${LOG_DIR}/${TASK_NAME}.out" \
               --error="${LOG_DIR}/${TASK_NAME}.err" \
               --cpus-per-task=48 \
               --time=02:00:00 \
               --wrap="bash 5.4batch_run.sh ${TASK_FILE}"
        
        echo -e "  ${GREEN}✅ 已提交: $TASK_NAME${NC}"
        submitted_count=$((submitted_count + 1))
        sleep 1
    done
    
    echo "批次 $BATCH_NUM 完成: 提交 $submitted_count 个，跳过 $skipped_count 个"
    
    # 如果是最后一批，或者所有剩余批次都已完成，结束循环
    if [ $END -ge $TOTAL_TASKS ]; then
        echo -e "${GREEN}✅ 所有批次处理完毕！${NC}"
        break
    fi
    
    # ============================================
    # 等待下一批：检查下一批（或后续批次）的完成情况
    # ============================================
    # 先检查下一批是否已完成
    NEXT_BATCH_NUM=$((BATCH_NUM + 1))
    
    # 如果下一批已完成，跳过等待
    if [ "${BATCH_COMPLETED[$NEXT_BATCH_NUM]}" = "true" ]; then
        echo -e "${GREEN}✅ 下一批（批次 $NEXT_BATCH_NUM）已完成，无需等待${NC}"
        continue
    fi
    
    # 获取下一批的完成率（用于决定等待时间）
    NEXT_RATE=${BATCH_COMPLETION_RATE[$NEXT_BATCH_NUM]:-0}
    
    # 如果下一批有部分完成，根据完成率决定等待时间
    if (( $(echo "$NEXT_RATE > 0" | bc -l) )); then
        WAIT_TIME=$(calculate_wait_time "$NEXT_RATE")
        echo -e "${YELLOW}⏱ 下一批已完成 ${NEXT_RATE}%，等待 $((WAIT_TIME/60)) 分钟${NC}"
    else
        WAIT_TIME=$MAX_INTERVAL
        echo -e "${YELLOW}⏱ 下一批尚未开始，等待 2 小时${NC}"
    fi
    
    echo "下次提交时间: $(date -d "+${WAIT_TIME} seconds" "+%Y-%m-%d %H:%M:%S")"
    echo "（可以按 Ctrl+C 中止，已提交的任务会继续运行）"
    
    # 分段睡眠，便于中途取消，同时动态检查任务是否提前完成
    remaining=$WAIT_TIME
    while [ $remaining -gt 0 ]; do
        sleep $CHECK_INTERVAL
        remaining=$((remaining - CHECK_INTERVAL))
        
        # 每5分钟重新检查下一批是否已完成
        if [ $((remaining % 300)) -eq 0 ]; then
            # 重新检查下一批的完成状态
            i_next=$((NEXT_BATCH_NUM - 1))
            i_next_start=$((i_next * BATCH_SIZE))
            i_next_end=$((i_next_start + BATCH_SIZE))
            [ $i_next_end -gt $TOTAL_TASKS ] && i_next_end=$TOTAL_TASKS
            
            total_next=0
            completed_next=0
            for ((k=i_next_start; k<i_next_end; k++)); do
                TASK_FILE="${TASK_FILES[$k]}"
                stats=$(get_completion_stats "$TASK_FILE")
                comp=$(echo $stats | cut -d: -f1)
                total=$(echo $stats | cut -d: -f2)
                completed_next=$((completed_next + comp))
                total_next=$((total_next + total))
            done
            
            if [ $total_next -gt 0 ] && [ $completed_next -eq $total_next ]; then
                echo -e "${GREEN}✅ 下一批已提前完成，退出等待！${NC}"
                # 更新批次状态
                BATCH_COMPLETED[$NEXT_BATCH_NUM]=true
                BATCH_COMPLETION_RATE[$NEXT_BATCH_NUM]=1.0
                remaining=0
                break
            fi
            
            echo "  剩余等待时间: $((remaining / 60)) 分钟"
        fi
    done
done

echo ""
echo "========================================="
echo -e "${GREEN}🎉 所有任务已提交/处理完成！${NC}"
echo "结束时间: $(date)"
echo "========================================="

# 显示最终统计
echo ""
echo "📊 最终状态摘要:"
total_all=0
completed_all=0
for ((i=0; i<TOTAL_TASKS; i++)); do
    TASK_FILE="${TASK_FILES[$i]}"
    stats=$(get_completion_stats "$TASK_FILE")
    comp=$(echo $stats | cut -d: -f1)
    total=$(echo $stats | cut -d: -f2)
    completed_all=$((completed_all + comp))
    total_all=$((total_all + total))
done

if [ $total_all -gt 0 ]; then
    rate=$(echo "scale=2; $completed_all * 100 / $total_all" | bc)
    echo "  总基因组数: $total_all"
    echo "  已完成: $completed_all"
    echo "  完成率: ${rate}%"
fi
