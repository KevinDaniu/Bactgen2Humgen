#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ==================== 配置文件参数 ====================
params.input_dir = "/home/liuzhh/bacteria_cluster/260811_allgenome_ani"
params.output_dir = "/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result"
params.query_genome = ""   // 必须指定，如：Alteromonas_226_GCF_000213655.1
params.query_cluster = ""  // 必须指定，如：Conserved_150
params.help = false
params.script_dir = "../../../script" // 脚本路径

// ==================== 帮助信息 ====================
if (params.help) {
    println """
    使用方法: nextflow run bactgen4eggnog.nf --query_genome <Genome ID> --query_cluster <Cluster ID> [选项]

    必需参数:
        --query_genome         查询Genome ID
        --query_cluster        查询Cluster ID 

    可选参数:
        --input_dir           输入序列目录 (默认: ${params.input_dir})
        --output_dir          输出目录 (默认: ${params.output_dir})
        --script_dir          脚本目录 (默认: ${params.script_dir})
        --help                显示帮助信息
    """
    exit 0
}

// 检查必需参数
if (!params.query_cluster) {
    error "必须指定 --query_cluster 参数！使用 --help 查看帮助"
}
if (!params.query_genome) {
    error "必须指定 --query_genome 参数！使用 --help 查看帮助"
}

// ==================== 定义路径 ====================
params.input_dir = params.input_dir + "/" + params.query_genome + "/gene_cluster/gene2human/sequence_file"
params.output_dir = params.output_dir + "/" + params.query_genome
params.eggnog_db = "/home/liuzhh/bacteria_cluster/bactgen2humgen/pipeline/eggnog-mapper-data"

// ==================== 调试信息 ====================
println "=== 调试信息 ==="
println "输入目录: ${params.input_dir}"
println "输出目录: ${params.output_dir}"
println "查询序列: ${params.query_cluster}"

// 检查输入文件是否存在
input_files = files("${params.input_dir}/*.fa")
println "找到FASTA文件数量: ${input_files.size()}"
if (input_files) {
    input_files.each { println "  - ${it.name}" }
} else {
    println "警告: 没有找到FASTA文件!"
}

// ==================== Channel定义 ====================
// 输入FASTA文件通道 - 使用ifEmpty防止空通道
Channel.fromPath("${params.input_dir}/${params.query_cluster}.fa")
    .ifEmpty { error "没有找到序列文件在: ${params.input_dir}" }
    .set { fa_file_ch }
// 创建fna目录通道 - 使用绝对路径
fa_dir_ch = Channel.value(file(params.input_dir))

// ==================== 流程模块 ====================

// eggNOG注释模块
process run_eggnog {
    tag "注释序列文件: ${params.query_cluster}"

    input:
    path cluster_dir
    path cluster_file

    output:
    path "eggnog_done.txt", emit: done_marker

    script:
    """
    set -e
    echo "========================================="
    echo "开始注释序列文件: \$(date)"
    echo "cluster_file: ${cluster_file}"
    echo "当前目录: \$(pwd)"

    # 创建输出目录
    mkdir -p ${params.output_dir}

    # 运行脚本
    echo "运行: bash ${params.script_dir}/2eggNOG4specific.sh ${params.query_cluster}"
    bash ${params.script_dir}/7eggNOG_anno/2eggNOG4specific.sh \
        ${cluster_dir} \
        ${cluster_file} \
        ${params.output_dir} \
        ${params.eggnog_db}

    # 创建完成标记
    echo "注释完成: \$(date)" > eggnog_done.txt
    echo "========================================="
    """
}

// ==================== 工作流定义 ====================

workflow {
    println "=== 开始执行流程 ==="
    println "查询序列: ${params.query_cluster}"

    // 注释基因
    eggnog_result = run_eggnog(fa_dir_ch, fa_file_ch)
    eggnog_result.done_marker.view { "✓ eggNOG 注释完成" }

    println "=== 所有进程已提交 ==="
}

// ==================== 完成通知 ====================

workflow.onComplete {
    println "========================================"
    println "流程执行完成!"
    println "完成时间: ${workflow.complete}"
    println "执行状态: ${workflow.success ? '成功' : '失败'}"
    println "输出目录: ${params.output_dir}"
    println "========================================"
}
