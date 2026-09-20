#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ==================== 配置文件参数 ====================
params.base_dir = "/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene"
params.query_genome = ""   // 必须指定，如：Alteromonas_226_GCF_000213655.1
params.help = false
params.script_dir = "../../../script" // 脚本路径

// ==================== 帮助信息 ====================
if (params.help) {
    println """
    使用方法: nextflow run bactgen4blast.nf --query_genome <Genome ID> --run_mmseqs true

    必需参数:
    --query_genome        查询基因组ID (例如: GCF_000009705.1)

    可选参数:
        --base_dir           输入GFF目录 (默认: ${params.input_dir})
        --run_mmseqs          是否运行mmseqs注释（默认：false）
        --script_dir          脚本目录 (默认: ${params.script_dir})
        --homo_pattern        同源基因匹配模式（默认：strict，可选strict, median, sensitive）
        --help                显示帮助信息
    """
    exit 0
 
}

// 检查必需参数
if (!params.query_genome) {
    error "必须指定 --query_genome 参数！使用 --help 查看帮助"
}

// ==================== 定义路径 ====================
params.input_dir = "${params.base_dir}/${params.query_genome}/gene_cluster"
params.output_dir = "${params.base_dir}/${params.query_genome}/gene_cluster"
params.sequence_output = "${params.output_dir}/gene2human/sequence_file"
params.blast_db = "/home/liuzhh/bacteria_cluster/bactgen2humgen/pipeline/blastdb_human/transcript_db"
params.blast_output = "${params.output_dir}/gene2human/blast_result"
params.blastp_db = "/home/liuzhh/bacteria_cluster/bactgen2humgen/pipeline/blastpdb_human/translation_db"
params.blastp_output = "${params.output_dir}/gene2human/blastp_result"
params.mmseqs_dir = "${params.output_dir}/gene2human/mmseqs"
params.human_prot = "/home/liuzhh/bacteria_cluster/bactgen2humgen/rawdata/gencode.v49.pc_translations.fa"
params.mmseqs_output = "${params.output_dir}/gene2human/mmseqs/human_clusters"
params.cluster_res = "${params.output_dir}/gene2human/mmseqs/mmseqs_result.csv"
params.cluster_sum = "${params.output_dir}/gene2human/mmseqs/mmseqs_summary.csv"
params.final_cluster = "${params.output_dir}/gene2human/final_result/final_homogroup_results.csv"
params.final_output = "${params.output_dir}/gene2human/final_result"

// ==================== 调试信息 ====================
println "=== 调试信息 ==="
println "输入目录: ${params.input_dir}"
println "输出目录: ${params.output_dir}"
println "查询基因组: ${params.query_genome}"

// ==================== Channel定义 ====================
// 创建fna目录通道 - 使用绝对路径
fa_dir_ch = Channel.value(file(params.input_dir))

// ==================== 流程模块 ====================

// BLAST序列准备
process prepare_blast_sequences {
    tag "准备BLAST序列: ${query_id}"
    
    input:
    val query_id

    output:
    val query_id, emit: query_id_out
    path "blast_prep_done.txt", emit: done_marker

    script:
    """
    set -e
    echo "========================================="
    echo "开始准备BLAST序列: \$(date)"

    # 创建输出目录
    mkdir -p ${params.sequence_output}

    echo "运行2extract_sequence.py"
    python ${params.script_dir}/6blast_human_gene/2extract_sequence.py \
        --input_dir ${params.input_dir} \
        --cds_file ${query_id} \
        --output_dir ${params.sequence_output}

    # 创建完成标记
    echo "BLAST序列准备完成: \$(date)" > blast_prep_done.txt
    echo "========================================="
    """
}

// BLAST比对
process blast_against_human {
        tag "BLAST比对: ${query_id}"

        input:
        val query_id
        path done_marker

        output:
        val query_id, emit: query_id_out
        path "blast_done.txt", emit: done_marker

        script:
        """
        set -e
        echo "========================================="
        echo "开始BLAST比对: \$(date)"

        # 创建输出目录
        mkdir -p ${params.blast_output}

        echo "运行3blast_sequence.sh"
        bash ${params.script_dir}/6blast_human_gene/3blast_sequence.sh \
            ${params.sequence_output} \
            ${params.blast_output} ${params.blast_db}

        echo "运行4extract_gene.py"
        python ${params.script_dir}/6blast_human_gene/4extract_gene.py \
            --blast_dir ${params.blast_output} \
            --output "${params.blast_output}/homology_summary.csv"

        echo "运行5blastp_sequence.sh"
        bash ${params.script_dir}/6blast_human_gene/5blastp_sequence.sh \
            ${params.sequence_output} \
            ${params.blastp_output} ${params.blastp_db}

        echo "运行6extract_prot.py"
        python ${params.script_dir}/6blast_human_gene/6extract_prot.py \
            --blast_dir ${params.blastp_output} \
            --output "${params.blastp_output}/homology_summary.csv"

        # 创建完成标记
        echo "BLAST比对完成: \$(date)" > blast_done.txt
        echo "========================================="
        """
    }

// MMseqs2聚类基因家族
if (params.run_mmseqs) {
    process run_mmseqs_gene {
        tag "进行聚类基因: ${query_id}"

        input:
        val query_id
        path done_marker

        output:
        path "mmseqs2_done.txt", emit: done_marker

        script:
        """
        set -e
        echo "========================================="
        echo "开始MMseqs2基因聚类: \$(date)"

        bash ${params.script_dir}/6blast_human_gene/7extract_human.sh \
            ${params.blastp_output} ${params.mmseqs_dir} ${params.human_prot}
        bash ${params.script_dir}/6blast_human_gene/8.1mmseqs_human.sh \
            ${params.mmseqs_dir}/human_sequences ${params.mmseqs_output}
        python ${params.script_dir}/6blast_human_gene/8.2extract_mmseqs.py \
            --cluster_dir ${params.mmseqs_output} \
            --output ${params.cluster_res} \
            --summary ${params.cluster_sum}
        python ${params.script_dir}/6blast_human_gene/9add_info.py \
            --blast_result ${params.final_output}/final_homology_results.csv \
            --cluster_map ${params.cluster_res} \
            --cluster_summary ${params.cluster_sum} \
            --output ${params.final_cluster}

        echo "基因聚类Group增添完成: \$(date)" > mmseqs2_done.txt
        echo "========================================="
        """
    }
}

// BLAST结果整合
process integrate_results_noeggnog_withblast {
    tag "整合最终结果: ${query_id}"

    input:
    val query_id
    path done_marker

    output:
    path "intergrate_done.txt", emit: final_result

    script:
    """
    set -e
    echo "========================================="
    echo "开始整合最终结果: \$(date)"

    # 创建输出目录
    mkdir -p ${params.final_output}

    python ${params.script_dir}/8product_result.py \
        --input_dir ${params.output_dir}

    echo "最终结果整合完成: \$(date)" > intergrate_done.txt
    echo "========================================="
    """
}

// ==================== 工作流定义 ====================

workflow {
    println "=== 开始执行流程 ==="
    println "查询基因组: ${params.query_genome}"
    println "运行mmseqs2: ${params.run_mmseqs}"

    // 1. 准备BLAST序列
    blast_prep_result = prepare_blast_sequences(
        params.query_genome,
    )
    blast_prep_result.done_marker.view { "✓ prepare_blast_sequences 完成" }

    // 2. BLAST比对
    blast_result = blast_against_human(
        params.query_genome,
        blast_prep_result.done_marker
    )
    blast_result.done_marker.view { "✓ blast_against_human 完成" }

    final_result = integrate_results_noeggnog_withblast(
        params.query_genome,
        blast_result.done_marker
    )
    // 3. 输出最终结果
    final_result.final_result
        .view { "✓ 最终结果文件: ${it}" }

    // 4. 添加聚类信息
    if (params.run_mmseqs) {
        // 运行MMseqs2
        mmseqs_result = run_mmseqs_gene(
            params.query_genome,
            final_result.final_result
        )
        mmseqs_result.done_marker.view { "✓ run_mmseqs_cluster 完成" }
    }

    println "=== 所有进程已提交 ==="
}

// ==================== 完成通知 ====================

workflow.onComplete {
    println "========================================"
    println "流程执行完成!"
    println "完成时间: ${workflow.complete}"
    println "执行状态: ${workflow.success ? '成功' : '失败'}"
    println "========================================"
}
