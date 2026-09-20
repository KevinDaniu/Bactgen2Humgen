#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ==================== 配置文件参数 ====================
params.query_genome = ""  // 必须指定，如：GCF_000009705.1
params.input_dir = ""
params.fna_dir = ""
params.gff_dir = ""
params.file_format = ""
params.output_dir = ""
params.run_blast = false  // 是否运行BLAST同源注释
params.run_eggnog = false  // 是否运行eggNOG注释
params.run_mmseqs = false // 是否运行mmseqs聚类
params.homo_pattern = "strict" // 同源基因匹配模式
params.max_parallel = 64    // 最大并行任务数
params.help = false
params.script_dir = "../../../script" // 脚本路径

// ==================== 帮助信息 ====================
if (params.help) {
    println """
    使用方法: nextflow run bactgen2humgen.nf --query_genome <基因组ID> [选项]
    
    参数:
        --query_genome        作为query的基因组ID
        --input_dir           输入工作目录
        --fna_dir             FNA文件目录
        --gff_dir             GFF文件目录
        --file_format         文件后缀
        --output_dir          输出目录
        --run_blast           是否运行BLAST同源注释（默认：false）
        --run_eggnog          是否运行eggNOG注释 (默认: false)
        --run_mmseqs          是否运行mmseqs注释（默认：false）
        --max_parallel        最大并行任务数 (默认: 64)
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

// ==================== 定义输出路径 ====================
params.bactdb = "${params.input_dir}/spacedust_db"
params.gff_res = "${params.input_dir}/all_cgff_files"
params.gff_meta = "${params.input_dir}/gff_meta"
params.spacedust_output = "${params.output_dir}/cluster_genes"
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
params.eggnog_db = "/home/liuzhh/bacteria_cluster/bactgen2humgen/pipeline/eggnog-mapper-data"
params.eggnog_output = "${params.output_dir}/gene2human/eggnog_result"
params.final_output = "${params.output_dir}/gene2human/final_result"

// ==================== 调试信息 ====================
println "=== 调试信息 ==="
println "输入目录: ${params.input_dir}"
println "FNA目录: ${params.fna_dir}"
println "GFF目录: ${params.gff_dir}"
println "输出目录: ${params.output_dir}"
println "查询基因组: ${params.query_genome}"

// ==================== Channel定义 ====================

// 创建fna目录通道 - 使用绝对路径
fna_dir_ch = Channel.value(file(params.fna_dir))
gff_dir_ch = Channel.value(file(params.gff_dir))

// ==================== 流程模块 ====================

// 1. Spacedust运行模块
process run_spacedust {
    tag "运行Spacedust分析: ${query_id}"
    
    input:
    val query_id
    
    output:
    val query_id, emit: query_id_out
    path "spacedust_done.txt", emit: done_marker
    
    script:
    """
    set -e
    echo "========================================="
    echo "开始Spacedust分析: \$(date)"
    echo "query_id: ${query_id}"
    
    mkdir -p ${params.spacedust_output}
    
    echo "运行1batchRun_spacedust.sh ${query_id} ${params.bactdb} ${params.spacedust_output}"
    bash ${params.script_dir}/3spacedust_result/1batchRun_spacedust.sh ${query_id} \
        ${params.bactdb} \
        ${params.spacedust_output} \
        ${params.max_parallel}
    
    # 创建完成标记
    echo "Spacedust分析完成: \$(date)" > spacedust_done.txt
    echo "========================================="
    """
}

// 2. 结果排序和统计模块
process order_and_stats {
    tag "排序和统计: ${query_id}"
    
    input:
    val query_id
    path done_marker
    
    output:
    val query_id, emit: query_id_out
    path "gene_count_simple.csv", emit: gene_count
    path "stats_done.txt", emit: done_marker
    
    script:
    """
    set -e
    echo "========================================="
    echo "开始排序和统计: \$(date)"
    
    
    # 排序结果
    echo "运行2.2batch_run.sh"
    bash ${params.script_dir}/3spacedust_result/2.2batch_run.sh \
        ${params.spacedust_output} \
        ${params.script_dir}/3spacedust_result \
        ${params.max_parallel}
    
    # 统计簇
    echo "运行3stat_cluster.py"
    python ${params.script_dir}/3spacedust_result/3stat_cluster.py \
        --input_dir ${params.spacedust_output}
    cp gene_count_simple.csv ${params.spacedust_output}/gene_count_simple.csv
    
    # 创建完成标记
    echo "排序和统计完成: \$(date)" > stats_done.txt
    echo "========================================="
    """
}

// 3. 提取保守基因簇
process extract_conserved_clusters {
    tag "提取保守基因簇: ${query_id}"
    
    input:
    val query_id
    path gene_count_file
    path done_marker
    
    output:
    val query_id, emit: query_id_out
    path "extract_done.txt", emit: done_marker
    
    script:
    """
    set -e
    echo "========================================="
    echo "开始提取保守基因簇: \$(date)"
    echo "输入文件: ${gene_count_file}"
    
    
    echo "运行4.1extract_cluster_***.py --input_file ${gene_count_file}"
    if [ "${params.homo_pattern}" == "strict" ]; then
        echo "运行4.1extract_cluster_strict.py --input_file ${gene_count_file}"
        python ${params.script_dir}/3spacedust_result/4.1extract_cluster_strict.py \
            --result_dir ${params.spacedust_output} \
            --input_file ${gene_count_file}
    elif [ "${params.homo_pattern}" == "median" ]; then
        echo "运行4.2extract_cluster_median.py --input_file ${gene_count_file}"
        python ${params.script_dir}/3spacedust_result/4.2extract_cluster_median.py \
            --result_dir ${params.spacedust_output} \
            --input_file ${gene_count_file}
    elif [ "${params.homo_pattern}" == "sensitive" ]; then
        echo "运行4.3extract_cluster_sensitive.py --input_file ${gene_count_file}"
        python ${params.script_dir}/3spacedust_result/4.3extract_cluster_sensitive.py \
            --result_dir ${params.spacedust_output} \
            --input_file ${gene_count_file}
    else
        echo "错误: 未知的 homo_pattern 参数: ${params.homo_pattern}"
        echo "请使用 strict, median 或 sensitive"
        exit 1
    fi
    
    # 创建完成标记
    echo "保守基因簇提取完成: \$(date)" > extract_done.txt
    echo "========================================="
    """
}

// 4. 用于BLAST
if (params.run_blast) {
    process prepare_blast_sequences {
        tag "准备BLAST序列: ${query_id}"
        
        input:
        val query_id
        path done_marker
        
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
            --input_dir ${params.output_dir} \
            --cds_file ${query_id} \
            --output_dir ${params.sequence_output}
        
        # 创建完成标记
        echo "BLAST序列准备完成: \$(date)" > blast_prep_done.txt
        echo "========================================="
        """
    }
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
}

// 5. eggNOG注释
if (params.run_eggnog) {
    process run_eggnog_annotation {
        tag "eggNOG注释: ${query_id}"
        
        input:
        val query_id
        path done_marker
        
        output:
        val query_id, emit: query_id_out
        path "eggnog_done.txt", emit: done_marker
        
        script:
        """
        set -e
        echo "========================================="
        echo "开始eggNOG注释: \$(date)"
        
        # 创建输出目录
        mkdir -p ${params.eggnog_output}
        
        # 运行eggNOG
        bash ${params.script_dir}/7eggNOG_anno/1.1eggNOG_diamond.sh ${params.sequence_output} ${params.eggnog_output} ${params.eggnog_db}
        
        # 创建完成标记
        echo "eggNOG注释完成: \$(date)" > eggnog_done.txt
        echo "========================================="
        """
    }
}

// 7. 最终结果整合（无eggNOG）
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

// 6. 最终结果整合（有eggNOG）
process integrate_results_witheggnog_withblast {
    tag "整合最终结果: ${query_id} (with eggNOG)"

    input:
    val query_id
    path eggnog_done_marker
    path blast_done_marker

    output:
    path "intergrate_done.txt", emit: final_result

    script:
    """
    set -e
    echo "========================================="
    echo "开始整合最终结果（含eggNOG）: \$(date)"
    
    # 创建输出目录
    mkdir -p ${params.final_output}
    
    
    echo "运行8product_result.py"
    python ${params.script_dir}/8product_result.py \
        --input_dir ${params.output_dir} \
        --eggnog_anno_dir ${params.eggnog_output} \
    
    echo "最终结果整合完成: \$(date)" > intergrate_done.txt
    echo "========================================="
    """
}

// 7. MMseqs2聚类基因家族
if (params.run_mmseqs && params.run_blast) {
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

// ==================== 工作流定义 ====================

workflow {
    println "=== 开始执行流程 ==="
    println "运行BLAST: ${params.run_blast}"
    println "运行eggNOG: ${params.run_eggnog}"
    println "运行mmseqs2: ${params.run_mmseqs}"
    
    // 1. 运行spacedust
    spacedust_result = run_spacedust(
        params.query_genome,
    )
    spacedust_result.done_marker.view { "✓ run_spacedust 完成" }
    
    // 2. 排序和统计
    order_result = order_and_stats(
        params.query_genome,
        spacedust_result.done_marker
    )
    order_result.done_marker.view { "✓ order_and_stats 完成" }
    
    // 3. 提取保守簇
    conserved_result = extract_conserved_clusters(
        params.query_genome,
        order_result.gene_count,
        order_result.done_marker
    )
    conserved_result.done_marker.view { "✓ extract_conserved_clusters 完成" }
    
    // 4. BLAST比对
    if (params.run_blast) {
        blast_prep_result = prepare_blast_sequences(
            params.query_genome,
            conserved_result.done_marker
        )
        blast_prep_result.done_marker.view { "✓ prepare_blast_sequences 完成" }
        
        blast_result = blast_against_human(
            params.query_genome,
            blast_prep_result.done_marker
        )
        blast_result.done_marker.view { "✓ blast_against_human 完成" }
    }
    
    // 5-6. eggNOG与整合结果
    if (params.run_eggnog) {
        if(params.run_blast){
            // 运行eggNOG
            eggnog_result = run_eggnog_annotation(
                params.query_genome,
                blast_prep_result.done_marker
            )
            eggnog_result.done_marker.view { "✓ run_eggnog_annotation 完成" }

            // 整合结果（含eggNOG）
            final_result = integrate_results_witheggnog_withblast(
                params.query_genome,
                eggnog_result.done_marker,
                blast_result.done_marker
            )
            // 输出最终结果
            final_result.final_result
                .view { "✓ 最终结果文件: ${it}" }
        }
    } else {
        if (params.run_blast) {
            // 整合结果（无eggNOG）
            final_result = integrate_results_noeggnog_withblast(
                params.query_genome,
                blast_result.done_marker
            )
            // 输出最终结果
            final_result.final_result
                .view { "✓ 最终结果文件: ${it}" }
         }
    }
    
    // 7. 添加聚类信息
    if (params.run_mmseqs && params.run_blast) {
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
