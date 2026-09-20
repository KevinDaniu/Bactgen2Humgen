import os
import sys
import argparse
from collections import defaultdict
import re
import glob
import logging

parser = argparse.ArgumentParser(description='Extract conserved gene cluster sequences')
parser.add_argument('--input_dir', type=str)
parser.add_argument('--gene_file', type=str,
        default="gene_cooccurrence_conserve.csv",
        help='Input gene cluster CSV file')
parser.add_argument('--cds_file', type=str, required=True,
        help='Input CDS sequences file')
parser.add_argument('--output_dir', type=str,
        default='result/gene2human/sequence_file',
        help='Output directory for cluster FASTA files')
parser.add_argument('--min_species', type=int, default=5,
        help='Minimum number of conserved species required (default: 5)')
args = parser.parse_args()

input_dir = args.input_dir
cds_files = glob.glob(input_dir + "/gff_meta/*_cds.fa")
total_species = len(cds_files)
thresholds = {
    'core': total_species - 1,      # 核心：几乎全部物种（严格保守）
    'soft_core': int(total_species * 0.95),  # 软核心：95%以上
    'habitat_specific': int(total_species * 0.7),  # 生境特异：70%以上
    'accessory': int(total_species * 0.5)   # 辅助基因：50%以上
}
args.min_species =  thresholds["accessory"]

os.makedirs(args.output_dir, exist_ok=True)

log_file = os.path.join(args.output_dir, "extract_cluster.log")
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file, encoding='utf-8'),
        logging.StreamHandler(sys.stdout)
    ]
)

def read_gene_clusters(gene_file, min_species=5):
    """
    读取基因簇CSV文件，筛选满足最小物种数的簇
    返回: {cluster_id: [gene_list], ...}
    """
    clusters = {}
    n = 1
    with open(gene_file, 'r') as f:
        for line in f:
            # 跳过首行
            line = line.strip()
            if n == 1:
                n = n + 1
                continue
            # 解析CSV行
            parts = line.split(',')
            if len(parts) < 5:
                logging.warning(f"跳过格式不正确的行: {line}")
                continue
            cluster_id = parts[0]
            gene_str = parts[1]
            gene_count = int(parts[2])
            species_count = int(parts[3])
            # 筛选满足物种数的簇
            if species_count >= min_species:
                # 解析基因列表
                genes = [g.strip() for g in gene_str.split('|')]
                clusters[cluster_id] = genes
                logging.info(f"筛选到簇 {cluster_id}: {len(genes)} 个基因, {species_count} 个物种")
    return clusters

def parse_cds_header(header):
    """
    解析CDS FASTA文件的头部
    格式: >cds-WP_158299063.1::NZ_CP073017.1:826971-826726(-)
    返回: (genome_id, start, end, strand) 返回原始坐标和链信息
    """
    # 去掉开头的 '>'
    header = header[1:] if header.startswith('>') else header

    # 提取位置信息部分 (::后面的内容)
    if '::' in header:
        loc_part = header.split('::')[1]
    else:
        loc_part = header

    # 解析 NZ_CP073017.1:826971-826726(-)
    if ':' in loc_part:
        genome_id, pos_part = loc_part.split(':')
        # 提取位置数字和链信息
        match = re.search(r'(\d+)-(\d+)(?:\(([+-])\))?', pos_part)
        if match:
            start = int(match.group(1))
            end = int(match.group(2))
            strand = match.group(3) if match.group(3) else '+'
            return genome_id, start, end, strand

    return None, None, None, None

def read_cds_sequences(cds_file):
    """
    读取CDS序列文件，返回字典，键为(genome_id, adjusted_pos1, adjusted_pos2)
    根据链信息调整坐标：
    - 正链 (+): 用 (start-1, end) 作为键
    - 负链 (-): 用 (start, end+1) 作为键
    """
    cds_dict = {}
    current_header = None
    current_seq = []

    logging.info("正在读取CDS序列文件...")

    with open(cds_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header and current_seq:
                    genome_id, start, end, strand = parse_cds_header(current_header)
                    if genome_id and start and end:
                        sequence = ''.join(current_seq)

                        if strand == '+':
                            # 正链: start-1, end
                            key1 = (genome_id, start+1, end)
                            key2 = (genome_id, end, start+1)  # 反向也存
                            cds_dict[key1] = sequence
                            cds_dict[key2] = sequence
                        else:
                            # 负链: start, end+1
                            key1 = (genome_id, start+1, end)
                            key2 = (genome_id, end, start+1)  # 反向也存
                            cds_dict[key1] = sequence
                            cds_dict[key2] = sequence

                current_header = line
                current_seq = []
            else:
                if line:
                    current_seq.append(line)

        # 保存最后一个基因
        if current_header and current_seq:
            genome_id, start, end, strand = parse_cds_header(current_header)
            if genome_id and start and end:
                sequence = ''.join(current_seq)

                if strand == '+':
                    key1 = (genome_id, start-1, end)
                    key2 = (genome_id, end, start-1)
                    cds_dict[key1] = sequence
                    cds_dict[key2] = sequence
                else:
                    key1 = (genome_id, start, end+1)
                    key2 = (genome_id, end+1, start)
                    cds_dict[key1] = sequence
                    cds_dict[key2] = sequence

    logging.info(f"读取了 {len(cds_dict)} 条CDS序列索引")
    return cds_dict

def parse_gene_info(gene_info):
    """
    解析基因信息字符串
    格式: NZ_CP073017.1_768_826971_826726
    返回: (genome_id, pos1, pos2)
    """
    parts = gene_info.split('_')
    if len(parts) >= 4:
        # 前两部分组成genome_id
        genome_id = '_'.join(parts[:2])
        # 最后两个是位置信息
        pos1 = int(parts[-2])
        pos2 = int(parts[-1])
        return genome_id, pos1, pos2
    return None, None, None

def extract_cluster_sequences(clusters, cds_dict, output_dir, min_species=5):
    """
    提取每个簇的序列并写入文件
    """
    os.makedirs(output_dir, exist_ok=True)

    for cluster_id, genes in clusters.items():
        output_file = os.path.join(output_dir, f"{cluster_id}.fa")
        cluster_seq_count = 0
        found_genes = set()

        logging.info(f"处理簇 {cluster_id}:")

        with open(output_file, 'w') as f_out:
            for gene_info in genes:
                gene_info = gene_info.strip()

                # 解析基因信息
                genome_id, pos1, pos2 = parse_gene_info(gene_info)

                if genome_id and pos1 and pos2:
                    # 尝试两种可能的键 (考虑位置顺序)
                    key1 = (genome_id, pos1, pos2)
                    key2 = (genome_id, pos2, pos1)
                    logging.debug(f"尝试键: {key1}")

                    if key1 in cds_dict:
                        sequence = cds_dict[key1]
                        header = f">{cluster_id}|{gene_info}"
                        f_out.write(f"{header}\n{sequence}\n")
                        cluster_seq_count += 1
                        found_genes.add(gene_info)
                        logging.info(f"  ✓ 匹配: {gene_info}")
                    elif key2 in cds_dict:
                        sequence = cds_dict[key2]
                        header = f">{cluster_id}|{gene_info}"
                        f_out.write(f"{header}\n{sequence}\n")
                        cluster_seq_count += 1
                        found_genes.add(gene_info)
                        logging.info(f"  ✓ 匹配 (反向): {gene_info}")
                    else:
                        logging.warning(f"  ✗ 未找到: {gene_info}")
                else:
                    logging.warning(f"  ✗ 解析失败: {gene_info}")

        if cluster_seq_count > 0:
            logging.info(f"簇 {cluster_id}: 提取了 {cluster_seq_count}/{len(genes)} 条序列")
            if cluster_seq_count < len(genes):
                missing = set(genes) - found_genes
                logging.warning(f"  缺失的基因: {missing}")
        else:
            logging.warning(f"警告: 簇 {cluster_id} 没有提取到任何序列")
            if os.path.exists(output_file) and os.path.getsize(output_file) == 0:
                os.remove(output_file)

# 构建完整的文件路径
gene_file = f"{input_dir}/cluster_genes/{args.gene_file}"
#args.cds_file = "_".join(args.cds_file.split("_")[:2])
cds_file = f"{input_dir}/gff_meta/{args.cds_file}_cds.fa"

logging.info(f"读取基因簇文件: {gene_file}")
logging.info(f"读取CDS序列文件: {cds_file}")
logging.info(f"输出目录: {args.output_dir}")
logging.info(f"最小物种数要求: {args.min_species}")

# 读取基因簇信息
clusters = read_gene_clusters(gene_file, args.min_species)
logging.info(f"筛选出 {len(clusters)} 个满足条件的基因簇")

if not clusters:
    logging.warning("没有找到满足条件的基因簇，程序退出")
    sys.exit(0)

# 读取CDS序列
cds_dict = read_cds_sequences(cds_file)

# 提取序列并写入文件
extract_cluster_sequences(clusters, cds_dict, args.output_dir, args.min_species)

logging.info(f"处理完成！结果保存在: {args.output_dir}")
logging.info(f"日志文件: {log_file}")

