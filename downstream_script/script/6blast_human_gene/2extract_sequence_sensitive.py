import os
import sys
import argparse
from collections import defaultdict
import re

parser = argparse.ArgumentParser(description='Extract conserved gene cluster sequences')
parser.add_argument('--gene_file', type=str,
        default="gene_cooccurrence_conserve.csv",
        help='Input gene cluster CSV file')
parser.add_argument('--cds_file', type=str, required=True,
        help='Input CDS sequences file')
parser.add_argument('--output_dir', type=str,
        default='result/gene2human/sequence_file',
        help='Output directory for cluster FASTA files')
args = parser.parse_args()

def read_gene_clusters(gene_file):
    """
    读取基因簇CSV文件
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
                print(f"警告: 跳过格式不正确的行: {line}")
                continue
            cluster_id = parts[0]
            gene_str = parts[1]
            gene_count = int(parts[2])
            species_count = int(parts[3])
            if species_count > 0:
                # 解析基因列表
                genes = [g.strip() for g in gene_str.split('|')]
                clusters[cluster_id] = genes
                print(f"筛选到簇 {cluster_id}: {len(genes)} 个基因, {species_count} 个物种")
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

    print("正在读取CDS序列文件...")

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

    print(f"读取了 {len(cds_dict)} 条CDS序列索引")
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

def extract_cluster_sequences(clusters, cds_dict, output_dir):
    """
    提取每个簇的序列并写入文件
    """
    os.makedirs(output_dir, exist_ok=True)
    
    for cluster_id, genes in clusters.items():
        output_file = os.path.join(output_dir, f"{cluster_id}.fa")
        cluster_seq_count = 0
        found_genes = set()
        
        print(f"\n处理簇 {cluster_id}:")
        
        with open(output_file, 'w') as f_out:
            for gene_info in genes:
                gene_info = gene_info.strip()
                
                # 解析基因信息
                genome_id, pos1, pos2 = parse_gene_info(gene_info)
                
                if genome_id and pos1 and pos2:
                    # 尝试两种可能的键 (考虑位置顺序)
                    key1 = (genome_id, pos1, pos2)
                    key2 = (genome_id, pos2, pos1)
                    print(key1)
                    
                    if key1 in cds_dict:
                        sequence = cds_dict[key1]
                        header = f">{cluster_id}|{gene_info}"
                        f_out.write(f"{header}\n{sequence}\n")
                        cluster_seq_count += 1
                        found_genes.add(gene_info)
                        print(f"  ✓ 匹配: {gene_info}")
                    elif key2 in cds_dict:
                        sequence = cds_dict[key2]
                        header = f">{cluster_id}|{gene_info}"
                        f_out.write(f"{header}\n{sequence}\n")
                        cluster_seq_count += 1
                        found_genes.add(gene_info)
                        print(f"  ✓ 匹配 (反向): {gene_info}")
                    else:
                        print(f"  ✗ 未找到: {gene_info}")
                        # 调试信息：显示可用的键
                        # available = [k for k in cds_dict.keys() if genome_id in k]
                        # if available:
                        #     print(f"    可用键示例: {available[:3]}")
                else:
                    print(f"  ✗ 解析失败: {gene_info}")
        
        if cluster_seq_count > 0:
            print(f"簇 {cluster_id}: 提取了 {cluster_seq_count}/{len(genes)} 条序列")
            if cluster_seq_count < len(genes):
                missing = set(genes) - found_genes
                print(f"  缺失的基因: {missing}")
        else:
            print(f"警告: 簇 {cluster_id} 没有提取到任何序列")
            if os.path.exists(output_file) and os.path.getsize(output_file) == 0:
                os.remove(output_file)

# 构建完整的文件路径
gene_file = f"../../../result/cluster_genes/{args.gene_file}"
cds_file = f"../../../result/gff_meta/{args.cds_file}_cds.fa"

print(f"读取基因簇文件: {gene_file}")
print(f"读取CDS序列文件: {cds_file}")
print(f"输出目录: {args.output_dir}")

# 读取基因簇信息
clusters = read_gene_clusters(gene_file)
print(f"\n筛选出 {len(clusters)} 个满足条件的基因簇")

if not clusters:
    print("没有找到满足条件的基因簇，程序退出")
    sys.exit(0)

# 读取CDS序列
cds_dict = read_cds_sequences(cds_file)

# 提取序列并写入文件
extract_cluster_sequences(clusters, cds_dict, args.output_dir)

print(f"\n处理完成！结果保存在: {args.output_dir}")
