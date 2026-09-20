import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description='Parse MMseqs2 cluster results')
parser.add_argument('--cluster_dir',
        default="../../result/gene2human/mmseqs/human_clusters",
        help='Directory containing cluster files')
parser.add_argument('--output',
        default="../../result/gene2human/mmseqs/mmseqs_result.csv",
        help='Output file for cluster mapping')
parser.add_argument('--summary',
        default="../../result/gene2human/mmseqs/mmseqs_summary.csv",
        help='Output file for cluster summary (optional)')
args = parser.parse_args()

def parse_cluster_file(cluster_file):
    """
    解析单个 MMseqs2 聚类文件
    返回: (gene_to_cluster, cluster_members)
    """
    gene_to_cluster = {}
    cluster_members = defaultdict(list)
    
    with open(cluster_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                rep = parts[0]
                member = parts[1]
                gene_to_cluster[member] = rep
                cluster_members[rep].append(member)
    
    return gene_to_cluster, cluster_members

def extract_gene_name(protein_id):
    """从蛋白ID中提取基因名称"""
    # 格式: ENSP00000387779.2|...|RGPD6-203|RGPD6|333
    parts = protein_id.split('|')
    if len(parts) >= 7:
        # 倒数第二个是基因名称
        return parts[-2]
    # 如果格式不对，返回原始ID
    return protein_id

cluster_dir = Path(args.cluster_dir)    
print(f"Scanning cluster files in: {cluster_dir}")
    
# 收集所有聚类结果
all_gene_to_cluster = {}
all_cluster_members = {}
all_cluster_summary = []
    
cluster_files = list(cluster_dir.glob("*_cluster_cluster.tsv"))
print(f"Found {len(cluster_files)} cluster files")
    
for cluster_file in cluster_files:
    conserved_id = cluster_file.stem.replace('_cluster_cluster', '')
    print(f"  Processing: {conserved_id}")

    gene_to_cluster, cluster_members = parse_cluster_file(cluster_file)        
    # 添加 conserved_id 前缀，避免不同文件中的相同ID冲突
    for gene, rep in gene_to_cluster.items():
        prefixed_gene = f"{conserved_id}|{gene}"
        prefixed_rep = f"{conserved_id}|{rep}"
        all_gene_to_cluster[prefixed_gene] = prefixed_rep

    for rep, members in cluster_members.items():
        prefixed_rep = f"{conserved_id}|{rep}"
        prefixed_members = [f"{conserved_id}|{m}" for m in members]
        all_cluster_members[prefixed_rep] = prefixed_members
        
    # 提取基因名称用于摘要
    member_genes = [extract_gene_name(m) for m in members]
    all_cluster_summary.append({
        'Conserved_ID': conserved_id,
        'Cluster_Representative': rep,
        'Cluster_Size': len(members),
        'Member_Genes': ','.join(sorted(set(member_genes))),
        'All_Members': ','.join(members)
        })
    
print(f"\nTotal genes in clusters: {len(all_gene_to_cluster)}")
print(f"Total clusters: {len(all_cluster_members)}")
    
# 保存基因到簇的映射文件
mapping_df = pd.DataFrame([
    {'Gene_ID': gene, 'Cluster_Representative': rep}
    for gene, rep in all_gene_to_cluster.items()
    ])
mapping_df.to_csv(args.output, sep='\t', index=False)
print(f"\nGene to cluster mapping saved to: {args.output}")
    
summary_df = pd.DataFrame(all_cluster_summary)
summary_df.to_csv(args.summary, sep='\t', index=False)
print(f"Cluster summary saved to: {args.summary}")

# 显示统计
print("\n" + "="*60)
print("Summary Statistics:")
print("="*60)
print(f"Total cluster files: {len(cluster_files)}")
print(f"Total clusters: {len(all_cluster_members)}")
print(f"Total genes clustered: {len(all_gene_to_cluster)}")
    
# 显示前10个最大的簇
if all_cluster_members:
    print("\nTop 10 largest clusters:")
    sorted_clusters = sorted(all_cluster_members.items(), 
            key=lambda x: len(x[1]), reverse=True)[:10]
    for rep, members in sorted_clusters:
        conserved_id = rep.split('|')[0]
        rep_gene = extract_gene_name(rep.split('|')[1])
        print(f"  {conserved_id} | {rep_gene}: {len(members)} members")
