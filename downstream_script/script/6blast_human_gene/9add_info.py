import argparse
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser(description='Add cluster group info for each criteria')
parser.add_argument('--blast_result',
        default="../../result/gene2human/final_result/final_homology_results.csv"
        )
parser.add_argument('--cluster_map',
        default="../../result/gene2human/mmseqs/mmseqs_result.csv"
        )
parser.add_argument('--cluster_summary',
        default="../../result/gene2human/mmseqs/mmseqs_summary.csv"
        )
parser.add_argument('--output',
        default="../../result/gene2human/final_result/final_homogroup_results.csv"
        )
args = parser.parse_args()

df = pd.read_csv(args.blast_result)
cluster_df = pd.read_csv(args.cluster_map, sep='\t')
summary_df = pd.read_csv(args.cluster_summary, sep='\t')

print(f"Loaded BLAST result: {len(df)} rows")
print(f"Loaded cluster map: {len(cluster_df)} mappings")
print(f"Loaded cluster summary: {len(summary_df)} clusters")

def extract_gene_name(protein_id):
    """从完整蛋白ID中提取基因名称"""
    parts = protein_id.split('|')
    # 基因名称通常在倒数第二个位置
    if len(parts) >= 2:
        return parts[-2]
    return protein_id

# 从 cluster_summary 直接获取每个簇的成员基因名称
# 建立基因名称 -> 簇内所有基因名称 的映射
gene_to_group_members = {}

for _, row in summary_df.iterrows():
    conserved_id = row['Conserved_ID']
    member_genes = row['Member_Genes']  # 已经是逗号分隔的基因名称，如 "RGPD6,RGPD8"
    all_members = row['All_Members']  # 完整的成员ID列表
    
    # 解析所有成员，提取基因名称
    member_names = set()
    for member in all_members.split(','):
        gene_name = extract_gene_name(member)
        member_names.add(gene_name)
    
    # 为每个成员基因建立映射（包括代表序列）
    for member in all_members.split(','):
        gene_name = extract_gene_name(member)
        if gene_name not in gene_to_group_members:
            gene_to_group_members[gene_name] = set()
        gene_to_group_members[gene_name].update(member_names)

print(f"\n建立了 {len(gene_to_group_members)} 个基因的组映射")

# 显示示例
print("\n=== 映射示例 ===")
for i, (gene, members) in enumerate(list(gene_to_group_members.items())[:10]):
    print(f"  {gene} -> {','.join(sorted(members))}")

# 定义四个指标
criteria_cols = [
    'Human Prot (by Length)',
    'Human Prot (by E-value)', 
    'Human Prot (by Identity)',
    'Human Prot (by Comprehensive)'
]

# 为每个指标添加 Group 列（簇内所有基因名称）
for col in criteria_cols:
    if col in df.columns:
        group_col = f'Group_{col.replace("Human Prot (by ", "").replace(")", "")}'
        df[group_col] = ''
        
        matched = 0
        for idx, row in df.iterrows():
            human_gene = row[col]
            if pd.notna(human_gene) and human_gene != 'No hit' and human_gene != '-':
                if human_gene in gene_to_group_members:
                    matched += 1
                    # 将簇内所有基因名称连接
                    group_genes = '|'.join(sorted(gene_to_group_members[human_gene]))
                    df.at[idx, group_col] = group_genes
        
        print(f"{col}: {matched}/{len(df)} matched ({matched/len(df)*100:.1f}%)")

# 添加一列综合的 Group（取第一个非空的）
def get_any_group(row):
    for col in criteria_cols:
        group_col = f'Group_{col.replace("Human Prot (by ", "").replace(")", "")}'
        if group_col in row and row[group_col]:
            return row[group_col]
    return ''

df['Group_All'] = df.apply(get_any_group, axis=1)

# 保存
df.to_csv(args.output, index=False)
print(f"\nSaved to: {args.output}")

# 显示示例
print("\n=== 示例数据（前5行）===")
display_cols = ['Conserved ID', 'Cluster Gene', 'Human Prot (by Length)', 'Group_by Length']
if 'Group_by Length' in df.columns:
    print(df[display_cols].head(5).to_string())
