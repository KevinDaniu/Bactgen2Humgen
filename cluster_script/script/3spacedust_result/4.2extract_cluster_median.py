import csv
from collections import defaultdict
from itertools import combinations
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", type=str)
parser.add_argument("--result_dir", type=str)
args = parser.parse_args()

input_file = args.input_file
result_dir = args.result_dir
output_file = result_dir + "/gene_cooccurrence_conserve.csv"

print("正在解析基因-簇关系...")
species_cluster_genes = defaultdict(lambda: defaultdict(set))
gene_species_clusters = defaultdict(lambda: defaultdict(list))

with open(input_file, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene_id = row['基因ID']
        clusters = row['所在簇列表'].split(' | ')
        for cluster in clusters:
            species = cluster.split(':')[0]
            cluster_num = cluster.split(':')[1]
            species_cluster_genes[species][cluster_num].add(gene_id)
            gene_species_clusters[gene_id][species].append(cluster_num)

print(f"共发现 {len(gene_species_clusters)} 个基因, {len(species_cluster_genes)} 个物种")

# 收集所有物种-簇（基因数>=2）
cluster_list = []  # [(species, cluster_num, gene_set), ...]
for species, clusters in species_cluster_genes.items():
    for cluster_num, genes in clusters.items():
        if len(genes) >= 2:
            cluster_list.append((species, cluster_num, genes))

print(f"共 {len(cluster_list)} 个候选簇")

# 找保守基因簇：允许包含关系，找"在多个物种中共同出现的基因组合"
conserved_groups = []

# 方法：对每一对物种-簇，找它们的共同基因，然后看其他物种是否有包含这些基因的簇
for i in range(len(cluster_list)):
    sp1, c1, genes1 = cluster_list[i]
    for j in range(i+1, len(cluster_list)):
        sp2, c2, genes2 = cluster_list[j]
        
        # 找共同基因（交集）
        common_genes = genes1 & genes2
        if len(common_genes) < 2:  # 至少2个共同基因
            continue
            
        # 检查其他物种是否有簇包含这些共同基因
        matching_species = [(sp1, c1), (sp2, c2)]
        
        for sp3, c3, genes3 in cluster_list:
            if sp3 in [sp1, sp2]:  # 跳过已检查的
                continue
            # 检查genes3是否包含common_genes
            if common_genes.issubset(genes3):
                matching_species.append((sp3, c3))
        
        # 如果至少在2个物种中找到（加上自己）
        if len(matching_species) >= 2:
            conserved_groups.append({
                'genes': sorted(list(common_genes)),
                'species': matching_species,
                'gene_count': len(common_genes),
                'species_count': len(matching_species)
            })

print(f"初步发现 {len(conserved_groups)} 个保守组")

# 去重：保留最大的（基因数最多，然后物种数最多）
conserved_groups.sort(key=lambda x: (-x['gene_count'], -x['species_count']))

unique_clusters = []
for group in conserved_groups:
    genes_set = set(group['genes'])
    # 检查是否已被更大的簇包含
    is_covered = False
    for existing in unique_clusters:
        if genes_set.issubset(set(existing['genes'])):
            # 还要检查物种是否也是子集
            existing_species = set(s for s, c in existing['species'])
            current_species = set(s for s, c in group['species'])
            if current_species.issubset(existing_species):
                is_covered = True
                break
    if not is_covered:
        unique_clusters.append(group)

# 最终排序
unique_clusters.sort(key=lambda x: (-x['gene_count'], -x['species_count']))

# 写入结果（同上）
with open(output_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(["保守簇ID", "基因列表", "基因数量", "保守物种数", "物种-簇对应关系"])
    for i, cluster in enumerate(unique_clusters, 1):
        species_info = " | ".join([f"{s}:{c}" for s, c in cluster['species']])
        writer.writerow([
            f"Conserved_{i:03d}",
            " | ".join(cluster['genes']),
            cluster['gene_count'],
            cluster['species_count'],
            species_info
        ])

print(f"\n结果保存到: {output_file}")
print(f"共发现 {len(unique_clusters)} 个保守基因簇")
