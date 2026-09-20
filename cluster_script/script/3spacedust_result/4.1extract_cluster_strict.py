import csv
from collections import defaultdict, Counter
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import sys

csv.field_size_limit(50 * 1024 * 1024)

parser = argparse.ArgumentParser()
parser.add_argument("--input_file", type=str)
parser.add_argument("--result_dir", type=str)
parser.add_argument("--workers", type=int, default=None, help="并行进程数，默认为CPU核心数")
args = parser.parse_args()

input_file = args.input_file
result_dir = args.result_dir
output_file = result_dir + "/gene_cooccurrence_conserve.csv"

max_workers = args.workers if args.workers else max(1, cpu_count() - 1)

print("=" * 80)
print(f"使用 {max_workers} 个进程并行处理")

# ============================================================
# 阶段1: 读取数据（单线程，I/O密集型）
# ============================================================
print("\n阶段1: 读取数据...")

species_cluster_genes = defaultdict(lambda: defaultdict(set))

with open(input_file, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene_id = row['基因ID']
        clusters = row['所在簇列表'].split(' | ')
        for cluster in clusters:
            species, cluster_num = cluster.split(':', 1)
            species_cluster_genes[species][cluster_num].add(gene_id)

gene_count = sum(len(clusters) for clusters in species_cluster_genes.values())
print(f"共发现 {gene_count} 个基因-簇关联")
print(f"共涉及 {len(species_cluster_genes)} 个物种")

# ============================================================
# 阶段2: 构建签名（并行）
# ============================================================
print("\n阶段2: 并行构建基因签名...")

def build_signatures_for_species(species_data_chunk):
    """
    处理一批物种数据，构建签名
    species_data_chunk: [(species, {cluster: set(genes)}), ...]
    """
    local_cluster_signatures = []
    
    for species, clusters in species_data_chunk:
        for cluster_num, genes in clusters.items():
            if len(genes) >= 2:
                local_cluster_signatures.append((species, cluster_num, frozenset(genes)))
    
    return local_cluster_signatures

# 将物种数据分块
species_items = list(species_cluster_genes.items())
chunk_size = max(1, len(species_items) // max_workers)
chunks = []
for i in range(0, len(species_items), chunk_size):
    chunks.append(species_items[i:i+chunk_size])

print(f"将 {len(species_items)} 个物种分为 {len(chunks)} 批处理")

cluster_signatures = []
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(build_signatures_for_species, chunk) for chunk in chunks]
    
    completed = 0
    for future in as_completed(futures):
        completed += 1
        print(f"  构建进度: {completed}/{len(chunks)}")
        cluster_signatures.extend(future.result())

print(f"共 {len(cluster_signatures)} 个候选簇签名")

# ============================================================
# 阶段3: 统计签名出现次数（并行）
# ============================================================
print("\n阶段3: 并行统计签名出现次数...")

def count_signatures_batch(signature_batch):
    """统计一批签名的出现次数"""
    local_counter = Counter()
    for sig in signature_batch:
        local_counter[sig] += 1
    return local_counter

# 提取所有签名（去掉物种信息）
signatures_only = [sig for _, _, sig in cluster_signatures]

# 分块
chunk_size = max(1, len(signatures_only) // max_workers)
sig_chunks = []
for i in range(0, len(signatures_only), chunk_size):
    sig_chunks.append(signatures_only[i:i+chunk_size])

print(f"将 {len(signatures_only)} 个签名分为 {len(sig_chunks)} 批统计")

signature_count = Counter()
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(count_signatures_batch, chunk) for chunk in sig_chunks]
    
    completed = 0
    for future in as_completed(futures):
        completed += 1
        print(f"  统计进度: {completed}/{len(sig_chunks)}")
        signature_count.update(future.result())

print(f"共 {len(signature_count)} 个唯一基因签名")

# ============================================================
# 阶段4: 构建签名到物种簇的映射（并行）
# ============================================================
print("\n阶段4: 并行构建签名-物种簇映射...")

def build_species_map_batch(signature_batch_with_info):
    """
    构建签名到物种簇的映射
    signature_batch_with_info: [(species, cluster_num, frozenset(genes)), ...]
    """
    local_map = defaultdict(list)
    for species, cluster_num, sig in signature_batch_with_info:
        local_map[sig].append((species, cluster_num))
    return local_map

# 分块处理cluster_signatures
chunk_size = max(1, len(cluster_signatures) // max_workers)
sig_with_info_chunks = []
for i in range(0, len(cluster_signatures), chunk_size):
    sig_with_info_chunks.append(cluster_signatures[i:i+chunk_size])

print(f"将 {len(cluster_signatures)} 个签名分为 {len(sig_with_info_chunks)} 批构建映射")

sig_to_species = defaultdict(list)
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(build_species_map_batch, chunk) 
               for chunk in sig_with_info_chunks]
    
    completed = 0
    for future in as_completed(futures):
        completed += 1
        print(f"  映射进度: {completed}/{len(sig_chunks)}")
        local_map = future.result()
        for sig, species_list in local_map.items():
            sig_to_species[sig].extend(species_list)

print(f"构建完成，共 {len(sig_to_species)} 个唯一签名映射")

# ============================================================
# 阶段5: 过滤和去重（并行）
# ============================================================
print("\n阶段5: 并行过滤和去重...")

# 获取所有物种
unique_species_in_data = {sig[0] for sig in cluster_signatures}
species_count = len(unique_species_in_data)
print(f"总物种数: {species_count}")

# 过滤条件
def filter_and_prepare_batch(batch_items):
    """
    过滤一批签名并准备数据
    batch_items: [(sig, species_list), ...]
    """
    results = []
    for sig, species_list in batch_items:
        if len(sig) < 2:
            continue
        
        unique_species = {s for s, _ in species_list}
        
        # 判断条件
        if species_count >= 3:
            if len(unique_species) < 2:
                continue
        
        results.append({
            'genes': set(sig),
            'species': species_list,
            'gene_count': len(sig),
            'species_count': len(unique_species),
        })
    
    return results

# 准备数据
sig_items = list(sig_to_species.items())
chunk_size = max(1, len(sig_items) // max_workers)
filter_chunks = []
for i in range(0, len(sig_items), chunk_size):
    filter_chunks.append(sig_items[i:i+chunk_size])

print(f"将 {len(sig_items)} 个签名分为 {len(filter_chunks)} 批过滤")

conserved_groups = []
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(filter_and_prepare_batch, chunk) for chunk in filter_chunks]
    
    completed = 0
    for future in as_completed(futures):
        completed += 1
        print(f"  过滤进度: {completed}/{len(filter_chunks)}")
        conserved_groups.extend(future.result())

print(f"初步发现 {len(conserved_groups)} 个保守组")

# ============================================================
# 阶段6: 去重（处理包含关系）- 这个必须串行
# ============================================================
print("\n阶段6: 去重（移除被包含的组）...")

# 按基因数量降序排序
conserved_groups.sort(key=lambda x: (-x['gene_count'], -x['species_count']))

unique_clusters = []
for group in conserved_groups:
    group_genes = group['genes']
    is_subset = False
    
    for other in unique_clusters:
        if group_genes.issubset(other['genes']):
            is_subset = True
            break
    
    if not is_subset:
        unique_clusters.append(group)

print(f"去重后剩余 {len(unique_clusters)} 个保守组")

# 转换set为list以便输出
for cluster in unique_clusters:
    cluster['genes'] = sorted(list(cluster['genes']))

# 按基因数量和物种数量排序
unique_clusters.sort(key=lambda x: (-x['gene_count'], -x['species_count']))

# ============================================================
# 阶段7: 输出结果
# ============================================================
print("\n阶段7: 输出结果...")

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

# 显示结果
print("\n保守基因簇 TOP 10:")
print("-" * 80)
for i, cluster in enumerate(unique_clusters[:10], 1):
    genes_str = ', '.join(cluster['genes'][:5]) + ('...' if len(cluster['genes'])>5 else '')
    species_str = ', '.join([f'{s}:{c}' for s, c in cluster['species'][:3]])
    if len(cluster['species']) > 3:
        species_str += '...'
    print(f"Conserved_{i:03d}: {cluster['gene_count']}个基因, 在{cluster['species_count']}个物种中保守")
    print(f"  基因: {genes_str}")
    print(f"  物种-簇: {species_str}")
    print()

print("=" * 80)
