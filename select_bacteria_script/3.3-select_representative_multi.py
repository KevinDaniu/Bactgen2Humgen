#!/usr/bin/env python3
"""
Mash距离文件聚类 - 多核并行版（选择最长序列作为代表）
"""

import os
import sys
import numpy as np
from collections import defaultdict
import time
import gc
import multiprocessing as mp
from functools import partial

# ===== 强制实时输出 =====
sys.stdout.reconfigure(line_buffering=True)

def log(msg):
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {msg}")
    sys.stdout.flush()

# ===== 配置 =====
DIST_FILE = "../select_bacteria_result/mash_work/aai_distances.txt"
REP_LIST = "../select_bacteria_result/representative_list_0.23.txt"
CLUSTER_INFO = "../select_bacteria_result/cluster_info.txt"
THRESHOLD = 0.23
TOTAL_LINES = 484000000
NUM_CORES = 48  # 使用48核
ASSEMBLY_SUMMARY = "/home/liuzhh/bacteria_cluster/bacteria_genomes/20260201/assembly_summary.tsv"

log("=" * 60)
log("Mash距离文件聚类 - 多核并行版（选择最长序列）")
log("=" * 60)
log(f"总行数: {TOTAL_LINES:,}")
log(f"使用核心数: {NUM_CORES}")
start_time = time.time()

# ===== 检查文件 =====
if not os.path.exists(DIST_FILE):
    log(f"错误: 文件不存在!")
    sys.exit(1)

file_size = os.path.getsize(DIST_FILE)
log(f"文件大小: {file_size/(1024**3):.2f} GB")

# ===== 第一步：获取基因组列表 =====
log("\n步骤1: 扫描基因组名称...")
genomes = set()
line_count = 0

with open(DIST_FILE, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        g1 = parts[0].split('/')[-1]
        g2 = parts[1].split('/')[-1]
        genomes.add(g1)
        genomes.add(g2)
        
        line_count += 1
        if line_count % 5000000 == 0:
            elapsed = time.time() - start_time
            speed = line_count / elapsed
            log(f"  进度: {line_count/1e6:.1f}M, 速度: {speed/1000:.1f}K/秒, 发现 {len(genomes)} 个基因组")

genome_list = sorted(genomes)
n = len(genome_list)
genome_to_idx = {g: i for i, g in enumerate(genome_list)}
log(f"  基因组总数: {n:,}")

# ===== 第一步补充：获取基因组大小信息 =====
log("\n步骤1.5: 获取基因组大小信息...")

# 从assembly_summary.tsv读取基因组大小
genome_size = {}
if os.path.exists(ASSEMBLY_SUMMARY):
    log(f"  从 {ASSEMBLY_SUMMARY} 读取基因组大小...")
    with open(ASSEMBLY_SUMMARY, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 10:
                continue
            # assembly_accession 在第1列 (索引0)
            # genome_size 在第14列 (索引13) (根据之前看到的列)
            acc = parts[0].split('.')[0]  # 去掉版本号
            # 查找对应的fna文件名
            ftp_path = parts[19] if len(parts) > 19 else ''
            if ftp_path:
                base_name = ftp_path.split('/')[-1]
                fna_file = f"{base_name}_genomic.fna.gz"
                # 提取大小
                if len(parts) > 13 and parts[13].strip():
                    try:
                        size = int(parts[13])
                        genome_size[fna_file] = size
                    except:
                        pass
    log(f"  获取到 {len(genome_size)} 个基因组的大小信息")
else:
    log(f"  警告: 找不到 {ASSEMBLY_SUMMARY}，将使用文件大小作为代理")

# 如果无法从assembly_summary获取，使用文件大小作为代理
if len(genome_size) == 0:
    log("  使用文件大小作为代理...")
    FNA_LINK_DIR = "../select_bacteria_result/linked_genomes/fna"
    for genome in genome_list:
        fna_path = os.path.join(FNA_LINK_DIR, genome)
        if os.path.exists(fna_path):
            genome_size[genome] = os.path.getsize(fna_path)
    log(f"  获取到 {len(genome_size)} 个基因组的大小信息")

# ===== 第二步：分块并行处理 =====
log("\n步骤2: 并行处理文件...")

# 将文件分成48块
chunk_size = file_size // NUM_CORES
chunks = []
for i in range(NUM_CORES):
    start = i * chunk_size
    end = (i + 1) * chunk_size if i < NUM_CORES - 1 else file_size
    chunks.append((DIST_FILE, start, end, genome_to_idx, n))

def process_chunk(args):
    """处理文件的一个块"""
    filename, start_byte, end_byte, genome_to_idx, n = args
    # 使用字典存储距离值（避免大数组）
    dist_dict = {}
    line_count = 0
    
    with open(filename, 'r') as f:
        f.seek(start_byte)
        if start_byte > 0:
            f.readline()  # 跳过不完整的行
        
        while f.tell() < end_byte:
            line = f.readline()
            if not line:
                break
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            
            g1 = parts[0].split('/')[-1]
            g2 = parts[1].split('/')[-1]
            dist = float(parts[2])
            
            i = genome_to_idx.get(g1)
            j = genome_to_idx.get(g2)
            if i is None or j is None or i == j:
                continue
            
            # 存储上三角
            if i > j:
                i, j = j, i
            key = i * n + j  # 用整数作为key
            dist_dict[key] = dist
            
            line_count += 1
            if line_count % 1000000 == 0:
                # 每100万行输出一次（但不会显示在终端）
                pass
    
    return dist_dict

log(f"  启动 {len(chunks)} 个进程并行处理...")

# 并行处理
with mp.Pool(processes=NUM_CORES) as pool:
    # 使用 imap 显示进度
    results = []
    for i, result in enumerate(pool.imap_unordered(process_chunk, chunks)):
        results.append(result)
        log(f"  已完成块 {i+1}/{len(chunks)}")

# 合并结果
log("  合并结果...")
dist_dict = {}
for result in results:
    dist_dict.update(result)

log(f"  距离值数量: {len(dist_dict):,}")
del results
gc.collect()

# ===== 第三步：构建上三角数组 =====
log("\n步骤3: 构建上三角数组...")
tri_size = n * (n - 1) // 2
log(f"  上三角元素数: {tri_size:,}")
log(f"  预计内存: {tri_size*4/(1024**3):.2f} GB (float32)")

dist_values = np.zeros(tri_size, dtype=np.float32)

def get_triu_idx(i, j, n):
    if i > j:
        i, j = j, i
    return i * n - i * (i + 1) // 2 + (j - i - 1)

log("  填充数组...")
count = 0
for key, dist in dist_dict.items():
    i = key // n
    j = key % n
    idx = get_triu_idx(i, j, n)
    dist_values[idx] = dist
    count += 1
    if count % 10000000 == 0:
        log(f"    已填充 {count/1e6:.1f}M 个值")

del dist_dict
gc.collect()
log(f"  数组构建完成")

# ===== 第四步：聚类 =====
log("\n步骤4: 连通分量聚类...")

log("  构建邻接矩阵...")
adjacency = np.zeros((n, n), dtype=np.int8)

for i in range(n):
    for j in range(i+1, n):
        idx = get_triu_idx(i, j, n)
        if dist_values[idx] <= THRESHOLD:
            adjacency[i, j] = 1
            adjacency[j, i] = 1
    if i % 1000 == 0:
        log(f"    处理行 {i}/{n}")

del dist_values
gc.collect()
log(f"  邻接矩阵构建完成")

# BFS找连通分量
def find_connected_components(adj_matrix):
    n = adj_matrix.shape[0]
    visited = np.zeros(n, dtype=bool)
    labels = np.zeros(n, dtype=np.int32) - 1
    cluster_id = 0
    
    for i in range(n):
        if not visited[i]:
            queue = [i]
            visited[i] = True
            labels[i] = cluster_id
            
            while queue:
                node = queue.pop(0)
                neighbors = np.where((adj_matrix[node] == 1) & (~visited))[0]
                for neighbor in neighbors:
                    visited[neighbor] = True
                    labels[neighbor] = cluster_id
                    queue.append(neighbor)
            
            cluster_id += 1
            if cluster_id % 100 == 0:
                log(f"    已找到 {cluster_id} 个聚类")
    
    return labels, cluster_id

log("  开始寻找连通分量...")
clusters, num_clusters = find_connected_components(adjacency)
del adjacency
gc.collect()
log(f"  找到 {num_clusters:,} 个聚类")

# ===== 第五步：选择代表性基因组（最长序列） =====
log("\n步骤5: 选择代表性基因组（最长序列）...")

cluster_dict = defaultdict(list)
for genome, cluster_id in zip(genome_list, clusters):
    cluster_dict[cluster_id].append(genome)

representatives = []
cluster_info_list = []

for cluster_id, members in cluster_dict.items():
    # 选择最长序列
    # 如果基因组大小信息可用，按大小排序
    if members and any(g in genome_size for g in members):
        # 按大小降序排列，取最大的
        members_sorted = sorted(members, key=lambda x: genome_size.get(x, 0), reverse=True)
        rep = members_sorted[0]
    else:
        # 如果没有大小信息，按名称排序取第一个
        rep = sorted(members)[0]
    
    representatives.append(rep)
    cluster_info_list.append({
        'cluster_id': cluster_id,
        'representative': rep,
        'size': len(members),
        'members': members
    })

# 统计
cluster_sizes = [info['size'] for info in cluster_info_list]
log(f"  聚类数: {num_clusters:,}")
log(f"  代表性基因组: {len(representatives):,}")
log(f"  压缩率: {(1 - len(representatives)/n)*100:.2f}%")
log(f"  平均聚类大小: {np.mean(cluster_sizes):.2f}")
log(f"  最大聚类大小: {max(cluster_sizes):,}")

# ===== 第六步：保存结果 =====
log("\n步骤6: 保存结果...")

# 保存代表性列表
with open(REP_LIST, 'w') as f:
    f.write("# 代表性基因组列表（按最长序列选择）\n")
    f.write(f"# 总计: {len(representatives)} 个\n")
    f.write(f"# 聚类数: {num_clusters}\n")
    f.write("# 列: 代表基因组文件名\n")
    for rep in sorted(representatives):
        f.write(f"{rep}\n")

# 保存详细的聚类信息
with open(CLUSTER_INFO, 'w') as f:
    f.write("Cluster_ID\tRepresentative\tSize\tRepresentative_Size\tMembers\n")
    for info in sorted(cluster_info_list, key=lambda x: x['cluster_id']):
        rep_size = genome_size.get(info['representative'], 0)
        members_str = ','.join(info['members'])
        f.write(f"{info['cluster_id']}\t{info['representative']}\t{info['size']}\t{rep_size}\t{members_str}\n")

log(f"  代表性列表: {REP_LIST}")
log(f"  聚类信息: {CLUSTER_INFO}")

# 显示代表性基因组的平均大小
rep_sizes = [genome_size.get(rep, 0) for rep in representatives if rep in genome_size]
if rep_sizes:
    log(f"  代表性基因组平均大小: {np.mean(rep_sizes)/1e6:.2f} MB")
    log(f"  代表性基因组最大大小: {max(rep_sizes)/1e6:.2f} MB")
    log(f"  代表性基因组最小大小: {min(rep_sizes)/1e6:.2f} MB")

elapsed = time.time() - start_time
log(f"\n总耗时: {elapsed/60:.2f} 分钟")
log("完成！")
