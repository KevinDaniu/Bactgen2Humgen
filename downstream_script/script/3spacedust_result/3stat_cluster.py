import os
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
import sys
import pandas as pd
import numpy as np

def process_folder(folder, base_dir):
    """处理单个文件夹，返回该文件夹的统计结果"""
    folder_path = os.path.join(base_dir, folder)
    if not os.path.isdir(folder_path):
        return None
    
    # 找结果文件
    result_file = os.path.join(folder_path, f"{folder}_sorted_clusters.txt")
    if not os.path.exists(result_file):
        result_file = os.path.join(folder_path, f"{folder}_result.tsv")
    if not os.path.exists(result_file):
        return None
    
    # 该文件夹的统计结果
    folder_stats = defaultdict(lambda: {"count": 0, "clusters": [], "folders": set()})
    
    try:
        with open(result_file, 'r') as f:
            current_cluster = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    current_cluster = line.split("\t")[0]
                elif line.startswith('>'):
                    parts = line.split('\t')
                    gene_id = parts[1]
                    folder_stats[gene_id]["count"] += 1
                    folder_stats[gene_id]["folders"].add(folder)
                    cluster_info = f"{folder}:{current_cluster}"
                    if cluster_info not in folder_stats[gene_id]["clusters"]:
                        folder_stats[gene_id]["clusters"].append(cluster_info)
        
        return (folder, {k: dict(v) for k, v in folder_stats.items()})
    except Exception as e:
        print(f"处理 {folder} 时出错: {e}", file=sys.stderr)
        return None

def merge_stats_fast(results):
    """最快合并方式：使用pandas完全向量化"""
    # 收集所有记录到列表（使用列表推导式提速）
    all_records = []
    for result in results:
        if result is None:
            continue
        folder, folder_stats = result
        for gene_id, stats in folder_stats.items():
            # 展开每个cluster记录
            for cluster in stats['clusters']:
                all_records.append((gene_id, folder, cluster))
    
    if not all_records:
        return {}
    
    # 转换为DataFrame（一次性创建，更快）
    df = pd.DataFrame(all_records, columns=['gene_id', 'folder', 'cluster'])
    
    # 使用pandas的groupby聚合（C++实现，极快）
    result = {}
    for gene_id, group in df.groupby('gene_id', observed=True):
        result[gene_id] = {
            'count': len(group),
            'folders': set(group['folder'].unique()),
            'clusters': group['cluster'].unique().tolist()
        }
    
    return result

def write_results_fastest(output_file, gene_stats):
    """最快写入方式：完全向量化 + numpy"""
    # 使用列表推导式提取数据（比for循环快）
    items = list(gene_stats.items())
    
    # 提取到numpy数组（更快）
    gene_ids = np.array([item[0] for item in items], dtype=object)
    counts = np.array([item[1]['count'] for item in items], dtype=np.int64)
    folder_counts = np.array([len(item[1]['folders']) for item in items], dtype=np.int64)
    cluster_lists = np.array([' | '.join(item[1]['clusters']) for item in items], dtype=object)
    
    # 创建DataFrame
    df = pd.DataFrame({
        '基因ID': gene_ids,
        '出现次数': counts,
        '出现文件夹数': folder_counts,
        '所在簇列表': cluster_lists
    })
    
    # 使用numpy的argsort排序（比pandas的sort_values快）
    sorted_idx = np.argsort(-counts)  # 负号实现降序
    df = df.iloc[sorted_idx].reset_index(drop=True)
    
    # 写入CSV（使用最快的引擎）
    df.to_csv(output_file, index=False, encoding='utf-8')
    
    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir", type=str, required=True)
    parser.add_argument("--output", type=str, default="gene_count_simple.csv")
    parser.add_argument("--processes", type=int, default=None,
                       help="并行进程数，默认为CPU核心数")
    args = parser.parse_args()
    
    base_dir = args.input_dir
    output_file = args.output
    
    # 获取所有文件夹
    folders = [f for f in os.listdir(base_dir) 
               if os.path.isdir(os.path.join(base_dir, f))]
    
    if not folders:
        print("警告：没有找到任何文件夹")
        return
    
    num_processes = args.processes or os.cpu_count()
    print(f"正在统计... (共 {len(folders)} 个文件夹)")
    print(f"使用 {num_processes} 个进程并行处理")
    print("使用最快处理模式（向量化计算）")
    
    # 并行处理所有文件夹
    results = []
    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        future_to_folder = {
            executor.submit(process_folder, folder, base_dir): folder 
            for folder in folders
        }
        
        completed = 0
        total = len(folders)
        for future in as_completed(future_to_folder):
            folder = future_to_folder[future]
            try:
                result = future.result()
                if result is not None:
                    results.append(result)
                    print(f"处理: {folder}")
                completed += 1
                if completed % max(1, total // 10) == 0 or completed == total:
                    print(f"进度: {completed}/{total} ({completed*100//total}%)")
            except Exception as e:
                print(f"处理 {folder} 失败: {e}", file=sys.stderr)
                completed += 1
    
    # 合并结果（使用最快方式）
    print("合并统计结果...")
    gene_stats = merge_stats_fast(results)
    
    # 写入结果（使用最快方式）
    print(f"写入结果文件 ({len(gene_stats)} 个基因)...")
    df = write_results_fastest(output_file, gene_stats)
    
    print(f"\n✅ 完成！共统计 {len(gene_stats)} 个基因")
    print(f"结果保存到: {output_file}")
    
    # 显示前20个最常见的基因
    print("\n出现次数最多的基因（前20）:")
    print("-" * 80)
    print(f"{'基因ID':<50} {'次数':<6} {'物种数':<6} {'簇数'}")
    print("-" * 80)
    
    # 直接从DataFrame取前20（已经排序好了）
    for _, row in df.head(20).iterrows():
        gene_id = row['基因ID']
        count = row['出现次数']
        folder_count = row['出现文件夹数']
        # 计算簇数
        cluster_count = len(row['所在簇列表'].split(' | ')) if row['所在簇列表'] else 0
        print(f"{gene_id[:50]:<50} {count:<6} {folder_count:<6} {cluster_count}")

if __name__ == "__main__":
    main()
