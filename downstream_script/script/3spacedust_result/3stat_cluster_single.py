import os
import csv
from collections import defaultdict
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", type=str)
args = parser.parse_args()

base_dir = args.input_dir
output_file = "gene_count_simple.csv"

# 存储结果：{基因ID: {"count": 出现次数, "clusters": [簇信息列表], "folders": set([文件夹名])}}
gene_stats = defaultdict(lambda: {"count": 0, "clusters": [], "folders": set()})

print("正在统计...")

# 遍历所有文件夹
for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    if not os.path.isdir(folder_path):
        continue
    # 找结果文件
    result_file = os.path.join(folder_path, f"{folder}_sorted_clusters.txt")
    if not os.path.exists(result_file):
        result_file = os.path.join(folder_path, f"{folder}_result.tsv")
    if not os.path.exists(result_file):
        continue
    print(f"处理: {folder}")
    # 解析文件
    with open(result_file, 'r') as f:
        current_cluster = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                # 保存当前簇头部
                current_cluster = line.split("\t")[0]
            elif line.startswith('>'):
                # 提取target基因ID
                parts = line.split('\t')
                gene_id = parts[1]  # NZ_CP073017.1_29_39817_42363
                # 更新统计
                gene_stats[gene_id]["count"] += 1
                gene_stats[gene_id]["folders"].add(folder)  # 记录出现在哪个文件夹
                # 记录这个基因出现在哪个簇（用文件夹名:簇头部表示）
                cluster_info = f"{folder}:{current_cluster}"
                if cluster_info not in gene_stats[gene_id]["clusters"]:
                    gene_stats[gene_id]["clusters"].append(cluster_info)

# 写入CSV
with open(output_file, 'w', newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    # 增加一列"出现文件夹数"
    writer.writerow(["基因ID", "出现次数", "出现文件夹数", "所在簇列表"])
    # 按出现次数从高到低排序
    for gene_id, info in sorted(gene_stats.items(), key=lambda x: -x[1]["count"]):
        writer.writerow([
            gene_id,
            info["count"],
            len(info["folders"]),  # 文件夹数 = 物种基因组数
            " | ".join(info["clusters"])
        ])

print(f"\n完成！共统计 {len(gene_stats)} 个基因")
print(f"结果保存到: {output_file}")

# 显示前20个最常见的基因
print("\n出现次数最多的基因（前20）:")
print("-" * 80)
print(f"{'基因ID':<50} {'次数':<6} {'物种数':<6} {'簇数'}")
print("-" * 80)

for gene_id, info in sorted(gene_stats.items(), key=lambda x: -x[1]["count"])[:20]:
    print(f"{gene_id[:50]:<50} {info['count']:<6} {len(info['folders']):<6} {len(info['clusters'])}")
