import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--name", type=str)
parser.add_argument("--input_dir", type=str)
args = parser.parse_args()

def sort_cluster_file(input_file, output_file):
    """
    对簇文件进行排序，保持每个簇的完整格式不变
    """
    clusters = []  # 存储每个簇的所有行
    current_cluster = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].rstrip('\n')
        
        # 找到簇的开始（以#开头）
        if line.startswith('#'):
            # 如果之前有簇，保存它
            if current_cluster:
                clusters.append(current_cluster)
            
            # 开始新簇
            current_cluster = [line]  # 保存簇头部
            
            # 读取接下来的基因行（以>开头）
            i += 1
            while i < len(lines) and lines[i].strip() and lines[i].startswith('>'):
                current_cluster.append(lines[i].rstrip('\n'))
                i += 1
            
            # 跳过可能的空行
            while i < len(lines) and not lines[i].strip():
                i += 1
            continue
        else:
            i += 1
    
    # 添加最后一个簇
    if current_cluster:
        clusters.append(current_cluster)
    
    print(f"找到 {len(clusters)} 个簇")
    
    # 定义排序键函数
    def get_sort_key(cluster):
        """从簇头部提取排序所需的值"""
        header = cluster[0]  # 第一行是#开头
        
        # 按制表符分割
        fields = header.split('\t')
        
        try:
            # 格式: #ID, query_gff, target_gff, cluster_p, multihit_p, gene_count
            cluster_p = float(fields[3])
            multihit_p = float(fields[4])
            gene_count = int(fields[5])
            
            # 返回排序键：cluster_p升序, multihit_p升序, gene_count降序
            return (cluster_p, multihit_p, -gene_count)
        except (IndexError, ValueError) as e:
            print(f"警告：解析簇头部时出错：{header[:100]}")
            print(f"字段数：{len(fields)}")
            print(f"错误信息：{e}")
            return (float('inf'), float('inf'), 0)
    
    # 排序前打印信息
    print("\n排序前：")
    for i, cluster in enumerate(clusters[:5]):
        header = cluster[0]
        fields = header.split('\t')
        print(f"簇{i+1}: {fields[0]} - p={fields[3]}, multi={fields[4]}, genes={fields[5]}")
    
    # 执行排序
    clusters.sort(key=get_sort_key)
    
    # 排序后打印信息
    print("\n排序后：")
    for i, cluster in enumerate(clusters[:5]):
        header = cluster[0]
        fields = header.split('\t')
        print(f"簇{i+1}: {fields[0]} - p={fields[3]}, multi={fields[4]}, genes={fields[5]}")
    
    # 写入文件，保持原始格式
    with open(output_file, 'w') as f:
        for i, cluster in enumerate(clusters):
            # 写入簇的所有行
            for line in cluster:
                f.write(line + '\n')
            # 在簇之间加空行（除了最后一个）
            if i < len(clusters) - 1:
                f.write('\n')
    
    print(f"\n排序完成！结果保存在 {output_file}")


INPUT = args.name
INPUT_DIR = args.input_dir
input_file = f"{INPUT_DIR}/{INPUT}/{INPUT}_result.tsv"
output_file = f"{INPUT_DIR}/{INPUT}/{INPUT}_sorted_clusters.txt"
if not os.path.exists(input_file):
    print(f"错误：输入文件不存在：{input_file}")
    sys.exit(1)    
sort_cluster_file(input_file, output_file)
