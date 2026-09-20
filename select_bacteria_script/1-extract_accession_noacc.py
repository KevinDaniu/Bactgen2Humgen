#!/usr/bin/env python3
"""
从 assembly_summary.tsv 匹配 bacteria_refseq.tsv 中的 GCF 编号
生成两个文件：
1. ncbi_refseq_index.txt - 不含版本号（原有格式，保留）
2. ncbi_refseq_index_with_version.txt - 包含完整版本号（新增）
3. matched_gcf_list.txt - 匹配到的GCF列表（新增）
"""

import os
import pandas as pd
from pathlib import Path

def extract_refseq_genomes(tsv_file, assembly_file, output_file):
    """
    从 assembly_summary.tsv 匹配 bacteria_refseq.tsv 中的 GCF 编号
    
    Parameters:
    -----------
    tsv_file : str
        bacteria_refseq.tsv 文件路径
    assembly_file : str
        assembly_summary.tsv 文件路径
    output_file : str
        输出文件路径（genome_refseq_index.txt）
    
    Returns:
    --------
    dict: 统计信息
    """
    
    print("=" * 70)
    print("从 assembly_summary.tsv 匹配 bacteria_refseq 的基因组")
    print("=" * 70)
    print(f"TSV文件: {tsv_file}")
    print(f"Assembly文件: {assembly_file}")
    print(f"输出文件: {output_file}")
    print()
    
    # 检查文件是否存在
    for file_path in [tsv_file, assembly_file]:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"文件不存在: {file_path}")
    
    # 1. 读取 bacteria_refseq.tsv
    print("步骤1: 读取 bacteria_refseq.tsv...")
    try:
        tsv_df = pd.read_csv(tsv_file, sep='\t', low_memory=False)
        print(f"  TSV文件总行数: {len(tsv_df)}")
        print(f"  TSV列名: {list(tsv_df.columns)}")
    except Exception as e:
        print(f"  读取TSV文件失败: {e}")
        raise
    
    # 提取GCF编号（第2列是Assembly Accession）
    if 'Assembly Accession' in tsv_df.columns:
        gcf_col = 'Assembly Accession'
    else:
        # 如果列名不对，尝试按位置获取（第2列）
        print("  警告: 未找到'Assembly Accession'列，使用第2列")
        gcf_col = tsv_df.columns[1]
    
    # 创建一个字典：基础GCF -> 完整版本号
    gcf_version_map = {}
    for acc in tsv_df[gcf_col].dropna():
        acc_str = str(acc).strip()
        base_gcf = acc_str.split('.')[0]
        # 如果有多个版本，保留最新的（通常版本号数字大的更新）
        if base_gcf not in gcf_version_map or acc_str > gcf_version_map[base_gcf]:
            gcf_version_map[base_gcf] = acc_str
    
    # 创建基础GCF集合用于匹配
    tsv_gcf_base = set(gcf_version_map.keys())
    
    print(f"  bacteria_refseq.tsv中唯一的GCF总数（基础编号）: {len(tsv_gcf_base)}")
    print(f"  bacteria_refseq.tsv中GCF版本号示例: {list(gcf_version_map.items())[:5]}")
    
    # 2. 读取 assembly_summary.tsv
    print("\n步骤2: 读取 assembly_summary.tsv...")
    try:
        assembly_df = pd.read_csv(assembly_file, sep='\t', comment='#', low_memory=False)
        print(f"  Assembly文件总行数: {len(assembly_df)}")
        print(f"  Assembly列名: {list(assembly_df.columns)}")
    except Exception as e:
        print(f"  读取Assembly文件失败: {e}")
        raise
    
    # 检查必要的列
    required_cols = ['assembly_accession', 'ftp_path']
    for col in required_cols:
        if col not in assembly_df.columns:
            raise ValueError(f"assembly_summary.tsv中缺少列: {col}")
    
    # 3. 匹配并生成输出文件
    print("\n步骤3: 匹配并生成输出文件...")
    
    # 创建输出目录
    os.makedirs(os.path.dirname(output_file) or '.', exist_ok=True)
    
    # 生成文件名
    output_file_base = output_file  # 不含版本号
    output_file_version = output_file.replace('.txt', '_with_version.txt')  # 含版本号
    output_file_matched_list = output_file.replace('ncbi_refseq_index.txt', 'matched_gcf_list.txt')  # 匹配到的GCF列表
    output_file_unmatched = output_file.replace('ncbi_refseq_index.txt', 'unmatched_gcf.txt')  # 未匹配的GCF列表
    
    matched_gcf = set()
    matched_gcf_with_version = set()
    output_lines_base = []  # 不含版本号
    output_lines_version = []  # 含版本号
    assembly_levels = {}
    
    # 遍历assembly_summary
    for idx, row in assembly_df.iterrows():
        acc = row['assembly_accession']
        # 提取基础编号（不含版本号）
        acc_base = str(acc).split('.')[0]
        
        # 检查是否在TSV的GCF列表中
        if acc_base in tsv_gcf_base:
            matched_gcf.add(acc_base)
            
            # 获取ftp_path
            ftp_path = row['ftp_path']
            if pd.isna(ftp_path) or ftp_path == '':
                print(f"  警告: {acc_base} 的ftp_path为空")
                continue
            
            # 提取基础文件名
            base_name = ftp_path.split('/')[-1]
            
            # 构建fna和gff文件路径
            fna_file = f"{base_name}_genomic.fna.gz"
            gff_file = f"{base_name}_genomic.gff.gz"
            
            # 获取完整的版本号（从TSV中获取）
            versioned_gcf = gcf_version_map.get(acc_base, acc_base)
            matched_gcf_with_version.add(versioned_gcf)
            
            # 输出格式1：不含版本号（原有格式）
            output_lines_base.append(f"{acc_base}\tfna\t{fna_file}")
            output_lines_base.append(f"{acc_base}\tgff\t{gff_file}")
            
            # 输出格式2：含版本号（新增格式）
            output_lines_version.append(f"{versioned_gcf}\tfna\t{fna_file}")
            output_lines_version.append(f"{versioned_gcf}\tgff\t{gff_file}")
            
            # 统计Assembly Level分布
            if 'assembly_level' in assembly_df.columns:
                level = row['assembly_level']
                if pd.isna(level):
                    level = 'Unknown'
                assembly_levels[level] = assembly_levels.get(level, 0) + 1
    
    # 写入文件1：不含版本号（原有格式）
    with open(output_file_base, 'w') as f:
        for line in output_lines_base:
            f.write(line + '\n')
    print(f"  已生成不含版本号文件: {output_file_base}")
    
    # 写入文件2：含版本号（新增格式）
    with open(output_file_version, 'w') as f:
        for line in output_lines_version:
            f.write(line + '\n')
    print(f"  已生成含版本号文件: {output_file_version}")
    
    # 写入文件3：匹配到的GCF列表（新增）
    with open(output_file_matched_list, 'w') as f:
        f.write("# 匹配到的GCF列表（基础编号）\n")
        f.write(f"# 总计: {len(matched_gcf)} 个\n")
        for gcf in sorted(matched_gcf):
            f.write(f"{gcf}\n")
    print(f"  已生成匹配GCF列表: {output_file_matched_list}")
    
    # 写入文件4：未匹配的GCF列表
    if len(tsv_gcf_base) - len(matched_gcf) > 0:
        unmatched_gcf = tsv_gcf_base - matched_gcf
        with open(output_file_unmatched, 'w') as f:
            f.write("# 未匹配的GCF列表（基础编号）\n")
            f.write(f"# 总计: {len(unmatched_gcf)} 个\n")
            for gcf in sorted(unmatched_gcf):
                f.write(f"{gcf}\n")
        print(f"  已生成未匹配GCF列表: {output_file_unmatched}")
    
    # 4. 统计信息
    total_lines_base = len(output_lines_base)
    total_lines_version = len(output_lines_version)
    unique_gcf = len(matched_gcf)
    unmatched = len(tsv_gcf_base) - unique_gcf
    
    print("\n" + "=" * 70)
    print("           统计报告")
    print("=" * 70)
    print(f"bacteria_refseq.tsv中的GCF总数: {len(tsv_gcf_base)}")
    print(f"匹配到的GCF数: {unique_gcf}")
    print(f"未匹配的GCF数: {unmatched}")
    print(f"匹配率: {unique_gcf/len(tsv_gcf_base)*100:.2f}%")
    print()
    print(f"生成的文件:")
    print(f"  1. {os.path.basename(output_file_base)}:")
    print(f"     - 总行数: {total_lines_base}")
    print(f"     - 唯一GCF数: {unique_gcf}")
    print(f"     - 文件大小: {os.path.getsize(output_file_base) / 1024:.2f} KB")
    print(f"     - 格式: 不含版本号（原有格式）")
    print(f"  2. {os.path.basename(output_file_version)}:")
    print(f"     - 总行数: {total_lines_version}")
    print(f"     - 唯一GCF数: {unique_gcf}")
    print(f"     - 文件大小: {os.path.getsize(output_file_version) / 1024:.2f} KB")
    print(f"     - 格式: 含完整版本号（新增格式）")
    print(f"  3. {os.path.basename(output_file_matched_list)}:")
    print(f"     - 唯一GCF数: {unique_gcf}")
    print(f"     - 文件大小: {os.path.getsize(output_file_matched_list) / 1024:.2f} KB")
    print(f"     - 格式: 匹配到的GCF列表（新增）")
    if unmatched > 0:
        print(f"  4. {os.path.basename(output_file_unmatched)}:")
        print(f"     - 未匹配GCF数: {unmatched}")
        print(f"     - 文件大小: {os.path.getsize(output_file_unmatched) / 1024:.2f} KB")
        print(f"     - 格式: 未匹配的GCF列表")
    print(f"  - 保存路径: {os.path.dirname(output_file_base)}")
    
    # 显示Assembly Level分布
    if assembly_levels:
        print("\n匹配基因组的Assembly Level分布:")
        for level, count in sorted(assembly_levels.items(), key=lambda x: x[1], reverse=True):
            print(f"  {level}: {count} ({count/unique_gcf*100:.1f}%)")
    
    # 显示版本号映射示例
    print("\n版本号映射示例（前10个）:")
    for i, (base, version) in enumerate(list(gcf_version_map.items())[:10]):
        print(f"  {base} -> {version}")
    if len(gcf_version_map) > 10:
        print(f"  ... (还有 {len(gcf_version_map) - 10} 个映射)")
    
    print("\n" + "=" * 70)
    print("处理完成！")
    print("=" * 70)
    
    return {
        'total_tsv': len(tsv_gcf_base),
        'matched': unique_gcf,
        'unmatched': unmatched,
        'total_lines_base': total_lines_base,
        'total_lines_version': total_lines_version,
        'assembly_levels': assembly_levels
    }


def main():
    # 文件路径
    BASE_DIR = "/home/liuzhh/bacteria_cluster/260811_allgenome_ani/select_bacteria_result"
    TSV_FILE = "/home/liuzhh/bacteria_cluster/bacteria_genomes/20260201/bacteria_refseq.tsv"
    ASSEMBLY_FILE = "/home/liuzhh/bacteria_cluster/bacteria_genomes/20260201/assembly_summary.tsv"
    OUTPUT_FILE = f"{BASE_DIR}/ncbi_refseq_index.txt"
    
    # 执行提取
    try:
        stats = extract_refseq_genomes(
            tsv_file=TSV_FILE,
            assembly_file=ASSEMBLY_FILE,
            output_file=OUTPUT_FILE
        )
        
        # 显示生成文件示例
        print("\n" + "=" * 70)
        print("生成文件示例:")
        print("=" * 70)
        
        print("\n1. 不含版本号文件 (原有格式) - 前10行:")
        with open(OUTPUT_FILE, 'r') as f:
            for i, line in enumerate(f):
                if i >= 10:
                    break
                print(f"  {line.strip()}")
        
        version_file = OUTPUT_FILE.replace('.txt', '_with_version.txt')
        print(f"\n2. 含版本号文件 (新增格式) - 前10行:")
        with open(version_file, 'r') as f:
            for i, line in enumerate(f):
                if i >= 10:
                    break
                print(f"  {line.strip()}")
        
        matched_list_file = OUTPUT_FILE.replace('ncbi_refseq_index.txt', 'matched_gcf_list.txt')
        print(f"\n3. 匹配到的GCF列表 (新增) - 前20行:")
        with open(matched_list_file, 'r') as f:
            for i, line in enumerate(f):
                if i >= 20:
                    break
                print(f"  {line.strip()}")
        
        print(f"\n文件位置: {os.path.dirname(OUTPUT_FILE)}")
        print("  - ncbi_refseq_index.txt (原有格式，不含版本号)")
        print("  - ncbi_refseq_index_with_version.txt (新增格式，含完整版本号)")
        print("  - matched_gcf_list.txt (新增，匹配到的22000个GCF列表)")
        
    except Exception as e:
        print(f"\n错误: {e}")
        raise


if __name__ == "__main__":
    main()
