import os
import pandas as pd
import sys
import glob
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input_dir", type=str)
parser.add_argument("--blast_file", type=str,
        default="/gene2human/blast_result/homology_summary.csv"
        )
parser.add_argument("--blastp_file", type=str,
        default="/gene2human/blastp_result/homology_summary.csv"
        )
parser.add_argument("--eggnog_file", type=str,
        default="/gene2human/eggnog_result/summary/eggnog_homology_summary.csv"
        )
parser.add_argument("--eggnog_anno_dir", type=str,
        default="/none/none"
        )
parser.add_argument("--output_dir", type=str,
        default="/gene2human/final_result"
        )
args = parser.parse_args()

# 拼接路径
input_dir = args.input_dir
args.blast_file = input_dir + args.blast_file
args.blastp_file = input_dir + args.blastp_file
args.eggnog_file = input_dir + args.eggnog_file
args.output_dir = input_dir + args.output_dir

# 设置路径
blast_file = args.blast_file
blastp_file = args.blastp_file
eggnog_file = args.eggnog_file
eggnog_anno_dir = args.eggnog_anno_dir
output_dir = args.output_dir

# 创建输出目录
Path(output_dir).mkdir(parents=True, exist_ok=True)
print("开始融合BLAST和eggNOG结果...")

# 读取BLAST结果
print("读取BLAST结果...")
blast_df = pd.read_csv(blast_file, sep=',')
blastp_df = pd.read_csv(blastp_file, sep=',')
# 过滤掉"No hit"的行
blast_df = blast_df[~blast_df['Cluster Gene'].str.contains('no_hit', na=False)]
blastp_df = blastp_df[~blastp_df['Cluster Gene'].str.contains('no_hit', na=False)]
print(f"BLAST结果: {len(blast_df)} 条记录")
print(f"BLASTp结果: {len(blastp_df)} 条记录")
blasts_df = pd.merge(
    blast_df,
    blastp_df,
    on=['Conserved ID', 'Cluster Gene'],
    how='outer'
)

# 检查eggNOG文件是否存在
if not os.path.exists(eggnog_file):
    print("警告: eggNOG结果文件不存在，只保存BLAST结果")
    blasts_df['Source'] = 'BLAST only'
    blasts_df['eggnog_Gene'] = 'NA'
    blasts_df['eggnog_Name'] = 'NA'
    blasts_df['eggnog_E-value'] = 'NA'
    blasts_df['eggnog_Score'] = 'NA'
    blasts_df['eggnog_Percent_Identity'] = 'NA'
    blasts_df['GOs'] = 'NA'
    blasts_df['Description'] = 'NA'
    blasts_df['COG_category'] = 'NA'
    blasts_df['KEGG_ko'] = 'NA'

    blasts_df.to_csv(os.path.join(output_dir, 'final_homology_results.csv'), index=False)
    blasts_df.to_csv(os.path.join(output_dir, 'final_homology_results.txt'), sep='\t', index=False)
    sys.exit(0)

# 读取eggNOG结果
print("读取eggNOG结果...")
eggnog_df = pd.read_csv(eggnog_file)
print(f"eggNOG结果: {len(eggnog_df)} 条记录")

# 重命名eggNOG的列，避免与BLAST列名冲突
eggnog_df = eggnog_df.rename(columns={
    'Human Gene': 'eggnog_Gene',
    'Preferred_name': 'eggnog_Name',
    'E-value': 'eggnog_E-value',
    'Score': 'eggnog_Score',
    'Percent Identity': 'eggnog_Percent_Identity'
})

# 读取所有eggNOG注释文件，提取GO信息
print("读取eggNOG注释文件，提取GO信息...")

# 获取所有annotations文件
anno_files = glob.glob(os.path.join(eggnog_anno_dir, "*.emapper.annotations"))
print(f"找到 {len(anno_files)} 个注释文件")

# 创建一个字典来存储每个基因的GO信息
name_info = {}
go_info = {}
description_info = {}
cog_info = {}
kegg_info = {}

for anno_file in anno_files:
    # 从文件名提取Conserved ID
    conserved_id = os.path.basename(anno_file).replace('.emapper.annotations', '')
    
    # 读取注释文件，跳过注释行（以#开头的行）
    try:
        # 先读取所有行，找到表头
        with open(anno_file, 'r') as f:
            lines = f.readlines()
        
        # 找到表头行（不以#开头的第一行）
        header_line = None
        data_lines = []
        for line in lines:
            if not line.startswith('##') and header_line is None:
                header_line = line.strip().split('\t')
            elif not line.startswith('##') and header_line is not None:
                data_lines.append(line.strip().split('\t'))
        
        if header_line and data_lines:
            # 创建DataFrame
            anno_df = pd.DataFrame(data_lines, columns=header_line)
            
            # 提取query名称（去掉Conserved ID前缀）
            anno_df['Cluster Gene'] = anno_df['#query'].apply(
                lambda x: x.replace(f'{conserved_id}|', '') if f'{conserved_id}|' in x else x
            )
            
            # 存储GO信息
            for _, row in anno_df.iterrows():
                gene_name = row['Cluster Gene']
                # 使用Conserved ID和Cluster Gene作为复合键
                key = f"{conserved_id}|{gene_name}"
                
                if 'Preferred_name' in anno_df.columns:
                    name_info[key] = row.get('Preferred_name', 'NA')
                if 'GOs' in anno_df.columns:
                    go_info[key] = row.get('GOs', 'NA')
                if 'Description' in anno_df.columns:
                    description_info[key] = row.get('Description', 'NA')
                if 'COG_category' in anno_df.columns:
                    cog_info[key] = row.get('COG_category', 'NA')
                if 'KEGG_ko' in anno_df.columns:
                    kegg_info[key] = row.get('KEGG_ko', 'NA')
                    
    except Exception as e:
        print(f"  警告: 读取文件 {anno_file} 时出错: {e}")

print(f"成功提取 {len(go_info)} 个基因的注释信息")

# 合并两个数据框
print("合并结果...")
merged_df = pd.merge(
    blasts_df,
    eggnog_df,
    on=['Conserved ID', 'Cluster Gene'],
    how='outer'
)

# 设置Source列
merged_df['Source'] = 'Both'
merged_df.loc[merged_df['Human Gene (by Comprehensive)'].isna() & merged_df['eggnog_Gene'].notna(), 'Source'] = 'eggNOG only'
merged_df.loc[merged_df['Human Gene (by Comprehensive)'].notna() & merged_df['eggnog_Gene'].isna(), 'Source'] = 'BLAST only'

# 创建复合键用于匹配注释信息
merged_df['gene_key'] = merged_df['Conserved ID'] + '|' + merged_df['Cluster Gene']

# 添加注释信息
merged_df['Name'] = merged_df['gene_key'].map(name_info).fillna('NA')
merged_df['GOs'] = merged_df['gene_key'].map(go_info).fillna('NA')
merged_df['Description'] = merged_df['gene_key'].map(description_info).fillna('NA')
merged_df['COG_category'] = merged_df['gene_key'].map(cog_info).fillna('NA')
merged_df['KEGG_ko'] = merged_df['gene_key'].map(kegg_info).fillna('NA')

# 删除辅助列
merged_df = merged_df.drop(columns=['gene_key'])

# 将所有空值替换为NA
merged_df = merged_df.fillna('NA')

# 检查所有需要的列是否存在
print("可用的列:", merged_df.columns.tolist())

# 重新排列列的顺序（只选择存在的列）
base_columns = ['Conserved ID', 'Cluster Gene', 'Source']
blast_columns1 = ['Human Gene (by Length)', 'Length (bp) (by Length)', 'Identity% (by Length)', 'E-value (by Length)']
blast_columns2 = ['Human Gene (by E-value)', 'Length (bp) (by E-value)', 'Identity% (by E-value)', 'E-value (by E-value)']
blast_columns3 = ['Human Gene (by Identity)', 'Length (bp) (by Identity)', 'Identity% (by Identity)', 'E-value (by Identity)']
blast_columns4 = ['Human Gene (by Comprehensive)', 'Length (bp) (by Comprehensive)', 'Identity% (by Comprehensive)', 'E-value (by Comprehensive)']
blastp_columns1 = ['Human Prot (by Length)', 'Prot Length (bp) (by Length)', 'Prot Identity% (by Length)', 'Prot E-value (by Length)']
blastp_columns2 = ['Human Prot (by E-value)', 'Prot Length (bp) (by E-value)', 'Prot Identity% (by E-value)', 'E-value (by E-value)']
blastp_columns3 = ['Human Prot (by Identity)', 'Prot Length (bp) (by Identity)', 'Prot Identity% (by Identity)', 'E-value (by Identity)']
blastp_columns4 = ['Human Prot (by Comprehensive)', 'Prot Length (bp) (by Comprehensive)', 'Prot Identity% (by Comprehensive)', 'Prot E-value (by Comprehensive)']
eggnog_columns = ['eggnog_Gene', 'eggnog_E-value', 'eggnog_Score', 'eggnog_Percent_Identity']
annotation_columns = ['Name', 'GOs', 'Description', 'COG_category', 'KEGG_ko']

# 只选择实际存在的列
all_columns = base_columns + blast_columns1 + blast_columns2 + blast_columns3 + blast_columns4
all_columns = all_columns + eggnog_columns + annotation_columns
all_columns = all_columns + blastp_columns1 + blastp_columns2 + blastp_columns3 + blastp_columns4
existing_columns = [col for col in all_columns if col in merged_df.columns]
merged_df = merged_df[existing_columns]

# 按Conserved ID排序
merged_df = merged_df.sort_values('Conserved ID')

# 保存结果
output_csv = os.path.join(output_dir, 'final_homology_results.csv')
output_txt = os.path.join(output_dir, 'final_homology_results.txt')

merged_df.to_csv(output_csv, index=False)
merged_df.to_csv(output_txt, sep='\t', index=False)

print(f"\n处理完成!")
print(f"结果已保存到: {output_dir}")
print(f"共处理 {len(merged_df)} 条记录")
print(f"统计:")
print(f"  - Both: {len(merged_df[merged_df['Source'] == 'Both'])}")
print(f"  - BLAST only: {len(merged_df[merged_df['Source'] == 'BLAST only'])}")
print(f"  - eggNOG only: {len(merged_df[merged_df['Source'] == 'eggNOG only'])}")
