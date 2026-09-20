#!/usr/bin/env python3
# 2.2-extract_gene_withGID.py - 使用CSV中的基因列表（按需加载 + 读表头修正 query/target）

import re
import csv
import sys
import gzip
from pathlib import Path
from collections import defaultdict
import argparse
import os

# ---------------- 工具函数 ----------------

def norm_genome(name):
    """统一基因组名：去掉常见后缀"""
    name = name.strip()
    for suf in ('.gff.gz', '.gff3.gz', '.gff', '.gff3'):
        if name.endswith(suf):
            name = name[:-len(suf)]
            break
    name = name.replace('_corrected', '')
    return name

def open_maybe_gzip(path):
    """根据后缀决定用 gzip 还是普通 open"""
    p = str(path)
    if p.endswith('.gz'):
        return gzip.open(p, 'rt', encoding='utf-8', errors='ignore')
    return open(p, 'r', encoding='utf-8', errors='ignore')

def find_gff_file(gff_dir, genome):
    """在 gff_dir 下找 genome 对应的 GFF，兼容多种后缀"""
    for suf in ('.gff.gz', '.gff3.gz', '.gff', '.gff3'):
        f = gff_dir / f"{genome}{suf}"
        if f.exists():
            return f
    g = genome.replace('_corrected', '')
    for suf in ('.gff.gz', '.gff3.gz', '.gff', '.gff3'):
        f = gff_dir / f"{g}{suf}"
        if f.exists():
            return f
    return None

# ---------------- GFF 解析 ----------------

def parse_gff_by_position(gff_file):
    """返回 (pos->gene, id->gene)"""
    genes_by_pos = {}
    genes_by_id = {}
    with open_maybe_gzip(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            try:
                start = int(fields[3]); end = int(fields[4])
            except ValueError:
                continue
            contig = fields[0]
            attr = fields[8]

            locus_match = re.search(r'locus_tag=([^;]+)', attr)
            gene_name = locus_match.group(1) if locus_match else "unknown"

            id_match = re.search(r'ID=([^;]+)', attr)
            gene_id = id_match.group(1) if id_match else None

            info = {
                'gene_name': gene_name,
                'contig': contig,
                'start': start,
                'end': end,
                'strand': fields[6]
            }
            genes_by_pos[(contig, start, end)] = info
            genes_by_pos[(contig, min(start, end), max(start, end))] = info
            if gene_id:
                genes_by_id[gene_id] = info
    return genes_by_pos, genes_by_id

# ---------------- 位置解析 ----------------

def parse_cluster_position(gene_id):
    """
    基因ID形如: NZ_OZ204251.1_1310_1427340_1426948
    格式: <contig>_<基因编号>_<坐标1>_<坐标2>
    末尾两段是坐标，再往前一段是基因编号，前面是 contig
    """
    gene_id = gene_id.lstrip('>')
    parts = gene_id.split('_')
    if len(parts) < 4:
        return None
    try:
        pos1 = int(parts[-2]); pos2 = int(parts[-1])
    except ValueError:
        return None
    contig = '_'.join(parts[:-3])
    return (contig, min(pos1, pos2), max(pos1, pos2))

def find_gene_in_gff(gff_data, genome, cluster_pos):
    if not cluster_pos or genome not in gff_data:
        return None
    return gff_data[genome][0].get(cluster_pos)

# ---------------- 主逻辑 ----------------

def extract_conserved_cluster(base_dir, gff_dir, refer_genome, conserved_id, output_file):
    print(f"开始处理保守簇: {conserved_id}")
    print(f"参考基因组为: {refer_genome}")
    print("-" * 60)

    reference_genome = norm_genome(refer_genome)

    # --- 读 conserve CSV ---
    csv_file = Path(base_dir) / 'gene_cooccurrence_conserve.csv'
    with open(csv_file, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['保守簇ID'] == conserved_id:
                target_relation = row['物种-簇对应关系']
                gene_list_str = row['基因列表']
                print(f"找到保守簇: {conserved_id}")
                print(f"基因数量: {row['基因数量']}, 保守物种数: {row['保守物种数']}")
                break
        else:
            print(f"错误: 未找到保守簇 {conserved_id}")
            return None

    ref_gene_ids = [gid.strip() for gid in gene_list_str.split('|') if gid.strip()]
    print(f"\n参考基因数量: {len(ref_gene_ids)}")

    # --- 收集需要的基因组 ---
    relations = [r.strip() for r in target_relation.split('|') if r.strip()]
    needed_genomes = {reference_genome}
    rel_parsed = []   # [(file_name, cluster_num)]
    for relation in relations:
        if ':#' not in relation:
            continue
        file_name, cluster_num = relation.rsplit(':#', 1)
        file_name = file_name.strip()
        cluster_num = cluster_num.strip()
        rel_parsed.append((file_name, cluster_num))
        if '_with_' in file_name:
            g1, g2 = file_name.split('_with_', 1)
            needed_genomes.add(norm_genome(g1))
            needed_genomes.add(norm_genome(g2))

    print(f"涉及 {len(rel_parsed)} 个比对文件, {len(needed_genomes)} 个基因组")

    # --- 按需加载 GFF ---
    print("\n按需读取 GFF ...")
    gff_dir_path = Path(gff_dir)
    gff_data = {}
    for genome in needed_genomes:
        gff_file = find_gff_file(gff_dir_path, genome)
        if gff_file is None:
            print(f"  警告: 找不到 {genome} 的 GFF")
            continue
        gff_data[genome] = parse_gff_by_position(gff_file)
        print(f"  加载 {genome}  ({gff_file.name})  解析到 {len(gff_data[genome][0])} 个位置")
    print(f"已加载 {len(gff_data)} 个基因组")

    if reference_genome not in gff_data:
        sys.exit(f"错误: 参考基因组 {reference_genome} 的 GFF 未加载")

    # --- 建立 ref_id -> 编号 映射 ---
    ref_id_to_num = {gid: i for i, gid in enumerate(ref_gene_ids, 1)}

    # --- 获取参考基因位置信息 ---
    ref_genes = []
    for gid in ref_gene_ids:
        gene_info = None
        pos = parse_cluster_position(gid)
        if pos:
            gene_info = find_gene_in_gff(gff_data, reference_genome, pos)
        if gene_info is None:
            for gid2, info in gff_data[reference_genome][1].items():
                if gid2 == gid or gid in gid2:
                    gene_info = info
                    break
        if gene_info is None:
            print(f"  警告: 参考基因组中未找到 {gid}  (pos={pos})")
        ref_genes.append({'ref_id': gid, 'gene_info': gene_info})

    n_found = sum(1 for r in ref_genes if r['gene_info'])
    print(f"参考基因匹配到位置: {n_found}/{len(ref_genes)}")

    # --- 存储结果 ---
    genome_genes = defaultdict(dict)
    for ref in ref_genes:
        if ref['gene_info'] is None:
            continue
        num = ref_id_to_num[ref['ref_id']]
        hid = f"Gene_{num:02d}"
        genome_genes[reference_genome][hid] = dict(ref['gene_info'], homologous_id=hid)

    # --- 处理每个比对关系 ---
    for file_name, cluster_num in rel_parsed:
        cluster_file = Path(base_dir) / file_name / f"{file_name}_result.tsv"
        if not cluster_file.exists():
            print(f"  警告: 文件不存在 {cluster_file}")
            continue

        # ★★★ 读表头拿 query/target ★★★
        with open(cluster_file, 'r', encoding='utf-8', errors='ignore') as f:
            header_line = None
            for line in f:
                if line.startswith(f'#{cluster_num}'):
                    header_line = line
                    break

        if not header_line:
            print(f"  警告: 未找到簇 #{cluster_num} 的表头")
            continue

        header_parts = header_line.strip().split('\t')
        if len(header_parts) < 3:
            print(f"  警告: 表头字段不足: {header_line.strip()}")
            continue

        query_genome = norm_genome(header_parts[1])
        target_genome = norm_genome(header_parts[2])

        if target_genome == reference_genome:
            other_genome = query_genome
            ref_is_target = True
        elif query_genome == reference_genome:
            other_genome = target_genome
            ref_is_target = False
        else:
            print(f"  警告: 参考基因组 {reference_genome} 不在表头: {query_genome} / {target_genome}")
            continue

        if other_genome not in gff_data:
            print(f"  警告: {other_genome} 的 GFF 未加载")
            continue

        # --- 读基因对 ---
        with open(cluster_file, 'r', encoding='utf-8', errors='ignore') as f:
            in_target = False
            n_hit = 0
            for line in f:
                line = line.rstrip('\n')
                if not line:
                    continue
                if line.startswith(f'#{cluster_num}'):
                    in_target = True
                    continue
                if in_target and line.startswith('#'):
                    nxt = line[1:].split()[0]
                    if nxt.isdigit():
                        break
                if in_target:
                    fields = line.split('\t')
                    if len(fields) < 2:
                        continue
                    query_id = fields[0].lstrip('>')
                    target_id = fields[1].lstrip('>')

                    if ref_is_target:
                        ref_id = target_id
                        other_pos = parse_cluster_position(query_id)
                    else:
                        ref_id = query_id
                        other_pos = parse_cluster_position(target_id)

                    if ref_id not in ref_id_to_num:
                        continue
                    hid = f"Gene_{ref_id_to_num[ref_id]:02d}"

                    gene_info = find_gene_in_gff(gff_data, other_genome, other_pos)
                    if gene_info and hid not in genome_genes[other_genome]:
                        genome_genes[other_genome][hid] = dict(gene_info, homologous_id=hid)
                        n_hit += 1
            print(f"  {file_name}:#{cluster_num}  命中 {n_hit} 个")

    # --- 输出 ---
    all_genes = []
    for genome, genes in genome_genes.items():
        for hid in sorted(genes.keys()):
            g = genes[hid]
            all_genes.append({
                'molecule': genome,
                'gene': g['gene_name'],
                'start': g['start'],
                'end': g['end'],
                'strand': g['strand'],
                'cluster': conserved_id,
                'homologous_id': hid
            })

    print("\n" + "=" * 60)
    print(f"总计提取 {len(all_genes)} 个基因")
    print(f"涉及 {len(genome_genes)} 个基因组")

    if not all_genes:
        print("错误: 没有提取到任何基因")
        return None

    import pandas as pd
    df = pd.DataFrame(all_genes).sort_values(['molecule', 'homologous_id'])
    out = Path(output_file)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    print(f"\n结果已保存到: {out}")
    print("\n各基因组基因数量统计:")
    for mol in sorted(df['molecule'].unique()):
        print(f"  {mol}: {len(df[df['molecule'] == mol])}")
    print("\n数据预览:")
    print(df[['molecule', 'homologous_id', 'gene', 'strand']].head(10))
    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_dir", type=str)
    parser.add_argument("--id", type=str)
    args = parser.parse_args()

    BASE_DIR = Path(f"/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/bactgene2humgene/{args.base_dir}/gene_cluster/cluster_genes")
    GFF_DIR = Path("/home/liuzhh/bacteria_cluster/260811_allgenome_ani/result/all_gff_files")

    conserved_id = args.id
    SAVE_DIR = f"/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result/{args.base_dir}"
    os.makedirs(SAVE_DIR, exist_ok=True)
    output_file = f"{SAVE_DIR}/{args.base_dir}_{conserved_id}_for_gggenes.csv"

    result = extract_conserved_cluster(BASE_DIR, GFF_DIR, args.base_dir, conserved_id, output_file)
    if result is not None:
        print(f"\n✓ 提取完成! 共 {len(result)} 个基因")
