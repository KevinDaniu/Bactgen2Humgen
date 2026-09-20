#!/bin/bash
set -euo pipefail

# ==================== 路径配置 ====================
SUMMARY=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result/high_conserved_summary.csv
ASM=/home/liuzhh/bacteria_cluster/bacteria_genomes/20260201/assembly_summary.tsv
OUTPUTDIR=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/downstream_result
OUTFILE="$OUTPUTDIR/high_conserved_summary_annotated.csv"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# ==================== 1. summary -> GCF + strain ====================
awk -F',' 'NR>1 {
    strain=$1
    if (match(strain, /GCF_[0-9]+/)) {
        gcf=substr(strain, RSTART, RLENGTH)
        print gcf "\t" strain
    }
}' "$SUMMARY" > "$TMPDIR/strain_gcf.tsv"

echo "[1] 含 GCF 的菌株: $(wc -l < "$TMPDIR/strain_gcf.tsv")"

# ==================== 2. GCF -> species_taxid + organism_name ====================
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR {keep[$1]=1; next}
    /^#/ {next}
    FNR==1 {next}
    {
        gcf=$1
        sub(/\.[0-9]+$/, "", gcf)
        if (gcf in keep) print gcf, $7, $8
    }
' "$TMPDIR/strain_gcf.tsv" "$ASM" > "$TMPDIR/gcf_taxid.tsv"

echo "[2] 匹配到 taxid 的 GCF: $(wc -l < "$TMPDIR/gcf_taxid.tsv")"
head -3 "$TMPDIR/gcf_taxid.tsv"

# ==================== 3. taxonkit 查询（带 rank 名） ====================
cut -f2 "$TMPDIR/gcf_taxid.tsv" | sort -u > "$TMPDIR/taxid_list.txt"
echo "[3] 唯一 taxid: $(wc -l < "$TMPDIR/taxid_list.txt")"

taxonkit lineage --show-lineage-ranks "$TMPDIR/taxid_list.txt" \
    > "$TMPDIR/taxid_lineage_rank.txt" 2> "$TMPDIR/taxonkit.err"

echo "[3] taxonkit 输出行数: $(wc -l < "$TMPDIR/taxid_lineage_rank.txt")"
head -2 "$TMPDIR/taxid_lineage_rank.txt"

# ==================== 4. taxid -> family, genus（按 rank 名提取） ====================
# taxonkit lineage --show-lineage-ranks 列结构：
#   $1 = taxid
#   $2 = lineage (名字, ; 分隔)
#   $3 = lineage_ranks (rank 名, ; 分隔)
awk -F'\t' 'BEGIN{OFS="\t"}
    {
        n=split($2, names, ";")
        m=split($3, ranks, ";")
        fam=""; gen=""
        for(i=1; i<=n && i<=m; i++){
            if(ranks[i]=="family") fam=names[i]
            if(ranks[i]=="genus")  gen=names[i]
        }
        print $1, fam, gen
    }
' "$TMPDIR/taxid_lineage_rank.txt" > "$TMPDIR/taxid_fg.tsv"

echo "[4] taxid->family/genus 行数: $(wc -l < "$TMPDIR/taxid_fg.tsv")"
head -3 "$TMPDIR/taxid_fg.tsv"

# ==================== 5. GCF -> organism_name + family + genus ====================
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR {fg[$1]=$2"\t"$3; next}
    {print $1, $2, $3, fg[$2]}
' "$TMPDIR/taxid_fg.tsv" "$TMPDIR/gcf_taxid.tsv" > "$TMPDIR/gcf_annot.tsv"
# 列: GCF, species_taxid, organism_name, family, genus

echo "[5] gcf_annot 行数: $(wc -l < "$TMPDIR/gcf_annot.tsv")"
head -3 "$TMPDIR/gcf_annot.tsv"

# ==================== 6. 合并回 summary ====================
awk -F'\t' 'BEGIN{OFS="\t"}
    # 读注释表：key=GCF, value="organism_name|family|genus"
    NR==FNR {
        annot[$1] = $3 "|" $4 "|" $5
        next
    }
    # 处理 summary（逗号分隔）
    FNR==1 { next }
    {
        strain = $1
        if (match(strain, /GCF_[0-9]+/)) {
            gcf = substr(strain, RSTART, RLENGTH)
            if (gcf in annot) {
                split(annot[gcf], v, "|")
                print $0 "," v[1] "," v[2] "," v[3]
            } else {
                print $0 ",,,"
            }
        } else {
            print $0 ",,,"
        }
    }
' "$TMPDIR/gcf_annot.tsv" "$SUMMARY" > "$TMPDIR/summary_body.csv"

# 加表头
{
    echo "菌株名称,保守簇数量,保守簇ID列表,species,family,genus"
    cat "$TMPDIR/summary_body.csv"
} > "$OUTFILE"

echo ""
echo "完成！输出文件: $OUTFILE"
echo ""
echo "预览:"
head -5 "$OUTFILE"
