#!/bin/bash
set -euo pipefail

WORKDIR=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/select_bacteria_script
MATCHED=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/select_bacteria_result/matched_gcf_list.txt
ASM=/home/liuzhh/bacteria_cluster/bacteria_genomes/20260201/assembly_summary.tsv
OUTPUTDIR=/home/liuzhh/bacteria_cluster/260811_allgenome_ani/select_bacteria_result/stat_level

mkdir -p "$OUTPUTDIR"

# ---------- 1. GCF -> taxid ----------
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR {if($0~/^#/)next; gsub(/^[ \t]+|[ \t]+$/,"",$0); if($0!="")m[$0]=1; next}
    /^#/{next} FNR==1{next}
    {gcf=$1; sub(/\.[0-9]+$/,"",gcf); if(gcf in m) print gcf,$6,$7,$8}
' "$MATCHED" "$ASM" > "$OUTPUTDIR/matched_gcf_taxid.tsv"

# ---------- 2. taxid list ----------
cut -f3 "$OUTPUTDIR/matched_gcf_taxid.tsv" | sort -u > "$OUTPUTDIR/taxid_list.txt"

# ---------- 3. taxonkit lineage + reformat ----------
cat "$OUTPUTDIR/taxid_list.txt" \
  | taxonkit lineage \
  | taxonkit reformat -f "{k};{p};{c};{o};{f};{g};{s}" \
  > "$OUTPUTDIR/taxid_lineage_reformat.txt"

# ---------- 4. 合并 taxid -> lineage 到 GCF 表（无表头，纯数据）----------
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==FNR {lin[$1]=$3; next}
    {print $1,$2,$3,$4,lin[$3]}
' "$OUTPUTDIR/taxid_lineage_reformat.txt" "$OUTPUTDIR/matched_gcf_taxid.tsv" \
  > "$OUTPUTDIR/matched_gcf_with_lineage.tsv"

# ---------- 5. 拆分 lineage 为 7 级（这里才补表头）----------
{
    printf "GCF\ttaxid\tspecies_taxid\torganism_name\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
    awk -F'\t' 'BEGIN{OFS="\t"}
        {
            n = split($5, a, ";")
            printf "%s\t%s\t%s\t%s", $1, $2, $3, $4
            for (i = 1; i <= 7; i++) printf "\t%s", (i <= n ? a[i] : "")
            printf "\n"
        }
    ' "$OUTPUTDIR/matched_gcf_with_lineage.tsv"
} > "$OUTPUTDIR/matched_gcf_taxonomy_full.tsv"

# ---------- 按科汇总：family_count.tsv ----------
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{next}
    {c[$9]++}
    END{for(k in c) print k, c[k]}
' "$OUTPUTDIR/matched_gcf_taxonomy_full.tsv" \
  | sort -t$'\t' -k2,2nr > "$OUTPUTDIR/family_count.tsv"

# ---------- 按属汇总：genus_count.tsv ----------
awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1{next}
    {c[$10]++}
    END{for(k in c) print k, c[k]}
' "$OUTPUTDIR/matched_gcf_taxonomy_full.tsv" \
  | sort -t$'\t' -k2,2nr > "$OUTPUTDIR/genus_count.tsv"

echo "Done. Outputs in $OUTPUTDIR"
