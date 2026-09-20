import os
import re
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser(description='Extract homologous gene information from BLAST results')
parser.add_argument('--blast_dir', type=str,
        default='result/gene2human/blastp_results',
        help='Directory containing BLASTx result files')
parser.add_argument('--output', type=str,
        default='result/gene2human/blastp_results/homology_summary.csv',
        help='Output CSV file name')
parser.add_argument('--min_length', type=int, default=1,
        help='Minimum alignment length (default: 1 - no filter)')
parser.add_argument('--max_evalue', type=float, default=9999,
        help='Maximum e-value (default: 9999 - sensitive)')
parser.add_argument('--min_identity', type=float, default=0,
        help='Minimum percent identity (default: 0 - no filter)')

# New arguments for weighted comprehensive score
parser.add_argument('--weight_length', type=float, default=1.0,
        help='Weight for alignment length in comprehensive score (default: 1.0)')
parser.add_argument('--weight_evalue', type=float, default=1.0,
        help='Weight for e-value in comprehensive score (default: 1.0)')
parser.add_argument('--weight_identity', type=float, default=1.0,
        help='Weight for percent identity in comprehensive score (default: 1.0)')

args = parser.parse_args()

def parse_blast_line(line, args):
    """Parse BLAST line with filtering"""
    fields = line.strip().split('\t')
    if len(fields) < 12:
        return None

    # Extract fields
    query_id = fields[0]
    subject_id = fields[1]
    perc_identity = float(fields[2])
    alignment_length = int(fields[3])
    mismatches = int(fields[4])
    gap_opens = int(fields[5])
    q_start = int(fields[6])
    q_end = int(fields[7])
    s_start = int(fields[8])
    s_end = int(fields[9])
    evalue = float(fields[10])
    bit_score = float(fields[11])

    # Apply filters
    if (alignment_length < args.min_length or
        evalue > args.max_evalue or
        perc_identity < args.min_identity):
        return None

    # Extract cluster gene (from query_id)
    if '|' in query_id:
        cluster_gene = query_id.split('|')[1]
    else:
        cluster_gene = query_id

    # Extract human gene name
    subject_parts = subject_id.split('|')
    if len(subject_parts) >= 6:
        human_gene = subject_parts[6]
        transcript = subject_parts[5]
    else:
        human_gene = "Unknown"
        transcript = "Unknown"

    # Extract conserved ID
    conserved_id = query_id.split('|')[0] if '|' in query_id else query_id

    return {
        'conserved_id': conserved_id,
        'cluster_gene': cluster_gene,
        'human_gene': human_gene,
        'transcript': transcript,
        'alignment_length': alignment_length,
        'perc_identity': perc_identity,
        'evalue': evalue,
        'bit_score': bit_score,
        'q_start': q_start,
        'q_end': q_end,
        's_start': s_start,
        's_end': s_end
    }

def calculate_comprehensive_score(hit, max_length, min_evalue, max_identity, weights):
    """Calculate weighted comprehensive score for a hit"""
    # Avoid division by zero
    max_length = max_length if max_length > 0 else 1
    min_evalue = min_evalue if min_evalue > 0 else 1e-300
    max_identity = max_identity if max_identity > 0 else 1
    
    # Normalize each metric to 0-1 scale
    length_score = hit['alignment_length'] / max_length
    evalue_score = -np.log10(hit['evalue'] + 1e-300) / -np.log10(min_evalue + 1e-300)
    identity_score = hit['perc_identity'] / max_identity
    
    # Calculate weighted average
    total_weight = weights['length'] + weights['evalue'] + weights['identity']
    weighted_score = (weights['length'] * length_score + 
                     weights['evalue'] * evalue_score + 
                     weights['identity'] * identity_score) / total_weight
    
    return weighted_score

def select_best_hits_by_all_criteria(hits, weights):
    """
    Select best hits using four different criteria:
    1. Best by length
    2. Best by e-value
    3. Best by identity
    4. Best by weighted comprehensive score
    Returns a dictionary with all four best hits including their metrics
    """
    if not hits:
        return None
    
    best_hits = {}
    
    # Best by length
    best_hits['by_length'] = sorted(hits, key=lambda x: -x['alignment_length'])[0]
    
    # Best by e-value
    best_hits['by_evalue'] = sorted(hits, key=lambda x: x['evalue'])[0]
    
    # Best by identity
    best_hits['by_identity'] = sorted(hits, key=lambda x: -x['perc_identity'])[0]
    
    # Best by weighted comprehensive score
    max_length = max(h['alignment_length'] for h in hits)
    min_evalue = min(h['evalue'] for h in hits)
    max_identity = max(h['perc_identity'] for h in hits)
    
    for hit in hits:
        hit['comprehensive_score'] = calculate_comprehensive_score(
            hit, max_length, min_evalue, max_identity, weights
        )
    
    best_hits['by_comprehensive'] = sorted(hits, key=lambda x: -x['comprehensive_score'])[0]
    
    return best_hits

print(f"Configuration:")
print(f"  BLAST directory: {args.blast_dir}")
print(f"  Min alignment length: {args.min_length}")
print(f"  Max e-value: {args.max_evalue}")
print(f"  Min identity: {args.min_identity}%")
print(f"\nComprehensive score weights:")
print(f"  Length weight: {args.weight_length}")
print(f"  E-value weight: {args.weight_evalue}")
print(f"  Identity weight: {args.weight_identity}")
print()

# Process files
blast_files = sorted(Path(args.blast_dir).glob("Conserved_*_blast.txt"))
print(f"Found {len(blast_files)} BLAST result files")

# Weights dictionary for comprehensive score
weights = {
    'length': args.weight_length,
    'evalue': args.weight_evalue,
    'identity': args.weight_identity
}

# Store results for each conserved_id and cluster_gene
results_dict = {}

for blast_file in blast_files:
    conserved_id = blast_file.stem.replace('_blast', '')
    print(f"\nProcessing {conserved_id}...")

    # Group hits by bacterial gene (cluster_gene)
    gene_hits = defaultdict(list)

    # Read all hits from this file
    with open(blast_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                result = parse_blast_line(line, args)
                if result:
                    gene_hits[result['cluster_gene']].append(result)

    print(f"  Found {len(gene_hits)} unique bacterial genes:")

    # Process each bacterial gene separately
    for cluster_gene, hits in gene_hits.items():
        print(f"    Gene: {cluster_gene}")
        print(f"      Total hits: {len(hits)}")

        # Get best hits by all four criteria
        best_hits_dict = select_best_hits_by_all_criteria(hits, weights)
        
        if best_hits_dict:
            # Store results with all four best hits
            if conserved_id not in results_dict:
                results_dict[conserved_id] = {}
            
            results_dict[conserved_id][cluster_gene] = {
                'length_hit': best_hits_dict['by_length'],
                'evalue_hit': best_hits_dict['by_evalue'],
                'identity_hit': best_hits_dict['by_identity'],
                'comprehensive_hit': best_hits_dict['by_comprehensive']
            }
            
            print(f"      Best by length: {best_hits_dict['by_length']['human_gene']} "
                  f"(len={best_hits_dict['by_length']['alignment_length']})")
            print(f"      Best by e-value: {best_hits_dict['by_evalue']['human_gene']} "
                  f"(e={best_hits_dict['by_evalue']['evalue']:.2e})")
            print(f"      Best by identity: {best_hits_dict['by_identity']['human_gene']} "
                  f"(id={best_hits_dict['by_identity']['perc_identity']:.1f}%)")
            print(f"      Best by comprehensive: {best_hits_dict['by_comprehensive']['human_gene']} "
                  f"(score={best_hits_dict['by_comprehensive']['comprehensive_score']:.3f})")

    # If no hits at all for this conserved_id
    if not gene_hits:
        print(f"  No hits passing filters for any bacterial gene")
        if conserved_id not in results_dict:
            results_dict[conserved_id] = {}
        
        # Add a placeholder for no hits
        no_hit_entry = {
            'conserved_id': conserved_id,
            'cluster_gene': f"No_hit_in_{conserved_id}",
            'human_gene': 'No hit',
            'transcript': 'No hit',
            'alignment_length': 0,
            'perc_identity': 0,
            'evalue': 1.0,
            'bit_score': 0,
            'q_start': 0,
            'q_end': 0,
            's_start': 0,
            's_end': 0
        }
        
        results_dict[conserved_id]['no_hit'] = {
            'length_hit': no_hit_entry,
            'evalue_hit': no_hit_entry,
            'identity_hit': no_hit_entry,
            'comprehensive_hit': no_hit_entry
        }

# Create the main summary DataFrame with four columns
summary_rows = []

for conserved_id, genes_dict in results_dict.items():
    for cluster_gene, hits_dict in genes_dict.items():
        row = {
            'Conserved ID': conserved_id,
            'Cluster Gene': cluster_gene,
            
            # Length-based best hit
            'Human Prot (by Length)': hits_dict['length_hit']['human_gene'],
            'Prot Length (bp) (by Length)': hits_dict['length_hit']['alignment_length'],
            'Prot Identity% (by Length)': hits_dict['length_hit']['perc_identity'],
            'Prot E-value (by Length)': hits_dict['length_hit']['evalue'],
            
            # E-value-based best hit
            'Human Prot (by E-value)': hits_dict['evalue_hit']['human_gene'],
            'Prot Length (bp) (by E-value)': hits_dict['evalue_hit']['alignment_length'],
            'Prot Identity% (by E-value)': hits_dict['evalue_hit']['perc_identity'],
            'Prot E-value (by E-value)': hits_dict['evalue_hit']['evalue'],
            
            # Identity-based best hit
            'Human Prot (by Identity)': hits_dict['identity_hit']['human_gene'],
            'Prot Length (bp) (by Identity)': hits_dict['identity_hit']['alignment_length'],
            'Prot Identity% (by Identity)': hits_dict['identity_hit']['perc_identity'],
            'Prot E-value (by Identity)': hits_dict['identity_hit']['evalue'],
            
            # Comprehensive-based best hit
            'Human Prot (by Comprehensive)': hits_dict['comprehensive_hit']['human_gene'],
            'Prot Length (bp) (by Comprehensive)': hits_dict['comprehensive_hit']['alignment_length'],
            'Prot Identity% (by Comprehensive)': hits_dict['comprehensive_hit']['perc_identity'],
            'Prot E-value (by Comprehensive)': hits_dict['comprehensive_hit']['evalue'],
        }
        
        # Add comprehensive score if available
        if 'comprehensive_score' in hits_dict['comprehensive_hit']:
            row['Prot Comprehensive Score'] = hits_dict['comprehensive_hit']['comprehensive_score']
        
        summary_rows.append(row)

# Create DataFrame
summary_df = pd.DataFrame(summary_rows)

# Sort by Conserved ID and Cluster Gene
summary_df = summary_df.sort_values(['Conserved ID', 'Cluster Gene'])

# Save the main summary
summary_df.to_csv(args.output, index=False)
print(f"\n{'='*50}")
print(f"Main summary saved to: {args.output}")

# Also create a simplified version with just the four gene columns
simple_df = summary_df[['Conserved ID', 'Cluster Gene', 
                        'Human Prot (by Length)', 
                        'Human Prot (by E-value)', 
                        'Human Prot (by Identity)', 
                        'Human Prot (by Comprehensive)']].copy()
simple_output = args.output.replace('.csv', '_simple.csv')
simple_df.to_csv(simple_output, index=False)
print(f"Simple summary (genes only) saved to: {simple_output}")

# Create detailed DataFrames for each criteria (for compatibility)
detailed_dfs = {}
for criteria, col_prefix in [('length', 'by Length'), 
                            ('evalue', 'by E-value'), 
                            ('identity', 'by Identity'), 
                            ('comprehensive', 'by Comprehensive')]:
    
    detailed_rows = []
    for _, row in summary_df.iterrows():
        detailed_rows.append({
            'conserved_id': row['Conserved ID'],
            'cluster_gene': row['Cluster Gene'],
            'human_prot': row[f'Human Prot ({col_prefix})'],
            'alignment_length': row[f'Prot Length (bp) ({col_prefix})'] if f'Prot Length (bp) ({col_prefix})' in row else 0,
            'perc_identity': row[f'Prot Identity% ({col_prefix})'] if f'Prot Identity% ({col_prefix})' in row else 0,
            'evalue': row[f'Prot E-value ({col_prefix})'] if f'Prot E-value ({col_prefix})' in row else 1.0
        })
    
    detailed_dfs[criteria] = pd.DataFrame(detailed_rows)

# Save detailed files for each criteria
for criteria in ['length', 'evalue', 'identity', 'comprehensive']:
    if criteria in detailed_dfs:
        detailed_output = args.output.replace('.csv', f'_{criteria}_detailed.csv')
        detailed_dfs[criteria].to_csv(detailed_output, index=False)
        print(f"  {criteria.capitalize()} detailed: {detailed_output}")

# Statistics
print(f"\n{'='*50}")
print("Summary Statistics:")
print(f"{'='*50}")
print(f"Total Conserved IDs processed: {len(blast_files)}")
print(f"Total bacterial gene entries: {len(summary_df)}")

# Count hits for each criteria
for criteria, col_name in [('Length', 'Human Prot (by Length)'),
                          ('E-value', 'Human Prot (by E-value)'),
                          ('Identity', 'Human Prot (by Identity)'),
                          ('Comprehensive', 'Human Prot (by Comprehensive)')]:
    hits_count = (summary_df[col_name] != 'No hit').sum()
    print(f"\n{criteria} criteria:")
    print(f"  Bacterial genes with hits: {hits_count}/{len(summary_df)} ({hits_count/len(summary_df)*100:.1f}%)")
    
    if hits_count > 0:
        hits_df = summary_df[summary_df[col_name] != 'No hit']
        if criteria == 'Length':
            print(f"  Average length: {hits_df['Prot Length (bp) (by Length)'].mean():.1f} bp")
        elif criteria == 'E-value':
            print(f"  Median e-value: {hits_df['Prot E-value (by E-value)'].median():.2e}")
        elif criteria == 'Identity':
            print(f"  Average identity: {hits_df['Prot Identity% (by Identity)'].mean():.1f}%")

# Show sample results
print(f"\n{'='*50}")
print("Sample of results (first 10 rows):")
print(f"{'='*50}")
display_cols = ['Conserved ID', 'Cluster Gene', 
                'Human Prot (by Length)', 
                'Human Prot (by E-value)', 
                'Human Prot (by Identity)', 
                'Human Prot (by Comprehensive)']
print(summary_df[display_cols].head(10).to_string(index=False))

# Also save a clean tab-separated version
clean_output = args.output.replace('.csv', '_clean.txt')
with open(clean_output, 'w') as f:
    # Write header
    f.write('\t'.join(display_cols) + '\n')
    # Write data
    for _, row in summary_df[display_cols].head(100).iterrows():  # Limit to first 100 for clean view
        f.write('\t'.join(str(val) for val in row) + '\n')
print(f"\nClean tab-separated version saved to: {clean_output}")

# Create comparison summary - show when different criteria select different genes
print(f"\n{'='*50}")
print("Criteria Agreement Analysis:")
print(f"{'='*50}")

# Check if all four criteria selected the same gene
summary_df['All_Same'] = (
    (summary_df['Human Prot (by Length)'] == summary_df['Human Prot (by E-value)']) &
    (summary_df['Human Prot (by Length)'] == summary_df['Human Prot (by Identity)']) &
    (summary_df['Human Prot (by Length)'] == summary_df['Human Prot (by Comprehensive)'])
)

all_same_count = summary_df['All_Same'].sum()
print(f"Entries where all criteria agree: {all_same_count}/{len(summary_df)} ({all_same_count/len(summary_df)*100:.1f}%)")

# Show examples where they disagree
if all_same_count < len(summary_df):
    print(f"\nExamples where criteria select different genes (first 5):")
    disagreement_df = summary_df[~summary_df['All_Same']].head(5)
    for _, row in disagreement_df.iterrows():
        print(f"\n  Conserved: {row['Conserved ID']}, Gene: {row['Cluster Gene']}")
        print(f"    Length: {row['Human Prot (by Length)']} (len={row['Prot Length (bp) (by Length)']})")
        print(f"    E-value: {row['Human Prot (by E-value)']} (e={row['Prot E-value (by E-value)']:.2e})")
        print(f"    Identity: {row['Human Prot (by Identity)']} (id={row['Prot Identity% (by Identity)']:.1f}%)")
        print(f"    Comprehensive: {row['Human Prot (by Comprehensive)']}")

print(f"\n{'='*50}")
print("Processing complete!")
print(f"{'='*50}")
