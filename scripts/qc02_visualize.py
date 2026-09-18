#!/usr/bin/env python3
"""
qc02_visualize.py
==================
Visualize the QC results from qc01_stats.sh as a four-panel comparison figure.
Reads stats_summary.tsv and produces a single overview PNG suitable for
reports or tutorial figures.

Usage:
    python3 qc02_visualize.py qc/stats_summary.tsv qc/qc_overview.png
"""
import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# Parse command-line arguments
if len(sys.argv) != 3:
    print("Usage: python3 qc02_visualize.py <input.tsv> <output.png>")
    sys.exit(1)

infile, outfile = sys.argv[1], sys.argv[2]

# Load data
df = pd.read_csv(infile, sep='\t')
df['label'] = df['sample'] + '_' + df['hap']

# ---------- Figure setup ----------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle(
    'Citrus Pangenome QC Summary: 3 cultivars x 2 haplotypes = 6 assemblies',
    fontsize=14, fontweight='bold'
)

# Assign one color per cultivar
color_map = {'CUN': '#e74c3c', 'CKI': '#3498db', 'CKU': '#27ae60'}
sample_colors = [color_map.get(s, '#95a5a6') for s in df['sample']]

# ============ Panel A: Total sequence length ============
ax = axes[0, 0]
bars = ax.bar(df['label'], df['sum_len_Mb'], color=sample_colors, alpha=0.85,
              edgecolor='black', linewidth=0.8)
# Hatch pattern to distinguish hap2 from hap1
for bar, hap in zip(bars, df['hap']):
    if hap == 'hap2':
        bar.set_hatch('///')

# Highlight the expected range for haploid Citrus (~290-370 Mb)
ax.axhspan(290, 370, alpha=0.15, color='green', label='Expected range (290-370 Mb)')
ax.set_ylabel('Total length (Mb)')
ax.set_title('A. Total sequence length per assembly', fontweight='bold', loc='left')
ax.tick_params(axis='x', rotation=45)
ax.legend(loc='upper right', fontsize=9)
ax.grid(axis='y', alpha=0.3)

# Value labels on top of bars
for bar, v in zip(bars, df['sum_len_Mb']):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
            f'{v:.0f}', ha='center', fontsize=8)

# ============ Panel B: Contig count (log scale) ============
ax = axes[0, 1]
bars = ax.bar(df['label'], df['num_seqs'], color=sample_colors, alpha=0.85,
              edgecolor='black', linewidth=0.8)
for bar, hap in zip(bars, df['hap']):
    if hap == 'hap2':
        bar.set_hatch('///')

# Reference lines for contiguity thresholds
ax.axhline(9, color='green', linestyle='--', label='Ideal: 9 (chromosome count)')
ax.axhline(20, color='orange', linestyle=':', label='CHR_SCALE threshold')
ax.axhline(200, color='red', linestyle=':', label='FRAGMENTED threshold')
ax.set_ylabel('num_seqs (contig count)')
ax.set_yscale('log')
ax.set_title('B. Assembly contiguity (contig count, log scale)', fontweight='bold', loc='left')
ax.tick_params(axis='x', rotation=45)
ax.legend(loc='upper left', fontsize=9)
ax.grid(axis='y', alpha=0.3, which='both')

# ============ Panel C: N50 ============
ax = axes[1, 0]
bars = ax.bar(df['label'], df['N50_Mb'], color=sample_colors, alpha=0.85,
              edgecolor='black', linewidth=0.8)
for bar, hap in zip(bars, df['hap']):
    if hap == 'hap2':
        bar.set_hatch('///')

# Reference lines for N50 quality
ax.axhline(20, color='green', linestyle='--', label='Chromosome-scale (>= 20 Mb)')
ax.axhline(1, color='orange', linestyle=':', label='Minimum acceptable (1 Mb)')
ax.set_ylabel('N50 (Mb)')
ax.set_title('C. N50: larger is better', fontweight='bold', loc='left')
ax.tick_params(axis='x', rotation=45)
ax.legend(loc='upper right', fontsize=9)
ax.grid(axis='y', alpha=0.3)

# ============ Panel D: GC content ============
ax = axes[1, 1]
bars = ax.bar(df['label'], df['GC_pct'], color=sample_colors, alpha=0.85,
              edgecolor='black', linewidth=0.8)
for bar, hap in zip(bars, df['hap']):
    if hap == 'hap2':
        bar.set_hatch('///')

# Expected GC range for Citrus (34-38%)
ax.axhspan(34, 38, alpha=0.15, color='green', label='Expected Citrus range (34-38%)')
ax.set_ylabel('GC (%)')
ax.set_title('D. GC content: reflects species characteristics', fontweight='bold', loc='left')
ax.tick_params(axis='x', rotation=45)
ax.legend(loc='upper right', fontsize=9)
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(30, 42)

# ---------- Combined legend (color = cultivar, pattern = haplotype) ----------
sample_patches = [
    mpatches.Patch(color=color_map['CUN'], label='Satsuma (CUN)'),
    mpatches.Patch(color=color_map['CKI'], label='Kishu (CKI)'),
    mpatches.Patch(color=color_map['CKU'], label='Kunenbo (CKU)'),
]
hap_patches = [
    mpatches.Patch(facecolor='white', edgecolor='black', label='hap1 (solid)'),
    mpatches.Patch(facecolor='white', edgecolor='black', hatch='///', label='hap2 (hatched)'),
]
fig.legend(handles=sample_patches + hap_patches, loc='lower center',
           ncol=5, bbox_to_anchor=(0.5, -0.02), fontsize=10)

plt.tight_layout(rect=(0, 0.03, 1, 0.96))
plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f"Figure saved: {outfile}")

# ---------- Print summary statistics ----------
print("\n--- Summary ---")
print(df[['sample', 'hap', 'num_seqs', 'sum_len_Mb', 'N50_Mb', 'GC_pct', 'verdict']].to_string(index=False))

# ---------- Print within-cultivar hap1 vs hap2 comparison ----------
print("\n--- Within-cultivar hap1 vs hap2 comparison ---")
for sid in df['sample'].unique():
    sub = df[df['sample'] == sid]
    if len(sub) == 2:
        h1, h2 = sub.iloc[0], sub.iloc[1]
        size_diff = abs(h1['sum_len_Mb'] - h2['sum_len_Mb']) / max(h1['sum_len_Mb'], h2['sum_len_Mb']) * 100
        gc_diff = abs(h1['GC_pct'] - h2['GC_pct'])
        print(f"  {sid}: size diff {size_diff:.1f}%, GC diff {gc_diff:.2f}pp")
