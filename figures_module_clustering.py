#!/usr/bin/env python3
"""
Fig S14: Module Gene Clustering Heatmap
Shows expression patterns of genes within top PD-positive modules,
validating that these are genuine co-expression clusters.
"""
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os, warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.weight': 'bold',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 5, 'ytick.labelsize': 5,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

C_RED  = '#D73027'
C_BLUE = '#4575B4'
C_BG   = '#FFFFFF'

os.makedirs('figures', exist_ok=True)

print("Loading data...")
muscle_expr = pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)
muscle_gm   = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')
muscle_mtc  = pd.read_csv('wgcna_output/muscle_module_trait_cor.csv', index_col=0)

def parse_sample(s):
    parts = s.split('_')
    return parts[0], int(parts[1].replace('kg', ''))

meta = pd.DataFrame({
    'Breed': [parse_sample(s)[0] for s in muscle_expr.index],
    'Stage': [parse_sample(s)[1] for s in muscle_expr.index],
}, index=muscle_expr.index)

# Top 4 PD-positive modules
pd_pos = [('green', +0.685), ('lightcyan', +0.678),
          ('lightgreen', +0.579), ('greenyellow', +0.529)]

# Select up to 50 genes per module (top by kME)
module_genes = {}
for mod, r_val in pd_pos:
    mod_df = muscle_gm[muscle_gm['Module'] == mod].dropna(subset=['kME_module'])
    top50 = mod_df.nlargest(min(50, len(mod_df)), 'kME_module')['Gene'].tolist()
    top50_in_expr = [g for g in top50 if g in muscle_expr.columns]
    module_genes[mod] = top50_in_expr
    print(f"  {mod}: {len(top50_in_expr)} genes selected for heatmap")

# ============================================================
# FIGURE: 4-panel module clustering heatmaps
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(16, 14))

# Sort samples: by breed then stage
sample_order = sorted(range(len(meta)), key=lambda i: (0 if meta.iloc[i]['Breed']=='DLY' else 1, meta.iloc[i]['Stage']))
sample_labels = [f"{'D' if meta.iloc[s]['Breed']=='DLY' else 'T'}{meta.iloc[s]['Stage']}" for s in sample_order]

# Breed color bar
breed_colors = [C_RED if meta.iloc[s]['Breed']=='DLY' else C_BLUE for s in sample_order]

for idx, (mod, r_val) in enumerate(pd_pos):
    ax = axes[idx // 2, idx % 2]
    genes = module_genes[mod]

    # Extract expression data
    data = muscle_expr[genes].iloc[sample_order].values.T

    # Cluster genes
    if data.shape[0] > 1:
        gene_order = leaves_list(linkage(pdist(data)))
        data = data[gene_order, :]
        gene_labels = [genes[i] for i in gene_order]
    else:
        gene_labels = genes

    # Z-score normalize by row (gene)
    data_z = (data - data.mean(axis=1, keepdims=True)) / (data.std(axis=1, keepdims=True) + 1e-10)
    vmax = min(3, np.percentile(np.abs(data_z), 99))
    im = ax.imshow(data_z, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax, interpolation='none')

    # Sample breed bar at top
    ax_bottom = ax.get_position()
    # Breed color bar as thin strip at top
    for i, c in enumerate(breed_colors):
        ax.add_patch(plt.Rectangle((i - 0.5, -0.8), 1, 0.6, facecolor=c, clip_on=False, alpha=0.8))

    # Highlight SRPX if in this module
    if 'SRPX' in gene_labels:
        srpx_idx = gene_labels.index('SRPX')
        ax.axhline(y=srpx_idx, color='black', linewidth=2, linestyle='-')
        ax.text(data.shape[1] + 1, srpx_idx, 'SRPX', fontsize=7, fontweight='bold',
                color='black', va='center')

    ax.set_yticks(range(len(gene_labels)))
    ax.set_yticklabels(gene_labels, fontsize=4.5)
    ax.set_xticks([])
    ax.set_xlabel('48 Samples (D=DLY red, T=TFB blue)', fontsize=7, fontweight='bold', labelpad=12)
    n_genes_total = (muscle_gm['Module'] == mod).sum()
    ax.set_title(f'{mod} Module (r_PD={r_val:+.3f})\n{n_genes_total} genes, showing top {len(genes)} by kME',
                 fontweight='bold', fontsize=9)

fig.suptitle('Module Gene Co-Expression Validation: PD-Positive Muscle Modules\n'
             '(Genes clustered by expression pattern | Red=DLY Blue=TFB | SRPX marked with black line)',
             fontweight='bold', fontsize=12, y=1.01)

plt.tight_layout()
fig.savefig('figures/FigS14_module_clustering.png', dpi=300, facecolor=C_BG)
plt.close(fig)
print("  -> figures/FigS14_module_clustering.png")

print("\nDone!")
