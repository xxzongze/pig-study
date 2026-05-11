#!/usr/bin/env python3
"""
Supplementary: Standard exploratory figures (PCA, Volcano, Clustering)
for pig liver-muscle transcriptome — JASB journal style.

48 samples: DLY/TFB × 15/45/75/105kg × 6 replicates
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, pearsonr
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os, warnings
warnings.filterwarnings('ignore')

# ============================================================
# Global style
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.weight': 'bold',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.linewidth': 0.6,
})

C_RED    = '#D73027'   # DLY
C_BLUE   = '#4575B4'   # TFB
C_PURPLE = '#762A83'   # combined/highlight
C_BG     = '#FFFFFF'
STAGE_COLORS = {15: '#FEE090', 45: '#FDAE61', 75: '#F46D43', 105: '#A50026'}
STAGE_SHAPES = {15: 'o', 45: 's', 75: '^', 105: 'D'}

os.makedirs('figures', exist_ok=True)

# ============================================================
# Load data
# ============================================================
print("Loading data...")
liver_expr = pd.read_csv('wgcna_output/liver_expr.csv', index_col=0)
muscle_expr= pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)
liver_gm   = pd.read_csv('wgcna_output/liver_gene_module_assignment.csv')
muscle_gm  = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')
liver_mtc  = pd.read_csv('wgcna_output/liver_module_trait_cor.csv', index_col=0)
muscle_mtc = pd.read_csv('wgcna_output/muscle_module_trait_cor.csv', index_col=0)

def parse_sample(sample_name):
    parts = sample_name.split('_')
    breed = parts[0]
    stage = int(parts[1].replace('kg', ''))
    rep   = parts[2] if len(parts) > 2 else '1'
    return breed, stage, rep

# Build metadata
meta_liver = pd.DataFrame({
    'Breed': [parse_sample(s)[0] for s in liver_expr.index],
    'Stage': [parse_sample(s)[1] for s in liver_expr.index],
    'Rep':   [parse_sample(s)[2] for s in liver_expr.index],
    'Group': [f"{parse_sample(s)[0]}_{parse_sample(s)[1]}kg" for s in liver_expr.index],
}, index=liver_expr.index)

meta_muscle = pd.DataFrame({
    'Breed': [parse_sample(s)[0] for s in muscle_expr.index],
    'Stage': [parse_sample(s)[1] for s in muscle_expr.index],
    'Rep':   [parse_sample(s)[2] for s in muscle_expr.index],
    'Group': [f"{parse_sample(s)[0]}_{parse_sample(s)[1]}kg" for s in muscle_expr.index],
}, index=muscle_expr.index)

# ============================================================
# Fig 1: PCA — Liver + Muscle (JASB style: Chen Fig 4A, Jia Fig style)
# ============================================================
print("Generating Fig 1: PCA...")

fig1, axes1 = plt.subplots(1, 2, figsize=(14, 6))

for idx, (tissue, expr, meta) in enumerate([
    ('Liver', liver_expr, meta_liver),
    ('Muscle', muscle_expr, meta_muscle)
]):
    ax = axes1[idx]

    # Use top 2000 most variable genes for PCA
    var_genes = expr.var().nlargest(2000).index
    expr_var = expr[var_genes]
    # Center
    expr_centered = expr_var - expr_var.mean()
    # SVD
    U, S, Vt = np.linalg.svd(expr_centered.values, full_matrices=False)
    pc1 = U[:, 0] * S[0]
    pc2 = U[:, 1] * S[1]
    var_pc1 = S[0]**2 / (S**2).sum() * 100
    var_pc2 = S[1]**2 / (S**2).sum() * 100

    # Scatter by stage (shape) and breed (color)
    for stage in [15, 45, 75, 105]:
        for breed, color in [('DLY', C_RED), ('TFB', C_BLUE)]:
            mask = (meta['Stage'] == stage) & (meta['Breed'] == breed)
            s = STAGE_SHAPES[stage]
            ax.scatter(pc1[mask], pc2[mask], c=color, marker=s, s=55,
                       edgecolors='black', linewidth=0.3, alpha=0.85, zorder=3,
                       label=f'{breed} {stage}kg' if stage == 15 else '')

    # Add 95% ellipses for breed
    for breed, color in [('DLY', C_RED), ('TFB', C_BLUE)]:
        mask = meta['Breed'] == breed
        ax.scatter(pc1[mask].mean(), pc2[mask].mean(), c=color, marker='X', s=120,
                   edgecolors='black', linewidth=1.0, zorder=5)

    ax.set_xlabel(f'PC1 ({var_pc1:.1f}%)', fontsize=9)
    ax.set_ylabel(f'PC2 ({var_pc2:.1f}%)', fontsize=9)
    ax.set_title(f'{tissue}', fontweight='bold', fontsize=10)

    # Custom legend: breeds (color) + stages (shape)
    legend_breeds = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_RED,
                   markersize=8, label='DLY', markeredgecolor='black', markeredgewidth=0.3),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_BLUE,
                   markersize=8, label='TFB', markeredgecolor='black', markeredgewidth=0.3),
    ]
    legend_stages = [
        plt.Line2D([0], [0], marker=s, color='w', markerfacecolor='#666666',
                   markersize=8, label=f'{st} kg', markeredgecolor='black', markeredgewidth=0.3)
        for st, s in STAGE_SHAPES.items()
    ]
    leg1 = ax.legend(handles=legend_breeds, loc='upper left', frameon=False, fontsize=7,
                     title='Breed', title_fontsize=7)
    ax.add_artist(leg1)
    ax.legend(handles=legend_stages, loc='upper right', frameon=False, fontsize=7,
              title='Stage', title_fontsize=7)

fig1.suptitle('Principal Component Analysis of Pig Liver and Skeletal Muscle Transcriptomes',
              fontweight='bold', fontsize=12)
plt.tight_layout()
fig1.savefig('figures/Fig1_pca.png', dpi=300, facecolor=C_BG)
fig1.savefig('figures/Fig1_pca.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig1)
print("  -> figures/Fig1_pca.png|tiff")

# ============================================================
# Fig 2: Volcano Plot — DLY vs TFB per tissue
# ============================================================
print("Generating Fig S4: Volcano Plot...")

def compute_degs(expr, meta, breed_a='DLY', breed_b='TFB', min_expr=1.0):
    """Welch's t-test for each gene between breeds."""
    mask_a = meta['Breed'] == breed_a
    mask_b = meta['Breed'] == breed_b
    results = []
    for gene in expr.columns:
        a = expr.loc[mask_a, gene].values
        b = expr.loc[mask_b, gene].values
        # Filter low expression
        if np.mean(a) < min_expr and np.mean(b) < min_expr:
            continue
        if np.std(a) == 0 and np.std(b) == 0:
            continue
        try:
            t_stat, p_val = ttest_ind(a, b, equal_var=False)
        except:
            continue
        fc = np.log2(np.mean(a) / np.mean(b)) if np.mean(b) > 0 else 0
        results.append({'Gene': gene, 'log2FC': fc, 'pvalue': p_val, 'mean_A': np.mean(a), 'mean_B': np.mean(b)})
    return pd.DataFrame(results)

liver_deg  = compute_degs(liver_expr, meta_liver)
muscle_deg = compute_degs(muscle_expr, meta_muscle)
liver_deg['-log10p'] = -np.log10(liver_deg['pvalue'].clip(lower=1e-50))
muscle_deg['-log10p'] = -np.log10(muscle_deg['pvalue'].clip(lower=1e-50))

fig2, axes2 = plt.subplots(1, 2, figsize=(14, 6))

key_genes_liver = ['AHSG', 'IGFBP3', 'IGF1', 'CPS1', 'SDS', 'ARG1', 'ASS1', 'ASL',
                   'GLUD1', 'STAT3', 'BCKDHA', 'OTC', 'C3', 'APOA1', 'TTR', 'ALB',
                   'G6PC1', 'PCK1', 'CYP1A2']
key_genes_muscle = ['CAV3', 'CCND2', 'ENHO', 'SRPX', 'NREP', 'MYOD1', 'MYOG',
                    'IGF1', 'FST', 'MSTN', 'FOXO1', 'AKT1', 'MTOR', 'DLK1',
                    'ANGPTL4', 'APOE', 'FGFR4']

for idx, (tissue, deg_df, key_genes) in enumerate([
    ('Liver', liver_deg, key_genes_liver),
    ('Muscle', muscle_deg, key_genes_muscle)
]):
    ax = axes2[idx]

    # All genes as background
    n_sig_up = (deg_df['log2FC'] > 1) & (deg_df['pvalue'] < 0.05)
    n_sig_dn = (deg_df['log2FC'] < -1) & (deg_df['pvalue'] < 0.05)

    ax.scatter(deg_df['log2FC'], deg_df['-log10p'], c='#CCCCCC', s=6, alpha=0.3,
               rasterized=True, zorder=1)
    # Sig genes
    ax.scatter(deg_df.loc[n_sig_up, 'log2FC'], deg_df.loc[n_sig_up, '-log10p'],
               c=C_RED, s=12, alpha=0.6, rasterized=True, zorder=2,
               label=f'DLY > TFB ({(n_sig_up).sum()})')
    ax.scatter(deg_df.loc[n_sig_dn, 'log2FC'], deg_df.loc[n_sig_dn, '-log10p'],
               c=C_BLUE, s=12, alpha=0.6, rasterized=True, zorder=2,
               label=f'TFB > DLY ({(n_sig_dn).sum()})')

    # Key genes
    for gene in key_genes:
        g_row = deg_df[deg_df['Gene'] == gene]
        if len(g_row) > 0:
            r = g_row.iloc[0]
            ax.scatter(r['log2FC'], r['-log10p'], c=C_PURPLE, s=50,
                       edgecolors='black', linewidth=0.5, zorder=4)
            ax.annotate(gene, (r['log2FC'], r['-log10p']),
                        fontsize=5.5, ha='center', va='bottom',
                        xytext=(0, 4), textcoords='offset points', fontweight='bold')

    ax.axhline(y=-np.log10(0.05), color='grey', linewidth=0.6, linestyle='--', alpha=0.5)
    ax.axvline(x=1, color='grey', linewidth=0.6, linestyle='--', alpha=0.4)
    ax.axvline(x=-1, color='grey', linewidth=0.6, linestyle='--', alpha=0.4)
    ax.axvline(x=0, color='grey', linewidth=0.4, alpha=0.3)

    ax.set_xlabel('log2 Fold-Change (DLY / TFB)', fontsize=9)
    ax.set_ylabel('-log10 P-value', fontsize=9)
    ax.set_title(f'{tissue}: DLY vs TFB', fontweight='bold', fontsize=10)
    ax.legend(loc='upper right', fontsize=6.5, frameon=False)

fig2.suptitle('Differential Gene Expression Between DLY and TFB Breeds',
              fontweight='bold', fontsize=12)
plt.tight_layout()
fig2.savefig('figures/FigS4_volcano.png', dpi=300, facecolor=C_BG)
fig2.savefig('figures/FigS4_volcano.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig2)
print("  -> figures/FigS4_volcano.png|tiff")

# ============================================================
# Fig 3: Sample-Sample Correlation Heatmap (Liver + Muscle)
# ============================================================
print("Generating Fig S3: Sample Correlation Heatmap...")

fig3, axes3 = plt.subplots(1, 2, figsize=(16, 7))

for idx, (tissue, expr, meta) in enumerate([
    ('Liver', liver_expr, meta_liver),
    ('Muscle', muscle_expr, meta_muscle)
]):
    ax = axes3[idx]

    # Top 1000 variable genes, log2 transform
    var_genes = expr.var().nlargest(1000).index
    expr_sub = np.log2(expr[var_genes] + 1)

    # Sample correlation
    cor_mat = expr_sub.T.corr()

    # Sort samples by breed then stage
    sample_order = meta.sort_values(['Breed', 'Stage']).index.tolist()
    cor_mat = cor_mat.loc[sample_order, sample_order]

    im = ax.imshow(cor_mat.values, aspect='auto', cmap='RdBu_r', vmin=0.7, vmax=1,
                   interpolation='nearest')

    # Breed/stage annotation bars
    meta_sorted = meta.loc[sample_order]
    for i, (s, row) in enumerate(meta_sorted.iterrows()):
        breed_color = C_RED if row['Breed'] == 'DLY' else C_BLUE
        stage_alpha = 0.3 + 0.7 * (row['Stage'] / 105)
        ax.add_patch(plt.Rectangle((-1.8, i - 0.5), 0.6, 1, facecolor=breed_color,
                                   alpha=stage_alpha, edgecolor='none'))
        ax.add_patch(plt.Rectangle((cor_mat.shape[0], i - 0.5), 0.6, 1, facecolor=breed_color,
                                   alpha=stage_alpha, edgecolor='none'))

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f'{tissue} (n={expr.shape[1]} genes)', fontweight='bold', fontsize=10)
    ax.set_xlabel(f'Samples (ordered by Breed × Stage)', fontsize=8)

    cbar = fig3.colorbar(im, ax=ax, shrink=0.8, aspect=25)
    cbar.set_label("Pearson r", fontsize=8)

# Add legend for annotation bars
legend_elements = [
    Patch(facecolor=C_RED, alpha=0.6, label='DLY'),
    Patch(facecolor=C_BLUE, alpha=0.6, label='TFB'),
]
fig3.legend(handles=legend_elements, loc='lower center', ncol=2, frameon=False,
            fontsize=8, bbox_to_anchor=(0.5, -0.02))

fig3.suptitle('Sample–Sample Expression Correlation (Top 1,000 Variable Genes)',
              fontweight='bold', fontsize=12)
plt.tight_layout()
fig3.savefig('figures/FigS3_sample_correlation.png', dpi=300, facecolor=C_BG)
fig3.savefig('figures/FigS3_sample_correlation.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig3)
print("  -> figures/FigS3_sample_correlation.png|tiff")

# ============================================================
# Fig S4: Per-Stage Volcano Grid (DLY vs TFB at each stage)
# ============================================================
print("Generating Fig S5: Per-Stage Volcano Grid...")

figS4, axesS4 = plt.subplots(2, 4, figsize=(18, 9))

for t_idx, (tissue, expr, meta, key_genes) in enumerate([
    ('Liver', liver_expr, meta_liver, key_genes_liver),
    ('Muscle', muscle_expr, meta_muscle, key_genes_muscle)
]):
    for s_idx, stage in enumerate([15, 45, 75, 105]):
        ax = axesS4[t_idx, s_idx]
        mask_stage = meta['Stage'] == stage
        expr_stage = expr.loc[mask_stage]
        meta_stage = meta.loc[mask_stage]
        deg_stage = compute_degs(expr_stage, meta_stage)

        if len(deg_stage) > 0:
            deg_stage['-log10p'] = -np.log10(deg_stage['pvalue'].clip(lower=1e-30))
            sig_up = (deg_stage['log2FC'] > 1) & (deg_stage['pvalue'] < 0.05)
            sig_dn = (deg_stage['log2FC'] < -1) & (deg_stage['pvalue'] < 0.05)

            ax.scatter(deg_stage['log2FC'], deg_stage['-log10p'], c='#CCCCCC', s=4, alpha=0.3)
            ax.scatter(deg_stage.loc[sig_up, 'log2FC'], deg_stage.loc[sig_up, '-log10p'],
                       c=C_RED, s=10, alpha=0.6, label=f'n={sig_up.sum()}')
            ax.scatter(deg_stage.loc[sig_dn, 'log2FC'], deg_stage.loc[sig_dn, '-log10p'],
                       c=C_BLUE, s=10, alpha=0.6, label=f'n={sig_dn.sum()}')

            for gene in key_genes[:6]:
                g = deg_stage[deg_stage['Gene'] == gene]
                if len(g) > 0:
                    r = g.iloc[0]
                    ax.scatter(r['log2FC'], r['-log10p'], c=C_PURPLE, s=30,
                               edgecolors='black', linewidth=0.3, zorder=4)
                    if abs(r['log2FC']) > 0.5 or r['pvalue'] < 0.05:
                        ax.annotate(gene, (r['log2FC'], r['-log10p']),
                                    fontsize=5, ha='center', va='bottom',
                                    xytext=(0, 3), textcoords='offset points')

        ax.axhline(y=-np.log10(0.05), color='grey', linewidth=0.5, linestyle='--', alpha=0.4)
        ax.axvline(x=0, color='grey', linewidth=0.4, alpha=0.3)
        ax.set_title(f'{tissue} {stage}kg', fontsize=9, fontweight='bold')
        if s_idx == 0:
            ax.set_ylabel('-log10 P', fontsize=8)
        if t_idx == 1:
            ax.set_xlabel('log2FC (DLY/TFB)', fontsize=8)
        if t_idx == 0:
            ax.legend(loc='upper right', fontsize=5, frameon=False)

figS4.suptitle('Per-Stage Differential Expression: DLY vs TFB',
               fontweight='bold', fontsize=12)
plt.tight_layout()
figS4.savefig('figures/FigS5_per_stage_volcano.png', dpi=300, facecolor=C_BG)
plt.close(figS4)
print("  -> figures/FigS5_per_stage_volcano.png")

# ============================================================
# Fig S5: Top DEG Heatmap (DLY vs TFB overall)
# ============================================================
print("Generating Fig S6: Top DEG Heatmap...")

figS5, axesS5 = plt.subplots(1, 2, figsize=(16, 8))

for idx, (tissue, expr, meta, deg_df) in enumerate([
    ('Liver', liver_expr, meta_liver, liver_deg),
    ('Muscle', muscle_expr, meta_muscle, muscle_deg)
]):
    ax = axesS5[idx]

    # Top 50 DEGs by |log2FC| * -log10p
    deg_df = deg_df.copy()
    deg_df['score'] = abs(deg_df['log2FC']) * deg_df['-log10p']
    top50 = deg_df.nlargest(50, 'score')

    # Build heatmap data
    heat_data = np.log2(expr[top50['Gene'].tolist()] + 1)
    heat_data = heat_data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)  # z-score per gene

    # Sort samples
    sample_order = meta.sort_values(['Breed', 'Stage']).index.tolist()
    heat_data = heat_data.loc[sample_order]

    im = ax.imshow(heat_data.T.values, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2,
                   interpolation='nearest')

    # Gene labels
    ax.set_yticks(range(len(top50)))
    ax.set_yticklabels(top50['Gene'].tolist(), fontsize=5.5)

    # Breed/Stage color bar at bottom
    meta_sorted = meta.loc[sample_order]
    for i, (s, row) in enumerate(meta_sorted.iterrows()):
        breed_color = C_RED if row['Breed'] == 'DLY' else C_BLUE
        ax.add_patch(plt.Rectangle((i - 0.5, -3.5), 1, 1.5, facecolor=breed_color,
                                   alpha=0.7, edgecolor='none', clip_on=False))

    ax.set_xticks([])
    # Stage divider lines
    prev_breed = None
    divider_positions = []
    for i, (s, row) in enumerate(meta_sorted.iterrows()):
        if row['Breed'] != prev_breed:
            ax.axvline(x=i - 0.5, color='black', linewidth=1.0)
            prev_breed = row['Breed']

    ax.set_title(f'{tissue}: Top 50 DEGs (DLY vs TFB)', fontweight='bold', fontsize=10)

    cbar = figS5.colorbar(im, ax=ax, shrink=0.8, aspect=25)
    cbar.set_label('Z-score', fontsize=7)

figS5.text(0.5, 0.01, 'DLY (red bar)                                           TFB (blue bar)',
           ha='center', fontsize=8, fontweight='bold',
           transform=figS5.transFigure)
figS5.suptitle('Expression Patterns of Top Differentially Expressed Genes',
               fontweight='bold', fontsize=12)
plt.tight_layout()
figS5.subplots_adjust(bottom=0.1)
figS5.savefig('figures/FigS6_deg_heatmap.png', dpi=300, facecolor=C_BG)
plt.close(figS5)
print("  -> figures/FigS6_deg_heatmap.png")

# ============================================================
# Fig S6: Number of DEGs per stage summary
# ============================================================
print("Generating Fig S7: DEG Summary...")

figS6, axS6 = plt.subplots(figsize=(10, 5))

deg_summary = []
for tissue, expr, meta in [
    ('Liver', liver_expr, meta_liver),
    ('Muscle', muscle_expr, meta_muscle)
]:
    for stage in [15, 45, 75, 105]:
        mask = meta['Stage'] == stage
        deg = compute_degs(expr.loc[mask], meta.loc[mask])
        n_up = ((deg['log2FC'] > 1) & (deg['pvalue'] < 0.05)).sum()
        n_dn = ((deg['log2FC'] < -1) & (deg['pvalue'] < 0.05)).sum()
        deg_summary.append({
            'Tissue': tissue, 'Stage': stage,
            'UP (DLY>T)': n_up, 'DOWN (TFB>D)': n_dn
        })
deg_sum_df = pd.DataFrame(deg_summary)

x_positions = np.arange(4)
bar_w = 0.3

for i, (tissue, offset, color_up, color_dn) in enumerate([
    ('Liver', -bar_w/2, '#D73027', '#FC8D59'),
    ('Muscle', +bar_w/2, '#4575B4', '#91BFDB')
]):
    sub = deg_sum_df[deg_sum_df['Tissue'] == tissue]
    axS6.bar(x_positions + offset - bar_w/3, sub['UP (DLY>T)'], bar_w*0.6,
             color=color_up, alpha=0.85, edgecolor='white', label=f'{tissue} DLY>T')
    axS6.bar(x_positions + offset + bar_w/3, -sub['DOWN (TFB>D)'], bar_w*0.6,
             color=color_dn, alpha=0.85, edgecolor='white', label=f'{tissue} TFB>D')

    # Annotate
    for j, (_, r) in enumerate(sub.iterrows()):
        axS6.text(x_positions[j] + offset - bar_w/3, r['UP (DLY>T)'] + 5,
                  str(r['UP (DLY>T)']), fontsize=6, ha='center', color=color_up, fontweight='bold')
        axS6.text(x_positions[j] + offset + bar_w/3, -r['DOWN (TFB>D)'] - 5,
                  str(r['DOWN (TFB>D)']), fontsize=6, ha='center', color=color_dn, fontweight='bold')

axS6.axhline(y=0, color='black', linewidth=0.5)
axS6.set_xticks(x_positions)
axS6.set_xticklabels(['15 kg', '45 kg', '75 kg', '105 kg'], fontsize=9)
axS6.set_ylabel('Number of DEGs (|log2FC| > 1, P < 0.05)', fontsize=9)
axS6.set_title('Differentially Expressed Genes: DLY vs TFB at Each Growth Stage',
               fontweight='bold', fontsize=11)
# Simplify legend (4 entries)
handles, labels = axS6.get_legend_handles_labels()
unique_labels = []
unique_handles = []
for h, l in zip(handles, labels):
    if l not in unique_labels:
        unique_labels.append(l)
        unique_handles.append(h)
axS6.legend(unique_handles, unique_labels, ncol=2, fontsize=7, frameon=False)
figS6.tight_layout()
figS6.savefig('figures/FigS7_deg_summary.png', dpi=300, facecolor=C_BG)
plt.close(figS6)
print("  -> figures/FigS7_deg_summary.png")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("EXPLORATORY FIGURES GENERATED")
print("=" * 60)
print("""
Main Figures:
  Fig 1: PCA (Liver + Muscle, breed × stage)

Supplementary:
  Fig S3: Sample-Sample Correlation Heatmap
  Fig S4: Volcano Plot (DLY vs TFB, overall)
  Fig S5: Per-Stage Volcano Grid (4 stages × 2 tissues)
  Fig S6: Top 50 DEG Heatmap
  Fig S7: DEG Count Summary per Stage
""")
print("Done!")
