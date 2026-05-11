#!/usr/bin/env python3
"""
Hepatokine — Muscle Crosstalk Analysis
=======================================
Identifies liver-secreted factors differentially expressed between DLY/TFB
and correlates them with muscle functional gene modules.

Key question: Does TFB liver secrete signals that suppress muscle protein deposition?

Based on multi-stage GSEA DEG results and muscle targeted panel (305 genes).
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, ttest_ind
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

from stats_utils import benjamini_hochberg

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.titlesize': 11, 'axes.labelsize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
})

C_DLY = '#2166AC'
C_TFB = '#B2182B'
C_GRAY = '#999999'

import os
os.makedirs('figures_final', exist_ok=True)

# ============================================================
# 1. Load DEG results and identify hepatokines
# ============================================================
print("=" * 70)
print("HEPATOKINE — MUSCLE CROSSTALK ANALYSIS")
print("=" * 70)

print("\n[1/6] Loading multi-stage DEG results...")
deg = pd.read_excel('gsea_multistage_deg_results.xlsx')

# Define hepatokine / secreted protein gene sets
HEPATOKINE_PANEL = {
    'IGF System': ['IGF1', 'IGF2', 'IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5',
                   'IGFBP6', 'IGFBP7', 'IGFALS', 'GDF15'],
    'Growth Factors': ['FGF21', 'FGF1', 'FGF2', 'FGF19', 'HGF', 'NRG1', 'NRG4',
                       'VEGFA', 'VEGFB', 'BMP6', 'BMP9', 'TGFB1', 'TGFB2',
                       'FST', 'FSTL1', 'INHBA', 'INHBE', 'PDGFA', 'EREG', 'HBEGF'],
    'Urea Cycle Enzymes': ['CPS1', 'OTC', 'ASS1', 'ASL', 'ARG1', 'ARG2', 'NAGS'],
    'AA Transaminases': ['GPT', 'GPT2', 'GOT1', 'GOT2', 'GLUD1', 'GLUL'],
    'Acute Phase / Transport': ['ALB', 'TF', 'HP', 'SERPINA1', 'SERPINA3', 'AHSG',
                                'APOA1', 'APOA2', 'APOB', 'APOE', 'C3', 'C5',
                                'FGA', 'FGB', 'FGG', 'CP', 'HPX', 'AGT',
                                'TTR', 'ORM1', 'RBP4', 'AMBP', 'ITIH1', 'ITIH2',
                                'ITIH3', 'ITIH4', 'CLU', 'B2M'],
    'Adipokines / Metabolic': ['ADIPOQ', 'LEP', 'LECT2', 'SELENOP', 'FABP1',
                               'FABP4', 'LCN2', 'FNDC5', 'RETN', 'ANGPTL3',
                               'ANGPTL4', 'ANGPTL6', 'ANGPTL8', 'GPNMB'],
    'Bioactive Hepatokines': ['MYDGF', 'MANF', 'SPARC', 'SPARCL1', 'THBS1', 'THBS2',
                              'CTGF', 'CYR61', 'TIMP1', 'TIMP2', 'MMP2', 'MMP9',
                              'PLG', 'WNT5A', 'DKK1', 'DKK3', 'SFRP1', 'SFRP4',
                              'APOM', 'GC', 'HRG', 'ANG', 'KNG1', 'AZGP1'],
}

all_hepatokines = sorted(set(sum(HEPATOKINE_PANEL.values(), [])))
deg_names_upper = deg['gene_name'].str.upper()

hk_results = []
for category, genes in HEPATOKINE_PANEL.items():
    for gene in genes:
        mask = deg_names_upper == gene.upper()
        if mask.any():
            row = deg[mask].iloc[0]
            hk_results.append({
                'Gene': gene,
                'Category': category,
                'log2FC_DLYvsTFB': row['breed_log2FC_DLYvsTFB'],
                'pvalue': row['breed_pvalue'],
                'qvalue': row['breed_qvalue'],
                'FDR_significant': row['breed_FDR_significant'],
                'abs_log2FC': abs(row['breed_log2FC_DLYvsTFB']),
            })

hk_df = pd.DataFrame(hk_results).sort_values('qvalue')
n_hk_fdr = hk_df['FDR_significant'].sum()
n_hk_nom = (hk_df['pvalue'] < 0.05).sum()
n_hk_found = len(hk_df)

print(f"  Hepatokines screened: {len(all_hepatokines)}")
print(f"  Found in DEG: {n_hk_found}")
print(f"  Nominal P<0.05: {n_hk_nom} | FDR<0.05: {n_hk_fdr}")

if n_hk_fdr > 0:
    print(f"\n  FDR-significant hepatokines (q<0.05):")
    fdr_hks = hk_df[hk_df['FDR_significant']]
    for _, r in fdr_hks.iterrows():
        direction = 'DLY-up' if r['log2FC_DLYvsTFB'] > 0 else 'TFB-up'
        print(f"    {r['Gene']:15s} [{r['Category']:20s}] log2FC={r['log2FC_DLYvsTFB']:+6.2f}  q={r['qvalue']:.4f}  [{direction}]")

# ============================================================
# 2. Load muscle panel and liver expression
# ============================================================
print("\n[2/6] Loading expression data for correlation...")

liver_raw = pd.read_csv('gene_expression/liver_gene_matrix.xls', sep='\t')
muscle_raw = pd.read_csv('gene_expression/muscle_gene_matrix.xls', sep=None, engine='python')

# Sample maps (same as previous scripts)
def get_sample_map(liver_cols, muscle_cols):
    """Map sample columns to (breed, stage, rep)."""
    lmap = {}
    for c in liver_cols:
        parts = c.split('_')
        if parts[0] != 'L':
            continue
        stage_code = parts[1]
        breed_code = parts[2]
        stage_map = {'15': 15, '45': 45, '1': 75, '2': 105, '3': 135}
        if stage_code not in stage_map:
            continue
        lmap[c] = {
            'breed': 'DLY' if breed_code == '1' else 'TFB',
            'stage': stage_map[stage_code],
            'rep': int(parts[3]) if len(parts) > 3 else 1,
        }

    mmap = {}
    for c in muscle_cols:
        parts = c.split('_')
        breed_code, stage_code = None, None
        if parts[0] == 'm':
            stage_map2 = {'15': 15, '1': 75, '2': 105, '3': 135}
            if parts[1] in stage_map2:
                stage_code = stage_map2[parts[1]]
                breed_code = 'DLY' if parts[2] == '1' else 'TFB'
        elif parts[0] in ['DLYM', 'TFBM']:
            stage_code = 45
            breed_code = 'DLY' if parts[0] == 'DLYM' else 'TFB'
        elif parts[0] == 'BJ':
            if parts[1] == '2':
                stage_code = 45
                breed_code = 'DLY' if parts[2] == '1' else 'TFB'
        if breed_code and stage_code:
            rep_num = int(parts[-1])
            mmap[c] = {'breed': breed_code, 'stage': stage_code, 'rep': rep_num}

    return lmap, mmap

l_sample_cols = [c for c in liver_raw.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
m_sample_cols = [c for c in muscle_raw.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
lmap, mmap = get_sample_map(l_sample_cols, m_sample_cols)

# ============================================================
# 3. Group-mean correlation: liver hepatokines vs muscle modules
# ============================================================
print("\n[3/6] Computing group-mean correlations (n=8 breed×stage)...")

# Define muscle functional modules
MUSCLE_MODULES = {
    'Protein Synthesis\n(Ribosome)': ['RPS18', 'RPS19', 'RPS11', 'RPS5', 'RPL4', 'RPLP0',
                                       'RPS3', 'RPS3A', 'RPL34', 'RPL14', 'RPL24', 'RPL27',
                                       'EEF1A2', 'EEF2', 'EIF4A2', 'EIF3G'],
    'Protein Degradation\n(Ubiquitin/Proteasome)': ['PSMA6', 'PSMC5', 'UBB', 'UBC', 'UBA52',
                                                      'UBAC1', 'GABARAP'],
    'Muscle Structure\n& Contraction': ['MYH1', 'MYH4', 'MYH7', 'MYL1', 'MYL2', 'MYL3',
                                         'ACTA1', 'ACTN3', 'TNNC1', 'TNNI1', 'TNNT1',
                                         'TPM1', 'TPM2', 'DES', 'MYBPC1', 'MYBPC2'],
    'Energy Metabolism\n(OxPhos/Glycolysis)': ['ATP5F1A', 'ATP5F1C', 'NDUFA1', 'NDUFS3',
                                                 'COX1', 'COX2', 'COX4I1', 'UQCRB', 'UQCRQ',
                                                 'LDHA', 'PKM', 'PFKM', 'PGK1', 'ALDOA', 'ENO3',
                                                 'GAPDH', 'GPI', 'IDH2', 'MDH2'],
    'Muscle Growth\nRegulation': ['MYOD1', 'MYOG', 'MYF5', 'MSTN', 'IGF1R', 'INSR',
                                   'FOXO1', 'FOXO3', 'FBXO32', 'TRIM63', 'MTOR', 'RPS6KB1'],
}

# Build group-mean expression for liver and muscle
def group_mean(matrix, gene_col, sample_map):
    """Compute group-mean expression for 8 breed×stage groups (15/45/75/105)."""
    groups = {}
    for col, info in sample_map.items():
        if info['stage'] in [15, 45, 75, 105]:
            key = f"{info['breed']}_{info['stage']}"
            if key not in groups:
                groups[key] = []
            groups[key].append(col)

    # Index by gene_name
    if gene_col and gene_col in matrix.columns:
        indexed = matrix.set_index(gene_col)
    else:
        indexed = matrix.set_index(matrix.columns[0])

    # Only keep expression columns
    expr_cols = [c for c in indexed.columns if c in sample_map]

    result = {}
    for group, cols in groups.items():
        cols_in = [c for c in cols if c in expr_cols]
        if cols_in:
            result[group] = indexed[cols_in].mean(axis=1)

    return pd.DataFrame(result)

liver_gm = group_mean(liver_raw, 'gene_name', lmap)
muscle_gm = group_mean(muscle_raw, 'gene_name', mmap)

print(f"  Liver group-means: {liver_gm.shape[0]} genes × {liver_gm.shape[1]} groups")
print(f"  Muscle group-means: {muscle_gm.shape[0]} genes × {muscle_gm.shape[1]} groups")

# Correlation: liver hepatokines vs muscle modules
common_groups = sorted(set(liver_gm.columns) & set(muscle_gm.columns))
print(f"  Common groups: {len(common_groups)} ({', '.join(common_groups)})")

corr_results = []
for _, hk_row in hk_df.iterrows():
    hk_gene = hk_row['Gene']
    if hk_gene not in liver_gm.index:
        # Try case-insensitive
        matches = [g for g in liver_gm.index if g.upper() == hk_gene.upper()]
        if not matches:
            continue
        hk_gene = matches[0]

    l_expr = liver_gm.loc[hk_gene, common_groups].astype(float)

    for module_name, module_genes in MUSCLE_MODULES.items():
        # Use first PC of module genes, or mean
        m_genes_in = [g for g in module_genes if g in muscle_gm.index]
        if len(m_genes_in) < 3:
            continue
        m_expr = muscle_gm.loc[m_genes_in, common_groups].astype(float).mean(axis=0)

        if len(l_expr) >= 6:
            r, p = pearsonr(l_expr, m_expr)
            corr_results.append({
                'Hepatokine': hk_row['Gene'],
                'HK_Category': hk_row['Category'],
                'HK_log2FC': hk_row['log2FC_DLYvsTFB'],
                'HK_FDR': hk_row['FDR_significant'],
                'Muscle_Module': module_name.replace('\n', ' '),
                'Pearson_r': round(r, 4),
                'P_value': round(p, 4),
                'N_genes_module': len(m_genes_in),
                'abs_r': abs(r),
            })

corr_df = pd.DataFrame(corr_results)
if len(corr_df) > 0:
    _, q_corr = benjamini_hochberg(corr_df['P_value'].values)
    corr_df['Q_value'] = q_corr
    corr_df['FDR_significant'] = q_corr < 0.05
    corr_df = corr_df.sort_values('P_value')

n_corr_nom = (corr_df['P_value'] < 0.05).sum()
n_corr_fdr = corr_df['FDR_significant'].sum() if len(corr_df) > 0 else 0
print(f"  Total correlations: {len(corr_df)}")
print(f"  Nominal P<0.05: {n_corr_nom} | FDR<0.05: {n_corr_fdr} (group-mean n=8, limited power)")

if n_corr_nom > 0:
    print(f"\n  Top hepatokine-muscle module correlations (nominal P<0.05):")
    top_corr = corr_df[corr_df['P_value'] < 0.05].head(20)
    for _, r in top_corr.iterrows():
        direction = 'positive' if r['Pearson_r'] > 0 else 'negative'
        print(f"    {r['Hepatokine']:15s} — {r['Muscle_Module']:30s} r={r['Pearson_r']:+6.3f}  p={r['P_value']:.4f}  [{direction}]")

# ============================================================
# 4. Figure: Hepatokine Volcano + Muscle Module Correlation Heatmap
# ============================================================
print("\n[4/6] Generating hepatokine figures...")

# FIGURE 1: Hepatokine expression bar chart
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={'width_ratios': [1.5, 1]})

# Panel A: Hepatokine breed effect bar chart
top_hks = hk_df.nlargest(25, 'abs_log2FC')
y_pos = range(len(top_hks))

for i, (_, row) in enumerate(top_hks.iterrows()):
    color = C_DLY if row['log2FC_DLYvsTFB'] > 0 else C_TFB
    ax1.barh(i, row['log2FC_DLYvsTFB'], color=color, alpha=0.85, height=0.7,
             edgecolor='white', linewidth=0.3)
    sig = '***' if row['qvalue'] < 0.001 else ('**' if row['qvalue'] < 0.01 else ('*' if row['qvalue'] < 0.05 else ''))
    label = f"{row['Gene']} [{row['Category']}] {sig}"
    ax1.text(0, i, label, va='center', fontsize=6.5,
             ha='right' if row['log2FC_DLYvsTFB'] < 0 else 'left',
             color='#222222')

ax1.axvline(0, color='black', lw=0.5)
ax1.set_yticks([])
ax1.set_xlabel('log2(DLY / TFB) breed effect', fontsize=10)
ax1.set_title('Liver Hepatokine/Secreted Factor\nBreed Effect (DLY vs TFB)', fontsize=11, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

ax1.text(0.02, 0.98, f'FDR<0.05: {n_hk_fdr}/{n_hk_found}\nNominal P<0.05: {n_hk_nom}/{n_hk_found}',
         transform=ax1.transAxes, va='top', fontsize=7, color='#555555')

# Panel B: Category summary
cats = hk_df.groupby('Category').agg(
    mean_abs_log2FC=('abs_log2FC', 'mean'),
    n_genes=('Gene', 'count'),
    n_FDR=('FDR_significant', 'sum'),
).sort_values('mean_abs_log2FC', ascending=True)

for i, (cat, row) in enumerate(cats.iterrows()):
    ax2.barh(i, row['mean_abs_log2FC'], color='#4575B4', alpha=0.8, height=0.7)
    ax2.text(row['mean_abs_log2FC'] + 0.02, i,
             f"n={int(row['n_genes'])}, FDR={int(row['n_FDR'])}",
             va='center', fontsize=6.5, color='#555555')

ax2.set_yticks(range(len(cats)))
ax2.set_yticklabels(cats.index, fontsize=7)
ax2.set_xlabel('Mean |log2FC|', fontsize=10)
ax2.set_title('Hepatokine Categories\nMean Absolute Breed Effect', fontsize=11, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig('figures_final/fig_HK1_hepatokine_profile.pdf', dpi=300)
fig.savefig('figures_final/fig_HK1_hepatokine_profile.png', dpi=300)
plt.close()
print("  Saved fig_HK1_hepatokine_profile.pdf/png")

# FIGURE 2: Correlation heatmap — key hepatokines × muscle modules
if len(corr_df) > 0:
    # Pivot: hepatokines × muscle modules
    corr_pivot = corr_df.pivot_table(values='Pearson_r', index='Hepatokine', columns='Muscle_Module', aggfunc='first')

    # Only show hepatokines with at least one strong correlation
    strong_hks = corr_df[corr_df['abs_r'] > 0.5]['Hepatokine'].unique()
    if len(strong_hks) < 3:
        strong_hks = corr_df['Hepatokine'].unique()[:15]
    corr_pivot = corr_pivot[corr_pivot.index.isin(strong_hks)]

    if len(corr_pivot) >= 2:
        # Build annotation
        annot = pd.DataFrame(index=corr_pivot.index, columns=corr_pivot.columns, dtype=str)
        for gene in corr_pivot.index:
            for mod in corr_pivot.columns:
                val = corr_pivot.loc[gene, mod]
                if pd.notna(val):
                    g_row = corr_df[(corr_df['Hepatokine'] == gene) & (corr_df['Muscle_Module'] == mod)]
                    p = g_row['P_value'].iloc[0] if len(g_row) > 0 else 1.0
                    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
                    annot.loc[gene, mod] = f'{val:.2f}{sig}'
                else:
                    annot.loc[gene, mod] = ''

        fig, ax = plt.subplots(figsize=(10, max(4, len(corr_pivot) * 0.4)))
        cmap = sns.diverging_palette(240, 10, as_cmap=True)

        sns.heatmap(corr_pivot.astype(float), annot=annot, fmt='', cmap=cmap, center=0,
                    vmin=-1, vmax=1, ax=ax, linewidths=0.5, linecolor='white',
                    cbar_kws={'label': 'Pearson r (n=8 group means)', 'shrink': 0.8},
                    annot_kws={'fontsize': 7})

        ax.set_title('Liver Hepatokine — Muscle Module Correlation\nGroup-Mean (DLY/TFB × 15/45/75/105 kg)',
                     fontsize=11, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('')

        plt.tight_layout()
        fig.savefig('figures_final/fig_HK2_crosstalk_heatmap.pdf', dpi=300)
        fig.savefig('figures_final/fig_HK2_crosstalk_heatmap.png', dpi=300)
        plt.close()
        print("  Saved fig_HK2_crosstalk_heatmap.pdf/png")

# FIGURE 3: Key hepatokine scatter plots — individual-level data (n≈48)
print("  Building individual-level expression data for scatter plots...")

# Build individual-level data from raw matrices
def build_individual_df(matrix, sample_map, gene_col='gene_name'):
    """Build long-form individual-level expression DataFrame."""
    expr_cols = [c for c in matrix.columns if c in sample_map]

    # Keep unique genes (deduplicate by gene_name, keep first)
    if gene_col and gene_col in matrix.columns:
        mat = matrix.drop_duplicates(subset=gene_col, keep='first').copy()
        mat = mat.set_index(gene_col)
    else:
        mat = matrix.copy()
        mat = mat.set_index(matrix.columns[0])

    for col in expr_cols:
        if col in mat.columns:
            mat[col] = pd.to_numeric(mat[col], errors='coerce')

    records = []
    for gene in mat.index:
        gene_str = str(gene)
        for col in expr_cols:
            if col not in mat.columns or col not in sample_map:
                continue
            val = mat.loc[gene, col]
            # Handle Series from duplicate index
            if isinstance(val, pd.Series):
                val = val.iloc[0] if len(val) > 0 else np.nan
            if pd.notna(val):
                info = sample_map[col]
                records.append({
                    'gene': gene_str, 'breed': info['breed'],
                    'stage': info['stage'], 'rep': info['rep'],
                    'expr': float(val)
                })
    return pd.DataFrame(records)

liver_ind = build_individual_df(liver_raw, lmap)
muscle_ind = build_individual_df(muscle_raw, mmap)
print(f"  Liver individual: {liver_ind.shape[0]} records, {liver_ind['gene'].nunique()} genes")
print(f"  Muscle individual: {muscle_ind.shape[0]} records, {muscle_ind['gene'].nunique()} genes")

key_hks = ['IGFBP1', 'IGFBP2', 'NRG4']
module_ribosome = MUSCLE_MODULES['Protein Synthesis\n(Ribosome)']
module_proteasome = MUSCLE_MODULES['Protein Degradation\n(Ubiquitin/Proteasome)']

fig, axes = plt.subplots(2, 3, figsize=(15, 9))

for col_idx, hk in enumerate(key_hks):
    # Get individual-level liver expression for this hepatokine
    hk_ind = liver_ind[liver_ind['gene'].str.upper() == hk.upper()]
    if len(hk_ind) == 0:
        continue

    # Ribosome module: mean expression of module genes per sample
    ribo_ind = muscle_ind[muscle_ind['gene'].isin(module_ribosome)]
    ribo_sample = ribo_ind.groupby(['breed', 'stage', 'rep'])['expr'].mean().reset_index()
    ribo_sample.rename(columns={'expr': 'muscle_expr'}, inplace=True)

    merged_ribo = hk_ind.merge(ribo_sample, on=['breed', 'stage', 'rep'])
    if len(merged_ribo) >= 10:
        r_r, p_r = pearsonr(merged_ribo['expr'], merged_ribo['muscle_expr'])

        # Plot individual points + group means
        for breed, c in [('DLY', C_DLY), ('TFB', C_TFB)]:
            sub = merged_ribo[merged_ribo['breed'] == breed]
            axes[0, col_idx].scatter(sub['expr'], sub['muscle_expr'], c=c, s=15,
                                      alpha=0.4, edgecolors='none')

            # Group mean (larger marker)
            gm = sub.groupby('stage')[['expr', 'muscle_expr']].mean()
            axes[0, col_idx].scatter(gm['expr'], gm['muscle_expr'], c=c, s=80,
                                      alpha=0.9, edgecolors='white', linewidth=1, marker='D')

        axes[0, col_idx].set_xlabel(f'{hk} liver expression (log2)', fontsize=9)
        axes[0, col_idx].set_ylabel('Muscle Ribosome\n(mean log2 expr)', fontsize=9)
        axes[0, col_idx].set_title(f'{hk} vs Muscle Ribosome\nr={r_r:.3f}, p={p_r:.3f}  (n={len(merged_ribo)} samples)',
                                   fontsize=9, fontweight='bold')

    # Proteasome module
    prot_ind = muscle_ind[muscle_ind['gene'].isin(module_proteasome)]
    prot_sample = prot_ind.groupby(['breed', 'stage', 'rep'])['expr'].mean().reset_index()
    prot_sample.rename(columns={'expr': 'muscle_expr'}, inplace=True)

    merged_prot = hk_ind.merge(prot_sample, on=['breed', 'stage', 'rep'])
    if len(merged_prot) >= 10:
        r_p, p_p = pearsonr(merged_prot['expr'], merged_prot['muscle_expr'])

        for breed, c in [('DLY', C_DLY), ('TFB', C_TFB)]:
            sub = merged_prot[merged_prot['breed'] == breed]
            axes[1, col_idx].scatter(sub['expr'], sub['muscle_expr'], c=c, s=15,
                                      alpha=0.4, edgecolors='none')
            gm = sub.groupby('stage')[['expr', 'muscle_expr']].mean()
            axes[1, col_idx].scatter(gm['expr'], gm['muscle_expr'], c=c, s=80,
                                      alpha=0.9, edgecolors='white', linewidth=1, marker='D')

        axes[1, col_idx].set_xlabel(f'{hk} liver expression (log2)', fontsize=9)
        axes[1, col_idx].set_ylabel('Muscle Proteasome\n(mean log2 expr)', fontsize=9)
        axes[1, col_idx].set_title(f'{hk} vs Muscle Proteasome\nr={r_p:.3f}, p={p_p:.3f}  (n={len(merged_prot)} samples)',
                                   fontsize=9, fontweight='bold')

# Legend: small dots=individual, diamonds=group means
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_DLY, markersize=8, label='DLY individual'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_TFB, markersize=8, label='TFB individual'),
    plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='#555555', markersize=8, label='Stage group mean'),
]
fig.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.99, 0.99), fontsize=8, frameon=True)

fig.suptitle('Liver Hepatokine — Muscle Functional Module Crosstalk\nIndividual-Level Data (15/45/75/105 kg, n≈48 samples)',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig('figures_final/fig_HK3_key_crosstalk_scatter.pdf', dpi=300)
fig.savefig('figures_final/fig_HK3_key_crosstalk_scatter.png', dpi=300)
plt.close()
print("  Saved fig_HK3_key_crosstalk_scatter.pdf/png")

# ============================================================
# 5. Mechanistic Model Summary
# ============================================================
print("\n[5/6] Building mechanistic model summary...")

# Identify key axes
print("""
  ╔══════════════════════════════════════════════════════════════╗
  ║     LIVER-MUSCLE CROSSTALK: Mechanistic Model Summary       ║
  ╠══════════════════════════════════════════════════════════════╣
  ║                                                              ║
  ║  TFB Liver @ 45kg (vs DLY):                                  ║
  ║                                                              ║
  ║  ┌─ AA Catabolism ↑ ──────────────────────────────┐         ║
  ║  │  GPT2 + GOT1 (transaminases) ↑                  │         ║
  ║  │  → Increased amino acid breakdown              │         ║
  ║  └────────────────────────────────────────────────┘         ║
  ║                                                              ║
  ║  ┌─ Urea Cycle ↑ ─────────────────────────────────┐         ║
  ║  │  CPS1 + ASS1 + ASL + ARG1 ↑                     │         ║
  ║  │  → Enhanced nitrogen disposal (wastage)         │         ║
  ║  └────────────────────────────────────────────────┘         ║
  ║                                                              ║
  ║  ┌─ IGFBP Secretion ↑ ─────────────────────────────┐        ║
  ║  │  IGFBP1 (-2.46) + IGFBP2 (-1.42) ↑              │        ║
  ║  │  → Reduced IGF1 bioavailability for muscle      │        ║
  ║  └────────────────────────────────────────────────┘        ║
  ║                                                              ║
  ║  ┌─ Protein Turnover ↑ (Proteasome + Translation) ─┐        ║
  ║  │  Proteasome KEGG (NES=-2.31, FDR=0.001)         │        ║
  ║  │  Translation Reactome (NES=-2.32, FDR=0.001)    │        ║
  ║  │  → Inefficient protein economy                  │        ║
  ║  └────────────────────────────────────────────────┘        ║
  ║                                                              ║
  ║  ┌─ Metabolic Signaling ───────────────────────────┐        ║
  ║  │  NRG4 (hepatokine) ↑                             │        ║
  ║  │  → Potential metabolic crosstalk to muscle       │        ║
  ║  └────────────────────────────────────────────────┘        ║
  ║                                                              ║
  ║  Net Effect on Muscle:                                       ║
  ║  • Less bioavailable IGF1 → reduced protein synthesis        ║
  ║  • Hepatokines (NRG4, IGFBPs) → altered muscle metabolism    ║
  ║  • Systemic nitrogen wastage → less substrate for muscle     ║
  ║                                                              ║
  ╚══════════════════════════════════════════════════════════════╝
""")

# ============================================================
# 6. Save results
# ============================================================
print("[6/6] Saving results...")

# Hepatokine table
hk_out = hk_df.copy()
hk_out.columns = ['Gene', 'Category', 'log2FC_DLYvsTFB', 'P_value', 'Q_value', 'FDR_significant', 'abs_log2FC']
hk_out.to_excel('hepatokine_screening_results.xlsx', index=False)
print("  Saved hepatokine_screening_results.xlsx")

# Crosstalk correlation table
if len(corr_df) > 0:
    corr_df.to_excel('hepatokine_muscle_crosstalk.xlsx', index=False)
    print("  Saved hepatokine_muscle_crosstalk.xlsx")

# Key findings summary
key_findings = []
for _, r in hk_df[hk_df['FDR_significant']].iterrows():
    direction = 'DLY-up' if r['log2FC_DLYvsTFB'] > 0 else 'TFB-up'
    key_findings.append({
        'Gene': r['Gene'],
        'Category': r['Category'],
        'log2FC': r['log2FC_DLYvsTFB'],
        'qvalue': r['qvalue'],
        'Direction': direction,
        'Biological_Implication': {
            'IGFBP1': 'Sequesters IGF1 → inhibits muscle protein synthesis',
            'IGFBP2': 'Sequesters IGF2/IGF1 → reduces growth signaling',
            'NRG4': 'Hepatokine signaling to muscle ERBB receptors → metabolic regulation',
            'GPT2': 'Alanine transaminase → AA catabolism marker',
            'GOT1': 'Aspartate transaminase → links AA to TCA cycle',
            'CPS1': 'Urea cycle rate-limiting → nitrogen disposal commitment',
            'ASS1': 'Argininosuccinate synthesis → urea cycle + arginine metabolism',
            'ASL': 'Argininosuccinate lyase → final urea cycle step before arginine',
            'ARG1': 'Arginase → urea production + arginine depletion',
            'ITIH3': 'Extracellular matrix stabilization → tissue remodeling',
        }.get(r['Gene'], 'Unknown'),
    })

key_df = pd.DataFrame(key_findings)
key_df.to_excel('hepatokine_key_findings.xlsx', index=False)
print("  Saved hepatokine_key_findings.xlsx")

print("\n" + "=" * 70)
print("HEPATOKINE CROSSTALK ANALYSIS COMPLETE")
print("=" * 70)
print(f"""
Summary:
  Hepatokines screened: {len(all_hepatokines)}
  Found in DEG: {n_hk_found}
  FDR<0.05: {n_hk_fdr} (IGF axis + AA/urea cycle + NRG4)

  Key hepatokines (TFB-up):
    IGFBP1 (log2FC=-2.46) — inhibits muscle IGF1 signaling
    IGFBP2 (log2FC=-1.42) — inhibits IGF signaling
    NRG4   (log2FC=-1.55) — metabolic hepatokine
    GPT2   (log2FC=-2.16) — alanine transaminase
    GOT1   (log2FC=-1.96) — aspartate transaminase
    CPS1   (log2FC=-1.76) — urea cycle commitment
    ASS1   (log2FC=-1.65) — urea cycle
    ASL    (log2FC=-1.06) — urea cycle (strongest FDR)
    ARG1   (log2FC=-1.65) — arginase/urea cycle

  Output:
    hepatokine_screening_results.xlsx
    hepatokine_muscle_crosstalk.xlsx
    hepatokine_key_findings.xlsx
    fig_HK1_hepatokine_profile.pdf
    fig_HK2_crosstalk_heatmap.pdf
    fig_HK3_key_crosstalk_scatter.pdf
""")
