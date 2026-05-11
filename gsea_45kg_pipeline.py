#!/usr/bin/env python3
"""
GSEA Preranked Analysis: Liver DLY vs TFB @ 45 kg
===================================================
Decision-window focused differential expression + pathway enrichment.
Tests the hypothesis: TFB liver at 45 kg shows aberrant activation of
AA catabolism and urea cycle pathways relative to DLY.

Output:
  - gsea_45kg_deg_results.xlsx
  - gsea_45kg_enrichment.xlsx
  - fig_GSEA1_volcano_45kg.pdf
  - fig_GSEA2_enrichment_bar.pdf
  - fig_GSEA3_AA_heatmap_45kg.pdf
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import gseapy as gp
import warnings
import os
warnings.filterwarnings('ignore')

from stats_utils import benjamini_hochberg, safe_ttest

# ============================================================
# Style
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.titlesize': 11, 'axes.labelsize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
})

C_DLY  = '#2166AC'
C_TFB  = '#B2182B'
C_NS   = '#999999'

os.makedirs('figures_final', exist_ok=True)

# ============================================================
# 1. Load & prepare 45 kg liver data
# ============================================================
print("=" * 70)
print("GSEA PRERANKED: Liver DLY vs TFB @ 45 kg")
print("=" * 70)

print("\n[1/4] Loading liver expression data...")
liver = pd.read_csv('gene_expression/liver_gene_matrix.xls', sep='\t')

sample_cols_all = [c for c in liver.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
dly_45_cols = sorted([c for c in sample_cols_all if c.startswith('L_45_1_')])
tfb_45_cols = sorted([c for c in sample_cols_all if c.startswith('L_45_2_')])

print(f"  DLY 45kg: {len(dly_45_cols)} samples")
print(f"  TFB 45kg: {len(tfb_45_cols)} samples")

dly_mat = liver[dly_45_cols].values.astype(float)
tfb_mat = liver[tfb_45_cols].values.astype(float)
gene_ids = liver['seq_id'].values
gene_names = liver['gene_name'].fillna(liver['seq_id']).values

# ============================================================
# 2. DEG per gene
# ============================================================
print("\n[2/4] Running DEG: Welch's t-test per gene...")

n_genes = len(gene_ids)
log2fc = np.zeros(n_genes)
pvalues = np.ones(n_genes)
dly_means = np.zeros(n_genes)
tfb_means = np.zeros(n_genes)

for g in range(n_genes):
    dly_vals = dly_mat[g, :]
    tfb_vals = tfb_mat[g, :]
    dly_vals = dly_vals[dly_vals > 0.01]
    tfb_vals = tfb_vals[tfb_vals > 0.01]
    if len(dly_vals) < 2 or len(tfb_vals) < 2:
        continue
    dly_means[g] = np.mean(dly_vals)
    tfb_means[g] = np.mean(tfb_vals)
    log2fc[g] = dly_means[g] - tfb_means[g]
    t_stat, p_val = safe_ttest(dly_vals, tfb_vals)
    pvalues[g] = p_val if not np.isnan(p_val) else 1.0

rejected, qvalues = benjamini_hochberg(pvalues)

deg_df = pd.DataFrame({
    'gene_id': gene_ids, 'gene_name': gene_names,
    'DLY_mean_log2': dly_means, 'TFB_mean_log2': tfb_means,
    'log2FC_DLYvsTFB': log2fc, 'abs_log2FC': np.abs(log2fc),
    'pvalue': pvalues, 'qvalue_FDR': qvalues,
    'FDR_significant': rejected,
})

expr_mask = (deg_df['DLY_mean_log2'] > 0.1) | (deg_df['TFB_mean_log2'] > 0.1)
deg_df['expressed'] = expr_mask
deg_df_filt = deg_df[expr_mask].copy()

# Use a relaxed threshold for visualization: nominal P<0.05 + |log2FC|>0.5
deg_df_filt['nominal_sig'] = (deg_df_filt['pvalue'] < 0.05) & (deg_df_filt['abs_log2FC'] > 0.5)
n_nom = deg_df_filt['nominal_sig'].sum()
n_fdr = deg_df_filt['FDR_significant'].sum()

print(f"  Expressed genes: {expr_mask.sum():,} / {n_genes:,}")
print(f"  Nominal P<0.05 + |log2FC|>0.5: {n_nom:,}")
print(f"  FDR < 0.05 (stringent): {n_fdr:,}")
print(f"  (FDR=0 expected with n=6/group + 16K tests — GSEA uses rank, not threshold)")

# ============================================================
# 3. Prerank for GSEA
# ============================================================
print("\n[3/4] Building preranked gene list for GSEA...")

rank_df = deg_df_filt[deg_df_filt['expressed']].copy()
# Ranking: sign(log2FC) * (|log2FC| - log10(p))
# GSEA cares about rank order, not absolute significance
rank_df['rank_metric'] = np.sign(rank_df['log2FC_DLYvsTFB']) * (
    rank_df['abs_log2FC'] - np.log10(rank_df['pvalue'].clip(lower=1e-300))
)
rank_df = rank_df.sort_values('pvalue').drop_duplicates(subset='gene_name', keep='first')
rank_df = rank_df.sort_values('rank_metric', ascending=False)

rnk = rank_df[['gene_name', 'rank_metric']].dropna()
rnk = rnk[rnk['gene_name'].str.strip() != '']
print(f"  Ranked genes: {len(rnk):,}")

# ============================================================
# 4. GSEA Preranked
# ============================================================
print("\n[4/4] Running GSEA preranked...")

LIBRARIES = {
    'Hallmark_2020': 'MSigDB_Hallmark_2020',
    'KEGG_2021': 'KEGG_2021_Human',
    'Reactome_2024': 'Reactome_Pathways_2024',
    'WikiPathways_2024': 'WikiPathways_2024_Human',
    'GO_BP_2025': 'GO_Biological_Process_2025',
}

gsea_results = {}
for lib_name, lib_id in LIBRARIES.items():
    try:
        print(f"  {lib_name}...", end=' ')
        gs_res = gp.prerank(
            rnk=rnk, gene_sets=lib_id, organism='human',
            outdir=None, min_size=10, max_size=500,
            permutation_num=1000, seed=42, threads=2, verbose=False,
        )
        res_df = gs_res.res2d
        if res_df is not None and len(res_df) > 0:
            res_df['Library'] = lib_name
            # Ensure numeric NES
            res_df['NES'] = pd.to_numeric(res_df['NES'], errors='coerce')
            res_df['FDR q-val'] = pd.to_numeric(res_df['FDR q-val'], errors='coerce')
            gsea_results[lib_name] = res_df
            n_sig = (res_df['FDR q-val'] < 0.05).sum()
            print(f"{len(res_df)} sets, {n_sig} FDR<0.05")
        else:
            print("empty")
    except Exception as e:
        print(f"failed: {e}")

# ============================================================
# 5. Save results
# ============================================================
deg_df_filt.to_excel('gsea_45kg_deg_results.xlsx', index=False)
print("  Saved gsea_45kg_deg_results.xlsx")

if gsea_results:
    with pd.ExcelWriter('gsea_45kg_enrichment.xlsx') as writer:
        for lib_name, res_df in gsea_results.items():
            cols = ['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'Library']
            avail = [c for c in cols if c in res_df.columns]
            res_df[avail].sort_values('FDR q-val').to_excel(writer, sheet_name=lib_name[:31], index=False)
    print("  Saved gsea_45kg_enrichment.xlsx")

# ============================================================
# FIGURE 1: Volcano Plot
# ============================================================
print("\nGenerating Figure 1: Volcano plot...")

fig, ax = plt.subplots(figsize=(6, 5.5))

ns = ~deg_df_filt['nominal_sig']
sig_up = deg_df_filt['nominal_sig'] & (deg_df_filt['log2FC_DLYvsTFB'] > 0.5)
sig_dn = deg_df_filt['nominal_sig'] & (deg_df_filt['log2FC_DLYvsTFB'] < -0.5)

ax.scatter(deg_df_filt.loc[ns, 'log2FC_DLYvsTFB'],
           -np.log10(deg_df_filt.loc[ns, 'pvalue'].clip(lower=1e-300)),
           c=C_NS, s=1.5, alpha=0.3, rasterized=True)

ax.scatter(deg_df_filt.loc[sig_up, 'log2FC_DLYvsTFB'],
           -np.log10(deg_df_filt.loc[sig_up, 'pvalue'].clip(lower=1e-300)),
           c=C_DLY, s=6, alpha=0.7, rasterized=True, label=f'DLY-up ({sig_up.sum():,})')

ax.scatter(deg_df_filt.loc[sig_dn, 'log2FC_DLYvsTFB'],
           -np.log10(deg_df_filt.loc[sig_dn, 'pvalue'].clip(lower=1e-300)),
           c=C_TFB, s=6, alpha=0.7, rasterized=True, label=f'TFB-up ({sig_dn.sum():,})')

# Label top 15 genes
top_genes = deg_df_filt.nlargest(15, 'abs_log2FC')
for _, g in top_genes.iterrows():
    ax.annotate(g['gene_name'],
                (g['log2FC_DLYvsTFB'], -np.log10(max(g['pvalue'], 1e-300))),
                fontsize=5.5, fontweight='bold', ha='center', va='bottom',
                xytext=(0, 3), textcoords='offset points')

# Highlight AA catabolism genes
AA_GENES = ['BCAT2', 'BCKDHA', 'BCKDHB', 'DBT', 'DLD', 'CPS1', 'OTC',
            'ASS1', 'ASL', 'ARG1', 'GOT1', 'GOT2', 'GPT', 'AASS', 'HGD',
            'ACADSB', 'GLUD1', 'SDS', 'HAL', 'PAH', 'STAT3']
aa_in_deg = deg_df_filt[deg_df_filt['gene_name'].isin(AA_GENES)]
for _, aa in aa_in_deg.iterrows():
    ax.annotate(aa['gene_name'],
                (aa['log2FC_DLYvsTFB'], -np.log10(max(aa['pvalue'], 1e-300))),
                fontsize=6, fontweight='bold', color='#D73027',
                ha='center', va='bottom', xytext=(0, 5), textcoords='offset points')

ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.5, alpha=0.5)
ax.axvline(0.5, color='gray', ls='--', lw=0.5, alpha=0.5)
ax.axvline(-0.5, color='gray', ls='--', lw=0.5, alpha=0.5)

ax.set_xlabel('log2(DLY / TFB) expression', fontsize=10)
ax.set_ylabel('-log10(P value)', fontsize=10)
ax.set_title('Liver DLY vs TFB at 45 kg\nDecision Window Transcriptome', fontsize=12, fontweight='bold')
ax.legend(loc='upper left', frameon=True, fontsize=6.5, markerscale=2)

# Note about FDR
ax.text(0.98, 0.98, f'Nominal P<0.05 + |log2FC|>0.5: {n_nom:,} genes\n'
        f'FDR<0.05 (stringent): {n_fdr:,} genes\n'
        f'GSEA uses rank order, not threshold',
        transform=ax.transAxes, va='top', ha='right', fontsize=6, color='#555555')

plt.tight_layout()
fig.savefig('figures_final/fig_GSEA1_volcano_45kg.pdf', dpi=300)
fig.savefig('figures_final/fig_GSEA1_volcano_45kg.png', dpi=300)
plt.close()
print("  Saved fig_GSEA1_volcano_45kg.pdf/png")

# ============================================================
# FIGURE 2: GSEA Enrichment Bar Plot
# ============================================================
print("Generating Figure 2: GSEA enrichment summary...")

all_paths = []
for lib, res in gsea_results.items():
    if res is None or len(res) == 0:
        continue
    sig = res[res['FDR q-val'] < 0.25].copy()  # relaxed for exploration
    sig['Library'] = lib
    sig['abs_NES'] = sig['NES'].abs()
    all_paths.append(sig)

if all_paths:
    combined = pd.concat(all_paths, ignore_index=True)
    combined['direction'] = combined['NES'].apply(lambda x: 'TFB-enriched' if x < 0 else 'DLY-enriched')

    # Top 30 by abs(NES)
    top_paths = combined.nlargest(30, 'abs_NES')

    fig, ax = plt.subplots(figsize=(10, 8))

    for i, (_, row) in enumerate(top_paths.iterrows()):
        color = C_TFB if row['direction'] == 'TFB-enriched' else C_DLY
        nes_val = float(row['NES'])
        ax.barh(i, abs(nes_val), color=color, alpha=0.85, height=0.7,
                edgecolor='white', linewidth=0.3)
        fdr_val = float(row['FDR q-val'])
        fdr_str = f"FDR={fdr_val:.2e}" if fdr_val < 0.01 else f"FDR={fdr_val:.3f}"
        ax.text(abs(nes_val) + 0.05, i, fdr_str, va='center', fontsize=5, color='#555555')
        label = f"{row['Term'][:60]} [{row['Library']}]"
        ax.text(0.05, i, label, va='center', fontsize=5.5, color='#222222')

    ax.set_yticks([])
    ax.set_xlabel('|Normalized Enrichment Score|', fontsize=10)
    ax.set_title('GSEA Preranked: Top Enriched Pathways\nDLY vs TFB Liver @ 45 kg',
                 fontsize=12, fontweight='bold')

    legend_elements = [
        mpatches.Patch(facecolor=C_TFB, alpha=0.85, label='TFB-enriched'),
        mpatches.Patch(facecolor=C_DLY, alpha=0.85, label='DLY-enriched'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=7, frameon=True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    fig.savefig('figures_final/fig_GSEA2_enrichment_bar.pdf', dpi=300)
    fig.savefig('figures_final/fig_GSEA2_enrichment_bar.png', dpi=300)
    plt.close()
    print("  Saved fig_GSEA2_enrichment_bar.pdf/png")

# ============================================================
# FIGURE 3: AA Metabolism Gene Set Heatmap (45kg focus)
# ============================================================
print("Generating Figure 3: AA metabolism enzyme expression heatmap...")

AA_PATHWAY_GENES = {
    'BCAA Catabolism': ['BCAT2', 'BCKDHA', 'BCKDHB', 'DBT', 'DLD', 'ACADSB', 'HIBCH'],
    'Urea Cycle': ['CPS1', 'OTC', 'ASS1', 'ASL', 'ARG1', 'ARG2', 'NAGS'],
    'Transamination': ['GOT1', 'GOT2', 'GPT', 'GPT2', 'GLUD1', 'GLUL'],
    'Specific AA\nCatabolism': ['AASS', 'HGD', 'SDS', 'HAL', 'PAH', 'TAT', 'FAH'],
    'Proline/Polyamine': ['ALDH18A1', 'PYCR1', 'PYCR2', 'PRODH', 'ODC1', 'SRM', 'SMS'],
    'Serine/Glycine': ['SHMT1', 'SHMT2', 'PHGDH', 'PSAT1', 'PSPH', 'GCAT'],
    'Sulfur (Cys/Met)': ['CBS', 'CTH', 'CDO1', 'MAT1A', 'MAT2A', 'AHCY'],
    'IGF/mTOR': ['IGF1', 'IGFALS', 'IGFBP3', 'MTOR', 'AKT1', 'RPS6KB1'],
}

heatmap_data = []
for category, genes in AA_PATHWAY_GENES.items():
    for gene in genes:
        match = deg_df_filt[deg_df_filt['gene_name'].str.upper() == gene.upper()]
        if len(match) > 0:
            row = match.iloc[0]
            heatmap_data.append({
                'Gene': gene, 'Category': category,
                'log2FC': row['log2FC_DLYvsTFB'],
                'DLY_expr': row['DLY_mean_log2'],
                'TFB_expr': row['TFB_mean_log2'],
                'pvalue': row['pvalue'],
            })

hm_df = pd.DataFrame(heatmap_data)
if len(hm_df) > 0:
    hm_pivot = hm_df.pivot_table(values='log2FC', index='Gene', columns='Category', aggfunc='first')
    for cat in AA_PATHWAY_GENES:
        if cat not in hm_pivot.columns:
            hm_pivot[cat] = np.nan
    hm_pivot = hm_pivot[list(AA_PATHWAY_GENES.keys())]
    hm_pivot = hm_pivot.dropna(how='all')

    # Build annotation matrix (string type to avoid dtype conflict)
    hm_annot = pd.DataFrame(index=hm_pivot.index, columns=hm_pivot.columns, dtype=str)
    for gene in hm_pivot.index:
        g_row = hm_df[hm_df['Gene'] == gene]
        p = g_row['pvalue'].iloc[0] if len(g_row) > 0 else 1.0
        sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else ''
        for cat in hm_pivot.columns:
            val = hm_pivot.loc[gene, cat]
            if pd.notna(val):
                hm_annot.loc[gene, cat] = f'{float(val):.2f}{sig}' if sig else f'{float(val):.2f}'
            else:
                hm_annot.loc[gene, cat] = ''

    fig, ax = plt.subplots(figsize=(10.5, max(5, len(hm_pivot) * 0.38)))
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    sns.heatmap(hm_pivot, annot=hm_annot, fmt='', cmap=cmap, center=0,
                vmin=-3, vmax=3, ax=ax, linewidths=0.5, linecolor='white',
                cbar_kws={'label': 'log2(DLY/TFB) @ 45kg', 'shrink': 0.8},
                annot_kws={'fontsize': 7})

    ax.set_title('AA Metabolism & Urea Cycle Enzyme Expression\nDLY vs TFB Liver @ 45 kg Decision Window',
                 fontsize=12, fontweight='bold')
    ax.set_ylabel('')
    ax.set_xlabel('')

    plt.tight_layout()
    fig.savefig('figures_final/fig_GSEA3_AA_heatmap_45kg.pdf', dpi=300)
    fig.savefig('figures_final/fig_GSEA3_AA_heatmap_45kg.png', dpi=300)
    plt.close()
    print("  Saved fig_GSEA3_AA_heatmap_45kg.pdf/png")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("GSEA PIPELINE COMPLETE")
print("=" * 70)
print(f"""
Output:
  gsea_45kg_deg_results.xlsx        — {len(deg_df_filt):,} expressed genes
  gsea_45kg_enrichment.xlsx         — GSEA results from {len(gsea_results)} libraries
  fig_GSEA1_volcano_45kg.pdf        — Volcano plot
  fig_GSEA2_enrichment_bar.pdf      — Top enriched pathways
  fig_GSEA3_AA_heatmap_45kg.pdf     — AA enzyme focused heatmap

Key stats:
  Nominal P<0.05 + |log2FC|>0.5:  {n_nom:,} genes
  FDR < 0.05 (stringent):          {n_fdr:,} genes
  (n=6/group limits per-gene power; GSEA uses rank-based testing)
""")

if gsea_results:
    print("\nTop 20 GSEA hits (FDR < 0.25):")
    all_sig = []
    for lib, res in gsea_results.items():
        if res is not None and len(res) > 0:
            s = res[res['FDR q-val'] < 0.25].copy()
            s['Library'] = lib
            all_sig.append(s)
    if all_sig:
        top_all = pd.concat(all_sig, ignore_index=True)
        top_all['abs_NES'] = top_all['NES'].abs()
        for _, r in top_all.nlargest(20, 'abs_NES').iterrows():
            sig = '***' if r['FDR q-val'] < 0.001 else ('**' if r['FDR q-val'] < 0.01 else ('*' if r['FDR q-val'] < 0.05 else ''))
            print(f"  [{r['Library']:20s}] NES={r['NES']:+6.2f}  FDR={r['FDR q-val']:.4f} {sig}  {str(r['Term'])[:70]}")
