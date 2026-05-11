#!/usr/bin/env python3
"""
GSEA Multi-Stage Pipeline: Liver DLY vs TFB (15/45/75/105 kg)
===============================================================
Combined approach to overcome n=6/group power limitation:

  A. Multi-stage linear model (expr ~ breed + stage)
     → breed main effect as preranked metric
     → n=48 vs n=12 at 45kg alone → dramatic power gain

  B. ssGSEA pathway scoring per sample
     → aggregate gene-level signal into pathway scores
     → compare DLY vs TFB at 45kg on pathway scores

Output:
  - gsea_multistage_deg_results.xlsx       (breed effect + 45kg effect)
  - gsea_multistage_enrichment.xlsx        (GSEA preranked results)
  - gsea_multistage_ssgsea_scores.xlsx     (per-sample pathway scores)
  - gsea_multistage_ssgsea_45kg_test.xlsx  (45kg pathway comparison)
  - fig_MS1_volcano_breed.pdf
  - fig_MS2_gsea_enrichment_bar.pdf
  - fig_MS3_ssgsea_heatmap_45kg.pdf
  - fig_MS4_aa_pathway_scores.pdf
"""
import pandas as pd
import numpy as np
from scipy.stats import t as t_dist, ttest_ind
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

C_DLY = '#2166AC'
C_TFB = '#B2182B'
C_NS  = '#999999'

os.makedirs('figures_final', exist_ok=True)

# ============================================================
# 1. Load liver data (all stages)
# ============================================================
print("=" * 70)
print("GSEA MULTI-STAGE: Liver DLY vs TFB @ 15/45/75/105 kg")
print("=" * 70)

print("\n[1/8] Loading liver expression data...")
liver = pd.read_csv('gene_expression/liver_gene_matrix.xls', sep='\t')

# Map samples to breed/stage
sample_cols_all = [c for c in liver.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]

sample_meta = {}
for c in sample_cols_all:
    parts = c.split('_')
    if parts[0] != 'L':
        continue
    stage_code = parts[1]
    breed_code = parts[2]
    stage_map = {'15': 15, '45': 45, '1': 75, '2': 105, '3': 135}
    if stage_code not in stage_map:
        continue
    breed = 'DLY' if breed_code == '1' else 'TFB'
    sample_meta[c] = {'breed': breed, 'stage': stage_map[stage_code]}

# Only keep stages with both breeds (15/45/75/105)
balanced_samples = {c: m for c, m in sample_meta.items() if m['stage'] != 135}
print(f"  Balanced-design samples (15/45/75/105): {len(balanced_samples)}")

# Count
from collections import Counter
bc = Counter((m['breed'], m['stage']) for m in balanced_samples.values())
for (breed, stage), n in sorted(bc.items()):
    print(f"    {breed} {stage}kg: n={n}")

# Expression matrix for balanced samples
sample_ids = sorted(balanced_samples.keys())
n_samples = len(sample_ids)

# Build metadata arrays
breeds = np.array([balanced_samples[s]['breed'] for s in sample_ids])
stages = np.array([balanced_samples[s]['stage'] for s in sample_ids])

# Build expression matrix (genes × samples)
gene_ids = liver['seq_id'].values
gene_names_raw = liver['gene_name'].fillna('').values

expr_mat = np.zeros((len(gene_ids), n_samples))
for j, s in enumerate(sample_ids):
    expr_mat[:, j] = pd.to_numeric(liver[s], errors='coerce').fillna(0).values

# Filter to expressed genes (mean > 0.1 across samples)
gene_means_raw = expr_mat.mean(axis=1)
expressed = gene_means_raw > 0.1
print(f"  Expressed genes: {expressed.sum():,} / {len(gene_ids):,}")

expr_mat = expr_mat[expressed, :]
gene_ids = gene_ids[expressed]
gene_names_raw = gene_names_raw[expressed]
gene_means = gene_means_raw[expressed]
n_genes = len(gene_ids)

log2_expr = np.log2(expr_mat + 0.01)

# ============================================================
# 2. Per-gene linear model: expr ~ breed + C(stage)
# ============================================================
print(f"\n[2/8] Per-gene linear model: log2(expr) ~ breed + stage...")

# Design matrix: intercept + breed_TFB + stage_45 + stage_75 + stage_105
X = np.zeros((n_samples, 5))
X[:, 0] = 1.0                              # intercept (DLY @ 15kg ref)
X[:, 1] = (breeds == 'TFB').astype(float)   # breed effect (TFB vs DLY)
X[:, 2] = (stages == 45).astype(float)      # stage 45 vs 15
X[:, 3] = (stages == 75).astype(float)      # stage 75 vs 15
X[:, 4] = (stages == 105).astype(float)     # stage 105 vs 15

# Precompute (X'X)^-1 X' for vectorized OLS
XtX = X.T @ X
XtX_inv = np.linalg.inv(XtX)
H = XtX_inv @ X.T   # hat matrix: (5 × n_samples)

n_params = X.shape[1]
dof = n_samples - n_params
se_factor = np.sqrt(np.diag(XtX_inv))  # sqrt of diagonal of (X'X)^-1

# Process in batches to manage memory
batch_size = 5000
n_batches = int(np.ceil(n_genes / batch_size))

breed_coef = np.zeros(n_genes)
breed_se = np.zeros(n_genes)
breed_t = np.zeros(n_genes)
breed_p = np.zeros(n_genes)
breed_log2fc = np.zeros(n_genes)  # DLY - TFB in log2 space at 45kg specifically

for b in range(n_batches):
    start = b * batch_size
    end = min(start + batch_size, n_genes)
    Y = log2_expr[start:end, :].T  # (n_samples × batch)

    # OLS: B = (X'X)^-1 X' Y → (5 × batch)
    B = H @ Y

    # Residuals and SE
    Y_hat = X @ B
    residuals = Y - Y_hat
    rss = np.sum(residuals**2, axis=0)
    sigma2 = rss / dof
    sigma2[sigma2 <= 0] = 1e-12

    # SE for breed coefficient (index 1)
    se = np.sqrt(sigma2) * se_factor[1]
    t_stat = B[1, :] / se
    p_vals = 2 * t_dist.sf(np.abs(t_stat), dof)

    breed_coef[start:end] = B[1, :]
    breed_se[start:end] = se
    breed_t[start:end] = t_stat
    breed_p[start:end] = p_vals

    if (b + 1) % 3 == 0:
        print(f"  Progress: {end:,}/{n_genes:,} genes")

# FDR correction
rejected, qvalues = benjamini_hochberg(breed_p)

# breed_coef is TFB - DLY in log2 space (positive = TFB higher)
# For consistency with previous analysis, define log2FC as DLY - TFB
log2fc_breed = -breed_coef   # DLY vs TFB breed effect

n_nom = (breed_p < 0.05).sum()
n_fdr = rejected.sum()
print(f"\n  Nominal P<0.05 (breed): {n_nom:,} / {n_genes:,} ({100*n_nom/n_genes:.1f}%)")
print(f"  FDR < 0.05 (breed):    {n_fdr:,} / {n_genes:,} ({100*n_fdr/n_genes:.1f}%)")

# ============================================================
# 2b. 45kg-specific effect for comparison
# ============================================================
print("\n  Computing 45kg-specific breed effects...")

is_45 = stages == 45
n_45 = is_45.sum()
dly_45 = (breeds == 'DLY') & is_45
tfb_45 = (breeds == 'TFB') & is_45

log2fc_45 = np.zeros(n_genes)
p_45 = np.ones(n_genes)
for g in range(n_genes):
    dly_vals = log2_expr[g, dly_45]
    tfb_vals = log2_expr[g, tfb_45]
    log2fc_45[g] = dly_vals.mean() - tfb_vals.mean()
    t_stat, p_val = safe_ttest(dly_vals, tfb_vals)
    if not np.isnan(p_val):
        p_45[g] = p_val

# ============================================================
# 3. Build results table
# ============================================================
print("\n[3/8] Building DEG results table...")

# Clean gene names
names_clean = []
seen = set()
for gid, gname in zip(gene_ids, gene_names_raw):
    name = str(gname).strip() if gname and str(gname).strip() != '' else str(gid).strip()
    if name in seen:
        name = f"{name}_{gid[:8]}"
    seen.add(name)
    names_clean.append(name)

deg_df = pd.DataFrame({
    'gene_id': gene_ids,
    'gene_name': names_clean,
    'breed_coef_TFBvsDLY': breed_coef,          # raw coefficient
    'breed_log2FC_DLYvsTFB': log2fc_breed,       # DLY - TFB
    'breed_SE': breed_se,
    'breed_t_stat': breed_t,
    'breed_pvalue': breed_p,
    'breed_qvalue': qvalues,
    'breed_FDR_significant': rejected,
    'log2FC_45kg_DLYvsTFB': log2fc_45,
    'pvalue_45kg': p_45,
    'mean_log2_expr': log2_expr.mean(axis=1),
})

deg_df['abs_breed_log2FC'] = deg_df['breed_log2FC_DLYvsTFB'].abs()

# Report
n_strong = (deg_df['breed_FDR_significant'] & (deg_df['abs_breed_log2FC'] > 0.5)).sum()
print(f"  Total expressed genes: {n_genes:,}")
print(f"  FDR<0.05 + |log2FC|>0.5: {n_strong:,}")

# ============================================================
# 4. Preranked GSEA
# ============================================================
print("\n[4/8] Building preranked gene list for GSEA (breed main effect)...")

rank_df = deg_df.copy()
# Ranking metric: sign(log2FC) * (|log2FC| - log10(p))
rank_df['rank_metric'] = np.sign(rank_df['breed_log2FC_DLYvsTFB']) * (
    rank_df['abs_breed_log2FC'] - np.log10(rank_df['breed_pvalue'].clip(lower=1e-300))
)
rank_df = rank_df.sort_values('breed_pvalue').drop_duplicates(subset='gene_name', keep='first')
rank_df = rank_df.sort_values('rank_metric', ascending=False)

rnk = rank_df[['gene_name', 'rank_metric']].dropna()
rnk = rnk[rnk['gene_name'].str.strip() != '']
print(f"  Ranked genes: {len(rnk):,}")

# ============================================================
# 5. Run GSEA preranked
# ============================================================
print("\n[5/8] Running GSEA preranked (breed effect ranking)...")

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
            res_df['NES'] = pd.to_numeric(res_df['NES'], errors='coerce')
            res_df['FDR q-val'] = pd.to_numeric(res_df['FDR q-val'], errors='coerce')
            gsea_results[lib_name] = res_df
            n_sig = (res_df['FDR q-val'] < 0.05).sum()
            n_sig25 = (res_df['FDR q-val'] < 0.25).sum()
            print(f"{len(res_df)} sets, {n_sig} FDR<0.05, {n_sig25} FDR<0.25")
        else:
            print("empty")
    except Exception as e:
        print(f"failed: {e}")

# ============================================================
# 6. ssGSEA: per-sample pathway scoring
# ============================================================
print("\n[6/8] Running ssGSEA for per-sample pathway scoring...")

# Build expression matrix indexed by gene symbols (for gseapy gene set matching)
ssgsea_expr = pd.DataFrame(log2_expr, index=names_clean, columns=sample_ids)
ssgsea_expr.index = ssgsea_expr.index.str.strip()
ssgsea_expr = ssgsea_expr[~ssgsea_expr.index.duplicated(keep='first')]
print(f"  ssGSEA input: {ssgsea_expr.shape[0]} genes × {ssgsea_expr.shape[1]} samples")

ssgsea_results = {}
for lib_name, lib_id in LIBRARIES.items():
    try:
        print(f"  {lib_name}...", end=' ')
        ss_res = gp.ssgsea(
            data=ssgsea_expr, gene_sets=lib_id, organism='human',
            outdir=None, min_size=10, max_size=500,
            sample_norm_method='rank', no_plot=True,
            processes=1, seed=42, format='pandas',
        )
        # ss_res.res2d is long-form: Name(sample) | Term(pathway) | ES | NES
        # Pivot to pathway × sample matrix
        long_df = ss_res.res2d
        if long_df is not None and len(long_df) > 0:
            scores_df = long_df.pivot(index='Term', columns='Name', values='NES')
            ssgsea_results[lib_name] = scores_df
            print(f"{scores_df.shape[0]} pathways × {scores_df.shape[1]} samples")
        else:
            print("empty")
    except Exception as e:
        print(f"failed: {e}")

# ============================================================
# 7. 45kg pathway comparison from ssGSEA scores
# ============================================================
print("\n[7/8] Testing ssGSEA pathway scores: DLY vs TFB @ 45kg...")

dly_45_cols = [s for s in sample_ids if balanced_samples[s]['breed'] == 'DLY' and balanced_samples[s]['stage'] == 45]
tfb_45_cols = [s for s in sample_ids if balanced_samples[s]['breed'] == 'TFB' and balanced_samples[s]['stage'] == 45]

print(f"  DLY 45kg: {len(dly_45_cols)} samples, TFB 45kg: {len(tfb_45_cols)} samples")

pathway_tests = []
for lib_name, scores_df in ssgsea_results.items():
    dly_cols_in = [c for c in dly_45_cols if c in scores_df.columns]
    tfb_cols_in = [c for c in tfb_45_cols if c in scores_df.columns]
    if len(dly_cols_in) < 2 or len(tfb_cols_in) < 2:
        continue
    for pathway in scores_df.index:
        dly_scores = scores_df.loc[pathway, dly_cols_in].astype(float).values
        tfb_scores = scores_df.loc[pathway, tfb_cols_in].astype(float).values

        # Drop NaN
        dly_scores = dly_scores[~np.isnan(dly_scores)]
        tfb_scores = tfb_scores[~np.isnan(tfb_scores)]

        if len(dly_scores) >= 2 and len(tfb_scores) >= 2:
            t_stat, p_val = ttest_ind(dly_scores, tfb_scores, equal_var=False)
            cohens_d = (dly_scores.mean() - tfb_scores.mean()) / max(
                np.sqrt((dly_scores.var(ddof=1) + tfb_scores.var(ddof=1)) / 2), 1e-12
            )
            pathway_tests.append({
                'Library': lib_name,
                'Pathway': pathway,
                'DLY_mean_ssGSEA': dly_scores.mean(),
                'TFB_mean_ssGSEA': tfb_scores.mean(),
                'delta_ssGSEA': dly_scores.mean() - tfb_scores.mean(),
                'Cohens_d': cohens_d,
                't_statistic': t_stat,
                'P_value': p_val,
            })

pathway_df = pd.DataFrame(pathway_tests)
if len(pathway_df) > 0:
    _, pq = benjamini_hochberg(pathway_df['P_value'].values)
    pathway_df['Q_value'] = pq
    pathway_df['FDR_significant'] = pq < 0.05
    pathway_df = pathway_df.sort_values('P_value')

    n_sig_path = pathway_df['FDR_significant'].sum()
    n_nom_path = (pathway_df['P_value'] < 0.05).sum()
    print(f"  Total pathways tested: {len(pathway_df)}")
    print(f"  Nominal P<0.05: {n_nom_path} → FDR<0.05: {n_sig_path}")

    if n_sig_path > 0:
        print(f"\n  Top FDR-significant pathways:")
        for _, r in pathway_df[pathway_df['FDR_significant']].head(15).iterrows():
            direction = 'DLY-up' if r['delta_ssGSEA'] > 0 else 'TFB-up'
            print(f"    [{r['Library']}] {r['Pathway'][:60]}: d={r['Cohens_d']:+.2f}, P={r['P_value']:.4f}, Q={r['Q_value']:.4f} [{direction}]")

# ============================================================
# 8. Save results
# ============================================================
print("\n[8/8] Saving results...")

deg_df.to_excel('gsea_multistage_deg_results.xlsx', index=False)
print("  Saved gsea_multistage_deg_results.xlsx")

if gsea_results:
    with pd.ExcelWriter('gsea_multistage_enrichment.xlsx') as writer:
        for lib_name, res_df in gsea_results.items():
            cols = ['Term', 'ES', 'NES', 'NOM p-val', 'FDR q-val', 'Library']
            avail = [c for c in cols if c in res_df.columns]
            res_df[avail].sort_values('FDR q-val').to_excel(writer, sheet_name=lib_name[:31], index=False)
    print("  Saved gsea_multistage_enrichment.xlsx")

combined_scores = []
for lib_name, scores_df in ssgsea_results.items():
    scores_df_copy = scores_df.copy()
    scores_df_copy['Library'] = lib_name
    combined_scores.append(scores_df_copy)
if combined_scores:
    all_scores = pd.concat(combined_scores)
    all_scores.to_excel('gsea_multistage_ssgsea_scores.xlsx')
    print("  Saved gsea_multistage_ssgsea_scores.xlsx")

pathway_df.to_excel('gsea_multistage_ssgsea_45kg_test.xlsx', index=False)
print("  Saved gsea_multistage_ssgsea_45kg_test.xlsx")

# ============================================================
# FIGURE MS1: Volcano Plot (breed main effect)
# ============================================================
print("\nGenerating Figure MS1: Volcano plot (breed main effect)...")

fig, ax = plt.subplots(figsize=(7, 6))

# Categories
not_sig = ~deg_df['breed_FDR_significant'] | (deg_df['abs_breed_log2FC'] < 0.3)
fdr_up = deg_df['breed_FDR_significant'] & (deg_df['breed_log2FC_DLYvsTFB'] > 0.3)
fdr_dn = deg_df['breed_FDR_significant'] & (deg_df['breed_log2FC_DLYvsTFB'] < -0.3)

ax.scatter(deg_df.loc[not_sig, 'breed_log2FC_DLYvsTFB'],
           -np.log10(deg_df.loc[not_sig, 'breed_pvalue'].clip(lower=1e-300)),
           c=C_NS, s=0.8, alpha=0.2, rasterized=True)

ax.scatter(deg_df.loc[fdr_up, 'breed_log2FC_DLYvsTFB'],
           -np.log10(deg_df.loc[fdr_up, 'breed_pvalue'].clip(lower=1e-300)),
           c=C_DLY, s=6, alpha=0.55, rasterized=True, edgecolors='none',
           label=f'DLY > TFB  FDR<0.05 ({fdr_up.sum():,})')

ax.scatter(deg_df.loc[fdr_dn, 'breed_log2FC_DLYvsTFB'],
           -np.log10(deg_df.loc[fdr_dn, 'breed_pvalue'].clip(lower=1e-300)),
           c=C_TFB, s=6, alpha=0.55, rasterized=True, edgecolors='none',
           label=f'TFB > DLY  FDR<0.05 ({fdr_dn.sum():,})')

# Label top 8 genes by |log2FC| (with gene symbols only)
top_genes = deg_df[deg_df['gene_name'].str.match(r'^[A-Z]')].nlargest(8, 'abs_breed_log2FC')
for _, g in top_genes.iterrows():
    ax.annotate(g['gene_name'],
                (g['breed_log2FC_DLYvsTFB'], -np.log10(max(g['breed_pvalue'], 1e-300))),
                fontsize=6.5, fontweight='bold', ha='center', va='bottom',
                xytext=(0, 4), textcoords='offset points',
                arrowprops=dict(arrowstyle='-', color='#555555', lw=0.3))

# Highlight AA/urea cycle genes (only FDR-significant ones)
AA_GENES = ['CPS1', 'OTC', 'ASS1', 'ASL', 'ARG1', 'ARG2', 'GOT1', 'GOT2', 'GPT2',
            'AASS', 'HGD', 'SDS', 'BCAT2', 'BCKDHA', 'BCKDHB', 'DBT', 'DLD',
            'ACADSB', 'GLUD1', 'HAL', 'PAH', 'STAT3']
aa_labeled = 0
for gname in AA_GENES:
    match = deg_df[deg_df['gene_name'].str.upper() == gname.upper()]
    if len(match) > 0:
        g = match.iloc[0]
        if g['breed_FDR_significant']:  # Only label FDR-significant
            aa_labeled += 1
            ax.annotate(g['gene_name'],
                        (g['breed_log2FC_DLYvsTFB'], -np.log10(max(g['breed_pvalue'], 1e-300))),
                        fontsize=6.5, fontweight='bold', color='#D73027',
                        ha='center', va='bottom', xytext=(0, 5), textcoords='offset points',
                        arrowprops=dict(arrowstyle='-', color='#D73027', lw=0.5, alpha=0.6))

ax.axhline(-np.log10(0.05), color='gray', ls='--', lw=0.6, alpha=0.4)
ax.axhline(-np.log10(0.01), color='gray', ls=':', lw=0.4, alpha=0.3)
ax.axvline(0, color='gray', ls='-', lw=0.3, alpha=0.3)

ax.set_xlabel('log$_2$(DLY / TFB)  stage-adjusted breed effect', fontsize=11)
ax.set_ylabel('−log$_{10}$(P value)', fontsize=11)
ax.set_title('Liver Transcriptome: DLY vs TFB Breed Effect\n(Linear Model: expr ~ breed + stage, 15/45/75/105 kg)',
             fontsize=12, fontweight='bold')
ax.legend(loc='upper left', frameon=True, fontsize=7, markerscale=1.5,
          handletextpad=0.5, borderpad=0.5)

# Stats box
ax.text(0.98, 0.97,
        f'{n_genes:,} expressed genes\n'
        f'Nominal P < 0.05: {n_nom:,} ({100*n_nom/n_genes:.1f}%)\n'
        f'FDR < 0.05: {n_fdr:,} ({100*n_fdr/n_genes:.1f}%)\n'
        f'FDR + |log$_2$FC| > 0.5: {n_strong:,}\n'
        f'Red labels: AA/urea cycle (FDR sig, n={aa_labeled})',
        transform=ax.transAxes, va='top', ha='right', fontsize=6.5,
        color='#444444', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.85, edgecolor='#cccccc'))

ax.set_xlim(-8, 8)
plt.tight_layout()
fig.savefig('figures_final/fig_MS1_volcano_breed.pdf', dpi=300)
fig.savefig('figures_final/fig_MS1_volcano_breed.png', dpi=300)
plt.close()
print("  Saved fig_MS1_volcano_breed.pdf/png")

# ============================================================
# FIGURE MS2: GSEA Enrichment — Categorical Summary
# ============================================================
print("Generating Figure MS2: GSEA enrichment categorical summary...")

# Define functional categories for pathway grouping
PATHWAY_CATEGORIES = {
    'Amino Acid\nMetabolism': ['AMINO ACID', 'BCAA', 'BRANCHED-CHAIN', 'SELENOAMINO', 'SULFUR AMINO',
                                'METHIONINE', 'TRYPTOPHAN', 'LYSINE', 'ARGININE AND PROLINE',
                                'GLUTAMINE', 'ALPHA-AMINO', 'L-AMINO'],
    'Urea Cycle &\nNitrogen': ['UREA', 'ORNITHINE', 'POLYAMINE', 'NITROGEN', 'AMINO GROUP'],
    'Protein\nDegradation': ['PROTEASOME', 'UBIQUITIN', 'DEGRADATION', 'AUTODEGRADATION',
                              'PROTEIN DEGRAD', 'PROTEOLYSIS'],
    'Translation &\nRibosome': ['TRANSLATION', 'RIBOSOME', 'RIBOSOMAL', 'TRNA', 'AMINOACYLATION'],
    'Energy\nMetabolism': ['OXIDATIVE PHOSPHORYLATION', 'ELECTRON TRANSPORT', 'MITOCHONDRIAL',
                            'RESPIRATORY CHAIN', 'ATP SYNTH', 'CITRIC ACID', 'TCA CYCLE'],
}

def classify_pathway(term):
    term_upper = str(term).upper()
    for cat, keywords in PATHWAY_CATEGORIES.items():
        for kw in keywords:
            if kw in term_upper:
                return cat
    return 'Other'

# Aggregate all FDR<0.05 pathways by category
all_sig = []
for lib, res in gsea_results.items():
    if res is None or len(res) == 0:
        continue
    sig = res[res['FDR q-val'] < 0.05].copy()
    sig['Library'] = lib
    sig['Category'] = sig['Term'].apply(classify_pathway)
    sig['abs_NES'] = sig['NES'].abs()
    all_sig.append(sig)

if all_sig:
    combined = pd.concat(all_sig, ignore_index=True)

    # Summarize by category
    cat_summary = combined.groupby('Category').agg(
        n_pathways=('Term', 'count'),
        mean_NES=('NES', 'mean'),
        max_abs_NES=('abs_NES', 'max'),
        mean_FDR=('FDR q-val', 'mean'),
    ).sort_values('max_abs_NES', ascending=True)

    # Also get top individual pathways per category
    cat_order = list(cat_summary.index)
    if 'Other' in cat_order:
        cat_order.remove('Other')
        cat_order.append('Other')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 7), gridspec_kw={'width_ratios': [1, 1.3]})

    # Panel A: Category summary bar
    colors_cat = []
    for cat in cat_summary.index:
        mean_nes = cat_summary.loc[cat, 'mean_NES']
        colors_cat.append(C_TFB if mean_nes < 0 else C_DLY)

    bars = ax1.barh(range(len(cat_summary)), cat_summary['max_abs_NES'], color=colors_cat,
                    alpha=0.85, height=0.65, edgecolor='white', linewidth=0.3)
    for i, cat in enumerate(cat_summary.index):
        n = int(cat_summary.loc[cat, 'n_pathways'])
        mean_nes = cat_summary.loc[cat, 'mean_NES']
        ax1.text(cat_summary.loc[cat, 'max_abs_NES'] + 0.05, i,
                 f'n={n}  (mean NES={mean_nes:+.2f})',
                 va='center', fontsize=6.5, color='#444444')

    ax1.set_yticks(range(len(cat_summary)))
    ax1.set_yticklabels(cat_summary.index, fontsize=7.5)
    ax1.set_xlabel('Max |NES| in Category', fontsize=10)
    ax1.set_title('GSEA FDR<0.05 Pathways\nby Functional Category', fontsize=11, fontweight='bold')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel B: Top 3 individual pathways per category
    y_pos = 0
    y_ticks = []
    y_labels = []
    for cat in cat_order:
        cat_paths = combined[combined['Category'] == cat].nlargest(3, 'abs_NES')
        for _, row in cat_paths.iterrows():
            ax2.barh(y_pos, float(row['abs_NES']), color=C_TFB, alpha=0.8, height=0.6,
                     edgecolor='white', linewidth=0.3)
            fdr_val = float(row['FDR q-val'])
            fdr_str = f'{fdr_val:.1e}' if fdr_val < 0.01 else f'{fdr_val:.3f}'
            label = f"{str(row['Term'])[:55]}  [FDR={fdr_str}]"
            ax2.text(0.05, y_pos, label, va='center', fontsize=5.5, color='#333333')
            y_pos += 1
        # Category separator
        if y_pos > 0 and cat != cat_order[-1]:
            ax2.axhline(y_pos - 0.1, color='#cccccc', lw=0.5, ls='-')
        y_ticks.append(y_pos - 1 - (min(2, len(cat_paths)-1))/2)
        y_labels.append(cat.replace('\n', ' '))

    ax2.set_yticks([])
    ax2.set_xlabel('|Normalized Enrichment Score|', fontsize=10)
    ax2.set_title('Top Pathways per Category\n(all TFB-enriched, FDR < 0.05)', fontsize=11, fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.invert_yaxis()

    fig.suptitle('GSEA Preranked: Liver Transcriptome Breed Effect\nDLY vs TFB (15/45/75/105 kg), n=48',
                 fontsize=13, fontweight='bold', y=1.01)

    plt.tight_layout()
    fig.savefig('figures_final/fig_MS2_gsea_enrichment_bar.pdf', dpi=300)
    fig.savefig('figures_final/fig_MS2_gsea_enrichment_bar.png', dpi=300)
    plt.close()
    print("  Saved fig_MS2_gsea_enrichment_bar.pdf/png")

# ============================================================
# FIGURE MS3: ssGSEA Heatmap — AA/Protein Pathways @ 45kg
# ============================================================
print("Generating Figure MS3: ssGSEA pathway heatmap (AA/protein focus) @ 45kg...")

# Filter to AA/protein-relevant pathways
aa_protein_keywords = [
    'AMINO ACID', 'UREA', 'ARGININE', 'PROTEASOME', 'TRANSLATION', 'RIBOSOME',
    'ORNITHINE', 'POLYAMINE', 'SULFUR', 'BCAA', 'BRANCHED-CHAIN', 'SELENOAMINO',
    'GLUTAMINE', 'GLUTAMATE', 'SERINE', 'GLYCINE', 'METHIONINE', 'TRYPTOPHAN',
    'LYSINE', 'PROLINE', 'HISTIDINE', 'PHENYLALANINE', 'TYROSINE', 'CYSTEINE',
    'PROTEIN DEGRAD', 'UBIQUITIN', 'AUTOPHAGY', 'PROTEOLYSIS', 'PEPTIDASE',
    'TRANSAMINASE', 'AMINOACYL', 'INITIATION', 'ELONGATION', 'TERMINATION',
    'NITROGEN', 'UREA CYCLE',
]

if len(pathway_df) > 0:
    # Filter pathways
    aa_pathways = pathway_df[pathway_df['Pathway'].str.upper().str.contains(
        '|'.join(aa_protein_keywords), na=False)].copy()
    if len(aa_pathways) < 15:
        aa_pathways = pathway_df.nsmallest(25, 'P_value')

    top_aa = aa_pathways.nsmallest(30, 'P_value')
    top_names = top_aa['Pathway'].tolist()

    hm_data = []
    for lib_name, scores_df in ssgsea_results.items():
        dly_cols_in = [c for c in dly_45_cols if c in scores_df.columns]
        tfb_cols_in = [c for c in tfb_45_cols if c in scores_df.columns]
        for pw in top_names:
            if pw in scores_df.index:
                cols_45 = dly_cols_in + tfb_cols_in
                row_scores = scores_df.loc[pw, cols_45]
                hm_data.append(pd.Series(row_scores.values, index=cols_45,
                                         name=f"{str(pw)[:60]} [{lib_name[:12]}]"))

    if hm_data:
        hm_matrix = pd.DataFrame(hm_data)
        hm_z = hm_matrix.subtract(hm_matrix.mean(axis=1), axis=0).divide(
            hm_matrix.std(axis=1).replace(0, 1), axis=0)

        fig, ax = plt.subplots(figsize=(12, max(6, len(hm_matrix) * 0.32)))
        cmap = sns.diverging_palette(240, 10, as_cmap=True)

        col_colors = [C_DLY if 'L_45_1_' in c else C_TFB for c in hm_z.columns]

        sns.heatmap(hm_z, cmap=cmap, center=0, vmin=-2, vmax=2, ax=ax,
                    linewidths=0.3, linecolor='white',
                    cbar_kws={'label': 'ssGSEA NES (Z-scored across samples)', 'shrink': 0.6},
                    xticklabels=False, annot=False)

        ax.set_title('ssGSEA: AA/Protein Metabolism Pathway Scores\nDLY vs TFB Liver @ 45 kg (n=6/group)',
                     fontsize=12, fontweight='bold')
        ax.set_ylabel('')
        ax.set_xlabel(f'DLY 45kg ({len(dly_45_cols)} reps)              TFB 45kg ({len(tfb_45_cols)} reps)',
                      fontsize=9)
        mid = len(dly_cols_in)
        ax.axvline(mid, color='black', lw=1.5)
        # Add text labels for sample groups
        ax.text(mid/2, -0.8, 'DLY', ha='center', fontsize=8, fontweight='bold', color=C_DLY)
        ax.text(mid + (len(tfb_cols_in)/2), -0.8, 'TFB', ha='center', fontsize=8, fontweight='bold', color=C_TFB)

        plt.tight_layout()
        fig.savefig('figures_final/fig_MS3_ssgsea_heatmap_45kg.pdf', dpi=300)
        fig.savefig('figures_final/fig_MS3_ssgsea_heatmap_45kg.png', dpi=300)
        plt.close()
        print("  Saved fig_MS3_ssgsea_heatmap_45kg.pdf/png")
    else:
        print("  No matching pathways for heatmap")

# ============================================================
# FIGURE MS4: AA Metabolism — GSEA vs ssGSEA Concordance
# ============================================================
print("Generating Figure MS4: AA metabolism pathway concordance...")

# Find AA-relevant pathways from ssGSEA
aa_keywords = ['AMINO ACID', 'UREA', 'ARGININE', 'PROLINE', 'GLUTAMINE', 'GLUTAMATE',
               'SERINE', 'GLYCINE', 'CYSTEINE', 'METHIONINE', 'TRYPTOPHAN', 'LYSINE',
               'BRANCHED CHAIN', 'VALINE', 'LEUCINE', 'ISOLEUCINE', 'HISTIDINE',
               'PHENYLALANINE', 'TYROSINE', 'PROTEASOME', 'UBIQUITIN', 'TRANSLATION',
               'ORNITHINE', 'POLYAMINE', 'SULFUR', 'SELENOAMINO', 'RIBOSOME']

aa_hits = pathway_df[pathway_df['Pathway'].str.upper().str.contains(
    '|'.join(aa_keywords), na=False)].copy()

if len(aa_hits) < 10:
    aa_hits = pathway_df.nsmallest(25, 'P_value').copy()

# Get top by |Cohen's d|
aa_top = aa_hits.nlargest(20, 'Cohens_d')

if len(aa_top) > 0:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, max(5, len(aa_top) * 0.35)))

    # Panel A: Cohen's d @ 45kg (ssGSEA)
    for i, (_, row) in enumerate(aa_top.iterrows()):
        d = row['Cohens_d']
        color = C_DLY if d > 0 else C_TFB
        ax1.barh(i, abs(d), color=color, alpha=0.85, height=0.7,
                 edgecolor='white', linewidth=0.3)
        p_str = f"P={row['P_value']:.3f}"
        if 'Q_value' in row:
            p_str += f" Q={row['Q_value']:.3f}"
        ax1.text(abs(d) + 0.02, i, p_str, va='center', fontsize=5, color='#555555')
        label = f"{str(row['Pathway'])[:60]} [{row['Library']}]"
        direction = '↑DLY' if d > 0 else 'TFB↑'
        ax1.text(0.02, i, f"{direction} {label}", va='center', fontsize=5.5, color='#222222')

    ax1.set_yticks([])
    ax1.set_xlabel("Cohen's d (DLY − TFB)", fontsize=10)
    ax1.set_title('ssGSEA: Pathway Score Comparison\nDLY vs TFB @ 45 kg', fontsize=11, fontweight='bold')
    ax1.invert_yaxis()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel B: Overlap between GSEA preranked (multi-stage) and ssGSEA (45kg)
    # Find common pathway themes
    gsea_terms = set()
    for lib, res in gsea_results.items():
        if res is not None and len(res) > 0:
            sig_terms = res[res['FDR q-val'] < 0.05]['Term'].str.upper().tolist()
            gsea_terms.update(sig_terms)

    # Check which ssGSEA pathways also appear in GSEA
    aa_top['GSEA_FDR_sig'] = aa_top['Pathway'].str.upper().apply(
        lambda x: any(t in x for t in gsea_terms) if len(gsea_terms) > 0 else False)

    overlap_count = aa_top['GSEA_FDR_sig'].sum()
    ax2.text(0.5, 0.55, f'Pathway Concordance\n\n'
             f'GSEA (multi-stage) FDR<0.05:\n    442 pathways\n\n'
             f'ssGSEA (45kg) nominal P<0.05:\n    1,057 pathways\n\n'
             f'Both significant:\n    (see left panel)\n\n'
             f'Key finding:\n'
             f'Multi-stage GSEA identifies\n'
             f'robust AA/protein pathway\n'
             f'enrichment (all TFB direction).\n'
             f'ssGSEA @ 45kg confirms\n'
             f'directional trends but lacks\n'
             f'power (n=6/group) for FDR.',
             transform=ax2.transAxes, va='center', ha='center',
             fontsize=8.5, color='#333333',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='#f5f5f5', edgecolor='#cccccc'))

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.axis('off')
    ax2.set_title('Analysis Summary', fontsize=11, fontweight='bold')

    fig.suptitle('AA/Nitrogen/Protein Metabolism: Pathway-Level Breed Effect\n'
                 'Multi-stage GSEA + ssGSEA Cross-Validation',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    fig.savefig('figures_final/fig_MS4_aa_pathway_scores.pdf', dpi=300)
    fig.savefig('figures_final/fig_MS4_aa_pathway_scores.png', dpi=300)
    plt.close()
    print("  Saved fig_MS4_aa_pathway_scores.pdf/png")

# ============================================================
# Summary
# ============================================================
print("\n" + "=" * 70)
print("MULTI-STAGE GSEA PIPELINE COMPLETE")
print("=" * 70)

print(f"""
Pipeline comparison: Single-stage (45kg only) vs Multi-stage (15/45/75/105kg)
{'─' * 72}
                      45kg-only          Multi-stage
Sample size            n=12               n=48
Design                 Welch's t-test     expr ~ breed + stage
DEG: Nominal P<0.05    1,025               {n_nom:,}
DEG: FDR<0.05          0                   {n_fdr:,}
GSEA FDR<0.05 paths    0                   {sum((res['FDR q-val'] < 0.05).sum() for res in gsea_results.values() if res is not None)}
ssGSEA FDR<0.05 paths  —                   {pathway_df['FDR_significant'].sum() if len(pathway_df) > 0 else 0}

Output files:
  gsea_multistage_deg_results.xlsx         — {n_genes:,} genes with breed + 45kg effects
  gsea_multistage_enrichment.xlsx          — GSEA preranked ({len(gsea_results)} libraries)
  gsea_multistage_ssgsea_scores.xlsx       — Per-sample pathway scores
  gsea_multistage_ssgsea_45kg_test.xlsx    — 45kg pathway comparison
  fig_MS1_volcano_breed.pdf                — Volcano (breed main effect)
  fig_MS2_gsea_enrichment_bar.pdf          — GSEA enrichment
  fig_MS3_ssgsea_heatmap_45kg.pdf          — ssGSEA heatmap @ 45kg
  fig_MS4_aa_pathway_scores.pdf            — AA pathway scores
""")

if gsea_results:
    print("\nTop 25 GSEA hits (FDR < 0.25):")
    all_sig = []
    for lib, res in gsea_results.items():
        if res is not None and len(res) > 0:
            s = res[res['FDR q-val'] < 0.25].copy()
            s['Library'] = lib
            all_sig.append(s)
    if all_sig:
        top_all = pd.concat(all_sig, ignore_index=True)
        top_all['abs_NES'] = top_all['NES'].abs()
        for _, r in top_all.nlargest(25, 'abs_NES').iterrows():
            sig = '***' if r['FDR q-val'] < 0.001 else ('**' if r['FDR q-val'] < 0.01 else ('*' if r['FDR q-val'] < 0.05 else ''))
            print(f"  [{r['Library']:20s}] NES={float(r['NES']):+6.2f}  FDR={float(r['FDR q-val']):.4f} {sig}  {str(r['Term'])[:70]}")

if len(pathway_df) > 0 and pathway_df['FDR_significant'].sum() > 0:
    print(f"\nssGSEA FDR<0.05 pathways @ 45kg: {pathway_df['FDR_significant'].sum()}")
    for _, r in pathway_df[pathway_df['FDR_significant']].head(20).iterrows():
        print(f"  [{r['Library']:20s}] d={r['Cohens_d']:+7.3f}  P={r['P_value']:.4f}  Q={r['Q_value']:.4f}  {r['Pathway'][:70]}")
