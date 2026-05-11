#!/usr/bin/env python3
"""
Diagnostic: Group-mean vs Individual-level correlation comparison.
Key question: When we correlate liver enzyme expression with serum urea,
do we lose power (and inflate r) by collapsing to n=8 group means?

Tests:
  A. Liver AA enzyme expr vs Serum Urea (master_analysis.py line 284)
  B. Cross-tissue correlation matrix (liver enzyme vs muscle ribosome)
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import sys
sys.path.insert(0, '.')
from stats_utils import benjamini_hochberg, safe_pearsonr

# ============================================================
# Load data
# ============================================================
serum_tidy = pd.read_csv('serum_all_tidy.csv')
muscle_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')
liver_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')

# Sample maps
smap_m = {
    'm_15_1_': ('DLY', 15), 'm_15_2_': ('TFB', 15),
    'BJ_2_1_': ('DLY', 45), 'BJ_2_2_': ('TFB', 45),
    'm_1_1_': ('DLY', 75), 'm_1_2_': ('TFB', 75),
    'm_2_1_': ('DLY', 105), 'm_2_2_': ('TFB', 105),
    'm_3_1_': ('DLY', 135),
}
smap_l = {
    'L_15_1_': ('DLY', 15), 'L_15_2_': ('TFB', 15),
    'L_45_1_': ('DLY', 45), 'L_45_2_': ('TFB', 45),
    'L_1_1_': ('DLY', 75), 'L_1_2_': ('TFB', 75),
    'L_2_1_': ('DLY', 105), 'L_2_2_': ('TFB', 105),
    'L_3_1_': ('DLY', 135),
}

def build_individual_df(mat, smap):
    val_cols = [c for c in mat.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    records = []
    for _, row in mat.iterrows():
        gene_name = str(row['gene_name']) if pd.notna(row['gene_name']) else row['seq_id']
        for col in val_cols:
            info = None
            for prefix, (breed, stage) in smap.items():
                if col.startswith(prefix):
                    rep_num = int(col.split('_')[-1])
                    info = (breed, stage, rep_num)
                    break
            if info and pd.notna(row[col]):
                records.append({
                    'gene_name': gene_name, 'breed': info[0],
                    'stage_kg': info[1], 'rep': info[2], 'expr': float(row[col])
                })
    return pd.DataFrame(records)

print("Building individual-level data...")
liver_ind = build_individual_df(liver_raw, smap_l)
muscle_ind = build_individual_df(muscle_raw, smap_m)

# AA catabolism genes
AA_GENES = ['BCAT2', 'BCKDHA', 'BCKDHB', 'DBT', 'DLD', 'CPS1', 'OTC',
            'ASS1', 'ASL', 'ARG1', 'GOT1', 'GOT2', 'GPT', 'AASS', 'HGD',
            'ACADSB', 'GLUD1', 'SDS', 'HAL', 'PAH']

# Find available genes
available = set(liver_ind['gene_name'].unique())
aa_found = [g for g in AA_GENES if g in available]
print(f"AA enzymes found in liver: {len(aa_found)}/{len(AA_GENES)}")

# Serum urea
serum_urea = serum_tidy[serum_tidy['metabolite'] == 'Urea'].copy()
urea_dedup = serum_urea.groupby(['breed', 'stage_kg', 'rep'])['value'].mean().reset_index()

# ============================================================
# TEST A: Liver AA enzyme vs Serum Urea
# ============================================================
print("\n" + "=" * 70)
print("TEST A: Liver AA Enzyme Expression vs Serum Urea")
print("=" * 70)

results = []
for gene in aa_found:
    gene_df = liver_ind[liver_ind['gene_name'] == gene]

    # --- Method 1: Group-mean (n=8: 2 breeds × 4 stages) ---
    gm = gene_df.groupby(['breed', 'stage_kg'])['expr'].mean().reset_index()
    um = urea_dedup.groupby(['breed', 'stage_kg'])['value'].mean().reset_index()
    merged_gm = gm.merge(um, on=['breed', 'stage_kg'])
    r_gm, p_gm = pearsonr(merged_gm['expr'], merged_gm['value']) if len(merged_gm) >= 6 else (np.nan, np.nan)

    # --- Method 2: Individual-level (n≈48: 6 reps × 2 breeds × 4 stages) ---
    merged_ind = gene_df.merge(urea_dedup, on=['breed', 'stage_kg', 'rep'])
    r_ind, p_ind = pearsonr(merged_ind['expr'], merged_ind['value']) if len(merged_ind) >= 10 else (np.nan, np.nan)

    results.append({
        'Gene': gene,
        'n_group_mean': len(merged_gm),
        'r_group_mean': round(r_gm, 3) if not np.isnan(r_gm) else np.nan,
        'p_group_mean': round(p_gm, 5) if not np.isnan(p_gm) else np.nan,
        'n_individual': len(merged_ind),
        'r_individual': round(r_ind, 3) if not np.isnan(r_ind) else np.nan,
        'p_individual': round(p_ind, 5) if not np.isnan(p_ind) else np.nan,
    })

res_df = pd.DataFrame(results)
# Sort by individual-level r
res_df = res_df.sort_values('r_individual', key=abs, ascending=False, na_position='last')

print(f"\n{'Gene':10s} {'n_grp':>5s} {'r_grp':>7s} {'p_grp':>9s}  |  {'n_ind':>5s} {'r_ind':>7s} {'p_ind':>9s}  {'Δr':>6s}")
print("-" * 70)
for _, r in res_df.iterrows():
    dr = (r['r_individual'] - r['r_group_mean']) if pd.notna(r['r_individual']) and pd.notna(r['r_group_mean']) else np.nan
    print(f"{r['Gene']:10s} {r['n_group_mean']:5.0f} {r['r_group_mean']:+7.3f} {r['p_group_mean']:9.5f}  |  "
          f"{r['n_individual']:5.0f} {r['r_individual']:+7.3f} {r['p_individual']:9.5f}  "
          f"{dr:+6.3f}" if not np.isnan(dr) else f"{'':6s}")

# Summary stats
n_both = res_df.dropna(subset=['r_group_mean', 'r_individual'])
print(f"\n--- Summary (n={len(n_both)} genes with both) ---")
print(f"Group-mean:   mean|r|={n_both['r_group_mean'].abs().mean():.3f},  P<0.05: {(n_both['p_group_mean']<0.05).sum()}/{len(n_both)}")
print(f"Individual:   mean|r|={n_both['r_individual'].abs().mean():.3f},  P<0.05: {(n_both['p_individual']<0.05).sum()}/{len(n_both)}")

# FDR on individual-level
if len(n_both) > 0:
    _, q_ind = benjamini_hochberg(n_both['p_individual'].values)
    n_both = n_both.copy()
    n_both['q_individual'] = q_ind
    n_fdr = (q_ind < 0.05).sum()
    print(f"Individual+FDR: q<0.05: {n_fdr}/{len(n_both)}")

# ============================================================
# TEST B: All liver genes vs Serum Urea at individual level
# ============================================================
print("\n" + "=" * 70)
print("TEST B: Genome-wide liver gene vs Serum Urea (individual-level)")
print("=" * 70)

all_genes = sorted(liver_ind['gene_name'].unique())
print(f"Total liver genes: {len(all_genes)}")

# Take top 5000 most variable genes for speed
gene_var = liver_ind.groupby('gene_name')['expr'].var().sort_values(ascending=False)
top_genes = gene_var.head(5000).index.tolist()

n_sig_nom = 0
n_sig_fdr = 0
pvals_all = []
for i, gene in enumerate(top_genes):
    gene_df = liver_ind[liver_ind['gene_name'] == gene]
    merged = gene_df.merge(urea_dedup, on=['breed', 'stage_kg', 'rep'])
    if len(merged) >= 10:
        r, p = pearsonr(merged['expr'], merged['value'])
        pvals_all.append(p)
        if p < 0.05:
            n_sig_nom += 1
    if (i + 1) % 1000 == 0:
        print(f"  Progress: {i+1}/{len(top_genes)} genes...")

print(f"\nGenes tested: {len(pvals_all)}")
print(f"Nominal P<0.05: {n_sig_nom} ({100*n_sig_nom/len(pvals_all):.1f}%)")
if len(pvals_all) > 0:
    _, q_all = benjamini_hochberg(np.array(pvals_all))
    n_fdr_all = (q_all < 0.05).sum()
    print(f"FDR<0.05: {n_fdr_all} ({100*n_fdr_all/len(pvals_all):.1f}%)")

# ============================================================
# CONCLUSION
# ============================================================
print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
Key observations:
1. Group-mean correlations (n=8) systematically INFLATE r-values because
   within-group variance is collapsed away.
2. Individual-level correlations (n≈48) have more degrees of freedom but
   smaller r values (due to biological + technical noise within groups).
3. With n=8, an r of ~0.71 is needed for nominal P<0.05 — the bar is
   so high that many real biological signals won't pass.
4. With n≈48, r of ~0.29 is sufficient for P<0.05 — more sensitive but
   these p-values MUST be FDR-corrected.

Recommendation:
  Report individual-level correlations with FDR-corrected p-values.
  Use mixed-effects models (gene ~ breed*stage + (1|stage)) for formal testing.
  The n=8 group-mean correlations are for visualization only, not inference.
""")
