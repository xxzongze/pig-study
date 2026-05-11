"""
Clean 3-tier liver-muscle axis gene screening.
Excludes DLY105 liver data (batch artifact confirmed by PCA).
Uses liver 15/45/75kg + muscle 15/45/75/105kg.
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, spearmanr
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 1. Load data
# ============================================================
liver = pd.read_csv('/Users/hezongze/pig_study/gene_expression/liver_gene_matrix.xls', sep='\t', index_col=0)
muscle = pd.read_csv('/Users/hezongze/pig_study/gene_expression/muscle_gene_matrix.xls', sep='\t', index_col=0)

# Expression columns only
liver_expr_cols = [c for c in liver.columns if c.startswith('L_')]
muscle_expr_cols = [c for c in muscle.columns if c.startswith('m_') or c.startswith('BJ_')]

liver_expr = liver[liver_expr_cols].apply(pd.to_numeric, errors='coerce')
muscle_expr = muscle[muscle_expr_cols].apply(pd.to_numeric, errors='coerce')

# Log2 transform
liver_log = np.log2(liver_expr.values + 1)
muscle_log = np.log2(muscle_expr.values + 1)

# Filter low expression (mean > 1 log2)
liver_keep = liver_log.mean(axis=1) > 1
muscle_keep = muscle_log.mean(axis=1) > 1
liver_log = liver_log[liver_keep, :]
muscle_log = muscle_log[muscle_keep, :]

print(f"Liver: {liver_log.shape[0]} genes x {liver_log.shape[1]} samples")
print(f"Muscle: {muscle_log.shape[0]} genes x {muscle_log.shape[1]} samples")

# ============================================================
# 2. Define groups
# ============================================================
# Liver group definitions (Sample Analysis Name based)
# L_15_1 = DLY 15kg, L_15_2 = TFB 15kg
# L_45_1 = DLY 45kg, L_45_2 = TFB 45kg
# L_1_1 = DLY 75kg, L_1_2 = TFB 75kg
# L_2_1 = DLY 105kg (EXCLUDED - batch artifact), L_2_2 = TFB 105kg
# L_3_1 = DLY 135kg (no TFB counterpart)

liver_groups = {
    '15kg_DLY': [i for i,c in enumerate(liver_expr_cols) if c.startswith('L_15_1')],
    '15kg_TFB': [i for i,c in enumerate(liver_expr_cols) if c.startswith('L_15_2')],
    '45kg_DLY': [i for i,c in enumerate(liver_expr_cols) if c.startswith('L_45_1')],
    '45kg_TFB': [i for i,c in enumerate(liver_expr_cols) if c.startswith('L_45_2')],
    '75kg_DLY': [i for i,c in enumerate(liver_expr_cols) if c.startswith('L_1_1')],
    '75kg_TFB': [i for i,c in enumerate(liver_expr_cols) if c.startswith('L_1_2')],
}

# Muscle group definitions
muscle_groups = {
    '15kg_DLY': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('m_15_1')],
    '15kg_TFB': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('m_15_2')],
    '45kg_DLY': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('BJ_2_1')],
    '45kg_TFB': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('BJ_2_2')],
    '75kg_DLY': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('m_1_1')],
    '75kg_TFB': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('m_1_2')],
    '105kg_DLY': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('m_2_1')],
    '105kg_TFB': [i for i,c in enumerate(muscle_expr_cols) if c.startswith('m_2_2')],
}

# ============================================================
# 3. Compute log2FC and p-values per stage
# ============================================================
def compute_stage_stats(expr_log, groups, stage_name):
    """Compute log2(DLY/TFB) and t-test p-value for a stage."""
    dly_idx = groups[f'{stage_name}_DLY']
    tfb_idx = groups[f'{stage_name}_TFB']

    dly_mean = expr_log[:, dly_idx].mean(axis=1)
    tfb_mean = expr_log[:, tfb_idx].mean(axis=1)
    log2fc = dly_mean - tfb_mean  # positive = DLY higher

    # Welch t-test per gene
    pvals = np.ones(expr_log.shape[0])
    for g in range(expr_log.shape[0]):
        if dly_mean[g] == tfb_mean[g]:
            continue
        try:
            t, p = ttest_ind(expr_log[g, dly_idx], expr_log[g, tfb_idx], equal_var=False)
            pvals[g] = p
        except:
            pass

    return log2fc, pvals, dly_mean, tfb_mean

# Liver stats for 3 reliable stages
liver_stages = ['15kg', '45kg', '75kg']
liver_results = {}
for stage in liver_stages:
    fc, pv, dm, tm = compute_stage_stats(liver_log, liver_groups, stage)
    liver_results[stage] = {'log2FC': fc, 'pvalue': pv, 'DLY_mean': dm, 'TFB_mean': tm}

# Muscle stats
muscle_stages = ['15kg', '45kg', '75kg', '105kg']
muscle_results = {}
for stage in muscle_stages:
    fc, pv, dm, tm = compute_stage_stats(muscle_log, muscle_groups, stage)
    muscle_results[stage] = {'log2FC': fc, 'pvalue': pv, 'DLY_mean': dm, 'TFB_mean': tm}

# ============================================================
# 4. Tier classification
# ============================================================
def classify_tiers(results, stages):
    """
    Tier 1 (Programming): |FC|>0.5 at 15kg AND direction consistent with later stages
    Tier 2 (Switch): |FC| at 45kg is max among all stages
    Tier 3 (Consequence): |FC| at 75kg is max among all stages
    """
    n_genes = len(results[stages[0]]['log2FC'])
    tier = np.full(n_genes, 'Unclassified', dtype=object)
    tier_score = np.zeros(n_genes)  # For ranking within tier

    fc_15 = results[stages[0]]['log2FC']
    fc_45 = results[stages[1]]['log2FC']
    fc_75 = results[stages[2]]['log2FC']

    abs_fc = {s: np.abs(results[s]['log2FC']) for s in stages}

    for g in range(n_genes):
        fcs = [results[s]['log2FC'][g] for s in stages]
        afcs = [abs_fc[s][g] for s in stages]

        # Tier 1: 15kg |FC| > 0.5 and direction agrees with majority of later stages
        if afcs[0] > 0.5:
            sign_15 = np.sign(fcs[0])
            later_signs = [np.sign(f) for f in fcs[1:] if abs(f) > 0.3]
            if not later_signs or all(s == sign_15 for s in later_signs):
                tier[g] = 'Tier1_Programming'
                tier_score[g] = afcs[0] + np.mean(afcs[1:])
                continue

        # Tier 2: 45kg has max |FC|
        if afcs[1] >= afcs[0] and afcs[1] >= afcs[2] and afcs[1] > 0.5:
            tier[g] = 'Tier2_Switch'
            tier_score[g] = afcs[1]
            continue

        # Tier 3: 75kg has max |FC|
        if afcs[2] >= afcs[0] and afcs[2] >= afcs[1] and afcs[2] > 0.5:
            tier[g] = 'Tier3_Consequence'
            tier_score[g] = afcs[2]
            continue

        # Low signal
        if max(afcs) <= 0.5:
            tier[g] = 'Low_Signal'
            tier_score[g] = max(afcs)
        else:
            tier[g] = 'Mixed'
            tier_score[g] = max(afcs)

    return tier, tier_score

liver_tier, liver_tier_score = classify_tiers(liver_results, liver_stages)
muscle_tier, muscle_tier_score = classify_tiers(muscle_results, muscle_stages)

# ============================================================
# 5. Build output tables
# ============================================================
# Gene IDs from filtered indices
liver_gene_ids = liver.index[liver_keep]
muscle_gene_ids = muscle.index[muscle_keep]

def build_gene_table(gene_ids, results, stages, tier, tier_score, tissue):
    """Build comprehensive gene table."""
    data = {'gene_id': gene_ids}

    for s in stages:
        data[f'{s}_log2FC'] = results[s]['log2FC']
        data[f'{s}_pvalue'] = results[s]['pvalue']
        data[f'{s}_DLY_mean'] = results[s]['DLY_mean']
        data[f'{s}_TFB_mean'] = results[s]['TFB_mean']

    data['Tier'] = tier
    data['Tier_Score'] = tier_score

    # Mean absolute FC across stages
    mean_abs_fc = np.mean([np.abs(results[s]['log2FC']) for s in stages], axis=0)
    data['Mean_abs_log2FC'] = mean_abs_fc

    df = pd.DataFrame(data)
    df['Tissue'] = tissue

    # Sort by Tier then Score
    tier_order = {'Tier1_Programming': 0, 'Tier2_Switch': 1, 'Tier3_Consequence': 2, 'Mixed': 3, 'Low_Signal': 4}
    df['_sort'] = df['Tier'].map(tier_order)
    df = df.sort_values(['_sort', 'Tier_Score'], ascending=[True, False])
    df = df.drop(columns=['_sort'])

    return df

liver_df = build_gene_table(liver_gene_ids, liver_results, liver_stages, liver_tier, liver_tier_score, 'Liver')
muscle_df = build_gene_table(muscle_gene_ids, muscle_results, muscle_stages, muscle_tier, muscle_tier_score, 'Muscle')

# ============================================================
# 6. Summary statistics
# ============================================================
print(f"\n=== Liver Gene Classification (15/45/75kg, DLY105 excluded) ===")
for t in ['Tier1_Programming', 'Tier2_Switch', 'Tier3_Consequence', 'Mixed', 'Low_Signal']:
    n = (liver_tier == t).sum()
    print(f"  {t}: {n} genes ({100*n/len(liver_tier):.1f}%)")

print(f"\n=== Muscle Gene Classification ===")
for t in ['Tier1_Programming', 'Tier2_Switch', 'Tier3_Consequence', 'Mixed', 'Low_Signal']:
    n = (muscle_tier == t).sum()
    print(f"  {t}: {n} genes ({100*n/len(muscle_tier):.1f}%)")

# ============================================================
# 7. Save results
# ============================================================
outdir = '/Users/hezongze/pig_study'

# Full tables
liver_df.to_excel(f'{outdir}/clean_liver_tier_genes.xlsx', index=False)
muscle_df.to_excel(f'{outdir}/clean_muscle_tier_genes.xlsx', index=False)

# Tier 1+2 only (the causal candidates)
liver_causal = liver_df[liver_df['Tier'].isin(['Tier1_Programming', 'Tier2_Switch'])]
muscle_causal = muscle_df[muscle_df['Tier'].isin(['Tier1_Programming', 'Tier2_Switch'])]

# Combined liver-muscle causal genes
combined = pd.concat([liver_causal, muscle_causal], ignore_index=True)
combined.to_excel(f'{outdir}/clean_combined_causal_genes.xlsx', index=False)

print(f"\n=== Output files ===")
print(f"  clean_liver_tier_genes.xlsx — {len(liver_df)} liver genes")
print(f"  clean_muscle_tier_genes.xlsx — {len(muscle_df)} muscle genes")
print(f"  clean_combined_causal_genes.xlsx — {len(combined)} causal genes (Tier1+2)")

# ============================================================
# 8. Key gene highlight
# ============================================================
# Search for AA metabolism genes in liver Tier 1
aa_genes_liver = ['CPS1','ASS1','ASL','ARG1','ARG2','GOT1','GOT2','GPT','BCAT1','BCAT2',
                  'BCKDHA','BCKDHB','DBT','DLD','AASS','HGD','SDS','GLUD1','OTC','PAH',
                  'HAL','IVD','MCCC1','MCCC2','AUH','HIBCH','HIBADH','ALDH6A1',
                  'IGF1','IGFBP1','IGFBP2','IGFBP3','FGF21','FST','MSTN','BDNF']

print(f"\n=== Key AA/Liver metabolism genes in Tier 1+2 ===")
for gene in aa_genes_liver:
    matches = liver_df[liver_df['gene_id'].str.contains(gene, case=False)]
    if len(matches) > 0:
        for _, row in matches.iterrows():
            fcs = f"15kg={row['15kg_log2FC']:.2f}, 45kg={row['45kg_log2FC']:.2f}, 75kg={row['75kg_log2FC']:.2f}"
            print(f"  {row['gene_id'][:40]:<42} {row['Tier']:<22} {fcs}")

print("\nDone.")
