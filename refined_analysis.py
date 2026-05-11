#!/usr/bin/env python3
"""
Refined tier classification: EXCLUDE DLY 105kg liver entirely.
Only use 15/45/75kg for breed comparison. Use TFB 105kg + DLY 135kg
as reference for longitudinal trends only.

Plus: STAT3 regulon deep-dive & prioritized validation gene list.
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import re

print("=" * 60)
print("Refined Analysis: Exclude DLY 105kg + STAT3 Regulon")
print("=" * 60)

# ---- Load ----
serum_tidy = pd.read_csv('serum_all_tidy.csv')
muscle = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')
liver = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')

sample_to_group_m = {
    'm_15_1_': ('DLY', 15), 'm_15_2_': ('TFB', 15),
    'BJ_2_1_': ('DLY', 45), 'BJ_2_2_': ('TFB', 45),
    'm_1_1_': ('DLY', 75), 'm_1_2_': ('TFB', 75),
    'm_2_1_': ('DLY', 105), 'm_2_2_': ('TFB', 105),
    'm_3_1_': ('DLY', 135),
}
sample_to_group_l = {
    'L_15_1_': ('DLY', 15), 'L_15_2_': ('TFB', 15),
    'L_45_1_': ('DLY', 45), 'L_45_2_': ('TFB', 45),
    'L_1_1_': ('DLY', 75), 'L_1_2_': ('TFB', 75),
    'L_2_1_': ('DLY', 105), 'L_2_2_': ('TFB', 105),
    'L_3_1_': ('DLY', 135),
}

def build_individual_df(mat, sample_map):
    val_cols = [c for c in mat.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    records = []
    for _, row in mat.iterrows():
        gene_name = str(row['gene_name']) if pd.notna(row['gene_name']) else row['seq_id']
        for col in val_cols:
            info = None
            for prefix, (breed, stage) in sample_map.items():
                if col.startswith(prefix):
                    info = (breed, stage, int(col.split('_')[-1]))
                    break
            if info and pd.notna(row[col]):
                records.append({
                    'gene_name': gene_name, 'breed': info[0],
                    'stage_kg': info[1], 'rep': info[2], 'expr': float(row[col])
                })
    return pd.DataFrame(records)

liver_ind = build_individual_df(liver, sample_to_group_l)
muscle_ind = build_individual_df(muscle, sample_to_group_m)

AA_GENES = [
    'AASS', 'ACADSB', 'ARG1', 'ARG2', 'ASL', 'ASS1', 'BCAT1', 'BCAT2',
    'BCKDHA', 'BCKDHB', 'CPS1', 'DBT', 'DLD', 'GLUD1', 'GOT1', 'GOT2',
    'HAL', 'HGD', 'OTC', 'PAH', 'SDS'
]

# ============================================================
# 1. Re-classify using ONLY 15/45/75kg (exclude DLY 105kg)
# ============================================================
print("\n" + "=" * 60)
print("1. Re-classification WITHOUT DLY 105kg Liver")
print("=" * 60)

# For reference: keep TFB 105kg and DLY 135kg as longitudinal context
# but DO NOT use for DLY vs TFB comparison

def classify_no105(gene_name, tissue_df):
    """Tier classification using ONLY 15/45/75kg for breed comparison."""
    gdf = tissue_df[tissue_df['gene_name'] == gene_name]
    if len(gdf) < 6:
        return None, None, {}

    # Per-stage log2FC (DLY/TFB) — only valid stages
    fcs = {}
    for s in [15, 45, 75]:
        dly_vals = gdf[(gdf['breed'] == 'DLY') & (gdf['stage_kg'] == s)]['expr']
        tfb_vals = gdf[(gdf['breed'] == 'TFB') & (gdf['stage_kg'] == s)]['expr']
        if len(dly_vals) > 0 and len(tfb_vals) > 0 and tfb_vals.mean() > 0:
            fcs[s] = np.log2(dly_vals.mean() / tfb_vals.mean())
        else:
            fcs[s] = np.nan

    # Also get TFB 105kg and DLY 135kg for longitudinal reference
    tfb_105 = gdf[(gdf['breed'] == 'TFB') & (gdf['stage_kg'] == 105)]['expr'].mean()
    dly_135 = gdf[(gdf['breed'] == 'DLY') & (gdf['stage_kg'] == 135)]['expr'].mean()

    valid = {s: fcs[s] for s in [15, 45, 75] if pd.notna(fcs[s])}
    if len(valid) < 2:
        return None, None, fcs

    mean_fc = np.mean(list(valid.values()))
    signs = [1 if v > 0 else -1 for v in valid.values()]
    consistent = len(set(signs)) == 1

    # Tier 1: consistent direction across 15+45+75, |FC15| > 0.5
    if 15 in valid and abs(valid[15]) > 0.5 and consistent:
        trend = 'TFB↑' if mean_fc < 0 else 'DLY↑'
        # Check if FC grows or shrinks
        fc_trend = 'growing' if abs(valid[75]) > abs(valid[15]) * 1.2 else \
                   ('shrinking' if abs(valid[75]) < abs(valid[15]) * 0.7 else 'stable')
        return 1, f'{trend}, {fc_trend}', fcs

    # Tier 2: peak at 45kg or divergence emerges at 45
    if 45 in valid and abs(valid[45]) > 0.5:
        fc15_abs = abs(valid.get(15, 0))
        fc45_abs = abs(valid[45])
        fc75_abs = abs(valid.get(75, 0))
        if fc45_abs >= fc15_abs and fc45_abs >= fc75_abs * 0.8:
            if fc15_abs < 0.3:
                return 2, f'45kg-emergent (TFB↑)' if valid[45] < 0 else '45kg-emergent (DLY↑)', fcs
            else:
                trend = 'TFB↑' if valid[45] < 0 else 'DLY↑'
                return 2, f'45kg-peak ({trend})', fcs

    # Tier 3: only at 75kg (since 105 excluded)
    if 75 in valid and abs(valid[75]) > 0.5:
        trend = 'TFB↑' if valid[75] < 0 else 'DLY↑'
        return 3, f'75kg-late ({trend})', fcs

    # Low signal
    if all(abs(v) < 0.5 for v in valid.values()):
        return 99, 'Low signal', fcs

    return None, None, fcs

# Apply to AA enzymes
print("\nAA Catabolism Enzymes (15/45/75kg only, NO DLY 105kg):")
print(f"{'Gene':10s} {'Tier':5s} {'Pattern':35s} {'FC15':>8s} {'FC45':>8s} {'FC75':>8s} {'TFB105':>8s} {'DLY135':>8s}")
print("-" * 100)
aa_refined = []
for gene in AA_GENES:
    tier, desc, fcs = classify_no105(gene, liver_ind)
    gdf = liver_ind[liver_ind['gene_name'] == gene]
    tfb105 = gdf[(gdf['breed'] == 'TFB') & (gdf['stage_kg'] == 105)]['expr'].mean()
    dly135 = gdf[(gdf['breed'] == 'DLY') & (gdf['stage_kg'] == 135)]['expr'].mean()
    tfb105_str = f'{tfb105:.2f}' if pd.notna(tfb105) else 'NA'
    dly135_str = f'{dly135:.2f}' if pd.notna(dly135) else 'NA'
    fc15 = f"{fcs.get(15, np.nan):+.3f}" if pd.notna(fcs.get(15, np.nan)) else 'NA'
    fc45 = f"{fcs.get(45, np.nan):+.3f}" if pd.notna(fcs.get(45, np.nan)) else 'NA'
    fc75 = f"{fcs.get(75, np.nan):+.3f}" if pd.notna(fcs.get(75, np.nan)) else 'NA'
    tier_str = str(tier) if tier else '?'
    desc_str = str(desc) if desc else ''
    print(f"{gene:10s} T{tier_str:3s} {desc_str:35s} {fc15:>8s} {fc45:>8s} {fc75:>8s} {tfb105_str:>8s} {dly135_str:>8s}")
    aa_refined.append({
        'Gene': gene, 'Tier_no105': tier if tier else 99,
        'Pattern': desc if desc else '',
        'log2FC_15': fcs.get(15, np.nan), 'log2FC_45': fcs.get(45, np.nan),
        'log2FC_75': fcs.get(75, np.nan),
        'TFB105_expr': tfb105, 'DLY135_expr': dly135,
    })

aa_refined_df = pd.DataFrame(aa_refined)

# Also re-classify cross-talk genes
CROSSTALK_GENES = [
    'IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5',
    'AHSG', 'RBP4', 'FGG', 'FGB', 'APOB', 'BDNF', 'SERPINC1',
    'XBP1', 'ATF4', 'FGF21', 'ANGPTL4', 'ANGPTL8',
    'MSTN', 'FST', 'FNDC5', 'SPARC', 'DCN', 'LIF', 'IL6', 'VEGFA',
]

ct_refined = []
for gene in CROSSTALK_GENES:
    l_tier, l_desc, l_fcs = classify_no105(gene, liver_ind)
    m_tier, m_desc, m_fcs = classify_no105(gene, muscle_ind)
    if l_tier or m_tier:
        ct_refined.append({
            'Gene': gene,
            'Liver_Tier': l_tier if l_tier else 99,
            'Liver_Pattern': l_desc if l_desc else '',
            'Liver_FC15': l_fcs.get(15, np.nan) if l_fcs else np.nan,
            'Liver_FC45': l_fcs.get(45, np.nan) if l_fcs else np.nan,
            'Liver_FC75': l_fcs.get(75, np.nan) if l_fcs else np.nan,
            'Muscle_Tier': m_tier if m_tier else 99,
            'Muscle_Pattern': m_desc if m_desc else '',
        })
ct_refined_df = pd.DataFrame(ct_refined)

print(f"\nCrosstalk genes re-classified:")
for tier in [1, 2, 3]:
    subset = ct_refined_df[ct_refined_df['Liver_Tier'] == tier]
    if len(subset) > 0:
        genes_str = ', '.join(subset['Gene'].tolist())
        print(f"  Liver T{tier} ({len(subset)}): {genes_str}")

# ============================================================
# 2. STAT3 Regulon Deep-Dive
# ============================================================
print("\n" + "=" * 60)
print("2. STAT3 Regulon Deep-Dive")
print("=" * 60)

# STAT3 expression pattern
stat3_df = liver_ind[liver_ind['gene_name'] == 'STAT3']
print("\nSTAT3 liver expression (mean ± SEM):")
for s in [15, 45, 75, 105]:
    for breed in ['DLY', 'TFB']:
        vals = stat3_df[(stat3_df['breed'] == breed) & (stat3_df['stage_kg'] == s)]['expr']
        if len(vals) > 0:
            print(f"  {breed} {s}kg: {vals.mean():.2f} ± {vals.std()/np.sqrt(len(vals)):.2f} (n={len(vals)})")

# STAT3 per-stage correlation with each AA enzyme at INDIVIDUAL level
print("\nSTAT3 vs AA Enzymes — Per-Stage Individual-Level Correlations:")
stat3_per_stage = {}
liver_bs_mean = liver_ind.groupby(['gene_name', 'breed', 'stage_kg'])['expr'].mean().reset_index()

stat3_corrs = []
for gene in AA_GENES:
    # Across all individuals
    stat3_all = liver_ind[liver_ind['gene_name'] == 'STAT3'][['breed', 'stage_kg', 'rep', 'expr']].rename(columns={'expr': 'stat3_expr'})
    gene_all = liver_ind[liver_ind['gene_name'] == gene][['breed', 'stage_kg', 'rep', 'expr']]

    merged = stat3_all.merge(gene_all, on=['breed', 'stage_kg', 'rep'])
    merged = merged.dropna()

    if len(merged) < 8:
        continue

    r_all, p_all = pearsonr(merged['stat3_expr'], merged['expr'])
    rho_all, ps_all = spearmanr(merged['stat3_expr'], merged['expr'])

    # Per-stage
    stage_data = {}
    for s in [15, 45, 75, 105]:
        sd = merged[merged['stage_kg'] == s]
        if len(sd) >= 4:
            rs, ps = pearsonr(sd['stat3_expr'], sd['expr'])
            stage_data[s] = (rs, ps)

    stat3_corrs.append({
        'Gene': gene,
        'r_all': round(r_all, 4), 'p_all': round(p_all, 8),
        'rho_all': round(rho_all, 4),
        'r_15': round(stage_data.get(15, (np.nan, np.nan))[0], 3),
        'r_45': round(stage_data.get(45, (np.nan, np.nan))[0], 3),
        'r_75': round(stage_data.get(75, (np.nan, np.nan))[0], 3),
        'r_105': round(stage_data.get(105, (np.nan, np.nan))[0], 3),
        'n_total': len(merged),
    })

stat3_corr_df = pd.DataFrame(stat3_corrs).sort_values('r_all', key=abs, ascending=False)
print(f"\n{'Gene':10s} {'r_all':>7s} {'p_all':>10s} {'r_15':>7s} {'r_45':>7s} {'r_75':>7s} {'r_105':>7s} {'n':>4s}")
print("-" * 65)
for _, r in stat3_corr_df.iterrows():
    sig = '***' if r['p_all'] < 0.001 else ('**' if r['p_all'] < 0.01 else ('*' if r['p_all'] < 0.05 else ''))
    print(f"{r['Gene']:10s} {r['r_all']:7.3f} {r['p_all']:10.6f}{sig:3s} {r['r_15']:7.3f} {r['r_45']:7.3f} {r['r_75']:7.3f} {r['r_105']:7.3f} {r['n_total']:4d}")

# STAT3 vs Serum Urea
print("\nSTAT3 vs Serum Urea:")
serum_urea_ind = serum_tidy[serum_tidy['metabolite'] == 'Urea'].copy()
parsed = serum_urea_ind['group'].apply(lambda g: ('DLY' if 'DLY' in g else 'TFB', int(re.search(r'(\d+)', g).group(1))))
serum_urea_ind['breed'] = [p[0] for p in parsed]
serum_urea_ind['stage_kg'] = [p[1] for p in parsed]
serum_urea_bs = serum_urea_ind.groupby(['breed', 'stage_kg'])['value'].mean().reset_index()
serum_urea_bs.rename(columns={'value': 'serum_urea'}, inplace=True)

stat3_merged = stat3_all.merge(serum_urea_bs, on=['breed', 'stage_kg'])
r_su, p_su = pearsonr(stat3_merged['stat3_expr'], stat3_merged['serum_urea'])
print(f"  STAT3 vs Serum Urea (breed×stage means): r={r_su:+.4f}, p={p_su:.6f}, n={len(stat3_merged)}")

# ============================================================
# 3. Generate prioritized validation gene list
# ============================================================
print("\n" + "=" * 60)
print("3. Prioritized Validation Gene List")
print("=" * 60)

# Merge tier info with STAT3 correlations and N balance correlations
validation_list = []
for _, row in aa_refined_df.iterrows():
    gene = row['Gene']
    tier = row['Tier_no105']

    # Get STAT3 correlation for this gene
    s3 = stat3_corr_df[stat3_corr_df['Gene'] == gene]
    s3_r = s3['r_all'].values[0] if len(s3) > 0 else np.nan
    s3_p = s3['p_all'].values[0] if len(s3) > 0 else np.nan

    # Get serum Urea correlation
    gene_df = liver_ind[liver_ind['gene_name'] == gene]
    gene_merged = gene_df.merge(serum_urea_bs, on=['breed', 'stage_kg'], how='left')
    gene_merged = gene_merged.dropna(subset=['expr', 'serum_urea'])
    if len(gene_merged) >= 8:
        r_urea, p_urea = pearsonr(gene_merged['expr'], gene_merged['serum_urea'])
    else:
        r_urea, p_urea = np.nan, np.nan

    # Consistency score: min FC across 15/45/75 if same direction, else 0
    fcs_valid = [row['log2FC_15'], row['log2FC_45'], row['log2FC_75']]
    fcs_valid = [v for v in fcs_valid if pd.notna(v)]
    if len(fcs_valid) >= 2:
        signs = [1 if v > 0 else -1 for v in fcs_valid]
        if len(set(signs)) == 1:
            consistency = min(abs(v) for v in fcs_valid)
        else:
            consistency = 0
    else:
        consistency = 0

    # Priority score:
    # P0: Tier 1 + consistent + strong STAT3 correlation
    # P1: Tier 1 or Tier 2 + moderate signal
    # P2: Tier 3 or low signal
    if tier == 1 and consistency > 0.5:
        priority = 'P0'
    elif tier == 1 or (tier == 2 and consistency > 0.3):
        priority = 'P1'
    elif tier == 2:
        priority = 'P1'
    else:
        priority = 'P2'

    validation_list.append({
        'Gene': gene,
        'Priority': priority,
        'Tier_no105': tier,
        'FC15': round(row['log2FC_15'], 3) if pd.notna(row['log2FC_15']) else '',
        'FC45': round(row['log2FC_45'], 3) if pd.notna(row['log2FC_45']) else '',
        'FC75': round(row['log2FC_75'], 3) if pd.notna(row['log2FC_75']) else '',
        'Consistency': round(consistency, 3),
        'STAT3_r': round(s3_r, 3) if pd.notna(s3_r) else '',
        'Urea_r': round(r_urea, 3) if pd.notna(r_urea) else '',
        'Pattern': row['Pattern'],
    })

val_df = pd.DataFrame(validation_list).sort_values(['Priority', 'Gene'])
print(f"\n{'Gene':10s} {'Pri':4s} {'Tier':5s} {'FC15':>7s} {'FC45':>7s} {'FC75':>7s} {'Cons':>6s} {'STAT3_r':>7s} {'Urea_r':>7s} {'Pattern'}")
print("-" * 100)
for _, r in val_df.iterrows():
    print(f"{r['Gene']:10s} {r['Priority']:4s} T{r['Tier_no105']:<4d} {str(r['FC15']):>7s} {str(r['FC45']):>7s} {str(r['FC75']):>7s} {str(r['Consistency']):>6s} {str(r['STAT3_r']):>7s} {str(r['Urea_r']):>7s} {r['Pattern']}")

# ============================================================
# 4. FIGURE: Old vs New Tier Comparison
# ============================================================
print("\nGenerating Figure: Tier Classification Comparison...")

# Load old tier data from the previous run
old_aa = pd.read_excel('advanced_analysis_results.xlsx', sheet_name='AA_Enzymes_Tier')
old_aa = old_aa.rename(columns={'Tier': 'Tier_old'})

merged = aa_refined_df.merge(old_aa[['Gene', 'Tier_old']], on='Gene', how='left')

fig, axes = plt.subplots(1, 2, figsize=(16, max(8, len(merged) * 0.45)))

# Panel A: Old classification (with DLY 105kg)
old_fc_cols = [c for c in old_aa.columns if c.startswith('log2FC_')]
old_heatmap = old_aa.set_index('Gene')[old_fc_cols].copy()
old_heatmap.columns = [c.replace('log2FC_', '').replace('kg', '') for c in old_fc_cols]
old_heatmap = old_heatmap.clip(-6, 6)
old_heatmap = old_heatmap.loc[[g for g in AA_GENES if g in old_heatmap.index]]

# Sort by old tier
old_tier_order = old_aa.set_index('Gene').loc[old_heatmap.index, 'Tier_old'].sort_values()
old_heatmap = old_heatmap.loc[old_tier_order.index]

sns.heatmap(old_heatmap, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
            vmin=-6, vmax=6, ax=axes[0], cbar_kws={'label': 'log2(DLY/TFB)'},
            linewidths=0.5, linecolor='white', annot_kws={'fontsize': 8})

# Add tier color bar
tier_colors_old = {1: '#4CAF50', 2: '#FF9800', 3: '#9E9E9E'}
for i, gene in enumerate(old_heatmap.index):
    tier = old_aa.set_index('Gene').loc[gene, 'Tier_old']
    axes[0].add_patch(plt.Rectangle((-0.10, i), 0.05, 1, facecolor=tier_colors_old.get(tier, 'white'),
                                     edgecolor='none', transform=axes[0].transData, clip_on=False))
axes[0].set_title('Old: With DLY 105kg\n(many genes → Tier 3 artifact)', fontsize=12, fontweight='bold')

# Panel B: New classification (15/45/75 only)
new_fc_data = merged.set_index('Gene')[['log2FC_15', 'log2FC_45', 'log2FC_75']].copy()
new_fc_data.columns = ['15', '45', '75']
new_fc_data = new_fc_data.clip(-6, 6)
new_fc_data = new_fc_data.loc[old_heatmap.index]  # Same gene order

sns.heatmap(new_fc_data, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
            vmin=-6, vmax=6, ax=axes[1], cbar_kws={'label': 'log2(DLY/TFB)'},
            linewidths=0.5, linecolor='white', annot_kws={'fontsize': 8})

tier_colors_new = {1: '#4CAF50', 2: '#FF9800', 3: '#9E9E9E', 99: '#BDBDBD'}
for i, gene in enumerate(new_fc_data.index):
    tier = merged.set_index('Gene').loc[gene, 'Tier_no105']
    axes[1].add_patch(plt.Rectangle((-0.10, i), 0.05, 1, facecolor=tier_colors_new.get(tier, 'white'),
                                     edgecolor='none', transform=axes[1].transData, clip_on=False))
axes[1].set_title('New: 15/45/75kg Only\n(DLY 105kg excluded)', fontsize=12, fontweight='bold')

fig.suptitle('Tier Classification Before vs After Excluding DLY 105kg Liver Data\n'
             'Green=Tier1 Orange=Tier2 Gray=Tier3 LightGray=LowSignal',
             fontsize=13, fontweight='bold')
plt.tight_layout()
fig.savefig('fig_tier_reclassification_comparison.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_tier_reclassification_comparison.pdf', bbox_inches='tight')
plt.close()
print("Saved fig_tier_reclassification_comparison.png/pdf")

# ============================================================
# 5. FIGURE: STAT3 Regulon — Per-Stage Correlation Profile
# ============================================================
print("Generating Figure: STAT3 Regulon...")

# STAT3 expression trajectory + top correlated enzymes
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Panel A: STAT3 liver expression
ax = axes[0]
for breed, color, marker in [('DLY', '#2196F3', 'o'), ('TFB', '#C62828', 's')]:
    bs = stat3_df[stat3_df['breed'] == breed].groupby('stage_kg')['expr']
    means = bs.mean()
    sems = bs.std() / np.sqrt(bs.count())
    ax.errorbar(means.index, means.values, yerr=sems.values, marker=marker,
               color=color, linewidth=2, markersize=8, capsize=4, label=breed)
ax.set_title('STAT3 Liver Expression', fontsize=12, fontweight='bold')
ax.set_xlabel('Stage (kg)'); ax.set_ylabel('Expression')
ax.legend(); ax.grid(alpha=0.3)
ax.set_xticks([15, 45, 75, 105])

# Panel B: STAT3 vs Top AA enzymes (scatter, all individuals)
ax = axes[1]
top_stat3_targets = stat3_corr_df.dropna(subset=['r_all']).head(5)['Gene'].tolist()
stat3_all_df = liver_ind[liver_ind['gene_name'] == 'STAT3'][['breed', 'stage_kg', 'rep', 'expr']].rename(columns={'expr': 'stat3_expr'})

x_pos = 0
yticks = []
for gene in top_stat3_targets:
    gene_all_df = liver_ind[liver_ind['gene_name'] == gene][['breed', 'stage_kg', 'rep', 'expr']]
    merged = stat3_all_df.merge(gene_all_df, on=['breed', 'stage_kg', 'rep'])
    r, p = pearsonr(merged['stat3_expr'], merged['expr'])
    ax.barh(x_pos, r, color='#E91E63' if r > 0 else '#2196F3', alpha=0.8, edgecolor='black')
    yticks.append(gene)
    x_pos += 1
ax.set_yticks(range(len(yticks)))
ax.set_yticklabels(yticks)
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_xlabel("Pearson's r (STAT3 vs Target)")
ax.set_title('Top 5 STAT3 ↔ AA Enzyme\n(Individual-Level)', fontsize=11, fontweight='bold')
ax.grid(axis='x', alpha=0.3)

# Panel C: STAT3 vs Serum Urea per stage
ax = axes[2]
serum_urea_bs_renamed = serum_urea_bs.rename(columns={'serum_urea': 'value'})
for s, color in zip([15, 45, 75, 105], ['#81C784', '#FFB74D', '#64B5F6', '#E57373']):
    sd = stat3_merged[stat3_merged['stage_kg'] == s]
    if len(sd) >= 4:
        ax.scatter(sd['stat3_expr'], sd['serum_urea'], c=color, s=60, alpha=0.8,
                  edgecolors='black', linewidth=0.5, label=f'{s} kg', zorder=5)

# Add trend line
r_all, p_all = pearsonr(stat3_merged['stat3_expr'], stat3_merged['serum_urea'])
x_range = np.linspace(stat3_merged['stat3_expr'].min(), stat3_merged['stat3_expr'].max(), 100)
z = np.polyfit(stat3_merged['stat3_expr'], stat3_merged['serum_urea'], 1)
ax.plot(x_range, np.polyval(z, x_range), 'k--', linewidth=1.5, alpha=0.6)
ax.set_xlabel('STAT3 Expression'); ax.set_ylabel('Serum Urea (mmol/L)')
ax.set_title(f'STAT3 vs Serum Urea\nr={r_all:+.3f}, p={p_all:.4f}', fontsize=11, fontweight='bold')
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.suptitle('STAT3 Regulon: A Master Regulator of Liver AA Catabolism',
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig('fig_STAT3_regulon.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_STAT3_regulon.pdf', bbox_inches='tight')
plt.close()
print("Saved fig_STAT3_regulon.png/pdf")

# ============================================================
# 6. FIGURE: Priority Gene Summary (Publication-Ready)
# ============================================================
print("Generating Figure: Priority Gene Summary...")

p0_genes = val_df[val_df['Priority'] == 'P0']['Gene'].tolist()
p1_genes = val_df[val_df['Priority'] == 'P1']['Gene'].tolist()

n_p0 = len(p0_genes)
n_p1 = len(p1_genes)
n_total = n_p0 + n_p1

fig, axes = plt.subplots(n_total, 1, figsize=(14, n_total * 2.2))

for i, (gene, pri) in enumerate([(g, 'P0') for g in p0_genes] + [(g, 'P1') for g in p1_genes]):
    ax = axes[i]
    gdf = liver_ind[liver_ind['gene_name'] == gene]

    # Per-stage expression
    stages = [15, 45, 75, 105]
    dly_means = []
    tfb_means = []
    dly_sems = []
    tfb_sems = []
    for s in stages:
        dly_v = gdf[(gdf['breed'] == 'DLY') & (gdf['stage_kg'] == s)]['expr']
        tfb_v = gdf[(gdf['breed'] == 'TFB') & (gdf['stage_kg'] == s)]['expr']
        dly_means.append(dly_v.mean() if len(dly_v) > 0 else np.nan)
        tfb_means.append(tfb_v.mean() if len(tfb_v) > 0 else np.nan)
        dly_sems.append(dly_v.std() / np.sqrt(len(dly_v)) if len(dly_v) > 0 else 0)
        tfb_sems.append(tfb_v.std() / np.sqrt(len(tfb_v)) if len(tfb_v) > 0 else 0)

    x = np.arange(len(stages))
    w = 0.35
    bars1 = ax.bar(x - w/2, dly_means, w, yerr=dly_sems, color='#2196F3', edgecolor='black',
                   linewidth=0.8, capsize=3, label='DLY')
    bars2 = ax.bar(x + w/2, tfb_means, w, yerr=tfb_sems, color='#C62828', edgecolor='black',
                   linewidth=0.8, capsize=3, label='TFB')

    # Add fold-change labels
    for j, s in enumerate(stages):
        if pd.notna(dly_means[j]) and pd.notna(tfb_means[j]) and tfb_means[j] > 0:
            fc = np.log2(dly_means[j] / tfb_means[j])
            color = '#E91E63' if fc > 0 else '#2196F3'
            ax.text(j, max(dly_means[j], tfb_means[j]) + max(dly_sems[j], tfb_sems[j]) * 1.5,
                   f'{fc:+.1f}', ha='center', fontsize=8, color=color, fontweight='bold')

    # Tier info
    tier_row = aa_refined_df[aa_refined_df['Gene'] == gene]
    tier = int(tier_row['Tier_no105'].values[0]) if len(tier_row) > 0 else 99
    pattern = tier_row['Pattern'].values[0] if len(tier_row) > 0 else ''

    # Color code the background
    tier_bg = {1: '#E8F5E9', 2: '#FFF3E0', 3: '#F5F5F5'}.get(tier, 'white')
    ax.set_facecolor(tier_bg)

    ax.set_ylabel(gene, fontsize=11, fontweight='bold', rotation=0, labelpad=40)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{s} kg' for s in stages])
    ax.set_title(f'[{pri}] Tier {tier}: {pattern}', fontsize=9, loc='right', color='#616161', fontstyle='italic')
    if i == 0:
        ax.legend(fontsize=9, loc='upper left')
    ax.grid(axis='y', alpha=0.2)

fig.suptitle('Prioritized Validation Genes: Expression Across Stages\n'
             'P0 = Genetic Programming Drivers | P1 = 45kg Switch Candidates',
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig('fig_priority_validation_genes.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_priority_validation_genes.pdf', bbox_inches='tight')
plt.close()
print("Saved fig_priority_validation_genes.png/pdf")

# ============================================================
# 7. Save results
# ============================================================
with pd.ExcelWriter('refined_tier_results.xlsx', engine='openpyxl') as writer:
    aa_refined_df.to_excel(writer, sheet_name='AA_Enzymes_Refined', index=False)
    ct_refined_df.to_excel(writer, sheet_name='Crosstalk_Refined', index=False)
    stat3_corr_df.to_excel(writer, sheet_name='STAT3_vs_AA_Enzymes', index=False)
    val_df.to_excel(writer, sheet_name='Validation_Priority_List', index=False)

print("\nSaved refined_tier_results.xlsx")
print("Sheets: AA_Enzymes_Refined | Crosstalk_Refined | STAT3_vs_AA_Enzymes | Validation_Priority_List")
print("\nDone!")
