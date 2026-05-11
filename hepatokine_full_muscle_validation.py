"""
Hepatokine/Myokine full muscle validation — use complete 24784-gene muscle matrix.
Validates liver→muscle signaling by:
1. Computing muscle FC for all hepatokines/myokines
2. Searching for known receptors in muscle
3. Computing liver-muscle cross-tissue correlations
4. Updated figures with real muscle data
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

out = '/Users/hezongze/pig_study'

# ============================================================
# 1. Load data
# ============================================================
print("Loading data...")

# Liver: corrected matrix
liver_corr = pd.read_csv(f'{out}/liver_expression_corrected.csv', index_col=0)
liver_raw = pd.read_csv(f'{out}/gene_expression/liver_gene_matrix.xls', sep='\t', index_col=0)
liver_cols_orig = [c for c in liver_raw.columns if c.startswith('L_')]
liver_mat = liver_corr.values  # already log2, filtered, corrected

# Full muscle matrix
muscle_full_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t', index_col=0)
# Expression columns: those starting with BJ_ or m_ (not gene_name, length, description)
muscle_expr_cols = [c for c in muscle_full_raw.columns if (c.startswith('BJ_') or c.startswith('m_'))]
muscle_expr_full = muscle_full_raw[muscle_expr_cols].apply(pd.to_numeric, errors='coerce')
muscle_log_full = np.log2(muscle_expr_full.values + 1)

# Filter low expression
muscle_keep = muscle_log_full.mean(axis=1) > 1
muscle_log_full = muscle_log_full[muscle_keep, :]
muscle_genes_full = muscle_full_raw.index[muscle_keep]
muscle_gene_name_full = muscle_full_raw['gene_name'].iloc[muscle_keep] if 'gene_name' in muscle_full_raw.columns else pd.Series(['']*muscle_keep.sum(), index=muscle_genes_full)
muscle_gene_desc_full = muscle_full_raw['description'].iloc[muscle_keep] if 'description' in muscle_full_raw.columns else pd.Series(['']*muscle_keep.sum(), index=muscle_genes_full)

print(f"  Liver: {liver_mat.shape[0]} genes x {liver_mat.shape[1]} samples")
print(f"  Muscle (full): {muscle_log_full.shape[0]} genes x {muscle_log_full.shape[1]} samples")

# ============================================================
# 2. Build stage-level FC matrices
# ============================================================
stages = ['15kg', '45kg', '75kg', '105kg']

# Liver stage specs: (DLY_prefix, TFB_prefix)
liver_stage_specs = [
    ('L_15_1', 'L_15_2'),
    ('L_45_1', 'L_45_2'),
    ('L_1_1', 'L_1_2'),
    ('L_2_1', 'L_2_2'),
]

# Muscle stage specs (full matrix has _X suffix for replicates)
muscle_stage_specs = [
    ('m_15_1_', 'm_15_2_'),
    ('BJ_2_1_', 'BJ_2_2_'),
    ('m_1_1_', 'm_1_2_'),
    ('m_2_1_', 'm_2_2_'),
]

# Compute liver FC per stage
liver_fc = {}
liver_expr_data = {}  # gene_symbol -> {stage: {DLY: mean, TFB: mean}}
for s, (dp, tp) in zip(stages, liver_stage_specs):
    di = [i for i, c in enumerate(liver_cols_orig) if c.startswith(dp)]
    ti = [i for i, c in enumerate(liver_cols_orig) if c.startswith(tp)]
    liver_fc[s] = liver_mat[:, di].mean(1) - liver_mat[:, ti].mean(1)

# Compute muscle FC per stage
muscle_fc = {}
for s, (dp, tp) in zip(stages, muscle_stage_specs):
    di = [i for i, c in enumerate(muscle_expr_cols) if c.startswith(dp)]
    ti = [i for i, c in enumerate(muscle_expr_cols) if c.startswith(tp)]
    muscle_fc[s] = muscle_log_full[:, di].mean(1) - muscle_log_full[:, ti].mean(1)
    print(f"  {s} muscle: {len(di)} DLY, {len(ti)} TFB samples")

# Build gene symbol maps
liver_gene_name = liver_raw['gene_name'] if 'gene_name' in liver_raw.columns else pd.Series(['']*len(liver_raw), index=liver_raw.index)
liver_gene_name_filt = pd.Series([str(liver_gene_name.get(gid, '')) for gid in liver_corr.index], index=liver_corr.index)

liver_symbols = {}
for i, gid in enumerate(liver_corr.index):
    sym = str(liver_gene_name_filt.iloc[i]).strip()
    if sym and sym != 'nan' and sym != '':
        liver_symbols[sym.upper()] = {
            'index': i,
            'gene_id': gid,
            'fcs': np.array([liver_fc[s][i] for s in stages]),
            'dly_means': {s: liver_mat[i, [j for j,c in enumerate(liver_cols_orig) if c.startswith(dp)]].mean()
                          for s, (dp, tp) in zip(stages, liver_stage_specs)},
            'tfb_means': {s: liver_mat[i, [j for j,c in enumerate(liver_cols_orig) if c.startswith(tp)]].mean()
                          for s, (dp, tp) in zip(stages, liver_stage_specs)},
        }

muscle_symbols = {}
for i, gid in enumerate(muscle_genes_full):
    sym = str(muscle_gene_name_full.iloc[i]).strip() if hasattr(muscle_gene_name_full, 'iloc') else ''
    if not sym or sym == 'nan':
        continue
    sym = sym.upper()
    muscle_symbols[sym] = {
        'index': i,
        'gene_id': gid,
        'fcs': np.array([muscle_fc[s][i] for s in stages]),
        'dly_means': {s: muscle_log_full[i, [j for j,c in enumerate(muscle_expr_cols) if c.startswith(dp)]].mean()
                      for s, (dp, tp) in zip(stages, muscle_stage_specs)},
        'tfb_means': {s: muscle_log_full[i, [j for j,c in enumerate(muscle_expr_cols) if c.startswith(tp)]].mean()
                      for s, (dp, tp) in zip(stages, muscle_stage_specs)},
    }

common_symbols = set(liver_symbols.keys()) & set(muscle_symbols.keys())
print(f"  Liver genes with symbols: {len(liver_symbols)}")
print(f"  Muscle genes with symbols: {len(muscle_symbols)}")
print(f"  Common symbols (both tissues): {len(common_symbols)}")

# ============================================================
# 3. Define hepatokines, myokines, and receptors
# ============================================================
hepatokines = [
    'FGF21', 'IGF1', 'IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFALS',
    'FST', 'ANGPTL4', 'ANGPTL6', 'ANGPTL8',
    'Fetuin-A', 'AHSG',
    'RBP4', 'APOB', 'APOA1', 'APOA2', 'APOC3', 'APOE',
    'TTR', 'ALB', 'TF', 'HP', 'SERPINA1', 'SERPINA3',
    'ORM1', 'LECT2', 'FGA', 'FGB', 'FGG',
    'HAMP', 'ERFE', 'GDF15', 'INHBE', 'MSTN',
    'CCN1', 'CCN2', 'CCN3', 'CCN4', 'CCN5',
    'SPARC', 'FSTL1', 'FSTL3',
]

myokines = [
    'MSTN', 'BDNF', 'FNDC5', 'IL6', 'IL15', 'MIF',
    'FGF2', 'FSTL1', 'SPARC', 'MUSK',
    'DCN', 'LUM', 'POSTN', 'CTGF',
]

# Hepatokine receptors (to search in muscle)
receptors = {
    'IGF1': ['IGF1R', 'INSR'],
    'IGFBP1': ['IGF1R', 'ITGB1'],
    'IGFBP2': ['IGF1R', 'ITGB1'],
    'IGFBP3': ['IGF1R', 'TMEM219'],
    'FGF21': ['FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'KLB'],
    'FST': ['ACVR2A', 'ACVR2B', 'ACVR1B'],
    'MSTN': ['ACVR2B', 'ACVR1B', 'TGFBR1'],
    'BDNF': ['NTRK2', 'NGFR'],
    'FNDC5': ['ITGA5', 'ITGB1'],
    'IL6': ['IL6R', 'IL6ST'],
    'IL15': ['IL2RA', 'IL2RB', 'IL2RG', 'IL15RA'],
    'ANGPTL4': ['ITGB1', 'ITGA5'],
    'ANGPTL6': ['ITGB1'],
    'ANGPTL8': ['ITGB1'],
    'GDF15': ['GFRAL', 'RET'],
    'CCN1': ['ITGA6', 'ITGB1', 'ITGAV', 'ITGB3'],
    'CCN2': ['ITGAV', 'ITGB1', 'LRP1'],
    'FSTL1': ['TLR4', 'DIP2A'],
    'FSTL3': ['ACVR2B'],
    'APOE': ['LRP1', 'LDLR'],
    'RBP4': ['STRA6'],
    'SPARC': ['ITGB1', 'ITGAV'],
    'ERFE': ['BMPR1A', 'BMPR1B'],
    'HAMP': ['SLC40A1'],
}

# Read tier data
tier_df = pd.read_excel(f'{out}/integrated_liver_4stage.xlsx')
tier_map = {}
for _, r in tier_df.iterrows():
    sym = str(r['Gene_Symbol']).upper()
    if sym and sym != 'nan':
        tier_map[sym] = r['Tier']

# ============================================================
# 4. Build comprehensive hepatokine table
# ============================================================
print("\nBuilding hepatokine signaling table with full muscle data...")

hk_rows = []
for gene in sorted(set(hepatokines + myokines)):
    gu = gene.upper()
    if gu not in liver_symbols:
        continue

    role = 'Hepatokine' if gene in hepatokines else 'Myokine'
    if gene in myokines and gene in hepatokines:
        role = 'Both'

    ld = liver_symbols[gu]
    row = {
        'Gene': gene,
        'Role': role,
        'Liver_Gene_ID': ld['gene_id'],
        'Liver_Tier': tier_map.get(gu, 'Unknown'),
    }

    # Liver FC per stage
    for si, s in enumerate(stages):
        row[f'Liver_{s}_log2FC'] = ld['fcs'][si]
        row[f'Liver_{s}_DLY_mean'] = ld['dly_means'][s]
        row[f'Liver_{s}_TFB_mean'] = ld['tfb_means'][s]

    # Liver summary
    row['Liver_MeanAbsFC'] = np.mean(np.abs(ld['fcs']))
    row['Liver_MaxAbsFC'] = np.max(np.abs(ld['fcs']))
    row['Liver_PeakStage'] = stages[np.argmax(np.abs(ld['fcs']))]

    # Muscle data if available
    if gu in muscle_symbols:
        md = muscle_symbols[gu]
        row['Muscle_Gene_ID'] = md['gene_id']
        for si, s in enumerate(stages):
            row[f'Muscle_{s}_log2FC'] = md['fcs'][si]
            row[f'Muscle_{s}_DLY_mean'] = md['dly_means'][s]
            row[f'Muscle_{s}_TFB_mean'] = md['tfb_means'][s]

        row['Muscle_MeanAbsFC'] = np.mean(np.abs(md['fcs']))
        row['Muscle_MaxAbsFC'] = np.max(np.abs(md['fcs']))

        # Liver-muscle FC correlation (Pearson)
        lfcs = ld['fcs']
        mfcs = md['fcs']
        if np.std(lfcs) > 0.05 and np.std(mfcs) > 0.05:
            row['LiverMuscle_FC_pearson'], _ = pearsonr(lfcs, mfcs)
            row['LiverMuscle_FC_spearman'], _ = spearmanr(lfcs, mfcs)
        else:
            row['LiverMuscle_FC_pearson'] = 0
            row['LiverMuscle_FC_spearman'] = 0

        # Direction agreement
        l_dir = np.sign(lfcs)
        m_dir = np.sign(mfcs)
        meaningful = sum(1 for i in range(4) if abs(lfcs[i]) > 0.2 and abs(mfcs[i]) > 0.2)
        same = sum(1 for i in range(4) if abs(lfcs[i]) > 0.2 and abs(mfcs[i]) > 0.2 and l_dir[i] == m_dir[i])
        row['Direction_Agreement'] = same / max(meaningful, 1) if meaningful > 0 else np.nan

        # 105kg-specific: liver-muscle FC difference
        row['105kg_LiverMuscle_Diff'] = abs(lfcs[3] - mfcs[3])
    else:
        for si, s in enumerate(stages):
            row[f'Muscle_{s}_log2FC'] = np.nan
            row[f'Muscle_{s}_DLY_mean'] = np.nan
            row[f'Muscle_{s}_TFB_mean'] = np.nan
        row['Muscle_MeanAbsFC'] = np.nan
        row['Muscle_MaxAbsFC'] = np.nan
        row['LiverMuscle_FC_pearson'] = np.nan
        row['LiverMuscle_FC_spearman'] = np.nan
        row['Direction_Agreement'] = np.nan
        row['105kg_LiverMuscle_Diff'] = np.nan

    hk_rows.append(row)

hk_df = pd.DataFrame(hk_rows)

# Classify signaling pattern
def classify_signaling(row):
    lfcs = np.array([row[f'Liver_{s}_log2FC'] for s in stages])
    mfcs = np.array([row[f'Muscle_{s}_log2FC'] for s in stages])

    l_max = np.max(np.abs(lfcs))
    m_max = np.max(np.abs(mfcs)) if not np.isnan(mfcs[0]) else 0

    if l_max <= 0.5 and (np.isnan(mfcs[0]) or m_max <= 0.5):
        return 'Low_Signal'

    if l_max > 0.5 and (np.isnan(mfcs[0]) or m_max <= 0.5):
        direction = 'DLY_UP' if np.mean(lfcs) > 0 else 'TFB_UP'
        return f'Liver_Only_{direction}'

    if l_max <= 0.5 and m_max > 0.5:
        direction = 'DLY_UP' if np.mean(mfcs) > 0 else 'TFB_UP'
        return f'Muscle_Only_{direction}'

    # Both have signal
    corr = row.get('LiverMuscle_FC_pearson', 0)
    agree = row.get('Direction_Agreement', 0)

    if not np.isnan(corr) and corr > 0.7 and isinstance(agree, (int, float)) and not np.isnan(agree) and agree > 0.75:
        return 'Coordinated_Positive'
    elif not np.isnan(corr) and corr < -0.5:
        return 'Opposing'
    elif isinstance(agree, (int, float)) and not np.isnan(agree) and agree < 0.5:
        return 'Discordant'
    else:
        return 'Complex'

hk_df['Signaling_Pattern'] = hk_df.apply(classify_signaling, axis=1)
hk_df = hk_df.sort_values('Liver_MeanAbsFC', ascending=False)

# ============================================================
# 5. Receptor expression in muscle
# ============================================================
print("\nSearching for hepatokine receptors in muscle...")

receptor_rows = []
all_receptor_genes = set()
for ligand, recs in receptors.items():
    for r in recs:
        all_receptor_genes.add(r.upper())

for ligand, rec_list in receptors.items():
    if ligand.upper() not in liver_symbols:
        continue

    for rec_name in rec_list:
        rec_upper = rec_name.upper()

        # Check in muscle
        muscle_info = None
        if rec_upper in muscle_symbols:
            md = muscle_symbols[rec_upper]
            muscle_info = {
                'gene_id': md['gene_id'],
                'fcs': md['fcs'],
                'dly_means': md['dly_means'],
                'tfb_means': md['tfb_means'],
                'mean_expr': np.mean([md['dly_means'][s] for s in stages] + [md['tfb_means'][s] for s in stages]),
            }

        # Check in liver (autocrine?)
        liver_info = None
        if rec_upper in liver_symbols:
            ld = liver_symbols[rec_upper]
            liver_info = {
                'gene_id': ld['gene_id'],
                'fcs': ld['fcs'],
                'mean_expr': np.mean([ld['dly_means'][s] for s in stages] + [ld['tfb_means'][s] for s in stages]),
            }

        row = {
            'Ligand': ligand,
            'Receptor': rec_name,
            'Ligand_Liver_Tier': tier_map.get(ligand.upper(), 'Unknown'),
        }

        # Ligand liver FC
        if ligand.upper() in liver_symbols:
            ld = liver_symbols[ligand.upper()]
            for si, s in enumerate(stages):
                row[f'Ligand_Liver_{s}_log2FC'] = ld['fcs'][si]
        else:
            for si, s in enumerate(stages):
                row[f'Ligand_Liver_{s}_log2FC'] = np.nan

        # Receptor in muscle
        if muscle_info:
            row['Receptor_In_Muscle'] = True
            row['Receptor_Muscle_ID'] = muscle_info['gene_id']
            row['Receptor_Muscle_MeanExpr'] = muscle_info['mean_expr']
            for si, s in enumerate(stages):
                row[f'Receptor_Muscle_{s}_log2FC'] = muscle_info['fcs'][si]
            # Ligand-receptor FC correlation
            if ligand.upper() in liver_symbols:
                ld_fcs = liver_symbols[ligand.upper()]['fcs']
                rec_fcs = muscle_info['fcs']
                if np.std(ld_fcs) > 0.05 and np.std(rec_fcs) > 0.05:
                    row['Ligand_Receptor_FC_corr'], _ = pearsonr(ld_fcs, rec_fcs)
                else:
                    row['Ligand_Receptor_FC_corr'] = 0
            else:
                row['Ligand_Receptor_FC_corr'] = np.nan
        else:
            row['Receptor_In_Muscle'] = False
            row['Receptor_Muscle_ID'] = ''
            row['Receptor_Muscle_MeanExpr'] = np.nan
            for si, s in enumerate(stages):
                row[f'Receptor_Muscle_{s}_log2FC'] = np.nan
            row['Ligand_Receptor_FC_corr'] = np.nan

        # Receptor in liver
        row['Receptor_In_Liver'] = liver_info is not None

        receptor_rows.append(row)

rec_df = pd.DataFrame(receptor_rows)

# Summary
rec_in_muscle = rec_df[rec_df['Receptor_In_Muscle']]
print(f"  Receptor-ligand pairs with receptor found in muscle: {len(rec_in_muscle)} / {len(rec_df)}")
print(f"  Unique receptors found: {len(set(r['Receptor'] for _, r in rec_in_muscle.iterrows()))}")

# ============================================================
# 6. Cross-tissue genome-wide correlation analysis
# ============================================================
print("\nGenome-wide cross-tissue correlation analysis...")

# For each stage, compute liver FC vs muscle FC for all common genes
stage_cross_corr = {}
for si, s in enumerate(stages):
    common_lfcs = []
    common_mfcs = []
    for sym in common_symbols:
        common_lfcs.append(liver_symbols[sym]['fcs'][si])
        common_mfcs.append(muscle_symbols[sym]['fcs'][si])
    common_lfcs = np.array(common_lfcs)
    common_mfcs = np.array(common_mfcs)

    # Remove zero-variance
    if np.std(common_lfcs) > 0 and np.std(common_mfcs) > 0:
        r, p = pearsonr(common_lfcs, common_mfcs)
        r_s, p_s = spearmanr(common_lfcs, common_mfcs)
    else:
        r, p, r_s, p_s = 0, 1, 0, 1

    stage_cross_corr[s] = {
        'n_genes': len(common_lfcs),
        'pearson_r': r,
        'pearson_p': p,
        'spearman_r': r_s,
        'spearman_p': p_s,
    }
    print(f"  {s}: n={len(common_lfcs)}, Pearson r={r:.4f} (p={p:.2e}), Spearman r={r_s:.4f} (p={p_s:.2e})")

# ============================================================
# 7. Identify top coordinated genes (high liver-muscle FC correlation)
# ============================================================
print("\nFinding top coordinated liver-muscle genes...")

gene_corrs = []
for sym in common_symbols:
    lfcs = liver_symbols[sym]['fcs']
    mfcs = muscle_symbols[sym]['fcs']

    if np.std(lfcs) > 0.05 and np.std(mfcs) > 0.05:
        r, p = pearsonr(lfcs, mfcs)
    else:
        r, p = 0, 1

    l_mean_abs = np.mean(np.abs(lfcs))
    m_mean_abs = np.mean(np.abs(mfcs))

    gene_corrs.append({
        'Gene': sym,
        'Liver_MeanAbsFC': l_mean_abs,
        'Muscle_MeanAbsFC': m_mean_abs,
        'CrossTissue_pearson_r': r,
        'CrossTissue_pearson_p': p,
        'Signal_Strength': l_mean_abs * m_mean_abs,
    })

gene_corr_df = pd.DataFrame(gene_corrs)
gene_corr_df = gene_corr_df.sort_values('CrossTissue_pearson_r', ascending=False)

# Filter to hepatokines/myokines
hk_upper = set(g.upper() for g in hepatokines + myokines)
hk_corr = gene_corr_df[gene_corr_df['Gene'].isin(hk_upper)]
print(f"  Hepatokines/myokines in cross-tissue correlation: {len(hk_corr)}")
print(f"  Top coordinated hepatokines:")
for _, r in hk_corr.sort_values('CrossTissue_pearson_r', ascending=False).head(10).iterrows():
    print(f"    {r['Gene']:<15} r={r['CrossTissue_pearson_r']:.3f}  Liver|FC|={r['Liver_MeanAbsFC']:.2f}  Muscle|FC|={r['Muscle_MeanAbsFC']:.2f}")

# ============================================================
# 8. FIGURE 1: Updated hepatokine crosstalk heatmap (with real muscle data)
# ============================================================
print("\nGenerating updated figures...")

cmap = sns.diverging_palette(250, 10, as_cmap=True)

# Select hepatokines with meaningful signal
hk_sig = hk_df[(hk_df['Liver_MeanAbsFC'] > 0.5) |
               (hk_df['Muscle_MeanAbsFC'].notna() & (hk_df['Muscle_MeanAbsFC'] > 0.5))]
hk_sig = hk_sig.sort_values('Liver_MeanAbsFC', ascending=False)

if len(hk_sig) > 0:
    n_hk = len(hk_sig)
    fig_h = max(7, n_hk * 0.35)
    fig, (ax_l, ax_m) = plt.subplots(1, 2, figsize=(14, fig_h),
                                      gridspec_kw={'width_ratios': [1, 1]})

    # Liver heatmap
    l_mat = np.array([[r[f'Liver_{s}_log2FC'] for s in stages] for _, r in hk_sig.iterrows()])
    vmax_hk = max(abs(l_mat[~np.isnan(l_mat)]).max(), 3) if l_mat.size > 0 else 3

    im_l = ax_l.imshow(l_mat, aspect='auto', cmap=cmap, vmin=-vmax_hk, vmax=vmax_hk)
    ax_l.set_xticks(range(4))
    ax_l.set_xticklabels(stages, fontsize=10)
    ax_l.set_yticks(range(n_hk))
    ax_l.set_yticklabels(hk_sig['Gene'].values, fontsize=9)
    ax_l.set_title('Liver: Hepatokine/Myokine log2FC\n(DLY vs TFB)', fontsize=11, fontweight='bold')

    role_colors = {'Hepatokine': '#d62728', 'Myokine': '#1f77b4', 'Both': '#9467bd'}
    for i, role in enumerate(hk_sig['Role'].values):
        ax_l.get_yticklabels()[i].set_color(role_colors.get(role, '#000000'))

    plt.colorbar(im_l, ax=ax_l, shrink=0.8)

    # Muscle heatmap
    m_mat = np.array([[r[f'Muscle_{s}_log2FC'] if not np.isnan(r[f'Muscle_{s}_log2FC']) else 0
                       for s in stages] for _, r in hk_sig.iterrows()])

    im_m = ax_m.imshow(m_mat, aspect='auto', cmap=cmap, vmin=-vmax_hk, vmax=vmax_hk)
    ax_m.set_xticks(range(4))
    ax_m.set_xticklabels(stages, fontsize=10)
    ax_m.set_yticks(range(n_hk))
    ax_m.set_yticklabels(hk_sig['Gene'].values, fontsize=9)
    ax_m.set_title('Muscle: Corresponding log2FC\n(Full 24784-gene matrix)', fontsize=11, fontweight='bold')

    for i, role in enumerate(hk_sig['Role'].values):
        ax_m.get_yticklabels()[i].set_color(role_colors.get(role, '#000000'))

    plt.colorbar(im_m, ax=ax_m, shrink=0.8)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=f'{r}') for r, c in role_colors.items()]
    fig.legend(handles=legend_patches, loc='lower center', ncol=3, fontsize=9)

    fig.suptitle('Hepatokine/Myokine Signaling: Liver → Muscle Crosstalk (Full Muscle Validation)',
                 fontsize=13, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(f'{out}/fig_hepatokine_crosstalk_full_muscle.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{out}/fig_hepatokine_crosstalk_full_muscle.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: fig_hepatokine_crosstalk_full_muscle.png/pdf")

# ============================================================
# 9. FIGURE 2: Receptor expression in muscle (heatmap)
# ============================================================
print("Generating receptor expression figure...")

rec_found = rec_df[rec_df['Receptor_In_Muscle']].drop_duplicates(subset=['Ligand', 'Receptor'])
# Pivot to get ligand x receptor matrix of correlations
rec_pivot_data = []
for _, r in rec_found.iterrows():
    rec_pivot_data.append({
        'Ligand': r['Ligand'],
        'Receptor': r['Receptor'],
        'Corr': r.get('Ligand_Receptor_FC_corr', np.nan),
        'Receptor_Muscle_MeanExpr': r.get('Receptor_Muscle_MeanExpr', np.nan),
    })

if rec_pivot_data:
    rec_pivot_df = pd.DataFrame(rec_pivot_data)
    rec_matrix = rec_pivot_df.pivot_table(values='Corr', index='Ligand', columns='Receptor', aggfunc='first')

    if rec_matrix.shape[0] > 0 and rec_matrix.shape[1] > 0:
        fig, ax = plt.subplots(figsize=(max(8, rec_matrix.shape[1]*0.8),
                                        max(5, rec_matrix.shape[0]*0.4)))

        mask = rec_matrix.isna()
        sns.heatmap(rec_matrix, annot=True, fmt='.2f', cmap='RdBu_r',
                    center=0, vmin=-1, vmax=1,
                    mask=mask, ax=ax,
                    cbar_kws={'label': 'Ligand-Receptor FC Correlation'},
                    linewidths=0.5, linecolor='#eeeeee')

        ax.set_title('Hepatokine → Muscle Receptor: Cross-Tissue FC Correlation\n(Liver ligand FC vs Muscle receptor FC across 4 stages)',
                     fontsize=12, fontweight='bold')
        ax.set_xlabel('Receptor (muscle expression)', fontsize=10)
        ax.set_ylabel('Ligand (liver expression)', fontsize=10)
        plt.tight_layout()
        plt.savefig(f'{out}/fig_hepatokine_receptor_correlation.png', dpi=200, bbox_inches='tight')
        plt.savefig(f'{out}/fig_hepatokine_receptor_correlation.pdf', bbox_inches='tight')
        plt.close()
        print(f"  Saved: fig_hepatokine_receptor_correlation.png/pdf")

# ============================================================
# 10. FIGURE 3: Cross-tissue correlation scatter per stage
# ============================================================
print("Generating cross-tissue scatter plots...")

fig, axes = plt.subplots(2, 2, figsize=(12, 11))
axes = axes.flatten()

for si, (s, ax) in enumerate(zip(stages, axes)):
    common_lfcs = np.array([liver_symbols[sym]['fcs'][si] for sym in common_symbols])
    common_mfcs = np.array([muscle_symbols[sym]['fcs'][si] for sym in common_symbols])

    # Filter extreme outliers for visualization
    l_q = np.percentile(common_lfcs, [1, 99])
    m_q = np.percentile(common_mfcs, [1, 99])
    keep = (common_lfcs >= l_q[0]) & (common_lfcs <= l_q[1]) & (common_mfcs >= m_q[0]) & (common_mfcs <= m_q[1])

    l_plot = common_lfcs[keep]
    m_plot = common_mfcs[keep]

    ax.scatter(l_plot, m_plot, s=1, alpha=0.15, color='#555555', rasterized=True)

    # Highlight hepatokines
    for sym in hk_upper:
        if sym in liver_symbols and sym in muscle_symbols:
            lfc = liver_symbols[sym]['fcs'][si]
            mfc = muscle_symbols[sym]['fcs'][si]
            if abs(lfc) > 0.5 or abs(mfc) > 0.5:
                ax.scatter([lfc], [mfc], s=40, alpha=0.9, color='#d62728', edgecolors='#333333', linewidth=0.5, zorder=5)
                ax.annotate(sym, (lfc, mfc), fontsize=6, xytext=(3, 3),
                           textcoords='offset points', alpha=0.8)

    r = stage_cross_corr[s]['pearson_r']
    r_s = stage_cross_corr[s]['spearman_r']
    ax.axhline(y=0, color='grey', linewidth=0.5, linestyle='--')
    ax.axvline(x=0, color='grey', linewidth=0.5, linestyle='--')
    ax.set_xlabel('Liver log2FC (DLY/TFB)')
    ax.set_ylabel('Muscle log2FC (DLY/TFB)')
    ax.set_title(f'{s} (n={len(common_symbols)}, r={r:.3f}, r_s={r_s:.3f})')
    ax.grid(True, alpha=0.2)

fig.suptitle('Liver-Muscle Cross-Tissue FC Correlation by Stage\n(Full muscle matrix, hepatokines highlighted in red)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{out}/fig_crosstissue_scatter_4stage.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{out}/fig_crosstissue_scatter_4stage.pdf', bbox_inches='tight')
plt.close()
print(f"  Saved: fig_crosstissue_scatter_4stage.png/pdf")

# ============================================================
# 11. FIGURE 4: Key hepatokine trajectories (liver + muscle)
# ============================================================
print("Generating hepatokine trajectory plots...")

# Select top hepatokines with both liver and muscle data
hk_both = hk_df[(hk_df['Liver_MeanAbsFC'] > 0.5) &
                hk_df['Muscle_MeanAbsFC'].notna()].sort_values('Liver_MeanAbsFC', ascending=False).head(12)

if len(hk_both) > 0:
    n_cols = 4
    n_rows = int(np.ceil(len(hk_both) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 4 * n_rows))
    axes = axes.flatten()

    x_pos = [15, 45, 75, 105]

    for idx, (_, row) in enumerate(hk_both.iterrows()):
        ax = axes[idx]
        gene = row['Gene']

        # Liver
        lfcs = [row[f'Liver_{s}_log2FC'] for s in stages]
        ax.plot(x_pos, lfcs, 'o-', color='#d62728', linewidth=2, markersize=6, label='Liver FC')

        # Muscle
        mfcs = [row[f'Muscle_{s}_log2FC'] for s in stages]
        ax.plot(x_pos, mfcs, 's--', color='#1f77b4', linewidth=2, markersize=6, label='Muscle FC')

        ax.axhline(y=0, color='grey', linewidth=0.5, linestyle='-')
        ax.set_title(f"{gene}\nr={row['LiverMuscle_FC_pearson']:.2f}",
                    fontsize=10, fontweight='bold')
        ax.set_xlabel('Weight (kg)')
        ax.set_ylabel('log2FC (DLY/TFB)')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    # Hide unused axes
    for idx in range(len(hk_both), len(axes)):
        axes[idx].set_visible(False)

    fig.suptitle('Key Hepatokine/Myokine FC Trajectories: Liver vs Muscle',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    plt.savefig(f'{out}/fig_hepatokine_trajectories_full_muscle.png', dpi=200, bbox_inches='tight')
    plt.savefig(f'{out}/fig_hepatokine_trajectories_full_muscle.pdf', bbox_inches='tight')
    plt.close()
    print(f"  Saved: fig_hepatokine_trajectories_full_muscle.png/pdf")

# ============================================================
# 12. Print key findings
# ============================================================
print(f"\n{'='*70}")
print("HEPATOKINE/MYOKINE FULL MUSCLE VALIDATION RESULTS")
print(f"{'='*70}")

print(f"\nGenome-wide cross-tissue correlations:")
for s in stages:
    sc = stage_cross_corr[s]
    print(f"  {s}: n={sc['n_genes']}, r={sc['pearson_r']:.4f} (p={sc['pearson_p']:.2e})")

print(f"\nHepatokine signaling classification ({len(hk_df)} genes):")
for pat, cnt in hk_df['Signaling_Pattern'].value_counts().items():
    print(f"  {pat}: {cnt}")

# Top hepatokines with muscle data
print(f"\nTop hepatokines with validated muscle correlation:")
hk_with_muscle = hk_df[hk_df['Muscle_MeanAbsFC'].notna()].sort_values('Liver_MeanAbsFC', ascending=False)
for _, r in hk_with_muscle.head(15).iterrows():
    lf = ' '.join([f"{r[f'Liver_{s}_log2FC']:>6.2f}" for s in stages])
    mf = ' '.join([f"{r[f'Muscle_{s}_log2FC']:>6.2f}" for s in stages])
    r_val = r['LiverMuscle_FC_pearson']
    print(f"  {r['Gene']:<12} L:[{lf}]  M:[{mf}]  r={r_val:.2f}  {r['Signaling_Pattern']}")

# Receptor findings
print(f"\nReceptor expression in muscle:")
for _, r in rec_found.sort_values('Ligand_Receptor_FC_corr', ascending=False).head(20).iterrows():
    lig = r['Ligand']
    rec = r['Receptor']
    corr = r.get('Ligand_Receptor_FC_corr', np.nan)
    expr = r.get('Receptor_Muscle_MeanExpr', np.nan)
    print(f"  {lig} → {rec}: corr={corr:.3f}" if not np.isnan(corr) else f"  {lig} → {rec}: no data")

# ============================================================
# 13. Save outputs
# ============================================================
print(f"\nSaving output files...")

# Updated hepatokine table with full muscle data
hk_df.to_excel(f'{out}/hepatokine_full_muscle_validation.xlsx', index=False)
print(f"  hepatokine_full_muscle_validation.xlsx — {len(hk_df)} genes")

# Receptor analysis
rec_df.to_excel(f'{out}/hepatokine_receptor_muscle_analysis.xlsx', index=False)
print(f"  hepatokine_receptor_muscle_analysis.xlsx — {len(rec_df)} receptor-ligand pairs")

# Cross-tissue correlation data
gene_corr_df.to_excel(f'{out}/crosstissue_correlation_genomewide.xlsx', index=False)
print(f"  crosstissue_correlation_genomewide.xlsx — {len(gene_corr_df)} genes")

# Stage correlation summary
pd.DataFrame(stage_cross_corr).T.to_excel(f'{out}/crosstissue_stage_correlation_summary.xlsx')
print(f"  crosstissue_stage_correlation_summary.xlsx")

# Master workbook
with pd.ExcelWriter(f'{out}/hepatokine_full_muscle_master.xlsx', engine='openpyxl') as writer:
    hk_df.to_excel(writer, sheet_name='Hepatokine_Signaling', index=False)
    rec_df.to_excel(writer, sheet_name='Receptor_Analysis', index=False)
    gene_corr_df.to_excel(writer, sheet_name='CrossTissue_Correlation', index=False)
    pd.DataFrame(stage_cross_corr).T.to_excel(writer, sheet_name='Stage_Correlation_Summary')

    # Add filtered coordinated hepatokines
    hk_corr_sorted = hk_corr.sort_values('CrossTissue_pearson_r', ascending=False)
    hk_corr_sorted.to_excel(writer, sheet_name='Hepatokine_Coordination', index=False)

print(f"  hepatokine_full_muscle_master.xlsx — Master workbook (5 sheets)")

print("\nDone. Full muscle validation complete.")
