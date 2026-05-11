"""
4-Stage Integrated Liver-Muscle Axis Analysis
Using corrected DLY105 liver data.
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, spearmanr, pearsonr
from stats_utils import benjamini_hochberg
import warnings
warnings.filterwarnings('ignore')

out = '/Users/hezongze/pig_study'

# ============================================================
# 1. Load corrected data
# ============================================================
print("Loading data...")

# Liver: corrected matrix
liver_corr_df = pd.read_csv(f'{out}/liver_expression_corrected.csv', index_col=0)

# Muscle: raw matrix
muscle_raw = pd.read_csv(f'{out}/gene_expression/muscle_gene_matrix.xls', sep='\t', index_col=0)
muscle_expr_cols = [c for c in muscle_raw.columns if c.startswith('m_') or c.startswith('BJ_')]
muscle_expr = muscle_raw[muscle_expr_cols].apply(pd.to_numeric, errors='coerce')
muscle_log = np.log2(muscle_expr.values + 1)
# Filter low expression
muscle_keep = muscle_log.mean(axis=1) > 1
muscle_log = muscle_log[muscle_keep, :]
muscle_genes = muscle_raw.index[muscle_keep]

# Original liver for reference (to extract gene annotations)
liver_raw = pd.read_csv(f'{out}/gene_expression/liver_gene_matrix.xls', sep='\t', index_col=0)
liver_expr_cols_orig = [c for c in liver_raw.columns if c.startswith('L_')]
liver_gene_name = liver_raw['gene_name'] if 'gene_name' in liver_raw.columns else pd.Series(['']*len(liver_raw), index=liver_raw.index)
liver_gene_desc = liver_raw['description'] if 'description' in liver_raw.columns else pd.Series(['']*len(liver_raw), index=liver_raw.index)

# Map corrected gene IDs to annotations
corr_gene_ids = liver_corr_df.index.tolist()
liver_gene_name_filt = pd.Series([str(liver_gene_name.get(gid, '')) for gid in corr_gene_ids], index=corr_gene_ids)
liver_gene_desc_filt = pd.Series([str(liver_gene_desc.get(gid, '')) for gid in corr_gene_ids], index=corr_gene_ids)

muscle_gene_name_raw = muscle_raw['gene_name'] if 'gene_name' in muscle_raw.columns else pd.Series(['']*len(muscle_raw), index=muscle_raw.index)
muscle_gene_desc_raw = muscle_raw['description'] if 'description' in muscle_raw.columns else pd.Series(['']*len(muscle_raw), index=muscle_raw.index)
muscle_gene_name_filt = pd.Series([str(muscle_gene_name_raw.get(gid, '')) for gid in muscle_genes], index=muscle_genes)
muscle_gene_desc_filt = pd.Series([str(muscle_gene_desc_raw.get(gid, '')) for gid in muscle_genes], index=muscle_genes)

# Expression values
liver_mat = liver_corr_df.values  # Already log2, filtered, corrected
muscle_mat = muscle_log

# Column indices
liver_cols = liver_expr_cols_orig
muscle_cols = muscle_expr_cols

print(f"  Liver: {liver_mat.shape[0]} genes x {liver_mat.shape[1]} samples")
print(f"  Muscle: {muscle_mat.shape[0]} genes x {muscle_mat.shape[1]} samples")

# ============================================================
# 2. Group definitions
# ============================================================
liver_stages = {
    '15kg': ([i for i,c in enumerate(liver_cols) if c.startswith('L_15_1')],
             [i for i,c in enumerate(liver_cols) if c.startswith('L_15_2')]),
    '45kg': ([i for i,c in enumerate(liver_cols) if c.startswith('L_45_1')],
             [i for i,c in enumerate(liver_cols) if c.startswith('L_45_2')]),
    '75kg': ([i for i,c in enumerate(liver_cols) if c.startswith('L_1_1')],
             [i for i,c in enumerate(liver_cols) if c.startswith('L_1_2')]),
    '105kg': ([i for i,c in enumerate(liver_cols) if c.startswith('L_2_1')],
              [i for i,c in enumerate(liver_cols) if c.startswith('L_2_2')]),
}

muscle_stages = {
    '15kg': ([i for i,c in enumerate(muscle_cols) if c.startswith('m_15_1')],
             [i for i,c in enumerate(muscle_cols) if c.startswith('m_15_2')]),
    '45kg': ([i for i,c in enumerate(muscle_cols) if c.startswith('BJ_2_1')],
             [i for i,c in enumerate(muscle_cols) if c.startswith('BJ_2_2')]),
    '75kg': ([i for i,c in enumerate(muscle_cols) if c.startswith('m_1_1')],
             [i for i,c in enumerate(muscle_cols) if c.startswith('m_1_2')]),
    '105kg': ([i for i,c in enumerate(muscle_cols) if c.startswith('m_2_1')],
              [i for i,c in enumerate(muscle_cols) if c.startswith('m_2_2')]),
}

# ============================================================
# 3. Differential expression (DLY vs TFB) for each tissue/stage
# ============================================================
def run_de(mat, stages_dict):
    results = {}
    for stage, (di, ti) in stages_dict.items():
        dm = mat[:, di].mean(axis=1)
        tm = mat[:, ti].mean(axis=1)
        fc = dm - tm
        pv = np.ones(mat.shape[0])
        for g in range(mat.shape[0]):
            if dm[g] == tm[g]: continue
            try:
                t, p = ttest_ind(mat[g, di], mat[g, ti], equal_var=False)
                pv[g] = p
            except: pass
        results[stage] = {'log2FC': fc, 'pvalue': pv, 'DLY_mean': dm, 'TFB_mean': tm}
    return results

liver_de = run_de(liver_mat, liver_stages)
muscle_de = run_de(muscle_mat, muscle_stages)

# FDR correction per tissue per stage
for de_dict, name in [(liver_de, 'Liver'), (muscle_de, 'Muscle')]:
    for stage_name, stage_results in de_dict.items():
        pvals = stage_results['pvalue']
        _, qvals = benjamini_hochberg(pvals)
        stage_results['qvalue'] = qvals
        stage_results['FDR_significant'] = qvals < 0.05
        n_nom = (pvals < 0.05).sum()
        n_fdr = (qvals < 0.05).sum()
        print(f"  {name} {stage_name}: nominal p<0.05={n_nom} → FDR<0.05={n_fdr}")

# ============================================================
# 4. Tier classification
# ============================================================
def classify_4tier(de_results):
    stages = ['15kg', '45kg', '75kg', '105kg']
    n = len(de_results[stages[0]]['log2FC'])
    tier = np.full(n, 'Low_Signal', dtype=object)
    score = np.zeros(n)

    fcs = {s: de_results[s]['log2FC'] for s in stages}

    for g in range(n):
        fc_list = [fcs[s][g] for s in stages]
        afc_list = [abs(f) for f in fc_list]
        m = max(afc_list)
        if m <= 0.5:
            tier[g] = 'Low_Signal'; score[g] = m; continue

        mi = np.argmax(afc_list)

        # Tier 1: 15kg consistent
        if afc_list[0] > 0.5:
            s0 = np.sign(fc_list[0])
            later_agree = sum(1 for i, f in enumerate(fc_list[1:]) if abs(f) > 0.3 and np.sign(f) == s0)
            if later_agree >= 2:
                tier[g] = 'Tier1_Programming'
                score[g] = afc_list[0] + np.mean(afc_list[1:])
                continue

        # Tier 4: 105kg-specific
        if mi == 3 and afc_list[3] > 0.5 and max(afc_list[:3]) <= 0.5:
            tier[g] = 'Tier4_LateSpecific'
            score[g] = afc_list[3]; continue

        # Tier 2: 45kg peak
        if mi == 1 and afc_list[1] > max(afc_list[0], afc_list[2]):
            tier[g] = 'Tier2_Switch'
            score[g] = afc_list[1]; continue

        # Tier 3: 75kg or 105kg peak
        if mi in [2, 3]:
            tier[g] = 'Tier3_Consequence'
            score[g] = m; continue

        tier[g] = 'Mixed'; score[g] = m

    return tier, score

liver_tier, liver_score = classify_4tier(liver_de)
muscle_tier, muscle_score = classify_4tier(muscle_de)

# ============================================================
# 5. Build gene tables with annotations
# ============================================================
def build_table(gene_ids, gene_names, gene_descs, de_results, tier, score, tissue):
    stages = ['15kg', '45kg', '75kg', '105kg']
    data = {
        'Gene_Symbol': gene_names,
        'Gene_ID': gene_ids,
        'Description': gene_descs,
        'Tissue': tissue,
        'Tier': tier,
        'Tier_Score': score,
    }
    for s in stages:
        data[f'{s}_log2FC_DLYvsTFB'] = de_results[s]['log2FC']
        data[f'{s}_pvalue'] = de_results[s]['pvalue']
        data[f'{s}_qvalue_FDR'] = de_results[s]['qvalue']
        data[f'{s}_DLY_log2Expr'] = de_results[s]['DLY_mean']
        data[f'{s}_TFB_log2Expr'] = de_results[s]['TFB_mean']

    # Summary metrics
    fc_mat = np.column_stack([de_results[s]['log2FC'] for s in stages])
    data['Mean_abs_log2FC'] = np.mean(np.abs(fc_mat), axis=1)
    data['CrossStage_FC_SD'] = np.std(fc_mat, axis=1)

    # Direction info
    signs = np.sign(fc_mat)
    data['Direction_15kg'] = ['DLY_UP' if s > 0 else ('TFB_UP' if s < 0 else 'NEUTRAL') for s in signs[:, 0]]
    data['Majority_Direction'] = ['DLY_UP' if np.mean(s) > 0 else ('TFB_UP' if np.mean(s) < 0 else 'MIXED') for s in signs]

    df = pd.DataFrame(data)
    tier_order = {'Tier1_Programming': 0, 'Tier2_Switch': 1, 'Tier3_Consequence': 2,
                   'Tier4_LateSpecific': 3, 'Mixed': 4, 'Low_Signal': 5}
    df['_s'] = df['Tier'].map(tier_order)
    df = df.sort_values(['_s', 'Tier_Score'], ascending=[True, False]).drop(columns=['_s'])
    return df

liver_df = build_table(corr_gene_ids, liver_gene_name_filt.values, liver_gene_desc_filt.values,
                       liver_de, liver_tier, liver_score, 'Liver')
muscle_df = build_table(muscle_genes.tolist(), muscle_gene_name_filt.values, muscle_gene_desc_filt.values,
                        muscle_de, muscle_tier, muscle_score, 'Muscle')

# ============================================================
# 6. Crosstalk genes (expressed in both tissues)
# ============================================================
print("\nIdentifying crosstalk genes...")

# Map gene symbols between tissues
liver_symbols = {str(s).upper(): (i, gid) for i, (s, gid) in enumerate(zip(liver_gene_name_filt.values, corr_gene_ids))
                 if str(s) != 'nan' and str(s) != ''}
muscle_symbols = {str(s).upper(): (i, gid) for i, (s, gid) in enumerate(zip(muscle_gene_name_filt.values, muscle_genes))
                  if str(s) != 'nan' and str(s) != ''}

common_symbols = set(liver_symbols.keys()) & set(muscle_symbols.keys())
print(f"  Genes with matching symbols in both tissues: {len(common_symbols)}")

# Build crosstalk table
crosstalk_rows = []
for sym in common_symbols:
    li, lgid = liver_symbols[sym]
    mi, mgid = muscle_symbols[sym]

    row = {
        'Gene_Symbol': sym,
        'Liver_Gene_ID': lgid,
        'Muscle_Gene_ID': mgid,
    }
    for s in ['15kg', '45kg', '75kg', '105kg']:
        row[f'Liver_{s}_log2FC'] = liver_de[s]['log2FC'][li]
        row[f'Muscle_{s}_log2FC'] = muscle_de[s]['log2FC'][mi]
        row[f'Liver_{s}_pvalue'] = liver_de[s]['pvalue'][li]
        row[f'Muscle_{s}_pvalue'] = muscle_de[s]['pvalue'][mi]

    # Concordance score
    l_fcs = np.array([row[f'Liver_{s}_log2FC'] for s in ['15kg','45kg','75kg','105kg']])
    m_fcs = np.array([row[f'Muscle_{s}_log2FC'] for s in ['15kg','45kg','75kg','105kg']])

    # Same direction across stages?
    l_dir = np.sign(l_fcs)
    m_dir = np.sign(m_fcs)
    # Correlation of FC patterns between tissues
    if np.std(l_fcs) > 0.1 and np.std(m_fcs) > 0.1:
        row['LiverMuscle_FC_corr'], _ = pearsonr(l_fcs, m_fcs)
    else:
        row['LiverMuscle_FC_corr'] = 0

    # Agreement: same direction
    same_dir = sum(1 for i in range(4) if abs(l_fcs[i]) > 0.3 and abs(m_fcs[i]) > 0.3 and l_dir[i] == m_dir[i])
    total_meaningful = sum(1 for i in range(4) if abs(l_fcs[i]) > 0.3 and abs(m_fcs[i]) > 0.3)
    row['Direction_Agreement'] = same_dir / max(total_meaningful, 1)

    # 105kg divergence
    row['105kg_LiverMuscle_FC_diff'] = abs(l_fcs[3] - m_fcs[3])

    crosstalk_rows.append(row)

ct_df = pd.DataFrame(crosstalk_rows)
ct_df = ct_df.sort_values('LiverMuscle_FC_corr', ascending=False)

# Categorize crosstalk genes
def categorize_crosstalk(row):
    """Classify crosstalk pattern."""
    l_fcs = np.array([row[f'Liver_{s}_log2FC'] for s in ['15kg','45kg','75kg','105kg']])
    m_fcs = np.array([row[f'Muscle_{s}_log2FC'] for s in ['15kg','45kg','75kg','105kg']])

    l_max = np.argmax(np.abs(l_fcs))
    m_max = np.argmax(np.abs(m_fcs))

    corr = row['LiverMuscle_FC_corr']
    agreement = row['Direction_Agreement']

    if np.max(np.abs(l_fcs)) <= 0.5 and np.max(np.abs(m_fcs)) <= 0.5:
        return 'Low_Signal_Both'
    if np.max(np.abs(l_fcs)) > 0.5 and np.max(np.abs(m_fcs)) <= 0.5:
        return 'Liver_Specific'
    if np.max(np.abs(l_fcs)) <= 0.5 and np.max(np.abs(m_fcs)) > 0.5:
        return 'Muscle_Specific'

    if corr > 0.7 and agreement > 0.75:
        return 'Coordinated'
    elif corr < -0.5:
        return 'Opposing'
    elif agreement < 0.5:
        return 'Discordant'
    else:
        return 'Complex'

ct_df['Pattern'] = ct_df.apply(categorize_crosstalk, axis=1)

# ============================================================
# 7. Stage-specific pattern analysis
# ============================================================
print("\nAnalyzing stage-specific patterns...")

# For each gene, identify its "expression trajectory pattern"
# Pattern 1: DLY consistently higher (all 4 stages positive)
# Pattern 2: TFB consistently higher (all 4 stages negative)
# Pattern 3: Divergence increasing (FC magnitude grows over stages)
# Pattern 4: Convergence (FC magnitude shrinks over stages)
# Pattern 5: U-shaped or flip

def classify_trajectory(fc_list):
    afc = np.abs(fc_list)
    signs = np.sign(fc_list)

    # Check consistency
    if all(s == signs[0] for s in signs[1:]):
        # Consistent direction
        if signs[0] > 0:
            if afc[3] >= afc[0] * 1.5:
                return 'DLY_Increasing'
            elif afc[3] <= afc[0] * 0.5:
                return 'DLY_Converging'
            else:
                return 'DLY_Stable'
        else:
            if afc[3] >= afc[0] * 1.5:
                return 'TFB_Increasing'
            elif afc[3] <= afc[0] * 0.5:
                return 'TFB_Converging'
            else:
                return 'TFB_Stable'
    else:
        # Direction change
        return 'Direction_Change'

liver_trajectory = np.array([classify_trajectory([liver_de[s]['log2FC'][g] for s in ['15kg','45kg','75kg','105kg']])
                             for g in range(liver_mat.shape[0])])
muscle_trajectory = np.array([classify_trajectory([muscle_de[s]['log2FC'][g] for s in ['15kg','45kg','75kg','105kg']])
                              for g in range(muscle_mat.shape[0])])

liver_df['Trajectory'] = liver_trajectory
muscle_df['Trajectory'] = muscle_trajectory

# ============================================================
# 8. Summary statistics
# ============================================================
print(f"\n{'='*60}")
print("LIVER CLASSIFICATION SUMMARY (4 stages, corrected)")
print(f"{'='*60}")
for t in ['Tier1_Programming','Tier2_Switch','Tier3_Consequence','Tier4_LateSpecific','Mixed','Low_Signal']:
    n = (liver_tier == t).sum()
    print(f"  {t:<25} {n:>6} ({100*n/len(liver_tier):.1f}%)")

print(f"\nLiver Trajectory Patterns:")
for pat in ['DLY_Stable','DLY_Increasing','DLY_Converging','TFB_Stable','TFB_Increasing','TFB_Converging','Direction_Change']:
    n = (liver_trajectory == pat).sum()
    if n > 0:
        print(f"  {pat:<25} {n:>6} ({100*n/len(liver_trajectory):.1f}%)")

print(f"\n{'='*60}")
print("MUSCLE CLASSIFICATION SUMMARY (4 stages)")
print(f"{'='*60}")
for t in ['Tier1_Programming','Tier2_Switch','Tier3_Consequence','Tier4_LateSpecific','Mixed','Low_Signal']:
    n = (muscle_tier == t).sum()
    print(f"  {t:<25} {n:>6} ({100*n/len(muscle_tier):.1f}%)")

print(f"\nCrosstalk Pattern Summary ({len(ct_df)} genes):")
for pat in ['Coordinated','Opposing','Discordant','Complex','Liver_Specific','Muscle_Specific','Low_Signal_Both']:
    n = (ct_df['Pattern'] == pat).sum()
    if n > 0:
        print(f"  {pat:<25} {n:>6} ({100*n/len(ct_df):.1f}%)")

# ============================================================
# 9. Key gene highlights
# ============================================================
# AA metabolism genes in liver
aa_genes = ['CPS1','ASS1','ASL','ARG1','ARG2','OTC','NAGS','GOT1','GOT2','GPT','GPT2',
            'BCAT1','BCAT2','BCKDHA','BCKDHB','DBT','DLD','AASS','HGD','SDS','GLUD1',
            'PAH','HAL','IVD','MCCC1','MCCC2','AUH','HIBCH','HIBADH','ALDH6A1']

print(f"\n{'='*60}")
print("KEY AA METABOLISM GENES — LIVER (Tier1 only)")
print(f"{'='*60}")
t1_aa = liver_df[(liver_df['Tier'] == 'Tier1_Programming') &
                 (liver_df['Gene_Symbol'].str.upper().isin([g.upper() for g in aa_genes]))]
for _, r in t1_aa.iterrows():
    fcs = ' | '.join([f"{s}={r[f'{s}_log2FC_DLYvsTFB']:.2f}" for s in ['15kg','45kg','75kg','105kg']])
    print(f"  {r['Gene_Symbol']:<12} [{r['Trajectory']:<20}] {fcs}")

print(f"\n{'='*60}")
print("KEY GROWTH/INSULIN AXIS GENES — LIVER (Tier1+2)")
print(f"{'='*60}")
growth_genes = ['IGF1','IGFBP1','IGFBP2','IGFBP3','IGFALS','FGF21','FST','MSTN','BDNF',
                'PPARGC1A','FOXO1','FOXO3','KLF15','XBP1','ATF4','ATF3','DDIT3']
t12_growth = liver_df[(liver_df['Tier'].isin(['Tier1_Programming','Tier2_Switch'])) &
                       (liver_df['Gene_Symbol'].str.upper().isin([g.upper() for g in growth_genes]))]
for _, r in t12_growth.iterrows():
    fcs = ' | '.join([f"{s}={r[f'{s}_log2FC_DLYvsTFB']:.2f}" for s in ['15kg','45kg','75kg','105kg']])
    print(f"  {r['Gene_Symbol']:<12} [{r['Tier']:<22}] {fcs}")

print(f"\n{'='*60}")
print("TOP COORDINATED CROSSTALK GENES (Liver-Muscle)")
print(f"{'='*60}")
coordinated = ct_df[ct_df['Pattern'] == 'Coordinated'].head(20)
for _, r in coordinated.iterrows():
    lf = ' | '.join([f"{r[f'Liver_{s}_log2FC']:.2f}" for s in ['15kg','45kg','75kg','105kg']])
    mf = ' | '.join([f"{r[f'Muscle_{s}_log2FC']:.2f}" for s in ['15kg','45kg','75kg','105kg']])
    print(f"  {r['Gene_Symbol']:<20} r={r['LiverMuscle_FC_corr']:.2f}")
    print(f"    Liver:  {lf}")
    print(f"    Muscle: {mf}")

# ============================================================
# 10. Save all outputs
# ============================================================
print(f"\n{'='*60}")
print("SAVING OUTPUTS")
print(f"{'='*60}")

liver_df.to_excel(f'{out}/integrated_liver_4stage.xlsx', index=False)
print(f"  integrated_liver_4stage.xlsx — {len(liver_df)} liver genes")

muscle_df.to_excel(f'{out}/integrated_muscle_4stage.xlsx', index=False)
print(f"  integrated_muscle_4stage.xlsx — {len(muscle_df)} muscle genes")

ct_df.to_excel(f'{out}/integrated_crosstalk_4stage.xlsx', index=False)
print(f"  integrated_crosstalk_4stage.xlsx — {len(ct_df)} crosstalk genes")

# Causal subset
causal_liver = liver_df[liver_df['Tier'].isin(['Tier1_Programming','Tier2_Switch'])]
causal_muscle = muscle_df[muscle_df['Tier'].isin(['Tier1_Programming','Tier2_Switch'])]
causal_all = pd.concat([causal_liver, causal_muscle], ignore_index=True)
causal_all.to_excel(f'{out}/integrated_causal_4stage.xlsx', index=False)
print(f"  integrated_causal_4stage.xlsx — {len(causal_all)} Tier1+2 causal genes")

# Master summary sheet
with pd.ExcelWriter(f'{out}/integrated_master_4stage.xlsx', engine='openpyxl') as writer:
    liver_df.to_excel(writer, sheet_name='Liver_All', index=False)
    muscle_df.to_excel(writer, sheet_name='Muscle_All', index=False)
    ct_df.to_excel(writer, sheet_name='Crosstalk', index=False)
    causal_all.to_excel(writer, sheet_name='Tier1_Tier2_Causal', index=False)

    # Add summary sheet
    summary_data = {
        'Category': [],
        'Count': [],
        'Percentage': []
    }
    for t in ['Tier1_Programming','Tier2_Switch','Tier3_Consequence','Tier4_LateSpecific','Mixed','Low_Signal']:
        summary_data['Category'].append(f'Liver_{t}')
        summary_data['Count'].append((liver_tier == t).sum())
        summary_data['Percentage'].append(f"{100*(liver_tier==t).sum()/len(liver_tier):.1f}%")
    for t in ['Tier1_Programming','Tier2_Switch','Tier3_Consequence','Tier4_LateSpecific','Mixed','Low_Signal']:
        summary_data['Category'].append(f'Muscle_{t}')
        summary_data['Count'].append((muscle_tier == t).sum())
        summary_data['Percentage'].append(f"{100*(muscle_tier==t).sum()/len(muscle_tier):.1f}%")
    pd.DataFrame(summary_data).to_excel(writer, sheet_name='Summary', index=False)

print(f"  integrated_master_4stage.xlsx — Master workbook (4 sheets)")

print("\nDone.")
