#!/usr/bin/env python3
"""
Step 1: WGCNA Module Detection + Module-Trait Correlation + Hub Gene Identification
Adapts Jia et al. (2026) framework for pig liver-muscle axis study.

Pipeline:
  1. Prepare expression matrices (filtered, normalized)
  2. Build trait matrix (protein deposition, serum urea, breed, stage)
  3. Call WGCNA (R) for each tissue
  4. Parse results: modules, hub genes, module-trait correlations
  5. Generate summary figures
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import subprocess
import os
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("WGCNA-BASED MULTI-TISSUE CO-EXPRESSION ANALYSIS")
print("Adapted from Jia et al. (2026) J Anim Sci Biotechnol")
print("=" * 70)

# ============================================================
# 0. LOAD DATA
# ============================================================
print("\n[0] Loading data...")

liver_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')
muscle_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')

# Phenotype anchor
PROTEIN_DEPOSITION = {
    ('DLY', 15): 1.58, ('TFB', 15): 1.26,
    ('DLY', 45): 1.59, ('TFB', 45): 1.12,
    ('DLY', 75): 1.11, ('TFB', 75): 0.68,
    ('DLY', 105): 0.87, ('TFB', 105): 0.49,
}
SERUM_UREA = {
    ('DLY', 15): 0.81, ('TFB', 15): 3.16,
    ('DLY', 45): 2.30, ('TFB', 45): 5.02,
    ('DLY', 75): 2.71, ('TFB', 75): 2.71,
    ('DLY', 105): 2.62, ('TFB', 105): 6.08,
}

sample_map_l = {
    'L_15_1_': ('DLY', 15), 'L_15_2_': ('TFB', 15),
    'L_45_1_': ('DLY', 45), 'L_45_2_': ('TFB', 45),
    'L_1_1_': ('DLY', 75), 'L_1_2_': ('TFB', 75),
    'L_2_1_': ('DLY', 105), 'L_2_2_': ('TFB', 105),
}
sample_map_m = {
    'm_15_1_': ('DLY', 15), 'm_15_2_': ('TFB', 15),
    'BJ_2_1_': ('DLY', 45), 'BJ_2_2_': ('TFB', 45),
    'm_1_1_': ('DLY', 75), 'm_1_2_': ('TFB', 75),
    'm_2_1_': ('DLY', 105), 'm_2_2_': ('TFB', 105),
}

print(f"  Liver genes: {liver_raw.shape[0]}, Muscle genes: {muscle_raw.shape[0]}")

# ============================================================
# 1. BUILD SAMPLE-LEVEL EXPRESSION MATRICES
# ============================================================
print("\n[1] Building sample-level expression matrices...")

def build_sample_matrix(mat, sample_map):
    """Build (genes × samples) expression matrix, excluding DLY 135kg."""
    meta_cols = ['seq_id', 'gene_name', 'length', 'description']
    val_cols = [c for c in mat.columns if c not in meta_cols]

    # Map each column to a sample
    sample_info = {}
    for col in val_cols:
        for prefix, (breed, stage) in sample_map.items():
            if col.startswith(prefix):
                sample_name = f"{breed}_{stage}kg_{col.split('_')[-1]}"
                sample_info[col] = (sample_name, breed, stage)
                break

    # Build matrix
    records = []
    for _, row in mat.iterrows():
        gn = str(row['gene_name']) if pd.notna(row['gene_name']) else row['seq_id']
        for col in val_cols:
            if col in sample_info:
                sname, breed, stage = sample_info[col]
                if pd.notna(row[col]):
                    records.append({'Gene': gn, 'Sample': sname, 'Breed': breed,
                                    'Weight': stage, 'Expr': float(row[col])})

    df = pd.DataFrame(records)
    # Pivot to gene × sample matrix
    mat_pivot = df.pivot_table(index='Gene', columns='Sample', values='Expr', aggfunc='mean')
    return mat_pivot, df

liver_mat, liver_long = build_sample_matrix(liver_raw, sample_map_l)
muscle_mat, muscle_long = build_sample_matrix(muscle_raw, sample_map_m)

print(f"  Liver matrix: {liver_mat.shape[0]} genes × {liver_mat.shape[1]} samples")
print(f"  Muscle matrix: {muscle_mat.shape[0]} genes × {muscle_mat.shape[1]} samples")

# ============================================================
# 2. FILTER LOW-EXPRESSION AND LOW-VARIANCE GENES
# ============================================================
print("\n[2] Filtering low-expression and low-variance genes...")

def filter_genes(mat, min_expr=0.1, top_n=8000):
    """Filter genes for WGCNA: remove noise, keep informative genes.

    min_expr=0.1: very lenient — removes only near-zero expression genes
    top_n=8000: keep top N most variable genes (WGCNA sweet spot)
    """
    mean_expr = mat.mean(axis=1)
    mat_f = mat.loc[mean_expr > min_expr]

    # Keep top N most variable genes by variance
    if mat_f.shape[0] > top_n:
        var_expr = mat_f.var(axis=1)
        top_idx = var_expr.nlargest(top_n).index
        mat_f = mat_f.loc[top_idx]

    return mat_f

liver_filt = filter_genes(liver_mat)
muscle_filt = filter_genes(muscle_mat)

print(f"  Liver filtered: {liver_filt.shape[0]} genes × {liver_filt.shape[1]} samples")
print(f"  Muscle filtered: {muscle_filt.shape[0]} genes × {muscle_filt.shape[1]} samples")

# ============================================================
# 3. BUILD TRAIT MATRIX
# ============================================================
print("\n[3] Building trait matrix...")

def build_trait_matrix(mat):
    """Build sample-level trait matrix."""
    samples = mat.columns.tolist()
    traits = {}
    for s in samples:
        parts = s.split('_')
        breed = parts[0]
        stage_str = parts[1].replace('kg', '')
        stage = int(stage_str)

        traits[s] = {
            'Breed': 1 if breed == 'DLY' else 0,
            'Weight': stage,
            'PD': PROTEIN_DEPOSITION.get((breed, stage), np.nan),
            'Urea': SERUM_UREA.get((breed, stage), np.nan),
            'Breed_x_Weight': (1 if breed == 'DLY' else 0) * stage,
        }

    return pd.DataFrame(traits).T

liver_traits = build_trait_matrix(liver_filt)
muscle_traits = build_trait_matrix(muscle_filt)

print(f"  Liver traits: {liver_traits.shape}")
print(f"  Muscle traits: {muscle_traits.shape}")
print(f"  Trait summary:\n{liver_traits.describe().to_string()}")

# ============================================================
# 4. EXPORT DATA FOR R WGCNA
# ============================================================
print("\n[4] Exporting data for WGCNA (R)...")

os.makedirs('wgcna_output', exist_ok=True)

# Transpose for WGCNA: genes as columns, samples as rows
liver_filt.T.to_csv('wgcna_output/liver_expr.csv')
muscle_filt.T.to_csv('wgcna_output/muscle_expr.csv')
liver_traits.to_csv('wgcna_output/liver_traits.csv')
muscle_traits.to_csv('wgcna_output/muscle_traits.csv')

print("  Exported to wgcna_output/")

# ============================================================
# 5. RUN WGCNA (R)
# ============================================================
print("\n[5] Running WGCNA for each tissue...")

r_script = 'wgcna_analysis_step1.R'

for tissue in ['Liver', 'Muscle']:
    tissue_lower = tissue.lower()
    print(f"\n  --- {tissue} ---")
    cmd = [
        'Rscript', r_script,
        tissue,
        f'wgcna_output/{tissue_lower}_expr.csv',
        f'wgcna_output/{tissue_lower}_traits.csv',
        f'wgcna_output/{tissue_lower}'
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    print(result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout)
    if result.stderr and 'Error' in result.stderr:
        print(f"  STDERR: {result.stderr[:1000]}")
    if result.returncode != 0:
        print(f"  WARNING: R script returned exit code {result.returncode}")

# ============================================================
# 6. PARSE WGCNA RESULTS
# ============================================================
print("\n[6] Parsing WGCNA results...")

def parse_wgcna_results(tissue):
    """Parse WGCNA output files."""
    prefix = f'wgcna_output/{tissue.lower()}'

    gene_mod = pd.read_csv(f'{prefix}_gene_module_assignment.csv')
    mod_trait_cor = pd.read_csv(f'{prefix}_module_trait_cor.csv', index_col=0)
    mod_trait_pval = pd.read_csv(f'{prefix}_module_trait_pvalue.csv', index_col=0)
    hub_genes = pd.read_csv(f'{prefix}_hub_genes.csv')
    mod_sizes = pd.read_csv(f'{prefix}_module_sizes.csv')

    return gene_mod, mod_trait_cor, mod_trait_pval, hub_genes, mod_sizes

liver_gm, liver_mtc, liver_mtp, liver_hub, liver_ms = parse_wgcna_results('Liver')
muscle_gm, muscle_mtc, muscle_mtp, muscle_hub, muscle_ms = parse_wgcna_results('Muscle')

# Print module summary
for tissue, gm, mtc, ms in [('Liver', liver_gm, liver_mtc, liver_ms),
                              ('Muscle', muscle_gm, muscle_mtc, muscle_ms)]:
    print(f"\n  {tissue} Modules:")
    for _, r in ms.iterrows():
        mod = r['Module']
        size = r['Size']
        if mod == 'grey':
            continue
        # Show correlation with key traits
        pd_cor = mtc.loc[mod, 'PD'] if mod in mtc.index else np.nan
        urea_cor = mtc.loc[mod, 'Urea'] if mod in mtc.index else np.nan
        print(f"    {mod:20s} size={size:5d}  r_PD={pd_cor:7.3f}  r_Urea={urea_cor:7.3f}")

# ============================================================
# 7. IDENTIFY PHENOTYPE-ASSOCIATED MODULES
# ============================================================
print("\n[7] Identifying phenotype-associated modules (like paper's liver M6/M7, muscle M16)...")

def find_significant_modules(mtc, mtp, tissue):
    """Find modules significantly associated with protein deposition or serum urea."""
    sig_mods = []
    for mod in mtc.index:
        if mod == 'grey':
            continue
        pd_r = mtc.loc[mod, 'PD']
        pd_p = mtp.loc[mod, 'PD']
        urea_r = mtc.loc[mod, 'Urea']
        urea_p = mtp.loc[mod, 'Urea']

        if pd_p < 0.05 or urea_p < 0.05:
            sig_mods.append({
                'Tissue': tissue,
                'Module': mod,
                'r_PD': round(pd_r, 3),
                'p_PD': round(pd_p, 5),
                'r_Urea': round(urea_r, 3),
                'p_Urea': round(urea_p, 5),
            })
    if len(sig_mods) == 0:
        return pd.DataFrame(columns=['Tissue', 'Module', 'r_PD',
                                     'p_PD', 'r_Urea', 'p_Urea'])
    return pd.DataFrame(sig_mods).sort_values('p_PD')

liver_sig = find_significant_modules(liver_mtc, liver_mtp, 'Liver')
muscle_sig = find_significant_modules(muscle_mtc, muscle_mtp, 'Muscle')

if len(liver_sig) > 0:
    print(f"\n  Liver phenotype-associated modules ({len(liver_sig)}):")
    for _, r in liver_sig.iterrows():
        print(f"    {r['Module']:20s} r_PD={r['r_PD']:.3f} p_PD={r['p_PD']:.5f}  r_Urea={r['r_Urea']:.3f} p_Urea={r['p_Urea']:.5f}")
else:
    print("\n  Liver: No modules with significant phenotype association (p<0.05)")
    print("  Showing all module-trait correlations instead:")
    for mod in liver_mtc.index:
        if mod == 'grey':
            continue
        pd_r = liver_mtc.loc[mod, 'PD'] if 'PD' in liver_mtc.columns else np.nan
        urea_r = liver_mtc.loc[mod, 'Urea'] if 'Urea' in liver_mtc.columns else np.nan
        print(f"    {mod:20s} r_PD={pd_r:.3f}  r_Urea={urea_r:.3f}")

if len(muscle_sig) > 0:
    print(f"\n  Muscle phenotype-associated modules ({len(muscle_sig)}):")
    for _, r in muscle_sig.iterrows():
        print(f"    {r['Module']:20s} r_PD={r['r_PD']:.3f} p_PD={r['p_PD']:.5f}  r_Urea={r['r_Urea']:.3f} p_Urea={r['p_Urea']:.5f}")
else:
    print("\n  Muscle: No modules with significant phenotype association (p<0.05)")
    print("  Showing all module-trait correlations instead:")
    for mod in muscle_mtc.index:
        if mod == 'grey':
            continue
        pd_r = muscle_mtc.loc[mod, 'PD'] if 'PD' in muscle_mtc.columns else np.nan
        urea_r = muscle_mtc.loc[mod, 'Urea'] if 'Urea' in muscle_mtc.columns else np.nan
        print(f"    {mod:20s} r_PD={pd_r:.3f}  r_Urea={urea_r:.3f}")

# ============================================================
# 8. CHECK KEY GENES IN MODULES
# ============================================================
print("\n[8] Key gene module assignments...")

key_genes_liver = ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'ASS1', 'ASL', 'HGD',
                   'GLUD1', 'AASS', 'HAL', 'BCKDHA', 'OTC', 'ARG2', 'NAGS']
key_genes_muscle = ['FOXO1', 'FOXO3', 'FBXO32', 'TRIM63', 'MYOG', 'MYOD1', 'IGF1',
                    'MTOR', 'AKT1', 'MSTN', 'RPS6KB1', 'EEF2', 'IRS2']

for gene in key_genes_liver:
    rows = liver_gm[liver_gm['Gene'] == gene]
    if len(rows) > 0:
        mod = rows.iloc[0]['Module']
        kme = rows.iloc[0].get('kME_module', np.nan)
        gs_pd = rows.iloc[0].get('GS_PD', np.nan)
        gs_urea = rows.iloc[0].get('GS_Urea', np.nan)
        # Check if module is phenotype-associated
        in_sig = mod in liver_sig['Module'].values if len(liver_sig) > 0 else False
        sig_mark = '*** SIG MODULE ***' if in_sig else ''
        print(f"  {gene:10s} Liver  Module={mod:15s} kME={kme:.3f} GS_PD={gs_pd:.3f} GS_Urea={gs_urea:.3f} {sig_mark}")
    else:
        print(f"  {gene:10s} Liver  NOT IN FILTERED SET")

for gene in key_genes_muscle:
    rows = muscle_gm[muscle_gm['Gene'] == gene]
    if len(rows) > 0:
        mod = rows.iloc[0]['Module']
        kme = rows.iloc[0].get('kME_module', np.nan)
        gs_pd = rows.iloc[0].get('GS_PD', np.nan)
        gs_urea = rows.iloc[0].get('GS_Urea', np.nan)
        in_sig = mod in muscle_sig['Module'].values if len(muscle_sig) > 0 else False
        sig_mark = '*** SIG MODULE ***' if in_sig else ''
        print(f"  {gene:10s} Muscle Module={mod:15s} kME={kme:.3f} GS_PD={gs_pd:.3f} GS_Urea={gs_urea:.3f} {sig_mark}")
    else:
        print(f"  {gene:10s} Muscle NOT IN FILTERED SET")

# ============================================================
# 9. HUB GENE SUMMARY
# ============================================================
print("\n[9] Hub gene summary (top 5 per significant module)...")

for tissue, hub_df, sig_mods in [('Liver', liver_hub, liver_sig), ('Muscle', muscle_hub, muscle_sig)]:
    print(f"\n  {tissue} Hub Genes:")
    if len(sig_mods) == 0:
        print("    No significant modules found")
        continue

    for _, mod_row in sig_mods.iterrows():
        mod = mod_row['Module']
        mod_hubs = hub_df[hub_df['Module'] == mod]
        if len(mod_hubs) > 0:
            top5 = mod_hubs.head(5)
            print(f"    Module {mod} (r_PD={mod_row['r_PD']:.3f}): "
                  f"{', '.join(top5['Gene'].tolist())}")

# ============================================================
# 10. GENERATE SUMMARY FIGURES
# ============================================================
print("\n[10] Generating summary figures...")

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

# Fig 1: Module-Trait Correlation Heatmap (mimicking paper's Fig 7A)
fig1, axes1 = plt.subplots(1, 2, figsize=(14, max(6, len(liver_mtc) * 0.35)))

for idx, (tissue, mtc, mtp, sig) in enumerate([
    ('Liver', liver_mtc, liver_mtp, liver_sig),
    ('Muscle', muscle_mtc, muscle_mtp, muscle_sig)
]):
    ax = axes1[idx]

    # Only show non-grey modules
    modules_show = [m for m in mtc.index if m != 'grey']
    traits_show = ['PD', 'Urea', 'Breed', 'Weight']

    cor_matrix = mtc.loc[modules_show, traits_show].values
    p_matrix = mtp.loc[modules_show, traits_show].values

    im = ax.imshow(cor_matrix, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)

    # Annotate
    for i in range(len(modules_show)):
        for j in range(len(traits_show)):
            r_val = cor_matrix[i, j]
            p_val = p_matrix[i, j]
            sig_str = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
            text = f'{r_val:.2f}{sig_str}'
            ax.text(j, i, text, ha='center', va='center', fontsize=7,
                   color='white' if abs(r_val) > 0.5 else 'black', fontweight='bold' if sig_str else 'normal')

    ax.set_xticks(range(len(traits_show)))
    ax.set_xticklabels(traits_show, rotation=30, ha='right', fontsize=8)
    ax.set_yticks(range(len(modules_show)))
    ax.set_yticklabels(modules_show, fontsize=8)
    ax.set_title(f'{tissue} Module-Trait Correlations\n(Adapted from Jia et al. 2026 Fig 7A)', fontweight='bold')

    # Highlight significant modules
    if len(sig) > 0:
        sig_mod_names = sig['Module'].tolist()
        for i, mod in enumerate(modules_show):
            if mod in sig_mod_names:
                ax.add_patch(plt.Rectangle((-0.5, i-0.5), len(traits_show), 1,
                                           fill=False, edgecolor='#D73027', linewidth=2, linestyle='--'))

    plt.colorbar(im, ax=ax, shrink=0.8, label='Pearson r')

fig1.suptitle('WGCNA Module-Trait Associations in Pig Liver and Muscle\n(Protein Deposition as Phenotype Anchor)',
             fontweight='bold', fontsize=12)
plt.tight_layout()
fig1.savefig('fig_wgcna_module_trait_heatmap.png')
print("  Saved fig_wgcna_module_trait_heatmap.png")

# Fig 2: Module size distribution
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))

for idx, (tissue, ms, sig) in enumerate([
    ('Liver', liver_ms, liver_sig),
    ('Muscle', muscle_ms, muscle_sig)
]):
    ax = axes2[idx]
    ms_plot = ms[ms['Module'] != 'grey'].sort_values('Size', ascending=False)
    colors = ['#D73027' if m in sig['Module'].values else '#91BFDB'
              for m in ms_plot['Module']] if len(sig) > 0 else '#91BFDB'

    ax.barh(range(len(ms_plot)), ms_plot['Size'], color=colors, edgecolor='white')
    ax.set_yticks(range(len(ms_plot)))
    ax.set_yticklabels(ms_plot['Module'], fontsize=7)
    ax.set_xlabel('Number of Genes')
    ax.set_title(f'{tissue} Module Sizes\n(Red = phenotype-associated, P<0.05)', fontweight='bold')
    ax.invert_yaxis()

fig2.suptitle('WGCNA Co-Expression Module Sizes', fontweight='bold', fontsize=12)
plt.tight_layout()
fig2.savefig('fig_wgcna_module_sizes.png')
print("  Saved fig_wgcna_module_sizes.png")

# ============================================================
# SAVE RESULTS
# ============================================================
print("\nSaving integrated results...")

with pd.ExcelWriter('wgcna_step1_results.xlsx', engine='openpyxl') as writer:
    # Module assignments
    liver_gm.to_excel(writer, sheet_name='Liver_Module_Assignment', index=False)
    muscle_gm.to_excel(writer, sheet_name='Muscle_Module_Assignment', index=False)

    # Module-trait correlations
    liver_mtc.to_excel(writer, sheet_name='Liver_ModuleTrait_Cor')
    muscle_mtc.to_excel(writer, sheet_name='Muscle_ModuleTrait_Cor')

    # Hub genes
    liver_hub.to_excel(writer, sheet_name='Liver_Hub_Genes', index=False)
    muscle_hub.to_excel(writer, sheet_name='Muscle_Hub_Genes', index=False)

    # Significant modules
    if len(liver_sig) > 0:
        liver_sig.to_excel(writer, sheet_name='Liver_SigModules', index=False)
    if len(muscle_sig) > 0:
        muscle_sig.to_excel(writer, sheet_name='Muscle_SigModules', index=False)

print("Saved wgcna_step1_results.xlsx")

print("\n" + "=" * 70)
print("WGCNA STEP 1 COMPLETE")
print("=" * 70)
print(f"""
Key Outputs:
  - wgcna_step1_results.xlsx (module assignments, hub genes, correlations)
  - fig_wgcna_module_trait_heatmap.png
  - fig_wgcna_module_sizes.png
  - wgcna_output/liver_gene_module_assignment.csv
  - wgcna_output/muscle_gene_module_assignment.csv
  - wgcna_output/liver_hub_genes.csv
  - wgcna_output/muscle_hub_genes.csv
  - wgcna_output/liver_module_trait_cor.csv
  - wgcna_output/muscle_module_trait_cor.csv

Next Step (Step 2):
  - GO/KEGG enrichment on phenotype-associated module genes
  - Cross-tissue module-module correlation (liver ME ~ muscle ME)
  - Serum metabolite integration as middle layer
  - Publication-quality integrated figures
""")

print("Done!")
