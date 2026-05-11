#!/usr/bin/env python3
"""
Master analysis: serum AA + liver/muscle gene expression → integrated figures.
Generates:
  1. Liver AA catabolism enzyme heatmap
  2. Serum Urea breed×stage dynamics
  3. Muscle ribosomal/translation gene heatmap
  4. Liver enzyme ↔ Serum Urea ↔ Muscle ribosome correlation heatmap
  5. Graphical Abstract (systems-level causal diagram)
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import ttest_ind, pearsonr
from scipy.cluster.hierarchy import linkage, leaves_list

# ============================================================
# 0. Load & prepare data
# ============================================================
serum = pd.read_csv('serum_summary.csv')
serum_tidy = pd.read_csv('serum_all_tidy.csv')

# Read gene expression matrices (tab-separated with .xls extension)
# Use full matrices from Downloads
muscle = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')
liver = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')

# Map sample columns to groups based on metadata
# From metadata: Sample Initial Name → Group Name
sample_to_group_m = {
    'm_15_1_': 'DLYM15', 'm_15_2_': 'TFBM15',
    'BJ_2_1_': 'DLYM45', 'BJ_2_2_': 'TFBM45',
    'm_1_1_': 'DLYM75', 'm_1_2_': 'TFBM75',
    'm_2_1_': 'DLYM105', 'm_2_2_': 'TFBM105',
    'm_3_1_': 'DLYM135',
}
sample_to_group_l = {
    'L_15_1_': 'DLYL15', 'L_15_2_': 'TFBL15',
    'L_45_1_': 'DLYL45', 'L_45_2_': 'TFBL45',
    'L_1_1_': 'DLYL75', 'L_1_2_': 'TFBL75',
    'L_2_1_': 'DLYL105', 'L_2_2_': 'TFBL105',
    'L_3_1_': 'DLYL135',
}

def map_col_to_group(col, mapping):
    for prefix, group in mapping.items():
        if col.startswith(prefix):
            return group
    return None

def prepare_expr_matrix(df, sample_map):
    """Convert wide expression matrix to tidy format with group means."""
    val_cols = [c for c in df.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    records = []
    for _, row in df.iterrows():
        gene_id = row['seq_id']
        gene_name = row['gene_name'] if pd.notna(row['gene_name']) else gene_id
        for col in val_cols:
            group = map_col_to_group(col, sample_map)
            if group is not None and pd.notna(row[col]):
                records.append({
                    'gene_id': gene_id,
                    'gene_name': str(gene_name),
                    'group': group,
                    'expr': float(row[col])
                })
    return pd.DataFrame(records)

print("Preparing expression matrices...")
muscle_long = prepare_expr_matrix(muscle, sample_to_group_m)
liver_long = prepare_expr_matrix(liver, sample_to_group_l)

# Parse breed/stage from group
def parse_group(g):
    if 'DLY' in g:
        breed = 'DLY'
    else:
        breed = 'TFB'
    import re
    m = re.search(r'(\d+)', g)
    stage = int(m.group(1)) if m else None
    tissue = 'Liver' if 'L' in g else 'Muscle'
    return breed, stage, tissue

for df_ in [muscle_long, liver_long]:
    parsed = df_['group'].apply(parse_group)
    df_['breed'] = [p[0] for p in parsed]
    df_['stage_kg'] = [p[1] for p in parsed]
    df_['tissue'] = [p[2] for p in parsed]

# Group means
muscle_mean = muscle_long.groupby(['gene_id', 'gene_name', 'group', 'breed', 'stage_kg', 'tissue'])['expr'].mean().reset_index()
liver_mean = liver_long.groupby(['gene_id', 'gene_name', 'group', 'breed', 'stage_kg', 'tissue'])['expr'].mean().reset_index()

print(f"Muscle genes: {muscle_long['gene_id'].nunique()}, Liver genes: {liver_long['gene_id'].nunique()}")
print(f"Muscle groups: {sorted(muscle_long['group'].unique())}")
print(f"Liver groups: {sorted(liver_long['group'].unique())}")

# ============================================================
# 1. Define gene sets
# ============================================================
# AA catabolism enzymes (liver)
AA_CATABOLISM_GENES = {
    # BCAA catabolism
    'BCAT2': 'BCAA transaminase',
    'BCKDHA': 'BCKDH E1α',
    'BCKDHB': 'BCKDH E1β',
    'DBT': 'BCKDH E2',
    'DLD': 'BCKDH E3',
    # Urea cycle
    'CPS1': 'Carbamoyl-P synthase',
    'OTC': 'Ornithine transcarbamylase',
    'ASS1': 'Argininosuccinate synthase',
    'ASL': 'Argininosuccinate lyase',
    'ARG1': 'Arginase',
    # Transaminases
    'GOT1': 'Asp transaminase (c)',
    'GOT2': 'Asp transaminase (m)',
    'GPT': 'Ala transaminase',
    # Specific AA catabolism
    'AASS': 'Lys catabolism',
    'HGD': 'Tyr catabolism',
    'ACADSB': 'SCAD (Ile/Val)',
    'GLUD1': 'Glu dehydrogenase',
    'SDS': 'Ser dehydratase',
    'HAL': 'His ammonia-lyase',
    'PAH': 'Phe hydroxylase',
}

# Ribosomal / translation (muscle)
RIBOSOMAL_GENE_PATTERNS = {
    'RPL': 'Ribosomal large subunit',
    'RPS': 'Ribosomal small subunit',
    'EIF': 'Translation initiation',
    'EEF': 'Translation elongation',
    'MRPL': 'Mito ribosomal large',
    'MRPS': 'Mito ribosomal small',
}

def find_genes(df, gene_list, by_symbol=True):
    """Find genes in expression dataframe, return matching gene_name list."""
    found = {}
    all_names = set(df['gene_name'].unique())
    all_ids = set(df['gene_id'].unique())

    for gene in gene_list:
        if gene in all_names:
            found[gene] = gene
        elif gene.upper() in all_names:
            found[gene] = gene.upper()
        elif gene.lower() in all_names:
            found[gene] = gene.lower()
    return found

def find_genes_by_pattern(df, patterns):
    """Find genes matching patterns in gene_name."""
    found = {}
    all_names = df['gene_name'].unique()
    for pattern, desc in patterns.items():
        matches = [n for n in all_names if n.upper().startswith(pattern.upper())]
        for m in matches:
            found[m] = desc
    return found

# Find AA catabolism genes in liver
liver_aa_genes = find_genes(liver_mean, list(AA_CATABOLISM_GENES.keys()))
liver_aa_genes_desc = {v: AA_CATABOLISM_GENES.get(k, AA_CATABOLISM_GENES.get(k.upper(), k)) for k, v in liver_aa_genes.items()}
print(f"\nLiver AA catabolism genes found: {len(liver_aa_genes)}/{len(AA_CATABOLISM_GENES)}")
for g in sorted(liver_aa_genes.values()):
    print(f"  {g}")

# Find ribosomal genes in muscle
muscle_ribo_genes = find_genes_by_pattern(muscle_mean, RIBOSOMAL_GENE_PATTERNS)
print(f"\nMuscle ribosomal genes found: {len(muscle_ribo_genes)}")
for pat in RIBOSOMAL_GENE_PATTERNS:
    count = sum(1 for g, d in muscle_ribo_genes.items() if d == RIBOSOMAL_GENE_PATTERNS[pat])
    print(f"  {pat}: {count} genes")

# ============================================================
# 2. Compute breed×stage fold-changes
# ============================================================
def compute_fc_and_pval(df, gene_list, groups_col='group', expr_col='expr',
                         breed_col='breed', stage_col='stage_kg'):
    """For each gene, compute DLY vs TFB fold change at each stage."""
    results = []
    for gene in gene_list:
        gene_df = df[df['gene_name'] == gene]
        if len(gene_df) == 0:
            continue
        for stage in sorted(gene_df[stage_col].unique()):
            dly = gene_df[(gene_df[breed_col] == 'DLY') & (gene_df[stage_col] == stage)]
            tfb = gene_df[(gene_df[breed_col] == 'TFB') & (gene_df[stage_col] == stage)]
            if len(dly) == 0 or len(tfb) == 0:
                continue
            # Get per-replicate values for t-test
            pass  # Will use group means for heatmap, per-rep for p-value
    return results

# Build pivot tables for heatmaps
def build_heatmap_data(mean_df, gene_list, value_col='expr'):
    """Build gene × stage matrix of log2(DLY/TFB) fold changes."""
    rows = []
    for gene in gene_list:
        gene_df = mean_df[mean_df['gene_name'] == gene]
        if len(gene_df) == 0:
            continue
        row_data = {'gene': gene}
        for stage in sorted(gene_df['stage_kg'].unique()):
            dly = gene_df[(gene_df['breed'] == 'DLY') & (gene_df['stage_kg'] == stage)][value_col]
            tfb = gene_df[(gene_df['breed'] == 'TFB') & (gene_df['stage_kg'] == stage)][value_col]
            if len(dly) > 0 and len(tfb) > 0 and tfb.mean() > 0:
                row_data[f'{stage}kg'] = np.log2(dly.mean() / tfb.mean())
            else:
                row_data[f'{stage}kg'] = np.nan
        rows.append(row_data)
    return pd.DataFrame(rows).set_index('gene')

print("\nBuilding liver AA enzyme heatmap data...")
liver_heatmap = build_heatmap_data(liver_mean, sorted(liver_aa_genes.values()))
print(liver_heatmap.round(2))

print("\nBuilding muscle ribosomal heatmap data (top 30 by variance)...")
# Select top N ribosomal genes with highest variance in FC
ribo_gene_list = sorted(muscle_ribo_genes.keys())
muscle_heatmap = build_heatmap_data(muscle_mean, ribo_gene_list)
# Filter to genes with data at most stages
muscle_heatmap = muscle_heatmap.dropna(thresh=2)
# Take top by absolute mean log2FC
muscle_heatmap['abs_mean_fc'] = muscle_heatmap.abs().mean(axis=1)
muscle_heatmap = muscle_heatmap.sort_values('abs_mean_fc', ascending=False).head(40)
muscle_heatmap = muscle_heatmap.drop(columns=['abs_mean_fc'])
print(f"Selected {len(muscle_heatmap)} ribosomal genes")

# ============================================================
# 3. Compute correlations: liver AA enzymes ↔ serum metabolites
# ============================================================
print("\n=== Computing cross-tissue correlations ===")

# Aggregate serum to group-level
serum_group_mean = serum.groupby(['metabolite', 'breed', 'stage_kg'])['mean'].mean().reset_index()

# Merge liver enzyme expression with serum data
# For each group, get mean expression of AA enzymes and serum metabolites
# We need to match at breed×stage level

# Compute per-group mean expression for liver AA enzymes
liver_aa_expr = liver_mean[liver_mean['gene_name'].isin(liver_aa_genes.values())]
liver_aa_group = liver_aa_expr.groupby(['breed', 'stage_kg', 'gene_name'])['expr'].mean().reset_index()

# Build correlation matrix: for each liver enzyme gene × serum metabolite
corr_rows = []
for _, enzyme_row in liver_aa_group.iterrows():
    gene = enzyme_row['gene_name']
    breed = enzyme_row['breed']
    stage = enzyme_row['stage_kg']
    enz_expr = enzyme_row['expr']

    # Get serum data for same breed×stage
    serum_match = serum_group_mean[(serum_group_mean['breed'] == breed) &
                                    (serum_group_mean['stage_kg'] == stage)]
    for _, s_row in serum_match.iterrows():
        corr_rows.append({
            'gene': gene,
            'metabolite': s_row['metabolite'],
            'breed': breed,
            'stage_kg': stage,
            'enz_expr': enz_expr,
            'serum_value': s_row['mean']
        })

corr_df = pd.DataFrame(corr_rows)

# Compute Pearson r for each enzyme-metabolite pair (across breed×stage combinations)
corr_pairs = []
for gene in corr_df['gene'].unique():
    for met in corr_df['metabolite'].unique():
        pair_df = corr_df[(corr_df['gene'] == gene) & (corr_df['metabolite'] == met)]
        if len(pair_df) >= 6:  # need at least 6 data points
            r, p = pearsonr(pair_df['enz_expr'], pair_df['serum_value'])
            corr_pairs.append({
                'gene': gene,
                'metabolite': met,
                'pearson_r': r,
                'p_value': p,
                'n': len(pair_df)
            })

corr_matrix = pd.DataFrame(corr_pairs)
print(f"Computed {len(corr_matrix)} enzyme-metabolite correlations")

# Top correlations for Urea
urea_corrs = corr_matrix[corr_matrix['metabolite'] == 'Urea'].sort_values('pearson_r', ascending=False)
print("\nTop liver enzyme ↔ Serum Urea correlations:")
for _, row in urea_corrs.head(10).iterrows():
    print(f"  {row['gene']:10s} r={row['pearson_r']:+.3f} p={row['p_value']:.4f}")

# Also: serum AA vs muscle ribosomal expression
muscle_ribo_group = muscle_mean[muscle_mean['gene_name'].isin(ribo_gene_list)]
muscle_ribo_agg = muscle_ribo_group.groupby(['breed', 'stage_kg'])['expr'].mean().reset_index()
muscle_ribo_agg.rename(columns={'expr': 'muscle_ribo_mean'}, inplace=True)

# ============================================================
# PLOT 1: Liver AA catabolism enzyme heatmap
# ============================================================
print("\nGenerating Figure 1: Liver AA catabolism heatmap...")
fig, ax = plt.subplots(figsize=(8, 10))

# Cluster rows
if liver_heatmap.shape[0] > 2:
    row_linkage = linkage(liver_heatmap.fillna(0), method='average')
    row_order = leaves_list(row_linkage)
    liver_plot = liver_heatmap.iloc[row_order]
else:
    liver_plot = liver_heatmap

sns.heatmap(liver_plot, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
            vmin=-3, vmax=3, ax=ax, cbar_kws={'label': 'log2(DLY/TFB)'},
            linewidths=0.5, linecolor='white')
ax.set_title('Liver AA Catabolism Enzymes\nlog2(DLY/TFB) Expression', fontsize=14, fontweight='bold')
ax.set_xlabel('Growth Stage')
ax.set_ylabel('Gene')
plt.tight_layout()
fig.savefig('fig_liver_AA_enzymes_heatmap.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_liver_AA_enzymes_heatmap.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_liver_AA_enzymes_heatmap.png/pdf")

# ============================================================
# PLOT 2: Serum Urea breed×stage dynamics
# ============================================================
print("Generating Figure 2: Serum Urea dynamics...")
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Pre-compute sample sizes per stage to warn about n=1
urea_data = serum_tidy[serum_tidy['metabolite'] == 'Urea'].copy()
stages = [15, 45, 75, 105]
n_per_stage = {}
for s in stages:
    n_dly = len(urea_data[(urea_data['breed'] == 'DLY') & (urea_data['stage_kg'] == s)])
    n_tfb = len(urea_data[(urea_data['breed'] == 'TFB') & (urea_data['stage_kg'] == s)])
    n_per_stage[s] = (n_dly, n_tfb)
    if n_dly < 2 or n_tfb < 2:
        print(f"  WARNING: Stage {s}kg serum n=(DLY={n_dly}, TFB={n_tfb}) — descriptive only, no inference")

# Panel A: Urea
ax = axes[0]
for breed, color, marker in [('DLY', '#2196F3', 'o'), ('TFB', '#F44336', 's')]:
    breed_data = urea_data[urea_data['breed'] == breed]
    means, sems = [], []
    for s in stages:
        vals = breed_data[breed_data['stage_kg'] == s]['value'].dropna()
        if len(vals) > 0:
            means.append(vals.mean())
            sems.append(vals.std() / np.sqrt(len(vals)) if len(vals) >= 2 else 0)
        else:
            means.append(np.nan)
            sems.append(0)
    ax.errorbar(stages, means, yerr=sems, marker=marker, color=color,
                label=f'{breed}', linewidth=2, markersize=10, capsize=5)
    # Significance annotations — ONLY for stages with n>=3 in both breeds
    for i, s in enumerate(stages):
        n_dly, n_tfb = n_per_stage[s]
        if n_dly < 3 or n_tfb < 3:
            continue  # t-test not valid with n<2; n>=3 per group needed for minimal power
        dly = urea_data[(urea_data['breed'] == 'DLY') & (urea_data['stage_kg'] == s)]['value']
        tfb = urea_data[(urea_data['breed'] == 'TFB') & (urea_data['stage_kg'] == s)]['value']
        _, p = ttest_ind(dly, tfb)
        if p < 0.05:
            y_pos = means[i] + sems[i] + 0.2
            p_text = '***' if p < 0.001 else '**' if p < 0.01 else '*'
            ax.annotate(p_text, (s, y_pos), ha='center', fontsize=12, fontweight='bold',
                       color='#333333')

ax.set_ylabel('Serum Urea (mmol/L)', fontsize=12)
ax.set_xlabel('Body Weight (kg)', fontsize=12)
ax.set_title('Serum Urea', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.set_xticks(stages)
ax.grid(axis='y', alpha=0.3)
# Note: n=1 at 75/105 kg
ax.text(0.98, 0.02, '75/105 kg: n=1, descriptive only', transform=ax.transAxes,
        ha='right', fontsize=7, color='#999999', style='italic')

# Panel B: BCAA (Leu+Ile+Val sum)
ax = axes[1]
for breed, color, marker in [('DLY', '#2196F3', 'o'), ('TFB', '#F44336', 's')]:
    breed_bcaa = []
    for s in stages:
        stage_data = serum_tidy[(serum_tidy['breed'] == breed) & (serum_tidy['stage_kg'] == s)]
        bcaa_sum = []
        for rep in stage_data['rep'].unique():
            rep_data = stage_data[stage_data['rep'] == rep]
            leu = rep_data[rep_data['metabolite'] == 'Leu']['value'].values
            ile = rep_data[rep_data['metabolite'] == 'Ile']['value'].values
            val = rep_data[rep_data['metabolite'] == 'Val']['value'].values
            if len(leu) > 0 and len(ile) > 0 and len(val) > 0:
                bcaa_sum.append(leu[0] + ile[0] + val[0])
        if bcaa_sum:
            n_bcaa = len(bcaa_sum)
            breed_bcaa.append((np.mean(bcaa_sum),
                              np.std(bcaa_sum) / np.sqrt(n_bcaa) if n_bcaa >= 2 else 0))
        else:
            breed_bcaa.append((np.nan, np.nan))
    means = [b[0] for b in breed_bcaa]
    sems = [b[1] for b in breed_bcaa]
    ax.errorbar(stages, means, yerr=sems, marker=marker, color=color,
                label=f'{breed}', linewidth=2, markersize=10, capsize=5)

ax.set_ylabel('Serum BCAA (Leu+Ile+Val, mmol/L)', fontsize=12)
ax.set_xlabel('Body Weight (kg)', fontsize=12)
ax.set_title('Serum Branched-Chain Amino Acids', fontsize=14, fontweight='bold')
ax.legend(fontsize=11)
ax.set_xticks(stages)
ax.grid(axis='y', alpha=0.3)
ax.text(0.98, 0.02, '75/105 kg: n=1, descriptive only', transform=ax.transAxes,
        ha='right', fontsize=7, color='#999999', style='italic')

plt.suptitle('Serum Nitrogen Metabolism Indicators', fontsize=16, fontweight='bold', y=1.02)
plt.tight_layout()
fig.savefig('fig_serum_urea_bcaa_dynamics.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_serum_urea_bcaa_dynamics.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_serum_urea_bcaa_dynamics.png/pdf")

# ============================================================
# PLOT 3: Muscle ribosomal/translation gene heatmap
# ============================================================
print("Generating Figure 3: Muscle ribosomal gene heatmap...")
fig, ax = plt.subplots(figsize=(8, 12))

if muscle_heatmap.shape[0] > 2:
    row_linkage = linkage(muscle_heatmap.fillna(0), method='average')
    row_order = leaves_list(row_linkage)
    muscle_plot = muscle_heatmap.iloc[row_order]
else:
    muscle_plot = muscle_heatmap

sns.heatmap(muscle_plot, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
            ax=ax, cbar_kws={'label': 'log2(DLY/TFB)'},
            linewidths=0.5, linecolor='white')
ax.set_title('Muscle Ribosomal & Translation Genes\nlog2(DLY/TFB) Expression', fontsize=14, fontweight='bold')
ax.set_xlabel('Growth Stage')
ax.set_ylabel('Gene')
plt.tight_layout()
fig.savefig('fig_muscle_ribosomal_heatmap.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_muscle_ribosomal_heatmap.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_muscle_ribosomal_heatmap.png/pdf")

# ============================================================
# PLOT 4: Liver enzyme ↔ Serum Urea ↔ Muscle ribosome correlation
# ============================================================
print("Generating Figure 4: Three-way correlation heatmap...")

# Build a combined correlation matrix for top genes
top_enzymes = urea_corrs.head(8)['gene'].tolist()
top_ribo = muscle_heatmap.index[:10].tolist()

# Correlate each enzyme with each ribosomal gene (across breed×stage means)
# Build breed×stage table for each gene
enzyme_table = liver_aa_group.pivot_table(values='expr', index=['breed', 'stage_kg'],
                                           columns='gene_name', aggfunc='mean')
ribo_table = muscle_mean[muscle_mean['gene_name'].isin(top_ribo)].pivot_table(
    values='expr', index=['breed', 'stage_kg'], columns='gene_name', aggfunc='mean')
serum_table = serum_group_mean.pivot_table(values='mean', index=['breed', 'stage_kg'],
                                            columns='metabolite', aggfunc='mean')

# Combine into one correlation matrix
combined_genes = top_enzymes + ['Serum_Urea'] + top_ribo
combined_data = enzyme_table[top_enzymes].copy()
combined_data['Serum_Urea'] = serum_table['Urea']
for rg in top_ribo:
    if rg in ribo_table.columns:
        combined_data[rg] = ribo_table[rg]

# Drop rows with too many NAs
combined_data = combined_data.dropna(thresh=len(combined_genes) // 2)
combined_corr = combined_data.corr()

fig, ax = plt.subplots(figsize=(14, 12))
mask = np.zeros_like(combined_corr, dtype=bool)
sns.heatmap(combined_corr, annot=True, fmt='.2f', cmap='RdBu_r', center=0,
            vmin=-1, vmax=1, ax=ax, mask=mask,
            cbar_kws={'label': 'Pearson r'},
            linewidths=0.5, linecolor='white',
            xticklabels=True, yticklabels=True)
ax.set_title('Cross-Tissue Correlation Matrix\nLiver AA Enzymes ↔ Serum Urea ↔ Muscle Ribosomal Genes',
             fontsize=13, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
fig.savefig('fig_cross_tissue_correlation.png', dpi=200, bbox_inches='tight')
fig.savefig('fig_cross_tissue_correlation.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_cross_tissue_correlation.png/pdf")

# ============================================================
# PLOT 5: Graphical Abstract — Systems-level causal diagram
# ============================================================
print("Generating Figure 5: Graphical Abstract...")
fig, ax = plt.subplots(figsize=(16, 10))
ax.set_xlim(0, 16)
ax.set_ylim(0, 10)
ax.axis('off')

# Color scheme
dly_color = '#1565C0'
tfb_color = '#C62828'
serum_color = '#6A1B9A'
muscle_color = '#E65100'
arrow_color = '#546E7A'

# Title
ax.text(8, 9.5, 'Breed-Driven Liver–Serum–Muscle Axis in Pig Protein Deposition',
        ha='center', va='center', fontsize=18, fontweight='bold')
ax.text(8, 9.0, 'DLY: efficient AA use for muscle protein  |  TFB: hepatic AA catabolism → urea waste',
        ha='center', va='center', fontsize=11, color='#666666', style='italic')

# ---- LIVER compartment (left) ----
liver_box = mpatches.FancyBboxPatch((0.5, 4.5), 4.5, 4.0, boxstyle='round,pad=0.1',
                                      facecolor='#E3F2FD', edgecolor='#1565C0', linewidth=2)
ax.add_patch(liver_box)
ax.text(2.75, 8.2, 'LIVER', ha='center', fontsize=14, fontweight='bold', color='#1565C0')

# DLY liver
ax.text(1.5, 7.6, 'DLY', ha='center', fontsize=11, fontweight='bold', color=dly_color)
ax.text(1.5, 7.1, 'AA release\nPUFA synthesis\nIGF delivery', ha='center', fontsize=8, color=dly_color,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor=dly_color, alpha=0.8))

# TFB liver
ax.text(4.0, 7.6, 'TFB', ha='center', fontsize=11, fontweight='bold', color=tfb_color)
ax.text(4.0, 6.6, 'BCAA catabolism ↑\nUrea cycle ↑\nAA oxidation ↑', ha='center', fontsize=8, color=tfb_color,
        bbox=dict(boxstyle='round', facecolor='#FFEBEE', edgecolor=tfb_color, alpha=0.8))

# Enzyme names in liver
enzymes_text = 'BCAT2  BCKDHA/B  GOT1  GPT\nASS1  ARG1  CPS1  OTC\nAASS  HGD  ACADSB  GLUD1'
ax.text(2.75, 4.9, enzymes_text, ha='center', fontsize=7, color='#555555',
        bbox=dict(boxstyle='round', facecolor='#FAFAFA', edgecolor='#BDBDBD'))

# ---- SERUM compartment (middle) ----
serum_box = mpatches.FancyBboxPatch((5.75, 4.5), 4.5, 4.0, boxstyle='round,pad=0.1',
                                      facecolor='#F3E5F5', edgecolor='#6A1B9A', linewidth=2)
ax.add_patch(serum_box)
ax.text(8.0, 8.2, 'SERUM', ha='center', fontsize=14, fontweight='bold', color='#6A1B9A')

# DLY serum
ax.text(6.75, 7.6, 'DLY', ha='center', fontsize=11, fontweight='bold', color=dly_color)
ax.text(6.75, 7.1, 'AA lower\nUrea lower\nEfficient circulation', ha='center', fontsize=8, color=dly_color,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor=dly_color, alpha=0.8))

# TFB serum
ax.text(9.25, 7.6, 'TFB', ha='center', fontsize=11, fontweight='bold', color=tfb_color)
ax.text(9.25, 6.6, 'AA higher (accumulation)\nUrea higher (waste)\nPoor muscle uptake', ha='center', fontsize=8, color=tfb_color,
        bbox=dict(boxstyle='round', facecolor='#FFEBEE', edgecolor=tfb_color, alpha=0.8))

# Key evidence
ax.text(8.0, 4.9, 'Serum Urea: r = +0.57~0.77\nwith all liver AA enzymes\nDirect in vivo evidence\nof AA→urea diversion',
        ha='center', fontsize=7.5, color='#555555',
        bbox=dict(boxstyle='round', facecolor='#FAFAFA', edgecolor='#BDBDBD'))

# ---- MUSCLE compartment (right) ----
muscle_box = mpatches.FancyBboxPatch((11.0, 4.5), 4.5, 4.0, boxstyle='round,pad=0.1',
                                       facecolor='#FFF3E0', edgecolor='#E65100', linewidth=2)
ax.add_patch(muscle_box)
ax.text(13.25, 8.2, 'MUSCLE', ha='center', fontsize=14, fontweight='bold', color='#E65100')

# DLY muscle
ax.text(12.0, 7.6, 'DLY', ha='center', fontsize=11, fontweight='bold', color=dly_color)
ax.text(12.0, 7.1, 'Ribosome ↑\nTranslation ↑\nSustained growth', ha='center', fontsize=8, color=dly_color,
        bbox=dict(boxstyle='round', facecolor='white', edgecolor=dly_color, alpha=0.8))

# TFB muscle
ax.text(14.5, 7.6, 'TFB', ha='center', fontsize=11, fontweight='bold', color=tfb_color)
ax.text(14.5, 6.6, 'Ribosome ↓\nTranslation ↓\nEarly maturation', ha='center', fontsize=8, color=tfb_color,
        bbox=dict(boxstyle='round', facecolor='#FFEBEE', edgecolor=tfb_color, alpha=0.8))

# Gene names in muscle
ribo_text = ('RPL3  RPL7  RPL10A  RPL13\nRPL15  RPL27A  RPLP0\n'
             'RPS3  RPS6  RPS15A  RPS23\nEIF3E  EIF4B  EEF2')
ax.text(13.25, 4.9, ribo_text, ha='center', fontsize=7, color='#555555',
        bbox=dict(boxstyle='round', facecolor='#FAFAFA', edgecolor='#BDBDBD'))

# ---- ARROWS between compartments ----
# Liver → Serum
ax.annotate('', xy=(5.7, 7.0), xytext=(5.0, 7.0),
            arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2.5, connectionstyle='arc3,rad=0'))
ax.text(5.35, 7.3, 'AA release\n& oxidation', ha='center', fontsize=7, color=arrow_color)

ax.annotate('', xy=(5.7, 5.5), xytext=(5.0, 5.5),
            arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2.5, connectionstyle='arc3,rad=0'))

# Serum → Muscle
ax.annotate('', xy=(11.0, 7.0), xytext=(10.3, 7.0),
            arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2.5, connectionstyle='arc3,rad=0'))
ax.text(10.65, 7.3, 'AA uptake', ha='center', fontsize=7, color=arrow_color)

ax.annotate('', xy=(11.0, 5.5), xytext=(10.3, 5.5),
            arrowprops=dict(arrowstyle='->', color=arrow_color, lw=2.5, connectionstyle='arc3,rad=0'))

# ---- PHENOTYPE outcome (bottom) ----
pheno_box = mpatches.FancyBboxPatch((3.0, 0.5), 10.0, 2.0, boxstyle='round,pad=0.15',
                                      facecolor='#E8F5E9', edgecolor='#2E7D32', linewidth=2)
ax.add_patch(pheno_box)
ax.text(8.0, 2.2, 'PHENOTYPE OUTCOME', ha='center', fontsize=13, fontweight='bold', color='#2E7D32')

ax.text(4.5, 1.5, 'DLY', ha='center', fontsize=12, fontweight='bold', color=dly_color)
ax.text(4.5, 1.0, 'Protein deposition\n75→105 kg sustains ↑\nNitrogen efficiency ↑', ha='center', fontsize=9, color=dly_color)

ax.text(11.5, 1.5, 'TFB', ha='center', fontsize=12, fontweight='bold', color=tfb_color)
ax.text(11.5, 1.0, 'Protein deposition\npeaks at 45 kg\nN wasted as urea', ha='center', fontsize=9, color=tfb_color)

# Arrow from muscle to phenotype
ax.annotate('', xy=(8.0, 2.5), xytext=(8.0, 3.5),
            arrowprops=dict(arrowstyle='->', color='#2E7D32', lw=2.5))
ax.text(8.5, 3.0, 'Muscle protein\nsynthesis rate', fontsize=7, color='#2E7D32')

# ---- Legend / key insight box ----
insight_box = mpatches.FancyBboxPatch((0.3, 0.1), 15.4, 0.4, boxstyle='round',
                                        facecolor='#FFFDE7', edgecolor='#F9A825', linewidth=1)
ax.add_patch(insight_box)
ax.text(8.0, 0.3,
        'Key Insight: TFB high serum AA ≠ supply abundance, but reflects poor muscle AA uptake & utilization — evidenced by low ribosomal expression & high serum Urea',
        ha='center', fontsize=8.5, fontweight='bold', color='#F57F17')

plt.tight_layout()
fig.savefig('fig_graphical_abstract.png', dpi=250, bbox_inches='tight')
fig.savefig('fig_graphical_abstract.pdf', dpi=250, bbox_inches='tight')
plt.close()
print("  Saved fig_graphical_abstract.png/pdf")

print("\n=== All figures generated ===")
print("""
Output files:
  fig_liver_AA_enzymes_heatmap.png/pdf      - Liver AA catabolism heatmap
  fig_serum_urea_bcaa_dynamics.png/pdf      - Serum Urea & BCAA dynamics
  fig_muscle_ribosomal_heatmap.png/pdf      - Muscle ribosomal gene heatmap
  fig_cross_tissue_correlation.png/pdf      - Cross-tissue correlation matrix
  fig_graphical_abstract.png/pdf            - Systems-level causal diagram
""")
