#!/usr/bin/env python3
"""
AA肝肌轴候选基因筛选：Venn图交集 + KEGG通路富集
逐层筛选逻辑：
  Layer 1: 81个肝脏AA通路基因 → 按PD窗口对齐模式筛选
  Layer 2: 与肌肉AA受体显著跨组织相关 (P<0.05)
  Layer 3: 与肌肉蛋白降解基因强相关 (|r|>0.7 with TRIM63/FBXO32)
  Layer 4: KEGG通路富集验证
  Final: 高置信度肝肌轴候选基因
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
from matplotlib_venn import venn3, venn2
from scipy.stats import fisher_exact
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 8,
})

# ============================================================
# 1. Load data
# ============================================================
aa_liver = pd.read_csv('/Users/hezongze/Downloads/liver_arachidonic_pathway_genes.csv')
corr = pd.read_excel('/Users/hezongze/pig_study/aa_crosstissue_full_results.xlsx')

print(f"Liver AA genes: {len(aa_liver)}")
print(f"Cross-tissue pairs: {len(corr)}")
print(f"Liver categories: {aa_liver['Category'].unique().tolist()}")
print(f"PD patterns: {aa_liver['Pattern'].unique().tolist()}")

# ============================================================
# 2. Define filtering sets for Venn
# ============================================================

# --- Set A: PD-window-aligned liver AA genes ---
# "TFB45与DLY75共同偏高" or "DLY75偏高，TFB45不突出"
pd_aligned = aa_liver[aa_liver['Pattern'].isin(
    ['TFB45与DLY75共同偏高', 'DLY75偏高，TFB45不突出']
)]
set_a = set(pd_aligned['Gene Name'].unique())
print(f"\nSet A (PD-aligned): {len(set_a)} genes")
for g in sorted(set_a):
    pat = aa_liver[aa_liver['Gene Name']==g]['Pattern'].values[0]
    cat = aa_liver[aa_liver['Gene Name']==g]['Category'].values[0]
    print(f"  {g} [{cat}] ({pat})")

# --- Set B: Liver genes with ANY significant cross-tissue correlation (P<0.05) ---
sig_corr = corr[corr['P_value'] < 0.05]
set_b = set(sig_corr['Liver_Gene'].unique())
print(f"\nSet B (P<0.05 cross-tissue): {len(set_b)} genes")

# --- Set C: Liver genes strongly correlated with muscle proteolysis genes ---
# TRIM63/MuRF1 and FBXO32/Atrogin-1 are the key muscle protein degradation markers
proteolysis_markers = ['TRIM63', 'FBXO32']
strong_proteolysis = corr[
    (corr['Muscle_Gene'].isin(proteolysis_markers)) &
    (corr['abs_r'] > 0.7) &
    (corr['P_value'] < 0.05)
]
set_c = set(strong_proteolysis['Liver_Gene'].unique())
print(f"\nSet C (|r|>0.7 with TRIM63/FBXO32): {len(set_c)} genes")
for _, row in strong_proteolysis.sort_values('abs_r', ascending=False).iterrows():
    print(f"  {row['Liver_Gene']} ↔ {row['Muscle_Gene']}: r={row['Pearson_r']:.3f}, P={row['P_value']:.4f}")

# --- Set D: Liver genes with strong muscle correlation to any muscle gene (|r|>0.75) ---
strong_any = sig_corr[sig_corr['abs_r'] > 0.75]
set_d = set(strong_any['Liver_Gene'].unique())
print(f"\nSet D (|r|>0.75, P<0.05): {len(set_d)} genes")

# ============================================================
# 3. Venn Diagram 1: A ∩ B ∩ C (PD-aligned × Significant × Proteolysis)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

# --- Panel a: A ∩ B ---
ax = axes[0]
only_a = set_a - set_b
only_b = set_b - set_a
both_ab = set_a & set_b

v = venn2([set_a, set_b], set_labels=('PD-Aligned\nLiver AA Genes', 'Significant Cross-Tissue\nCorrelation (P<0.05)'),
          set_colors=('#E64B35', '#4DBBD5'), alpha=0.6, ax=ax)
ax.set_title('A. PD Window × Cross-Tissue Sig', fontsize=9, fontweight='bold', loc='left')

# Annotate intersection
if v.get_label_by_id('11'):
    v.get_label_by_id('11').set_fontsize(9)
    v.get_label_by_id('11').set_fontweight('bold')

# --- Panel b: (A ∩ B) ∩ C ---
ax = axes[1]
ab_union = set_a & set_b
only_ab = ab_union - set_c
only_c = set_c - ab_union
both_abc = ab_union & set_c

v2 = venn2([ab_union, set_c],
           set_labels=('PD-Aligned × Sig\n(A ∩ B)', 'Strong Proteolysis\n(|r|>0.7, TRIM63/FBXO32)'),
           set_colors=('#F39B7F', '#DC0000'), alpha=0.6, ax=ax)
ax.set_title('B. (A ∩ B) × Proteolysis Link', fontsize=9, fontweight='bold', loc='left')

if v2.get_label_by_id('11'):
    v2.get_label_by_id('11').set_fontsize(9)
    v2.get_label_by_id('11').set_fontweight('bold')

# --- Panel c: 3-way Venn: A ∩ B ∩ C ---
ax = axes[2]
v3 = venn3([set_a, set_b, set_c],
           set_labels=('PD-Aligned', 'Sig\n(P<0.05)', 'Proteolysis\n(|r|>0.7)'),
           set_colors=('#E64B35', '#4DBBD5', '#DC0000'), alpha=0.5, ax=ax)
ax.set_title('C. Three-Way Intersection', fontsize=9, fontweight='bold', loc='left')

# Highlight the triple intersection
if v3.get_label_by_id('111'):
    v3.get_label_by_id('111').set_fontsize(10)
    v3.get_label_by_id('111').set_fontweight('bold')
    v3.get_label_by_id('111').set_color('#DC0000')

plt.tight_layout()
plt.savefig('/Users/hezongze/pig_study/fig_Venn1_gene_screening.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.savefig('/Users/hezongze/pig_study/fig_Venn1_gene_screening.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("\nVenn Fig 1 saved.")

# ============================================================
# 4. Tiered candidate gene list
# ============================================================

# Tier 1: A ∩ B ∩ C (PD-aligned + sig + strong proteolysis link)
tier1 = set_a & set_b & set_c
print(f"\n{'='*60}")
print(f"TIER 1 (A ∩ B ∩ C): {len(tier1)} genes — Highest Priority")
print(f"{'='*60}")
tier1_details = []
for g in sorted(tier1):
    cat = aa_liver[aa_liver['Gene Name']==g]['Category'].values[0]
    pat = aa_liver[aa_liver['Gene Name']==g]['Pattern'].values[0]
    # Get best correlation
    g_corr = corr[(corr['Liver_Gene']==g) & (corr['Muscle_Gene'].isin(proteolysis_markers))]
    best = g_corr.sort_values('abs_r', ascending=False).iloc[0]
    tier1_details.append({
        'Gene': g, 'Category': cat, 'Pattern': pat,
        'Best_Muscle': best['Muscle_Gene'], 'r': best['Pearson_r'], 'P': best['P_value']
    })
    print(f"  {g} [{cat}] ↔ {best['Muscle_Gene']}: r={best['Pearson_r']:.3f}, P={best['P_value']:.4f}")

tier1_df = pd.DataFrame(tier1_details)

# Tier 2: A ∩ B (PD-aligned + sig, but without strong proteolysis link)
tier2 = (set_a & set_b) - set_c
print(f"\n{'='*60}")
print(f"TIER 2 (A ∩ B, not C): {len(tier2)} genes — Secondary Priority")
print(f"{'='*60}")
for g in sorted(tier2):
    cat = aa_liver[aa_liver['Gene Name']==g]['Category'].values[0]
    pat = aa_liver[aa_liver['Gene Name']==g]['Pattern'].values[0]
    g_corr = sig_corr[sig_corr['Liver_Gene']==g].sort_values('abs_r', ascending=False)
    best = g_corr.iloc[0]
    print(f"  {g} [{cat}] ↔ {best['Muscle_Gene']}: r={best['Pearson_r']:.3f}, P={best['P_value']:.4f}")

# ============================================================
# 5. Venn Diagram 2: Category-level intersections
# ============================================================
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5.5))

# Group liver AA genes by category
categories = aa_liver['Category'].unique()

# --- Panel a: Category composition of each tier ---
ax = axes2[0]
cat_tier1 = aa_liver[aa_liver['Gene Name'].isin(tier1)]['Category'].value_counts()
cat_tier2 = aa_liver[aa_liver['Gene Name'].isin(tier2)]['Category'].value_counts()
cat_all = aa_liver['Category'].value_counts()

x = np.arange(len(categories))
w = 0.25
ax.bar(x - w, [cat_all.get(c, 0) for c in categories], w, label='All 81 AA genes', color='#CCCCCC', edgecolor='white')
ax.bar(x, [cat_tier2.get(c, 0) for c in categories], w, label='Tier 2 (A∩B)', color='#4DBBD5', edgecolor='white')
ax.bar(x + w, [cat_tier1.get(c, 0) for c in categories], w, label='Tier 1 (A∩B∩C)', color='#DC0000', edgecolor='white')

ax.set_xticks(x)
ax.set_xticklabels([c.replace('：','\n') for c in categories], rotation=45, ha='right', fontsize=6)
ax.set_ylabel('Gene Count', fontsize=7)
ax.set_title('D. Category Enrichment Across Tiers', fontsize=9, fontweight='bold', loc='left')
ax.legend(fontsize=6, frameon=False)
ax.spines[['top','right']].set_visible(False)

# --- Panel b: Heatmap of Tier1 genes correlations with key muscle genes ---
ax = axes2[1]
key_muscle = ['TRIM63', 'FBXO32', 'FOXO3', 'PPARA', 'LTB4R', 'IL6R', 'PTGER4']

tier1_genes = sorted(tier1)
heatmap_data = np.zeros((len(tier1_genes), len(key_muscle)))
annot_data = []
for i, lg in enumerate(tier1_genes):
    row_annot = []
    for j, mg in enumerate(key_muscle):
        row_match = corr[(corr['Liver_Gene']==lg) & (corr['Muscle_Gene']==mg)]
        if len(row_match) > 0:
            heatmap_data[i, j] = row_match.iloc[0]['Pearson_r']
            p = row_match.iloc[0]['P_value']
            row_annot.append(f"{row_match.iloc[0]['Pearson_r']:.2f}\n{p:.1e}" if p < 0.01 else f"{row_match.iloc[0]['Pearson_r']:.2f}\n{p:.2f}")
        else:
            heatmap_data[i, j] = np.nan
            row_annot.append('')
    annot_data.append(row_annot)

im = ax.imshow(heatmap_data, cmap='RdBu_r', aspect='auto', vmin=-1, vmax=1)
for i in range(len(tier1_genes)):
    for j in range(len(key_muscle)):
        if not np.isnan(heatmap_data[i, j]):
            val = heatmap_data[i, j]
            tc = 'white' if abs(val) > 0.6 else 'black'
            ax.text(j, i, f'{val:.2f}', ha='center', va='center', fontsize=5.5, color=tc)

ax.set_xticks(range(len(key_muscle)))
ax.set_xticklabels(key_muscle, rotation=45, ha='right', fontsize=6.5)
ax.set_yticks(range(len(tier1_genes)))
ax.set_yticklabels(tier1_genes, fontsize=6.5)
ax.set_title('E. Tier 1 Liver Genes × Key Muscle Genes', fontsize=9, fontweight='bold', loc='left')

cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
cbar.set_label('Pearson r', fontsize=6, labelpad=1)
cbar.ax.tick_params(labelsize=5, width=0.5)

plt.tight_layout()
fig2.savefig('/Users/hezongze/pig_study/fig_Venn2_category_heatmap.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig2.savefig('/Users/hezongze/pig_study/fig_Venn2_category_heatmap.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Venn Fig 2 saved.")

# ============================================================
# 6. KEGG Pathway Analysis
# ============================================================

print(f"\n{'='*60}")
print("KEGG Pathway Enrichment Analysis")
print(f"{'='*60}")

# Gene-to-KEGG pathway mapping for AA metabolism relevant pathways
# Using known pig KEGG pathway annotations
kegg_pathways = {
    'Arachidonic acid metabolism': {
        'ssc00590': ['PLA2G1B', 'PLA2G2A', 'PLA2G4A', 'PLA2G6', 'PLA2G7',
                      'PTGS1', 'PTGS2', 'PTGDS', 'PTGES', 'PTGES2', 'PTGES3',
                      'CBR1', 'CBR2', 'CBR3', 'CBR4',
                      'ALOX5', 'ALOX5AP', 'ALOX12', 'ALOX12B', 'ALOX15', 'ALOX15B',
                      'LTC4S', 'GPX1', 'GPX2', 'GPX3', 'GPX4', 'GPX7', 'GPX8',
                      'CYP2C', 'CYP2E1', 'CYP2J2', 'CYP4A', 'CYP4F',
                      'EPHX1', 'EPHX2', 'PTGER1', 'PTGER2', 'PTGER3', 'PTGER4',
                      'PTGIR', 'TBXA2R', 'CYSLTR1', 'CYSLTR2', 'LTB4R', 'PTGFR'],
    },
    'Linoleic acid metabolism': {
        'ssc00591': ['FADS1', 'FADS2', 'FADS6', 'ELOVL2', 'ELOVL5',
                      'PLA2G1B', 'PLA2G2A', 'PLA2G4A', 'PLA2G6', 'CYP2E1', 'CYP2C'],
    },
    'PPAR signaling pathway': {
        'ssc03320': ['PPARA', 'PPARD', 'PPARG', 'PPARGC1A', 'FABP1', 'FABP2', 'FABP3',
                      'FABP4', 'FABP5', 'FABP6', 'CD36', 'ACSL1', 'ACSL3', 'ACSL4',
                      'ACSL5', 'ACSL6', 'SCD', 'CPT1A', 'CPT1B', 'CPT2',
                      'ACOX1', 'ACOX2', 'ACOX3', 'EHHADH', 'HMGCS2', 'CYP4A',
                      'CYP8B1', 'ME1', 'UCP1', 'ADIPOQ', 'ILK', 'PDK1', 'PCK1', 'GK'],
    },
    'Fatty acid metabolism': {
        'ssc01212': ['ACSL1', 'ACSL3', 'ACSL4', 'ACSL5', 'ACSL6',
                      'CPT1A', 'CPT1B', 'CPT2', 'ACOX1', 'ACOX2', 'ACOX3',
                      'EHHADH', 'HADHA', 'HADHB', 'ACAA1', 'ACAA2',
                      'FADS1', 'FADS2', 'FADS6', 'ELOVL2', 'ELOVL5', 'ELOVL6',
                      'SCD', 'SCD5', 'FASN', 'ACACA', 'ACACB'],
    },
    'Fatty acid degradation': {
        'ssc00071': ['ACSL1', 'ACSL3', 'ACSL4', 'ACSL5', 'ACSL6',
                      'CPT1A', 'CPT1B', 'CPT2', 'ACOX1', 'ACOX3',
                      'EHHADH', 'HADHA', 'HADHB', 'ACAA1', 'ACAA2'],
    },
    'Biosynthesis of unsaturated fatty acids': {
        'ssc01040': ['FADS1', 'FADS2', 'FADS6', 'ELOVL2', 'ELOVL5', 'ELOVL6',
                      'SCD', 'SCD5', 'HACD1', 'HACD2', 'ACOT1', 'ACOT2', 'TECR'],
    },
    'Glycerophospholipid metabolism': {
        'ssc00564': ['PLA2G1B', 'PLA2G2A', 'PLA2G4A', 'PLA2G6', 'PLA2G7',
                      'PLD1', 'PLD2', 'PLD3', 'PLD4',
                      'LPCAT1', 'LPCAT2', 'LPCAT3', 'LPCAT4',
                      'CEPT1', 'CHPT1', 'PEMT', 'SELENOI'],
    },
    'Ubiquitin mediated proteolysis': {
        'ssc04120': ['FBXO32', 'TRIM63', 'MURF1', 'FBXO30', 'FBXO8',
                      'UBE2D1', 'UBE2D2', 'UBE2D3', 'CUL3', 'RBX1'],
    },
    'PI3K-Akt signaling pathway': {
        'ssc04151': ['FOXO1', 'FOXO3', 'MTOR', 'RPS6KB1', 'MYOD1', 'MYOG',
                      'IL6R', 'IGF1R', 'INSR', 'IRS1', 'PIK3CA', 'AKT1', 'AKT2'],
    },
    'cAMP signaling pathway': {
        'ssc04024': ['PTGER2', 'PTGER4', 'PTGIR', 'PTGER1', 'PTGER3',
                      'PTGFR', 'ADORA1', 'ADORA2A',
                      'CREB1', 'CREB3', 'CREB5', 'ATF2', 'ATF4', 'ATF6B'],
    },
    'Calcium signaling pathway': {
        'ssc04020': ['PTGER1', 'PTGER3', 'TBXA2R', 'CYSLTR1', 'CYSLTR2',
                      'LTB4R', 'PTGFR', 'EGFR', 'PDGFRA', 'PDGFRB',
                      'CALM1', 'CALM2', 'CALM3', 'CAMK2A', 'CAMK2D', 'CAMK2G'],
    },
    'cGMP-PKG signaling pathway': {
        'ssc04022': ['ADORA1', 'ADORA2A', 'ADORA2B', 'PTGIR',
                      'CREB1', 'CREB3', 'CREB5', 'ATF2', 'ATF4'],
    },
    'Insulin signaling pathway': {
        'ssc04910': ['FOXO1', 'FOXO3', 'MTOR', 'RPS6KB1', 'INSR', 'IGF1R',
                      'IRS1', 'PIK3CA', 'AKT1', 'AKT2', 'SREBF1', 'PPARGC1A'],
    },
}

# Map our gene symbols to KEGG pathways
def map_genes_to_kegg(gene_set, kegg_dict):
    """Map a set of gene symbols to KEGG pathways"""
    result = {}
    for pathway_name, pathway_info in kegg_dict.items():
        for kegg_id, kegg_genes in pathway_info.items():
            overlap = gene_set & set(kegg_genes)
            if overlap:
                if pathway_name not in result:
                    result[pathway_name] = {'kegg_id': kegg_id, 'genes': set(), 'all_kegg_genes': set(kegg_genes)}
                result[pathway_name]['genes'].update(overlap)
    return result

# Map all 81 liver AA genes
all_liver_genes = set(aa_liver['Gene Name'].unique())
all_muscle_genes = set(corr['Muscle_Gene'].unique())

# Tier 1 genes (our key candidates)
tier1_kegg = map_genes_to_kegg(tier1, kegg_pathways)

print(f"\n--- KEGG Pathways for Tier 1 Genes ({len(tier1)} genes: {', '.join(sorted(tier1))}) ---")
for pname, pinfo in sorted(tier1_kegg.items()):
    print(f"  {pname} ({pinfo['kegg_id']}): {', '.join(sorted(pinfo['genes']))}")

# Full PD-aligned gene set (A ∩ B)
ab_genes = set_a & set_b
ab_kegg = map_genes_to_kegg(ab_genes, kegg_pathways)
print(f"\n--- KEGG Pathways for A ∩ B Genes ({len(ab_genes)} genes) ---")
for pname, pinfo in sorted(ab_kegg.items()):
    print(f"  {pname} ({pinfo['kegg_id']}): {len(pinfo['genes'])} genes — {', '.join(sorted(pinfo['genes']))}")

# Liver AA pathway specific genes (from our 81)
liver_aa_kegg = map_genes_to_kegg(all_liver_genes, kegg_pathways)
print(f"\n--- KEGG Pathways for All 81 Liver AA Genes ---")
for pname, pinfo in sorted(liver_aa_kegg.items()):
    print(f"  {pname} ({pinfo['kegg_id']}): {len(pinfo['genes'])} genes — {', '.join(sorted(list(pinfo['genes'])[:8]))}{'...' if len(pinfo['genes'])>8 else ''}")

# ============================================================
# 7. KEGG Enrichment Bar Plot
# ============================================================

# Count genes per KEGG pathway for each tier
def calc_enrichment(gene_set, kegg_dict, bg_size=20000):
    """Calculate enrichment stats for a gene set against KEGG pathways"""
    results = []
    set_size = len(gene_set)
    for pathway_name, pathway_info in kegg_dict.items():
        for kegg_id, kegg_genes in pathway_info.items():
            overlap = gene_set & set(kegg_genes)
            k = len(overlap)
            if k == 0:
                continue
            K = len(kegg_genes)  # Total genes in KEGG pathway
            # Fisher's exact test
            # Contingency table:
            #               In pathway   Not in pathway
            # In gene set       k           set_size - k
            # Not in set      K - k        bg - set_size - K + k
            bg = bg_size
            table = [[k, K - k], [set_size - k, bg - set_size - K + k]]
            odds_ratio, p_value = fisher_exact(table, alternative='greater')
            results.append({
                'Pathway': pathway_name,
                'KEGG_ID': kegg_id,
                'Overlap': k,
                'Pathway_Size': K,
                'Set_Size': set_size,
                'P_value': p_value,
                'Genes': ', '.join(sorted(overlap)),
            })
    return pd.DataFrame(results).sort_values('P_value')

# Run enrichment for each tier
tier1_enrich = calc_enrichment(tier1, kegg_pathways) if len(tier1) > 0 else pd.DataFrame()
ab_enrich = calc_enrichment(ab_genes, kegg_pathways)
all_enrich = calc_enrichment(all_liver_genes, kegg_pathways)

print(f"\n--- KEGG Enrichment (Fisher's exact test) ---")
if len(tier1_enrich) > 0:
    print(f"\nTier 1 significant pathways:")
    for _, row in tier1_enrich.iterrows():
        print(f"  {row['Pathway']}: {row['Overlap']}/{row['Pathway_Size']} genes, P={row['P_value']:.4f}")
        print(f"    Genes: {row['Genes']}")

print(f"\nA ∩ B significant pathways:")
for _, row in ab_enrich.head(10).iterrows():
    print(f"  {row['Pathway']}: {row['Overlap']}/{row['Pathway_Size']} genes, P={row['P_value']:.4f}")

# ============================================================
# 8. KEGG Enrichment Figure
# ============================================================
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 5.5))

# --- Panel a: Top KEGG pathways for A∩B genes ---
ax = axes3[0]
top_paths = ab_enrich.head(10)
colors = plt.cm.Reds_r(np.linspace(0.3, 0.9, len(top_paths)))
bars = ax.barh(range(len(top_paths)), -np.log10(top_paths['P_value'].values), color=colors, edgecolor='white', height=0.7)

ax.set_yticks(range(len(top_paths)))
ax.set_yticklabels([f"{r['Pathway']}\n({r['Overlap']}/{r['Pathway_Size']} genes)" for _, r in top_paths.iterrows()],
                    fontsize=6)
ax.set_xlabel('-log10(P-value)', fontsize=7)
ax.axvline(x=-np.log10(0.05), color='#CCCCCC', linestyle='--', lw=0.8, label='P=0.05')
ax.set_title('F. KEGG Pathway Enrichment (A ∩ B Genes)', fontsize=9, fontweight='bold', loc='left')
ax.spines[['top','right']].set_visible(False)

# Add gene labels
for i, (_, row) in enumerate(top_paths.iterrows()):
    ax.text(0.1, i, row['Genes'][:60], fontsize=5, va='center', ha='left', color='#333333')

# --- Panel b: Final Candidate Summary ---
ax = axes3[1]
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# Title
ax.text(0, 9.8, 'FINAL CANDIDATE GENES — AA Liver-Muscle Axis',
        fontsize=11, fontweight='bold', color='#B71C1C')

# Tier 1 genes with details
y = 9.0
for _, row in tier1_df.iterrows():
    # Gene symbol as header
    ax.text(0.2, y, f"{row['Gene']}", fontsize=10, fontweight='bold', color='#DC0000')
    ax.text(3.0, y, f"↔ {row['Best_Muscle']}", fontsize=9, color='#333333')
    ax.text(5.5, y, f"r={row['r']:.3f}", fontsize=8, fontweight='bold', color='#E64B35')
    ax.text(7.0, y, f"P={row['P']:.4f}", fontsize=8, color='#888888')
    ax.text(0.4, y - 0.35, f"[{row['Category']}]", fontsize=6.5, color='#666666', style='italic')
    ax.text(0.4, y - 0.65, f"Pattern: {row['Pattern']}", fontsize=6, color='#999999')
    y -= 1.2

# Screening summary
y -= 0.3
ax.text(0, y, f'Screening Pipeline: 81 AA genes → {len(set_a)} PD-aligned → {len(set_a & set_b)} sig cross-tissue → {len(tier1)} final',
        fontsize=7, color='#555555', style='italic')

# KEGG summary
y -= 0.6
if len(tier1_enrich) > 0:
    top_kegg_names = ', '.join(tier1_enrich.head(3)['Pathway'].tolist())
    ax.text(0, y, f'Top KEGG Pathways: {top_kegg_names}', fontsize=7, color='#555555')

plt.tight_layout()
fig3.savefig('/Users/hezongze/pig_study/fig_Venn3_KEGG_enrichment.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig3.savefig('/Users/hezongze/pig_study/fig_Venn3_KEGG_enrichment.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("\nKEGG Figure saved.")

# ============================================================
# 9. Save final candidate list
# ============================================================
output = []
# Tier 1
for _, row in tier1_df.iterrows():
    output.append({
        'Tier': 1,
        'Gene': row['Gene'],
        'Category': row['Category'],
        'Pattern': row['Pattern'],
        'Best_Muscle_Partner': row['Best_Muscle'],
        'Pearson_r': row['r'],
        'P_value': row['P'],
        'Rationale': 'PD-aligned + significant cross-tissue + strong proteolysis link'
    })

# Tier 2
for g in sorted(tier2):
    cat = aa_liver[aa_liver['Gene Name']==g]['Category'].values[0]
    pat = aa_liver[aa_liver['Gene Name']==g]['Pattern'].values[0]
    g_corr = sig_corr[sig_corr['Liver_Gene']==g].sort_values('abs_r', ascending=False)
    best = g_corr.iloc[0]
    output.append({
        'Tier': 2,
        'Gene': g,
        'Category': cat,
        'Pattern': pat,
        'Best_Muscle_Partner': best['Muscle_Gene'],
        'Pearson_r': best['Pearson_r'],
        'P_value': best['P_value'],
        'Rationale': 'PD-aligned + significant cross-tissue (not proteolysis-specific)'
    })

output_df = pd.DataFrame(output)
output_df.to_csv('/Users/hezongze/pig_study/aa_candidate_genes_screened.csv', index=False)
print(f"\nFinal candidate list saved: aa_candidate_genes_screened.csv ({len(output_df)} genes)")

print("\n=== Screening Complete ===")
print(f"Pipeline: 81 → {len(set_a)} (PD-aligned) → {len(set_a & set_b)} (sig cross-tissue) → {len(tier1)} final")
