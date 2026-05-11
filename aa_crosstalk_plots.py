#!/usr/bin/env python3
"""AA代谢肝肌轴——可视化"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 11,
    'axes.labelsize': 9,
    'figure.dpi': 150,
})

# Load detail data
detail = pd.read_csv('/Users/hezongze/pig_study/aa_crosstalk_detail_v2.csv')
corr = pd.read_excel('/Users/hezongze/pig_study/aa_crosstissue_full_results.xlsx')

# ============================================================
# Figure 1: Muscle AA receptor expression heatmap
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={'width_ratios': [1, 1.3]})

# 1A: Muscle AA receptor Z-score heatmap
muscle_genes = ['PTGER2', 'PTGER4', 'PTGIR', 'TBXA2R', 'CYSLTR1', 'CYSLTR2', 'LTB4R',
                'PPARA', 'PPARD', 'FABP3', 'IL6R', 'FOXO1', 'FOXO3', 'FBXO32', 'TRIM63']
detail_m = detail[detail['Muscle_Gene'].isin(muscle_genes)].copy()

# Pivot to matrix
mp = detail_m.pivot_table(index='Muscle_Gene', columns=['Breed', 'Stage'], values='Muscle_Mean_Expr', aggfunc='mean')
# Z-score within gene
mp_z = mp.subtract(mp.mean(axis=1), axis=0).divide(mp.std(axis=1), axis=0)

# Reorder
gene_order = ['PTGER2', 'PTGER4', 'PTGIR', 'TBXA2R', 'CYSLTR1', 'CYSLTR2', 'LTB4R',
              'PPARA', 'PPARD', 'FABP3', 'IL6R', 'FOXO1', 'FOXO3', 'FBXO32', 'TRIM63']
mp_z = mp_z.reindex([g for g in gene_order if g in mp_z.index])

sns.heatmap(mp_z, ax=ax1, cmap='RdBu_r', center=0, annot=True, fmt='.1f',
            linewidths=0.5, cbar_kws={'label': 'Z-score', 'shrink': 0.7})
ax1.set_title('A. Muscle AA Receptor / Effector Expression Z-score', fontweight='bold', loc='left')

# 1B: Top cross-tissue correlations
sig_corr = corr[corr['P_value'] < 0.05].copy()
top_corr = sig_corr.nlargest(15, 'abs_r')

# Group by category
categories = top_corr['Liver_Category'].unique()
colors = plt.cm.tab10(np.linspace(0, 1, len(categories)))
cat_colors = dict(zip(categories, colors))

y_positions = range(len(top_corr))
bars = ax2.barh(y_positions, top_corr['Pearson_r'].values,
                color=[cat_colors[c] for c in top_corr['Liver_Category']])

# Label
labels = [f"{r['Liver_Gene']} → {r['Muscle_Gene']}" for _, r in top_corr.iterrows()]
ax2.set_yticks(y_positions)
ax2.set_yticklabels(labels, fontsize=6.5)
ax2.set_xlabel('Pearson r')
ax2.set_title('B. Top 15 Liver AA Enzyme — Muscle Gene Cross-Tissue Correlations', fontweight='bold', loc='left')
ax2.axvline(x=0, color='black', linewidth=0.5)
ax2.set_xlim(-1, 1)

# Legend
from matplotlib.patches import Patch
legend_patches = [Patch(color=c, label=cat) for cat, c in cat_colors.items()]
ax2.legend(handles=legend_patches, fontsize=5, loc='lower right', ncol=1)

plt.tight_layout()
plt.savefig('/Users/hezongze/pig_study/fig_AA_muscle_liver_crosstalk.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/hezongze/pig_study/fig_AA_muscle_liver_crosstalk.pdf', bbox_inches='tight')
plt.close()
print("Figure 1 saved.")

# ============================================================
# Figure 2: Scatter plots for key AA-Muscle pairs
# ============================================================
key_pairs = [
    ('CBR2', 'TRIM63', 'COX→Prostaglandin'),
    ('PLA2G6', 'TRIM63', 'Membrane AA Release'),
    ('FADS6', 'TRIM63', 'LA→AA Synthesis'),
    ('CYP2E1', 'PPARA', 'CYP450→EET'),
    ('GPX3', 'TRIM63', 'LOX→Leukotriene'),
    ('PLA2G7', 'LTB4R', 'Membrane AA Release'),
]

fig, axes = plt.subplots(2, 3, figsize=(14, 9))
axes = axes.flatten()

for i, (l_gene, m_gene, cat) in enumerate(key_pairs):
    ax = axes[i]
    pair_data = detail[(detail['Liver_Gene'] == l_gene) & (detail['Muscle_Gene'] == m_gene)]

    for breed, color, marker in [('DLY', '#2196F3', 'o'), ('TFB', '#FF5722', 's')]:
        breed_data = pair_data[pair_data['Breed'] == breed]
        ax.scatter(breed_data['Liver_Mean_Expr'], breed_data['Muscle_Mean_Expr'],
                  c=color, marker=marker, s=80, label=breed, edgecolors='white', linewidth=0.5, zorder=3)
        # Connect stages with arrow
        stages_order = ['15', '45', '75', '105']
        points = breed_data[breed_data['Stage'].isin(stages_order)].copy()
        points['stage_num'] = points['Stage'].astype(int)
        points = points.sort_values('stage_num')
        if len(points) > 1:
            ax.plot(points['Liver_Mean_Expr'], points['Muscle_Mean_Expr'],
                   color=color, alpha=0.3, linewidth=1, zorder=2)
            # Annotate stages
            for _, pt in points.iterrows():
                ax.annotate(pt['Stage'], (pt['Liver_Mean_Expr'], pt['Muscle_Mean_Expr']),
                           textcoords="offset points", xytext=(5, 5), fontsize=6, color=color, alpha=0.8)

    # Correlation
    r_val = pair_data['Pearson_r'].iloc[0]
    p_val = pair_data['P_value'].iloc[0]
    ax.set_title(f'{l_gene} ↔ {m_gene}\n({cat})', fontsize=9, fontweight='bold')
    ax.text(0.95, 0.05, f'r = {r_val:.3f}\nP = {p_val:.4f}',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=7,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    ax.set_xlabel(f'Liver {l_gene} (mean TPM)', fontsize=7)
    ax.set_ylabel(f'Muscle {m_gene} (mean TPM)', fontsize=7)
    ax.legend(fontsize=6)

plt.suptitle('Liver AA Enzyme vs Muscle Gene Expression Across 4 Stages (15-105 kg)', fontweight='bold', fontsize=12)
plt.tight_layout()
plt.savefig('/Users/hezongze/pig_study/fig_AA_key_pairs_scatter.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/hezongze/pig_study/fig_AA_key_pairs_scatter.pdf', bbox_inches='tight')
plt.close()
print("Figure 2 saved.")

# ============================================================
# Figure 3: Pathway summary diagram data
# ============================================================
fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)
ax.axis('off')

# Title
ax.text(5, 9.5, 'Arachidonic Acid Metabolism — Liver-Muscle Axis Crosstalk Model',
        ha='center', fontsize=13, fontweight='bold')

# Liver side (left)
ax.text(1.5, 8.5, 'LIVER', fontsize=11, fontweight='bold', color='#1565C0',
        ha='center', bbox=dict(boxstyle='round', facecolor='#E3F2FD', edgecolor='#1565C0'))

# AA pathway categories with key genes
categories_liver = [
    ('LA→AA Synthesis', ['FADS1', 'FADS6', 'FADS2', 'ELOVL2'], 8.0),
    ('Membrane AA Release', ['PLA2G6', 'PLA2G7', 'PLD1/2'], 7.0),
    ('COX→Prostaglandins', ['PTGDS', 'CBR2', 'PTGES', 'PTGS1'], 6.0),
    ('LOX→Leukotrienes', ['GPX3', 'GPX7', 'LTC4S', 'ALOX5'], 5.0),
    ('CYP450→EET/HETE', ['CYP2E1', 'CYP4V2', 'CYP2U1', 'EPHX2'], 4.0),
]

for cat, genes, y in categories_liver:
    ax.text(0.5, y, cat, fontsize=8, fontweight='bold', ha='right', va='center')
    ax.text(1.0, y, ', '.join(genes), fontsize=7, ha='left', va='center', style='italic')

# Blood circulation (middle)
ax.annotate('', xy=(8.5, 5.5), xytext=(1.5, 5.5),
           arrowprops=dict(arrowstyle='->', color='#D32F2F', lw=2, connectionstyle='arc3,rad=0'))
ax.text(5, 6.2, 'Blood Circulation', ha='center', fontsize=9, fontweight='bold', color='#D32F2F')
ax.text(5, 5.7, 'Eicosanoids: PGs, LTs, EETs/HETEs', ha='center', fontsize=8, color='#555555')

# Muscle side (right)
ax.text(8.5, 8.5, 'SKELETAL MUSCLE', fontsize=11, fontweight='bold', color='#C62828',
        ha='center', bbox=dict(boxstyle='round', facecolor='#FFEBEE', edgecolor='#C62828'))

# Muscle receptors with strongest correlations
muscle_targets = [
    ('Membrane Receptors', ['PTGER4 (r=0.745)', 'TBXA2R (r=0.715)', 'LTB4R (r=0.931)'], 8.0),
    ('Nuclear Receptors', ['PPARA (r=0.884)', 'PPARD'], 7.0),
    ('Protein Degradation', ['TRIM63/MuRF1 (r=0.927)', 'FBXO32/Atrogin-1 (r=0.852)', 'FOXO3 (r=0.841)'], 6.0),
    ('Inflammatory', ['IL6R (r=0.784)', 'FABP3 (r=0.802)'], 5.0),
]

for cat, genes, y in muscle_targets:
    ax.text(7.5, y, cat, fontsize=8, fontweight='bold', ha='right', va='center')
    ax.text(8.0, y, ', '.join(genes), fontsize=7, ha='left', va='center', style='italic')

# Bottom summary
summary_text = (
    "Key Finding: Liver AA pathway enzymes show strong positive cross-tissue correlations (r>0.7, P<0.05) "
    "with muscle protein degradation genes (TRIM63, FBXO32) and transcription factors (FOXO3, PPARA).\n"
    "The 'COX→Prostaglandin' and 'CYP450→EET' categories show the strongest coupling (mean r=0.77-0.82).\n"
    "This suggests a liver→muscle signaling axis where AA-derived eicosanoids influence muscle proteolysis programs."
)
ax.text(5, 1.5, summary_text, ha='center', fontsize=8,
        bbox=dict(boxstyle='round', facecolor='#FFF9C4', edgecolor='#F9A825'), wrap=True)

# Bottom stats box
stats_text = f"1458 cross-tissue gene pairs tested | 41 significant (P<0.05) | 39 strong (|r|>0.7)"
ax.text(5, 0.5, stats_text, ha='center', fontsize=7, style='italic', color='#888888')

plt.tight_layout()
plt.savefig('/Users/hezongze/pig_study/fig_AA_axis_model.png', dpi=200, bbox_inches='tight')
plt.savefig('/Users/hezongze/pig_study/fig_AA_axis_model.pdf', bbox_inches='tight')
plt.close()
print("Figure 3 saved.")

print("\nAll figures generated successfully!")
