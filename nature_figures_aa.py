#!/usr/bin/env python3
"""
Nature-style Figures for AA Metabolism Liver-Muscle Axis
Fig 1: Muscle AA receptor expression heatmap
Fig 2: Key cross-tissue scatter plots
Fig 3: Mechanistic model diagram
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, Rectangle
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Nature Global Style
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica'],
    'font.size': 7,
    'axes.titlesize': 8,
    'axes.labelsize': 7,
    'xtick.labelsize': 6.5,
    'ytick.labelsize': 6.5,
    'legend.fontsize': 6,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    'xtick.major.pad': 2,
    'ytick.major.pad': 2,
    'lines.linewidth': 1.0,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

# Nature color palette
NATURE_COLORS = ['#E64B35', '#4DBBD5', '#00A087', '#3C5488', '#F39B7F', '#8491B4', '#91D1C2', '#DC0000']
DLY_COLOR = '#4DBBD5'   # Blue
TFB_COLOR = '#E64B35'   # Red/Coral

# Load data
detail = pd.read_csv('/Users/hezongze/pig_study/aa_crosstalk_detail_v2.csv')
corr = pd.read_excel('/Users/hezongze/pig_study/aa_crosstissue_full_results.xlsx')

# ============================================================
# FIGURE 1: Muscle AA Receptor/Effector Z-score Heatmap
# ============================================================
# Single column width: 89mm → slightly wider for heatmap: 120mm
fig_width_mm = 135
fig_height_mm = 120
fig1, ax1 = plt.subplots(figsize=(fig_width_mm/25.4, fig_height_mm/25.4))

# Gene order with functional groups
gene_groups = {
    'GPCR\nReceptors': ['PTGER2', 'PTGER4', 'PTGIR', 'TBXA2R', 'CYSLTR1', 'CYSLTR2', 'LTB4R'],
    'Nuclear\nReceptors': ['PPARA', 'PPARD'],
    'Transport': ['FABP3'],
    'Signaling': ['IL6R'],
    'Proteolysis': ['FOXO1', 'FOXO3', 'FBXO32', 'TRIM63'],
}

all_genes = []
group_boundaries = []
for group, genes in gene_groups.items():
    start = len(all_genes)
    all_genes.extend(genes)
    group_boundaries.append((start, len(all_genes), group))

# Build z-score matrix
muscle_genes_in_detail = detail['Muscle_Gene'].unique()
available = [g for g in all_genes if g in muscle_genes_in_detail]
# Also add back genes not in detail but needed
missing = set(all_genes) - set(available)

# Pivot
mp = detail[detail['Muscle_Gene'].isin(available)].pivot_table(
    index='Muscle_Gene', columns=['Breed', 'Stage'],
    values='Muscle_Mean_Expr', aggfunc='mean')

# Z-score within gene
mp_z = mp.subtract(mp.mean(axis=1), axis=0).divide(mp.std(axis=1).replace(0, np.nan), axis=0)

# Reorder columns: DLY 15-45-75-105, TFB 15-45-75-105
col_order = []
for breed in ['DLY', 'TFB']:
    for s in ['15', '45', '75', '105']:
        col_key = (breed, s)
        if col_key in mp_z.columns:
            col_order.append(col_key)
mp_z = mp_z[col_order]

# Reorder rows
row_order = [g for g in all_genes if g in mp_z.index]
mp_z = mp_z.reindex(row_order)

# Column labels
col_labels = [f'{b} {s}kg' for b, s in mp_z.columns]
# Row labels
row_labels = list(mp_z.index)

# Plot heatmap
im = ax1.imshow(mp_z.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)

# Annotate cells
for i in range(len(mp_z.index)):
    for j in range(len(mp_z.columns)):
        val = mp_z.values[i, j]
        if not np.isnan(val):
            text_color = 'white' if abs(val) > 1.2 else 'black'
            ax1.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=5.5,
                    color=text_color, fontweight='bold' if abs(val) > 1.5 else 'normal')

# Axis labels
ax1.set_xticks(range(len(col_labels)))
ax1.set_xticklabels(col_labels, rotation=45, ha='right', fontsize=6)
ax1.set_yticks(range(len(row_labels)))
ax1.set_yticklabels(row_labels, fontsize=6.5)

# Group separators (horizontal lines between functional groups)
for start, end, group_name in group_boundaries:
    visible_start = max(0, row_order.index(all_genes[start]) if all_genes[start] in row_order else -1)
    if visible_start >= 0 and visible_start > 0:
        ax1.axhline(y=visible_start - 0.5, color='black', linewidth=0.8, linestyle='-')
    # Group labels on right side
    mid = (start + end) / 2
    mid_visible = None
    for g in all_genes[start:end]:
        if g in row_order:
            mid_visible = row_order.index(g) + (end - start - 1) / 2
            break
    if mid_visible is not None:
        # Adjust to center of visible group
        visible_genes = [g for g in all_genes[start:end] if g in row_order]
        if visible_genes:
            first_idx = row_order.index(visible_genes[0])
            last_idx = row_order.index(visible_genes[-1])
            mid_visible = (first_idx + last_idx) / 2
            ax1.text(len(col_labels) + 0.8, mid_visible, group_name.replace('\n', ' '),
                    fontsize=6, va='center', ha='left', fontweight='bold',
                    color='#333333')

# Breed separator
ax1.axvline(x=3.5, color='black', linewidth=1.0, linestyle='-')

# Colorbar
cbar = plt.colorbar(im, ax=ax1, fraction=0.025, pad=0.02)
cbar.set_label('Z-score', fontsize=6, labelpad=1)
cbar.ax.tick_params(labelsize=5, width=0.5)

# Title
ax1.set_title('Muscle AA Receptor / Effector Expression\nZ-score Across Breeds and Growth Stages',
              fontsize=8, fontweight='bold', pad=8, loc='left')

# Breed labels on top
ax1.text(1.5, -1.2, 'DLY', fontsize=7, fontweight='bold', ha='center', color=DLY_COLOR)
ax1.text(5.5, -1.2, 'TFB', fontsize=7, fontweight='bold', ha='center', color=TFB_COLOR)

ax1.spines[['top', 'right']].set_visible(False)

plt.tight_layout()
fig1.savefig('/Users/hezongze/pig_study/fig_Nature1_AA_heatmap.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig1.savefig('/Users/hezongze/pig_study/fig_Nature1_AA_heatmap.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Nature Figure 1 saved.")

# ============================================================
# FIGURE 2: Key Cross-Tissue Scatter Panels
# ============================================================
key_pairs = [
    ('CBR2', 'TRIM63', 'COX→PGS → MuRF1'),
    ('CYP2E1', 'PPARA', 'CYP450→EET → PPARα'),
    ('GPX3', 'TRIM63', 'LOX→LTs → MuRF1'),
    ('FADS6', 'TRIM63', 'LA→AA → MuRF1'),
    ('PLA2G6', 'TRIM63', 'PLA2→AA → MuRF1'),
    ('PLA2G7', 'LTB4R', 'PLA2→AA → BLT'),
]

fig_width_mm = 183  # Double column
fig_height_mm = 130
fig2 = plt.figure(figsize=(fig_width_mm/25.4, fig_height_mm/25.4))
gs = gridspec.GridSpec(2, 3, figure=fig2, hspace=0.45, wspace=0.40,
                        left=0.08, right=0.98, top=0.92, bottom=0.12)

for i, (l_gene, m_gene, subtitle) in enumerate(key_pairs):
    ax = fig2.add_subplot(gs[i // 3, i % 3])

    pair_data = detail[(detail['Liver_Gene'] == l_gene) & (detail['Muscle_Gene'] == m_gene)]

    for breed, color, marker in [('DLY', DLY_COLOR, 'o'), ('TFB', TFB_COLOR, 's')]:
        breed_data = pair_data[pair_data['Breed'] == breed].copy()
        # Sort by stage
        breed_data['stage_num'] = breed_data['Stage'].astype(int)
        breed_data = breed_data.sort_values('stage_num')

        ax.scatter(breed_data['Liver_Mean_Expr'], breed_data['Muscle_Mean_Expr'],
                  c=color, marker=marker, s=35, label=breed, edgecolors='white',
                  linewidth=0.3, zorder=3, alpha=0.9)

        if len(breed_data) > 1:
            ax.plot(breed_data['Liver_Mean_Expr'], breed_data['Muscle_Mean_Expr'],
                   color=color, alpha=0.25, linewidth=0.7, zorder=2)
            for _, pt in breed_data.iterrows():
                offset_x = 3 if breed == 'DLY' else -8
                offset_y = 3 if pt['Stage'] != '75' else -8
                ax.annotate(f"{pt['Stage']}kg", (pt['Liver_Mean_Expr'], pt['Muscle_Mean_Expr']),
                           textcoords="offset points", xytext=(offset_x, offset_y),
                           fontsize=5, color=color, alpha=0.7)

    # Stats
    r_val = pair_data['Pearson_r'].iloc[0]
    p_val = pair_data['P_value'].iloc[0]
    p_str = f'{p_val:.4f}' if p_val >= 0.001 else f'{p_val:.2e}'

    ax.text(0.97, 0.03, f'r = {r_val:.3f}\nP = {p_str}',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=5.5,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='#CCCCCC',
                     alpha=0.85, linewidth=0.5))

    ax.set_title(subtitle, fontsize=7, fontweight='bold', pad=4)
    ax.set_xlabel(f'{l_gene} (Liver TPM)', fontsize=6, labelpad=2)
    ax.set_ylabel(f'{m_gene} (Muscle TPM)', fontsize=6, labelpad=2)

    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelsize=5.5)

    if i == 0:
        ax.legend(fontsize=5.5, frameon=False, loc='upper left',
                 markerscale=0.7, handletextpad=0.5)

# Panel labels
for i, label in enumerate(['a', 'b', 'c', 'd', 'e', 'f']):
    ax = fig2.add_subplot(gs[i // 3, i % 3])
    ax.text(-0.15, 1.08, label, transform=ax.transAxes, fontsize=9, fontweight='bold')

fig2.suptitle('Liver AA Enzyme — Muscle Gene Cross-Tissue Correlation Across Growth Stages',
              fontsize=9, fontweight='bold', x=0.08, ha='left')

fig2.savefig('/Users/hezongze/pig_study/fig_Nature2_AA_scatter.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig2.savefig('/Users/hezongze/pig_study/fig_Nature2_AA_scatter.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Nature Figure 2 saved.")

# ============================================================
# FIGURE 3: Mechanistic Model Diagram
# ============================================================
fig_width_mm = 183
fig_height_mm = 110
fig3, ax3 = plt.subplots(figsize=(fig_width_mm/25.4, fig_height_mm/25.4))
ax3.set_xlim(0, 12)
ax3.set_ylim(0, 8)
ax3.axis('off')

# --- Color scheme ---
LIVER_BG = '#E3F2FD'
LIVER_BORDER = '#1565C0'
BLOOD_COLOR = '#D32F2F'
MUSCLE_BG = '#FFEBEE'
MUSCLE_BORDER = '#C62828'
ARROW_COLOR = '#555555'
HIGHLIGHT = '#FF6F00'

# --- LIVER compartment (left) ---
liver_rect = FancyBboxPatch((0.2, 0.3), 4.2, 7.2, boxstyle="round,pad=0.3",
                             facecolor=LIVER_BG, edgecolor=LIVER_BORDER, linewidth=1.5, alpha=0.5)
ax3.add_patch(liver_rect)
ax3.text(2.3, 7.7, 'LIVER', fontsize=11, fontweight='bold', ha='center', color=LIVER_BORDER)

# Liver pathways (top to bottom)
liver_paths = [
    (0.5, 6.8, 'LA (18:2n-6)', '#4DBBD5'),
    (0.5, 5.8, 'FADS1/6, ELOVL2/5', '#666666'),
    (0.5, 5.0, 'Arachidonic Acid (20:4n-6)', '#E64B35'),
    (0.5, 4.2, 'PLA2G6/7, PLD1/2', '#666666'),
    (0.5, 3.4, 'Free AA', '#E64B35'),
]

for x, y, text, color in liver_paths:
    ax3.text(x, y, text, fontsize=7.5, ha='left', color=color, fontweight='bold' if 'AA' in text or 'LA' in text else 'normal')

# Three branches from Free AA
branch_start = (2.8, 3.2)

# COX branch
ax3.annotate('', xy=(1.2, 2.0), xytext=branch_start,
            arrowprops=dict(arrowstyle='->', color='#555555', lw=0.8))
ax3.text(1.5, 2.6, 'COX / PTGDS\nCBR2', fontsize=6, ha='center', color='#555555')
ax3.text(0.6, 1.6, 'PGE₂ / PGD₂ / TXA₂', fontsize=7, ha='left', fontweight='bold', color='#E64B35')

# LOX branch
ax3.annotate('', xy=(2.8, 2.0), xytext=branch_start,
            arrowprops=dict(arrowstyle='->', color='#555555', lw=0.8))
ax3.text(2.8, 2.6, 'LOX / GPX3\nLTC4S', fontsize=6, ha='center', color='#555555')
ax3.text(2.3, 1.6, 'LTs / HETEs', fontsize=7, ha='left', fontweight='bold', color='#DC0000')

# CYP branch
ax3.annotate('', xy=(4.0, 2.0), xytext=branch_start,
            arrowprops=dict(arrowstyle='->', color='#555555', lw=0.8))
ax3.text(4.0, 2.6, 'CYP450\nCYP2E1/4V2', fontsize=6, ha='center', color='#555555')
ax3.text(3.5, 1.6, 'EETs / 20-HETE', fontsize=7, ha='left', fontweight='bold', color='#00A087')

# --- CIRCULATION arrow ---
ax3.annotate('', xy=(10.0, 5.5), xytext=(5.0, 5.5),
            arrowprops=dict(arrowstyle='->', color=BLOOD_COLOR, lw=2.5,
                          connectionstyle='arc3,rad=0', linestyle='-'))
ax3.text(7.5, 6.5, 'PORTAL / SYSTEMIC\nCIRCULATION', fontsize=8, fontweight='bold',
         ha='center', color=BLOOD_COLOR)
ax3.text(7.5, 6.0, 'Eicosanoids as endocrine signals', fontsize=6.5, ha='center',
         color='#888888', style='italic')

# --- MUSCLE compartment (right) ---
muscle_rect = FancyBboxPatch((7.6, 0.3), 4.2, 7.2, boxstyle="round,pad=0.3",
                              facecolor=MUSCLE_BG, edgecolor=MUSCLE_BORDER, linewidth=1.5, alpha=0.5)
ax3.add_patch(muscle_rect)
ax3.text(9.7, 7.7, 'SKELETAL MUSCLE', fontsize=11, fontweight='bold', ha='center', color=MUSCLE_BORDER)

# Muscle signaling cascade (top to bottom)
muscle_cascade = [
    (8.0, 6.8, 'EP4 (PTGER4) / TP (TBXA2R)\nBLT (LTB4R) — Membrane GPCRs', '#666666'),
    (8.0, 5.8, 'cAMP-PKA / Ca²⁺-PKC\nSecond Messengers', '#888888'),
    (8.0, 4.8, 'PPARα / PPARδ — Nuclear Receptors\n(EET as ligand)', '#888888'),
    (8.0, 3.8, 'FOXO3 — Transcription Factor\n(r = 0.841 with CBR1)', '#E64B35'),
    (8.0, 2.8, 'TRIM63 / FBXO32\nMuRF1 / Atrogin-1', '#DC0000'),
    (8.0, 1.8, 'Ubiquitin-Proteasome\nSystem Activation', '#888888'),
    (8.0, 0.8, '↑ Muscle Protein Degradation\n↓ Protein Deposition', '#DC0000'),
]

for x, y, text, color in muscle_cascade:
    fontsize = 7.5 if any(kw in text for kw in ['FOXO3', 'TRIM63', 'Degradation']) else 6.5
    fontweight = 'bold' if any(kw in text for kw in ['FOXO3', 'TRIM63', 'Degradation']) else 'normal'
    ax3.text(x, y, text, fontsize=fontsize, ha='left', color=color, fontweight=fontweight)

# Arrows between cascade levels
for y_from, y_to in [(6.5, 6.1), (5.5, 5.1), (4.5, 4.1), (3.5, 3.1), (2.5, 2.1), (1.5, 1.1)]:
    ax3.annotate('', xy=(9.0, y_to), xytext=(9.0, y_from),
                arrowprops=dict(arrowstyle='->', color='#AAAAAA', lw=0.6))

# --- Correlation callouts ---
callouts = [
    (None, None, 'Key Cross-Tissue Correlations'),
    (0.5, 0.1, 'CBR2 ↔ TRIM63: r = 0.927, P = 0.0009'),
    (0.5, -0.3, 'CYP2E1 ↔ PPARA: r = 0.884, P = 0.0036'),
    (0.5, -0.7, 'GPX3 ↔ TRIM63: r = 0.839, P = 0.0092'),
]

for i, (x, y, text) in enumerate(callouts):
    fs = 8 if i == 0 else 6.5
    fw = 'bold' if i == 0 else 'normal'
    ax3.text(6.0, y if y else -0.5 + i*0.4, text, fontsize=fs, ha='center',
            fontweight=fw, color='#333333' if i > 0 else HIGHLIGHT)

# --- Bottom explanatory note ---
note = ("Model: Liver AA metabolism → eicosanoid production (PGE₂, LTs, EETs) → "
        "circulation → muscle GPCR/PPAR activation → FOXO3 → MuRF1/Atrogin-1 → "
        "ubiquitin-proteasome → protein degradation")
ax3.text(6.0, -1.5, note, fontsize=6, ha='center', va='top', style='italic',
        color='#888888',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#F5F5F5', edgecolor='#DDDDDD', linewidth=0.5))

fig3.savefig('/Users/hezongze/pig_study/fig_Nature3_AA_model.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig3.savefig('/Users/hezongze/pig_study/fig_Nature3_AA_model.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Nature Figure 3 saved.")

print("\n=== All Nature-style figures generated ===")
print("fig_Nature1_AA_heatmap.pdf/png")
print("fig_Nature2_AA_scatter.pdf/png")
print("fig_Nature3_AA_model.pdf/png")
