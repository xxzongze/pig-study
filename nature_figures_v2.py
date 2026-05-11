#!/usr/bin/env python3
"""
Nature-style Figures v2 — 修正版
- 从原始TPM直接计算所有肌肉AA基因的Z-score
- 确保热图数据完整
- 散点图使用所有8点(4阶段x2品种)
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from scipy.stats import pearsonr
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Nature Style
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
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

DLY_COLOR = '#4DBBD5'
TFB_COLOR = '#E64B35'

# ============================================================
# Load raw TPM data
# ============================================================
tpm = pd.read_excel('/Users/hezongze/Downloads/gene.tpm.matrix.xlsx')
tpm = tpm.set_index('seq_id')
print(f"TPM: {tpm.shape}")

# Sample grouping
liver_map = {}
muscle_map = {}
for col in tpm.columns:
    parts = col.split('_')
    if parts[0] == 'L':
        sc, bc = parts[1], parts[2]
        stage_map = {'15':'15','45':'45','1':'75','2':'105','3':'135'}
        if sc in stage_map:
            stage = stage_map[sc]
            breed = 'DLY' if bc == '1' else 'TFB'
            liver_map[col] = f"{breed}_L_{stage}"
    elif parts[0] == 'm':
        sc, bc = parts[1], parts[2]
        stage_map = {'15':'15','1':'75','2':'105','3':'135'}
        if sc in stage_map:
            stage = stage_map[sc]
            breed = 'DLY' if bc == '1' else 'TFB'
            muscle_map[col] = f"{breed}_M_{stage}"
    elif parts[0] in ['DLYM', 'TFBM']:
        breed = 'DLY' if parts[0] == 'DLYM' else 'TFB'
        muscle_map[col] = f"{breed}_M_45"

def group_stats(tpm, gmap):
    groups = sorted(set(gmap.values()))
    means = pd.DataFrame(index=tpm.index)
    for g in groups:
        cols = [c for c in gmap if gmap[c] == g]
        if cols:
            means[g] = tpm[cols].mean(axis=1)
    z = means.subtract(means.mean(axis=1), axis=0).divide(means.std(axis=1).replace(0, np.nan), axis=0)
    return means, z

liver_mean, liver_z = group_stats(tpm, liver_map)
muscle_mean, muscle_z = group_stats(tpm, muscle_map)

# ============================================================
# Muscle AA genes (Ensembl IDs)
# ============================================================
muscle_aa = {
    'PTGER2': 'ENSSSCG00000021862',
    'PTGER4': 'ENSSSCG00000024439',
    'PTGIR': 'ENSSSCG00000026602',
    'TBXA2R': 'ENSSSCG00000033759',
    'CYSLTR1': 'ENSSSCG00000039203',
    'CYSLTR2': 'ENSSSCG00000009399',
    'LTB4R': 'ENSSSCG00000020941',
    'PPARA': 'ENSSSCG00000010316',
    'PPARD': 'ENSSSCG00000015482',
    'FABP3': 'ENSSSCG00000016497',
    'IL6R': 'ENSSSCG00000013829',
    'FOXO1': 'ENSSSCG00000001779',
    'FOXO3': 'ENSSSCG00000009668',
    'FBXO32': 'ENSSSCG00000017676',
    'TRIM63': 'ENSSSCG00000008312',
}

# Verify availability
available_muscle = {}
for sym, eid in muscle_aa.items():
    if eid in muscle_mean.index:
        available_muscle[sym] = eid
    else:
        print(f"  MISSING: {sym} ({eid})")

# Build muscle Z-score matrix for ALL available genes
muscle_z_data = {}
for sym, eid in available_muscle.items():
    cols = [f'{b}_M_{s}' for b in ['DLY','TFB'] for s in ['15','45','75','105']]
    vals = [muscle_mean.loc[eid, c] if c in muscle_mean.columns else np.nan for c in cols]
    muscle_z_data[sym] = vals

mz_df = pd.DataFrame(muscle_z_data, index=cols).T
# Z-score within gene
mz_df = mz_df.subtract(mz_df.mean(axis=1), axis=0).divide(mz_df.std(axis=1).replace(0, np.nan), axis=0)

# ============================================================
# Liver AA genes
# ============================================================
aa_liver = pd.read_csv('/Users/hezongze/Downloads/liver_arachidonic_pathway_genes.csv')
liver_aa_data = {}
for _, row in aa_liver.iterrows():
    gid = row['Gene ID']
    if gid in liver_mean.index:
        liver_aa_data[row['Gene Name']] = {
            'eid': gid,
            'category': row['Category'],
            'pattern': row['Pattern'],
        }

# Compute cross-tissue correlations (all pairs, all 8 points)
print("\nComputing cross-tissue correlations...")
corr_results = []
for l_gene, l_info in liver_aa_data.items():
    for m_gene, m_eid in available_muscle.items():
        l_vals, m_vals, labels = [], [], []
        for breed in ['DLY', 'TFB']:
            for s in ['15','45','75','105']:
                lk = f"{breed}_L_{s}"
                mk = f"{breed}_M_{s}"
                if lk in liver_mean.columns and mk in muscle_mean.columns:
                    l_vals.append(liver_mean.loc[l_info['eid'], lk])
                    m_vals.append(muscle_mean.loc[m_eid, mk])
                    labels.append(f"{breed}{s}")
        if len(l_vals) >= 6:
            r, p = pearsonr(l_vals, m_vals)
            corr_results.append({
                'Liver_Gene': l_gene,
                'Liver_Category': l_info['category'],
                'Liver_Pattern': l_info['pattern'],
                'Muscle_Gene': m_gene,
                'Pearson_r': r, 'P_value': p,
                'N': len(l_vals),
                'Liver_vals': l_vals,
                'Muscle_vals': m_vals,
                'Labels': labels,
            })

corr_df = pd.DataFrame(corr_results)
sig = corr_df[corr_df['P_value'] < 0.05]
print(f"Total pairs: {len(corr_df)}, Significant: {len(sig)}")

# ============================================================
# FIGURE 1: Heatmap
# ============================================================
fig1, ax1 = plt.subplots(figsize=(135/25.4, 110/25.4))

gene_groups = [
    ('GPCR Receptors', ['PTGER2', 'PTGER4', 'PTGIR', 'TBXA2R', 'CYSLTR1', 'CYSLTR2', 'LTB4R']),
    ('Nuclear Receptors', ['PPARA', 'PPARD']),
    ('FA Transporter', ['FABP3']),
    ('Cytokine Receptor', ['IL6R']),
    ('Proteolysis / TF', ['FOXO1', 'FOXO3', 'FBXO32', 'TRIM63']),
]

row_genes = []
group_positions = []
for gname, genes in gene_groups:
    start = len(row_genes)
    for g in genes:
        if g in mz_df.index:
            row_genes.append(g)
    end = len(row_genes)
    if end > start:
        group_positions.append((start, end, gname))

mz_plot = mz_df.reindex(row_genes)
# Column order
col_order = []
for b in ['DLY', 'TFB']:
    for s in ['15','45','75','105']:
        c = f'{b}_M_{s}'
        if c in mz_plot.columns:
            col_order.append(c)
mz_plot = mz_plot[col_order]

# Plot
im = ax1.imshow(mz_plot.values, cmap='RdBu_r', aspect='auto', vmin=-2, vmax=2)
for i in range(len(mz_plot)):
    for j in range(len(col_order)):
        val = mz_plot.values[i, j]
        if not np.isnan(val):
            tc = 'white' if abs(val) > 1.2 else 'black'
            ax1.text(j, i, f'{val:.1f}', ha='center', va='center', fontsize=5.5, color=tc)

# Labels
col_lbl = [c.replace('_M_',' ').replace('_',' ') + 'kg' for c in col_order]
ax1.set_xticks(range(len(col_lbl)))
ax1.set_xticklabels(col_lbl, rotation=45, ha='right', fontsize=6)
ax1.set_yticks(range(len(row_genes)))
ax1.set_yticklabels(row_genes, fontsize=6.5)

# Group lines
for start, end, gname in group_positions:
    if start > 0:
        ax1.axhline(y=start-0.5, color='black', lw=0.8)
    mid = (start+end-1)/2
    ax1.text(len(col_order)+0.6, mid, gname, fontsize=5.5, va='center', ha='left',
            fontweight='bold', color='#444444', fontstyle='italic')

# Breed separator
ax1.axvline(x=3.5, color='black', lw=1.0)
ax1.text(1.5, -1.0, 'DLY', fontsize=7, fontweight='bold', ha='center', color=DLY_COLOR)
ax1.text(5.5, -1.0, 'TFB', fontsize=7, fontweight='bold', ha='center', color=TFB_COLOR)

cbar = plt.colorbar(im, ax=ax1, fraction=0.025, pad=0.02)
cbar.set_label('Z-score', fontsize=6, labelpad=1)
cbar.ax.tick_params(labelsize=5, width=0.5)

ax1.set_title('Muscle AA Receptor / Effector Expression\nZ-score Across Breeds and Growth Stages',
              fontsize=8, fontweight='bold', pad=8, loc='left')
ax1.spines[['top','right']].set_visible(False)

plt.tight_layout()
fig1.savefig('/Users/hezongze/pig_study/fig_Nature1_AA_heatmap_v2.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig1.savefig('/Users/hezongze/pig_study/fig_Nature1_AA_heatmap_v2.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Fig 1 saved.")

# ============================================================
# FIGURE 2: Scatter panels
# ============================================================
key_pairs = [
    ('CBR2', 'TRIM63', 'COX pathway → MuRF1'),
    ('CYP2E1', 'PPARA', 'CYP450→EET → PPARα'),
    ('GPX3', 'TRIM63', 'LOX pathway → MuRF1'),
    ('FADS6', 'TRIM63', 'LA→AA synthesis → MuRF1'),
    ('PLA2G6', 'TRIM63', 'PLA2→AA release → MuRF1'),
    ('PLA2G7', 'LTB4R', 'PLA2→AA → BLT receptor'),
]

fig2 = plt.figure(figsize=(183/25.4, 125/25.4))
gs = gridspec.GridSpec(2, 3, figure=fig2, hspace=0.50, wspace=0.42,
                        left=0.08, right=0.97, top=0.90, bottom=0.10)

for i, (lg, mg, sub) in enumerate(key_pairs):
    ax = fig2.add_subplot(gs[i//3, i%3])

    # Find data
    row = corr_df[(corr_df['Liver_Gene']==lg) & (corr_df['Muscle_Gene']==mg)]
    if len(row) == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center')
        continue
    row = row.iloc[0]

    lv = row['Liver_vals']
    mv = row['Muscle_vals']
    labels = row['Labels']

    # Separate by breed
    for j, (lx, mx, lb) in enumerate(zip(lv, mv, labels)):
        breed = 'DLY' if 'DLY' in lb else 'TFB'
        stage = lb.replace('DLY','').replace('TFB','')
        color = DLY_COLOR if breed == 'DLY' else TFB_COLOR
        marker = 'o' if breed == 'DLY' else 's'
        ax.scatter(lx, mx, c=color, marker=marker, s=40,
                  edgecolors='white', linewidth=0.5, zorder=3)
        ax.annotate(stage, (lx, mx), textcoords="offset points",
                   xytext=(6, 4), fontsize=5.5, color=color, alpha=0.8)

    # Connect each breed
    for breed, color in [('DLY', DLY_COLOR), ('TFB', TFB_COLOR)]:
        breed_pts = [(lx, mx) for lx, mx, lb in zip(lv, mv, labels) if breed in lb]
        breed_stages = [int(lb.replace('DLY','').replace('TFB','')) for lb in labels if breed in lb]
        if len(breed_pts) > 1:
            sorted_pts = [p for _, p in sorted(zip(breed_stages, breed_pts))]
            ax.plot([p[0] for p in sorted_pts], [p[1] for p in sorted_pts],
                   color=color, alpha=0.25, lw=0.8, zorder=2)

    r, p = row['Pearson_r'], row['P_value']
    p_str = f'{p:.4f}' if p >= 0.001 else f'{p:.2e}'
    ax.text(0.97, 0.03, f'r = {r:.3f}\nP = {p_str}',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=5.5,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='#CCCCCC', alpha=0.85, lw=0.5))

    ax.set_title(sub, fontsize=7, fontweight='bold', pad=4)
    ax.set_xlabel(f'{lg} (Liver)', fontsize=6, labelpad=2)
    ax.set_ylabel(f'{mg} (Muscle)', fontsize=6, labelpad=2)
    ax.spines[['top','right']].set_visible(False)
    ax.tick_params(labelsize=5.5)

    if i == 0:
        ax.legend([plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=DLY_COLOR, markersize=6),
                   plt.Line2D([0],[0], marker='s', color='w', markerfacecolor=TFB_COLOR, markersize=6)],
                  ['DLY', 'TFB'], fontsize=5.5, frameon=False, loc='upper left',
                  handletextpad=0.3, markerscale=0.8)

# Panel labels
for i, label in enumerate(['a', 'b', 'c', 'd', 'e', 'f']):
    ax = fig2.add_subplot(gs[i//3, i%3])
    ax.text(-0.18, 1.12, label, transform=ax.transAxes, fontsize=9, fontweight='bold')

fig2.suptitle('Liver AA Enzyme — Muscle Gene Cross-Tissue Correlations (4 Stages × 2 Breeds)',
              fontsize=9, fontweight='bold', x=0.08, ha='left')

fig2.savefig('/Users/hezongze/pig_study/fig_Nature2_AA_scatter_v2.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig2.savefig('/Users/hezongze/pig_study/fig_Nature2_AA_scatter_v2.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Fig 2 saved.")

# ============================================================
# FIGURE 3: Model diagram
# ============================================================
fig3, ax = plt.subplots(figsize=(183/25.4, 95/25.4))
ax.set_xlim(0, 12)
ax.set_ylim(0, 7)
ax.axis('off')

# Colors
L_BG, L_BD = '#E3F2FD', '#1565C0'
M_BG, M_BD = '#FFEBEE', '#C62828'
BLD = '#B71C1C'

# --- Liver box ---
liver_box = FancyBboxPatch((0.2, 0.5), 4.3, 6.0, boxstyle="round,pad=0.3",
                           fc=L_BG, ec=L_BD, lw=1.5, alpha=0.4, zorder=0)
ax.add_patch(liver_box)
ax.text(2.35, 6.7, 'LIVER', fontsize=10, fontweight='bold', ha='center', color=L_BD)

liver_items = [
    (0.5, 6.0, 'LA (18:2n-6)', 'normal', '#4DBBD5'),
    (0.5, 5.2, 'FADS1/6, ELOVL2/5', 'small', '#666666'),
    (0.5, 4.5, 'Arachidonic Acid', 'bold', '#E64B35'),
    (0.5, 3.8, 'PLA2G6/7, PLD1/2', 'small', '#666666'),
    (0.5, 3.1, 'Free AA', 'bold', '#E64B35'),
]
for x, y, txt, wt, clr in liver_items:
    fs = 7.5 if wt == 'bold' else 6.5
    fw = 'bold' if wt == 'bold' else 'normal'
    ax.text(x, y, txt, fontsize=fs, ha='left', color=clr, fontweight=fw)

# Three branches
branches = [
    (1.3, 2.3, 'COX/PTGDS\nCBR1/2 → PGs', '#E64B35'),
    (2.5, 2.3, 'LOX/GPX3\nLTC4S → LTs', '#DC0000'),
    (3.7, 2.3, 'CYP450\nCYP2E1/4V2 → EETs', '#00A087'),
]
for x, y, txt, clr in branches:
    ax.text(x, y, txt, fontsize=6.5, ha='center', color=clr, fontweight='bold')
    ax.annotate('', xy=(x, y-0.4), xytext=(x-0.3, 3.0),
               arrowprops=dict(arrowstyle='->', color='#888888', lw=0.6))

# --- Circulation ---
ax.annotate('', xy=(10.0, 3.5), xytext=(4.8, 1.5),
           arrowprops=dict(arrowstyle='->', color=BLD, lw=2.0, connectionstyle='arc3,rad=0'))
ax.text(7.4, 4.5, 'PORTAL / SYSTEMIC', fontsize=8, fontweight='bold', ha='center', color=BLD)
ax.text(7.4, 4.0, 'CIRCULATION', fontsize=8, fontweight='bold', ha='center', color=BLD)
ax.text(7.4, 3.5, 'Eicosanoids as\nendocrine signals', fontsize=6, ha='center', color='#888888', style='italic')

# --- Muscle box ---
muscle_box = FancyBboxPatch((7.5, 0.5), 4.3, 6.0, boxstyle="round,pad=0.3",
                            fc=M_BG, ec=M_BD, lw=1.5, alpha=0.4, zorder=0)
ax.add_patch(muscle_box)
ax.text(9.65, 6.7, 'SKELETAL MUSCLE', fontsize=10, fontweight='bold', ha='center', color=M_BD)

muscle_items = [
    (7.8, 6.0, 'EP4 / TP / BLT', 'bold', '#1565C0'),
    (7.8, 5.3, '(PTGER4/TBXA2R/LTB4R)', 'small', '#888888'),
    (7.8, 4.6, 'cAMP-PKA / Ca²⁺-PKC', 'small', '#666666'),
    (7.8, 3.9, 'PPARα (EET sensor)', 'bold', '#00A087'),
    (7.8, 3.2, 'FOXO3', 'bold', '#E64B35'),
    (7.8, 2.5, 'TRIM63 / FBXO32', 'bold', '#DC0000'),
    (7.8, 1.8, 'Ubiquitin-Proteasome', 'small', '#666666'),
    (7.8, 1.1, '↑ Protein Degradation', 'bold', '#DC0000'),
]
for x, y, txt, wt, clr in muscle_items:
    fs = 7.5 if wt == 'bold' else 6
    fw = 'bold' if wt == 'bold' else 'normal'
    ax.text(x, y, txt, fontsize=fs, ha='left', color=clr, fontweight=fw)

# Arrows in muscle cascade
for yf, yt in [(5.8,5.5), (5.1,4.8), (4.4,4.1), (3.7,3.4), (3.0,2.7), (2.3,2.0)]:
    ax.annotate('', xy=(9.0, yt), xytext=(9.0, yf),
               arrowprops=dict(arrowstyle='->', color='#BBBBBB', lw=0.5))

# --- Key correlations callout ---
callout_box = FancyBboxPatch((0.5, -1.5), 11.0, 1.2, boxstyle="round,pad=0.3",
                             fc='#FFF9C4', ec='#F9A825', lw=1.0, alpha=0.5, zorder=0)
ax.add_patch(callout_box)

corr_text = (
    "Key Cross-Tissue Correlations:  "
    "CBR2 ↔ TRIM63 (r=0.927, P=0.0009)  |  "
    "CYP2E1 ↔ PPARA (r=0.884, P=0.0036)  |  "
    "GPX3 ↔ TRIM63 (r=0.839, P=0.0092)  |  "
    "FADS6 ↔ TRIM63 (r=0.843, P=0.0086)  |  "
    "PLA2G6 ↔ TRIM63 (r=0.808, P=0.0152)"
)
ax.text(6.0, -0.9, corr_text, fontsize=6.5, ha='center', va='center', fontweight='bold')

fig3.savefig('/Users/hezongze/pig_study/fig_Nature3_AA_model_v2.pdf', dpi=300, bbox_inches='tight', pad_inches=0.05)
fig3.savefig('/Users/hezongze/pig_study/fig_Nature3_AA_model_v2.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
plt.close()
print("Fig 3 saved.")

print("\n=== All Nature v2 figures saved ===")
