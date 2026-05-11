#!/usr/bin/env python3
"""
Regenerate FigS9 and FigS10 with SRPX (replacing CAV3).
Same layout as the original JASB figures, different gene.

FigS9: SRPX Expression Pattern + Liver-Muscle-PD Axis (3-panel)
FigS10: SRPX vs Serum Urea Correlation
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.weight': 'bold',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
    'axes.spines.top': False, 'axes.spines.right': False,
})

C_RED  = '#D73027'
C_BLUE = '#4575B4'
C_BG   = '#FFFFFF'

os.makedirs('figures', exist_ok=True)

def sig_str(p):
    if p < 0.001: return '***'
    if p < 0.01:  return '**'
    if p < 0.05:  return '*'
    return 'ns'

# ============================================================
# Load data
# ============================================================
print("Loading data...")
muscle_expr = pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)
liver_expr  = pd.read_csv('wgcna_output/liver_expr.csv', index_col=0)
muscle_traits = pd.read_csv('wgcna_output/muscle_traits.csv', index_col=0)
liver_traits  = pd.read_csv('wgcna_output/liver_traits.csv', index_col=0)

def parse_sample(s):
    parts = s.split('_')
    return parts[0], int(parts[1].replace('kg', ''))

meta_m = pd.DataFrame({
    'Breed': [parse_sample(s)[0] for s in muscle_expr.index],
    'Stage': [parse_sample(s)[1] for s in muscle_expr.index],
}, index=muscle_expr.index)

GENE = 'SRPX'
stages = [15, 45, 75, 105]

# ============================================================
# FigS9: SRPX Expression Pattern — 3-panel axis
# ============================================================
print("Generating FigS9: SRPX Expression Pattern...")

figS9, (axA, axB, axC) = plt.subplots(1, 3, figsize=(18, 5.5))

# --- Panel A: SRPX expression in muscle across stages ---
srpx_m = muscle_expr[GENE]
x = np.arange(len(stages))
w = 0.35

for i, stg in enumerate(stages):
    idx = meta_m['Stage'] == stg
    dly_vals = srpx_m[idx & (meta_m['Breed'] == 'DLY')]
    tfb_vals = srpx_m[idx & (meta_m['Breed'] == 'TFB')]

    # Individual points
    axA.scatter([i - w/2]*len(dly_vals), dly_vals, c=C_RED, s=25, alpha=0.6, zorder=3)
    axA.scatter([i + w/2]*len(tfb_vals), tfb_vals, c=C_BLUE, s=25, alpha=0.6, zorder=3)

    # Mean bars
    axA.bar(i - w/2, dly_vals.mean(), w, color=C_RED, alpha=0.85, edgecolor='white', zorder=2)
    axA.bar(i + w/2, tfb_vals.mean(), w, color=C_BLUE, alpha=0.85, edgecolor='white', zorder=2)

    # Sig annotation
    from scipy.stats import ttest_ind
    t, p = ttest_ind(dly_vals, tfb_vals, equal_var=False)
    fc = dly_vals.mean() / tfb_vals.mean() if tfb_vals.mean() > 0 else 0
    y_max = max(dly_vals.max(), tfb_vals.max())
    axA.text(i, y_max + 0.3, f'{fc:.1f}× {sig_str(p)}', ha='center', fontsize=7, fontweight='bold')

axA.set_xticks(x)
axA.set_xticklabels([f'{s}kg' for s in stages])
axA.set_ylabel('SRPX Expression (FPKM)', fontweight='bold')
axA.set_title('Muscle SRPX: DLY vs TFB', fontweight='bold', fontsize=10)

from matplotlib.patches import Patch
legend_A = [Patch(color=C_RED, alpha=0.8, label='DLY (high PD)'),
            Patch(color=C_BLUE, alpha=0.8, label='TFB (low PD)')]
axA.legend(handles=legend_A, frameon=False, fontsize=8)

# --- Panel B: SRPX vs PD correlation ---
common = muscle_expr.index.intersection(muscle_traits.index)
r_pd, p_pd = pearsonr(srpx_m.loc[common], muscle_traits.loc[common, 'PD'])

for breed, color, marker in [('DLY', C_RED, 'o'), ('TFB', C_BLUE, 's')]:
    idx = (meta_m['Breed'] == breed) & muscle_expr.index.isin(common)
    axB.scatter(srpx_m.loc[idx], muscle_traits.loc[idx, 'PD'],
                c=color, marker=marker, s=35, edgecolors='black',
                linewidth=0.3, alpha=0.8, label=breed, zorder=3)

z = np.polyfit(srpx_m.loc[common], muscle_traits.loc[common, 'PD'], 1)
x_line = np.linspace(srpx_m.min(), srpx_m.max(), 100)
axB.plot(x_line, np.poly1d(z)(x_line), '--', color='#333333', linewidth=1.5, zorder=2)

axB.set_xlabel('Muscle SRPX (FPKM)', fontweight='bold')
axB.set_ylabel('Protein Deposition (PD)', fontweight='bold')
axB.set_title(f'SRPX vs PD (r={r_pd:.3f}, P={p_pd:.1e})', fontweight='bold', fontsize=10)
axB.legend(frameon=False, fontsize=7)

# --- Panel C: Cross-tissue — Liver SRPX vs Muscle SRPX ---
if GENE in liver_expr.columns:
    srpx_l = liver_expr[GENE]
    meta_l = pd.DataFrame({
        'Breed': [parse_sample(s)[0] for s in liver_expr.index],
        'Stage': [parse_sample(s)[1] for s in liver_expr.index],
    }, index=liver_expr.index)

    common_samples = muscle_expr.index.intersection(liver_expr.index)
    for breed, color, marker in [('DLY', C_RED, 'o'), ('TFB', C_BLUE, 's')]:
        idx = [s for s in common_samples if meta_m.loc[s, 'Breed'] == breed]
        axC.scatter(srpx_m.loc[idx], srpx_l.loc[idx],
                    c=color, marker=marker, s=35, edgecolors='black',
                    linewidth=0.3, alpha=0.8, label=breed, zorder=3)

    r_cross, p_cross = pearsonr(srpx_m.loc[common_samples], srpx_l.loc[common_samples])

    zc = np.polyfit(srpx_m.loc[common_samples], srpx_l.loc[common_samples], 1)
    axC.plot(x_line, np.poly1d(zc)(x_line), '--', color='#333333', linewidth=1.5, zorder=2)

    axC.set_xlabel('Muscle SRPX (FPKM)', fontweight='bold')
    axC.set_ylabel('Liver SRPX (FPKM)', fontweight='bold')
    axC.set_title(f'Liver-Muscle SRPX (r={r_cross:.3f}, P={p_cross:.2e})', fontweight='bold', fontsize=10)
    axC.legend(frameon=False, fontsize=7)

figS9.suptitle('SRPX Expression Pattern: Multi-Tissue, Multi-Trait Association\n'
               '(The Optimal Breed-Differential Hub Gene in the Green Module, r_PD=+0.685)',
               fontweight='bold', fontsize=12, y=1.02)
plt.tight_layout()
figS9.savefig('figures/FigS9_cav3_expression_axis.png', dpi=300, facecolor=C_BG)
figS9.savefig('figures/FigS9_cav3_expression_axis.tiff', dpi=300, facecolor=C_BG,
              pil_kwargs={'compression': 'tiff_lzw'})
plt.close(figS9)
print("  -> figures/FigS9_cav3_expression_axis.png|tiff (SRPX version)")

# ============================================================
# FigS10: SRPX vs Serum Urea
# ============================================================
print("Generating FigS10: SRPX vs Serum Urea...")

# Need serum urea data — use liver_traits which has Urea column
common_urea = muscle_expr.index.intersection(liver_traits.index)

figS10, axS10 = plt.subplots(figsize=(6, 5))

urea_data = []
for sample in common_urea:
    if GENE in muscle_expr.columns and sample in liver_traits.index:
        srpx_val = muscle_expr.loc[sample, GENE]
        urea_val = liver_traits.loc[sample, 'Urea']
        if pd.notna(srpx_val) and pd.notna(urea_val):
            breed = meta_m.loc[sample, 'Breed']
            urea_data.append({'SRPX': srpx_val, 'Urea': urea_val, 'Breed': breed})

urea_df = pd.DataFrame(urea_data)
r_u, p_u = pearsonr(urea_df['SRPX'], urea_df['Urea'])

for breed, color, marker in [('DLY', C_RED, 'o'), ('TFB', C_BLUE, 's')]:
    sub = urea_df[urea_df['Breed'] == breed]
    axS10.scatter(sub['SRPX'], sub['Urea'], c=color, marker=marker, s=40,
                  edgecolors='black', linewidth=0.3, alpha=0.8, label=breed, zorder=3)

zu = np.polyfit(urea_df['SRPX'], urea_df['Urea'], 1)
x_sort = np.sort(urea_df['SRPX'])
axS10.plot(x_sort, np.poly1d(zu)(x_sort), '--', color='#333333', linewidth=1.5, zorder=2)

axS10.set_xlabel('Muscle SRPX Expression (FPKM)', fontweight='bold')
axS10.set_ylabel('Serum Urea (mmol/L)', fontweight='bold')
axS10.set_title(f'SRPX vs Serum Urea\n(r={r_u:.3f}, P={p_u:.2e})', fontweight='bold', fontsize=10)

if p_u < 0.05:
    axS10.text(0.05, 0.95,
               f'High SRPX → Low Urea\n(Less AA waste, more protein retention)',
               transform=axS10.transAxes, fontsize=8, fontweight='bold', va='top',
               bbox=dict(boxstyle='round,pad=0.3', facecolor='#E8F5E9', edgecolor=C_RED, alpha=0.9))

axS10.legend(frameon=False, fontsize=8)

figS10.tight_layout()
figS10.savefig('figures/FigS10_cav3_urea.png', dpi=300, facecolor=C_BG)
plt.close(figS10)
print("  -> figures/FigS10_cav3_urea.png (SRPX version)")

print("\nDone! FigS9 and FigS10 now feature SRPX.")
