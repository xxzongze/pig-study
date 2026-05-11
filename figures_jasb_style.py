#!/usr/bin/env python3
"""
Publication-quality figures in JASB journal style.
Adapts figure logic from:
  - Jia et al. (2026) JASB 17:19 (WGCNA multi-tissue, Hu sheep)
  - Chen et al. (2026) JASB 17:24 (liver-muscle axis, broiler)

Color/styling rules:
  - White background, no unnecessary gridlines
  - RdBu_r for correlation heatmaps
  - Red (#D73027) = PD-positive / DLY, Blue (#4575B4) = PD-negative / TFB
  - Arial/Helvetica, 300 dpi, .tiff or .png output
  - Significance: *** P<0.001, ** P<0.01, * P<0.05
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, zscore
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
import matplotlib.ticker as ticker
import os
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# Global style — JASB journal compatible
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.weight': 'bold',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
})

# JASB color palette
C_RED    = '#D73027'   # positive / DLY / significance
C_BLUE   = '#4575B4'   # negative / TFB
C_PURPLE = '#762A83'   # cross-tissue highlight
C_GREEN  = '#1B7837'   # muscle modules
C_ORANGE = '#E66101'   # liver modules
C_GREY   = '#999999'   # grey module / non-sig
C_LIGHT  = '#F7F7F7'   # background
C_BG     = '#FFFFFF'

os.makedirs('figures', exist_ok=True)

# ============================================================
# LOAD ALL DATA
# ============================================================
print("Loading data...")
liver_gm   = pd.read_csv('wgcna_output/liver_gene_module_assignment.csv')
muscle_gm  = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')
liver_mtc  = pd.read_csv('wgcna_output/liver_module_trait_cor.csv', index_col=0)
muscle_mtc = pd.read_csv('wgcna_output/muscle_module_trait_cor.csv', index_col=0)
liver_mtp  = pd.read_csv('wgcna_output/liver_module_trait_pvalue.csv', index_col=0)
muscle_mtp = pd.read_csv('wgcna_output/muscle_module_trait_pvalue.csv', index_col=0)
liver_ms   = pd.read_csv('wgcna_output/liver_module_sizes.csv')
muscle_ms  = pd.read_csv('wgcna_output/muscle_module_sizes.csv')
liver_expr = pd.read_csv('wgcna_output/liver_expr.csv', index_col=0)
muscle_expr= pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)
liver_traits = pd.read_csv('wgcna_output/liver_traits.csv', index_col=0)

# Helper: significance string
def sig_str(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    return ''

def sig_str_long(p):
    if p < 0.0001: return 'P < 0.0001'
    elif p < 0.001: return f'P = {p:.4f}'
    elif p < 0.01: return f'P = {p:.4f}'
    elif p < 0.05: return f'P = {p:.4f}'
    return f'P = {p:.4f} (ns)'

# ============================================================
# Fig 2: Module-Trait Correlation Heatmap (mirrors Jia Fig 7A)
# ============================================================
print("Generating Fig 2: Module-Trait Correlation Heatmap...")

fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5.5),
                           gridspec_kw={'width_ratios': [0.35, 1.0]})

for idx, (tissue, mtc, mtp, ms) in enumerate([
    ('Liver', liver_mtc, liver_mtp, liver_ms),
    ('Muscle', muscle_mtc, muscle_mtp, muscle_ms)
]):
    ax = axes2[idx]

    modules_show = [m for m in mtc.index if m != 'grey']
    traits_show  = ['PD', 'Urea', 'Breed', 'Weight']

    cor_mat = mtc.loc[modules_show, traits_show].values.astype(float)
    p_mat   = mtp.loc[modules_show, traits_show].values.astype(float)

    im = ax.imshow(cor_mat, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1,
                   interpolation='nearest')

    # Annotate each cell
    for i in range(len(modules_show)):
        for j in range(len(traits_show)):
            r_val = cor_mat[i, j]
            p_val = p_mat[i, j]
            ss = sig_str(p_val)
            text = f'{r_val:.2f}'
            if ss:
                text += f'\n{ss}'
            color = 'white' if abs(r_val) > 0.55 else 'black'
            ax.text(j, i, text, ha='center', va='center', fontsize=6.5,
                    color=color, fontweight='bold' if ss else 'normal',
                    linespacing=1.1)

    ax.set_xticks(range(len(traits_show)))
    ax.set_xticklabels(traits_show, rotation=0, ha='center', fontsize=8)
    ax.xaxis.set_ticks_position('bottom')
    ax.set_yticks(range(len(modules_show)))
    ax.set_yticklabels(modules_show, fontsize=7.5)
    ax.set_title(tissue, fontweight='bold', fontsize=10, pad=8)

    # Highlight rows with |r_PD| > 0.5
    for i, mod in enumerate(modules_show):
        if abs(mtc.loc[mod, 'PD']) > 0.5:
            ax.add_patch(plt.Rectangle((-0.5, i - 0.5), len(traits_show), 1,
                                       fill=False, edgecolor=C_RED, linewidth=1.5,
                                       linestyle='-'))

# Unified colorbar
cbar_ax = fig2.add_axes([0.94, 0.15, 0.012, 0.7])
cbar = fig2.colorbar(im, cax=cbar_ax)
cbar.set_label('Pearson r', fontsize=8)
cbar.ax.tick_params(labelsize=6)

fig2.suptitle('Module–Trait Correlations in Pig Liver and Skeletal Muscle',
              fontweight='bold', fontsize=12, y=1.02)
fig2.text(0.5, -0.02, 'Red outline: |r_PD| > 0.5  |  *** P < 0.001  ** P < 0.01  * P < 0.05',
          ha='center', fontsize=7, style='italic')
plt.subplots_adjust(wspace=0.25, top=0.88, bottom=0.10, right=0.92)
fig2.savefig('figures/Fig2_module_trait_heatmap.png', dpi=300, facecolor=C_BG)
fig2.savefig('figures/Fig2_module_trait_heatmap.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig2)
print("  -> figures/Fig2_module_trait_heatmap.png|tiff")

# ============================================================
# Fig 3: Module-PD Association + Module Size (combined info)
# ============================================================
print("Generating Fig 3: Module-PD Association Barplot...")

fig3, axes3 = plt.subplots(1, 2, figsize=(12, 5))

for idx, (tissue, mtc, mtp, ms) in enumerate([
    ('Liver', liver_mtc, liver_mtp, liver_ms),
    ('Muscle', muscle_mtc, muscle_mtp, muscle_ms)
]):
    ax = axes3[idx]

    # Sort by |r_PD|
    modules_show = [m for m in mtc.index if m != 'grey']
    r_pd_vals = {m: mtc.loc[m, 'PD'] for m in modules_show}
    sorted_mods = sorted(modules_show, key=lambda m: abs(r_pd_vals[m]), reverse=True)

    r_vals  = [r_pd_vals[m] for m in sorted_mods]
    p_vals  = [mtp.loc[m, 'PD'] for m in sorted_mods]
    sizes   = [ms[ms['Module'] == m]['Size'].values[0] for m in sorted_mods]
    colors  = [C_RED if r > 0 else C_BLUE for r in r_vals]
    alphas  = [0.9 if p < 0.05 else 0.35 for p in p_vals]

    bars = ax.barh(range(len(sorted_mods)), r_vals, color=colors,
                   edgecolor='white', height=0.7)
    for bar, alpha in zip(bars, alphas):
        bar.set_alpha(alpha)  # set alpha per bar
    ax.axvline(x=0, color='black', linewidth=0.5)

    # Add size annotation on right side
    for i, (r, s, p) in enumerate(zip(r_vals, sizes, p_vals)):
        ss = sig_str(p) if p < 0.05 else ''
        label = f'n={s} {ss}'
        ax.text(max(r, 0) + 0.03 if r >= 0 else min(r, 0) - 0.03, i,
                label, va='center', fontsize=6, color='#333333')

    ax.set_yticks(range(len(sorted_mods)))
    ax.set_yticklabels(sorted_mods, fontsize=7.5)
    ax.set_xlabel("Pearson r (Module ME ~ PD)", fontsize=8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_title(tissue, fontweight='bold', fontsize=10)
    ax.invert_yaxis()

    # Legend
    legend_elements = [
        Patch(facecolor=C_RED, alpha=0.9, label='r_PD > 0 (positive)'),
        Patch(facecolor=C_BLUE, alpha=0.9, label='r_PD < 0 (negative)'),
        Patch(facecolor=C_GREY, alpha=0.35, label='P ≥ 0.05 (ns)'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=6.5,
              frameon=False)

fig3.suptitle('Co-Expression Module Associations with Protein Deposition',
              fontweight='bold', fontsize=12)
plt.tight_layout()
fig3.savefig('figures/Fig3_module_pd_barplot.png', dpi=300, facecolor=C_BG)
fig3.savefig('figures/Fig3_module_pd_barplot.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig3)
print("  -> figures/Fig3_module_pd_barplot.png|tiff")

# ============================================================
# Fig 4: Cross-Tissue Module Eigengene Correlation (mirrors Jia style)
# ============================================================
print("Generating Fig 4: Cross-Tissue Module Correlation...")

# Compute module eigengenes
def compute_module_eigengenes(expr, gene_module):
    mes = {}
    for mod in gene_module['Module'].unique():
        if mod == 'grey':
            continue
        mod_genes = gene_module[gene_module['Module'] == mod]['Gene'].tolist()
        mod_genes_in_expr = [g for g in mod_genes if g in expr.columns]
        if len(mod_genes_in_expr) < 5:
            continue
        mod_expr = expr[mod_genes_in_expr]
        mod_expr_scaled = mod_expr.apply(zscore, nan_policy='omit').fillna(0)
        if mod_expr_scaled.shape[1] > 1:
            U, S, Vt = np.linalg.svd(mod_expr_scaled.values, full_matrices=False)
            pc1 = U[:, 0] * S[0]
            mean_expr = mod_expr.mean(axis=1)
            if np.corrcoef(pc1, mean_expr)[0, 1] < 0:
                pc1 = -pc1
            mes[mod] = pd.Series(pc1, index=mod_expr.index)
    return pd.DataFrame(mes)

liver_MEs  = compute_module_eigengenes(liver_expr, liver_gm)
muscle_MEs = compute_module_eigengenes(muscle_expr, muscle_gm)
common_samples = liver_MEs.index.intersection(muscle_MEs.index)
liver_MEs_aligned  = liver_MEs.loc[common_samples]
muscle_MEs_aligned = muscle_MEs.loc[common_samples]

cross_cor  = pd.DataFrame(index=liver_MEs_aligned.columns,
                          columns=muscle_MEs_aligned.columns, dtype=float)
cross_pval = pd.DataFrame(index=liver_MEs_aligned.columns,
                          columns=muscle_MEs_aligned.columns, dtype=float)
for lm in liver_MEs_aligned.columns:
    for mm in muscle_MEs_aligned.columns:
        r, p = pearsonr(liver_MEs_aligned[lm], muscle_MEs_aligned[mm])
        cross_cor.loc[lm, mm] = r
        cross_pval.loc[lm, mm] = p

# Sort rows/cols by mean |r|
row_order = cross_cor.abs().mean(axis=1).sort_values(ascending=False).index.tolist()
col_order = cross_cor.abs().mean(axis=0).sort_values(ascending=False).index.tolist()
cross_cor = cross_cor.loc[row_order, col_order]
cross_pval = cross_pval.loc[row_order, col_order]

fig4, ax4 = plt.subplots(figsize=(max(7, len(col_order)*0.45),
                                   max(5.5, len(row_order)*0.38)))

im4 = ax4.imshow(cross_cor.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1,
                 interpolation='nearest')

# Annotate cells with |r| > 0.3
for i in range(cross_cor.shape[0]):
    for j in range(cross_cor.shape[1]):
        r_val = cross_cor.values[i, j]
        p_val = cross_pval.values[i, j]
        if abs(r_val) > 0.3:
            ss = sig_str(p_val)
            text = f'{r_val:.2f}'
            if ss:
                text += f'\n{ss}'
            color = 'white' if abs(r_val) > 0.6 else 'black'
            ax4.text(j, i, text, ha='center', va='center', fontsize=5.5,
                     color=color, fontweight='bold' if ss else 'normal')

ax4.set_xticks(range(len(col_order)))
ax4.set_xticklabels(col_order, rotation=45, ha='right', fontsize=7)
ax4.set_yticks(range(len(row_order)))
ax4.set_yticklabels(row_order, fontsize=7)
ax4.set_xlabel('Muscle Modules', fontweight='bold', fontsize=9)
ax4.set_ylabel('Liver Modules', fontweight='bold', fontsize=9)
ax4.set_title('Cross-Tissue Module Eigengene Correlation\n(Liver ↔ Skeletal Muscle)',
              fontweight='bold', fontsize=11)

cbar4 = fig4.colorbar(im4, ax=ax4, shrink=0.75, aspect=25)
cbar4.set_label('Pearson r', fontsize=8)
cbar4.ax.tick_params(labelsize=6)

fig4.text(0.5, -0.03, 'Annotated cells: |r| > 0.3  |  *** P < 0.001  ** P < 0.01  * P < 0.05',
          ha='center', fontsize=7, style='italic')
fig4.tight_layout()
fig4.savefig('figures/Fig4_cross_tissue_correlation.png', dpi=300, facecolor=C_BG)
fig4.savefig('figures/Fig4_cross_tissue_correlation.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig4)
print("  -> figures/Fig4_cross_tissue_correlation.png|tiff")

# ============================================================
# Fig 5: CAV3 Expression Pattern (breed x stage) + AHSG-CAV3-PD axis
# ============================================================
print("Generating Fig S9: CAV3 Expression and AHSG-CAV3-PD Axis...")

fig5 = plt.figure(figsize=(14, 5))

# --- Panel A: CAV3 muscle expression across breed x stage ---
axA = fig5.add_subplot(1, 3, 1)

# Build sample-level data for CAV3
muscle_long_for_cav3 = []
for sample in muscle_expr.index:
    if 'CAV3' not in muscle_expr.columns:
        continue
    parts = sample.split('_')
    breed = parts[0]
    stage_str = parts[1].replace('kg', '')
    stage = int(stage_str)
    cav3_val = muscle_expr.loc[sample, 'CAV3']
    if pd.notna(cav3_val):
        muscle_long_for_cav3.append({'Sample': sample, 'Breed': breed,
                                     'Stage': stage, 'CAV3': cav3_val})
cav3_long = pd.DataFrame(muscle_long_for_cav3)

stages = [15, 45, 75, 105]
positions = np.arange(len(stages))
bar_width = 0.35

for i, (breed, color, hatch) in enumerate([
    ('DLY', C_RED, ''),
    ('TFB', C_BLUE, '//')
]):
    means, sems = [], []
    for st in stages:
        vals = cav3_long[(cav3_long['Breed'] == breed) & (cav3_long['Stage'] == st)]['CAV3']
        means.append(vals.mean())
        sems.append(vals.sem())
    x_pos = positions + (i - 0.5) * bar_width
    bars = axA.bar(x_pos, means, bar_width, yerr=sems,
                   color=color, edgecolor='black', linewidth=0.4,
                   capsize=3, error_kw={'linewidth': 0.8}, label=breed)

axA.set_xticks(positions)
axA.set_xticklabels([f'{s} kg' for s in stages], fontsize=8)
axA.set_ylabel('CAV3 Expression (FPKM)', fontsize=9)
axA.set_xlabel('Growth Stage', fontsize=9)
axA.set_title('CAV3 in Skeletal Muscle', fontweight='bold', fontsize=10)
axA.legend(frameon=False, fontsize=8)

# --- Panel B: CAV3 vs PD scatter ---
axB = fig5.add_subplot(1, 3, 2)

# Build per-sample CAV3 vs PD
cav3_pd_data = []
for sample in cav3_long['Sample'].unique():
    if sample in muscle_expr.index and sample in liver_traits.index:
        c = cav3_long[cav3_long['Sample'] == sample]['CAV3'].values[0]
        pd_val = liver_traits.loc[sample, 'PD']
        breed = cav3_long[cav3_long['Sample'] == sample]['Breed'].values[0]
        stage = cav3_long[cav3_long['Sample'] == sample]['Stage'].values[0]
        cav3_pd_data.append({'CAV3': c, 'PD': pd_val, 'Breed': breed, 'Stage': stage})
cav3_pd_df = pd.DataFrame(cav3_pd_data)

r_cav3_pd, p_cav3_pd = pearsonr(cav3_pd_df['CAV3'], cav3_pd_df['PD'])
for breed, color, marker in [('DLY', C_RED, 'o'), ('TFB', C_BLUE, 's')]:
    sub = cav3_pd_df[cav3_pd_df['Breed'] == breed]
    axB.scatter(sub['CAV3'], sub['PD'], c=color, marker=marker, s=35,
                edgecolors='black', linewidth=0.3, alpha=0.8, label=breed, zorder=3)

# Trend line
z = np.polyfit(cav3_pd_df['CAV3'], cav3_pd_df['PD'], 1)
p_trend = np.poly1d(z)
x_trend = np.linspace(cav3_pd_df['CAV3'].min(), cav3_pd_df['CAV3'].max(), 100)
axB.plot(x_trend, p_trend(x_trend), '--', color='#333333', linewidth=1, zorder=2)

axB.set_xlabel('CAV3 Expression (FPKM)', fontsize=9)
axB.set_ylabel('Protein Deposition (N g/kg·⁵/d)', fontsize=9)
axB.set_title('CAV3 ↔ Protein Deposition', fontweight='bold', fontsize=10)
axB.text(0.05, 0.95, f'r = {r_cav3_pd:.3f}\n{sig_str_long(p_cav3_pd)}',
         transform=axB.transAxes, fontsize=7, va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
axB.legend(frameon=False, fontsize=7)

# --- Panel C: AHSG(liver) vs CAV3(muscle) scatter ---
axC = fig5.add_subplot(1, 3, 3)

common_ahsg = []
for sample in common_samples:
    if 'AHSG' in liver_expr.columns and 'CAV3' in muscle_expr.columns:
        a = liver_expr.loc[sample, 'AHSG']
        c = muscle_expr.loc[sample, 'CAV3']
        if pd.notna(a) and pd.notna(c):
            parts = sample.split('_')
            breed = parts[0]
            common_ahsg.append({'AHSG_liver': a, 'CAV3_muscle': c, 'Breed': breed})
ahsg_df = pd.DataFrame(common_ahsg)

r_ahsg, p_ahsg = pearsonr(ahsg_df['AHSG_liver'], ahsg_df['CAV3_muscle'])

for breed, color, marker in [('DLY', C_RED, 'o'), ('TFB', C_BLUE, 's')]:
    sub = ahsg_df[ahsg_df['Breed'] == breed]
    axC.scatter(sub['AHSG_liver'], sub['CAV3_muscle'], c=color, marker=marker, s=35,
                edgecolors='black', linewidth=0.3, alpha=0.8, label=breed, zorder=3)

z2 = np.polyfit(ahsg_df['AHSG_liver'], ahsg_df['CAV3_muscle'], 1)
p2 = np.poly1d(z2)
x2 = np.linspace(ahsg_df['AHSG_liver'].min(), ahsg_df['AHSG_liver'].max(), 100)
axC.plot(x2, p2(x2), '--', color='#333333', linewidth=1, zorder=2)

axC.set_xlabel('Liver AHSG Expression (FPKM)', fontsize=9)
axC.set_ylabel('Muscle CAV3 Expression (FPKM)', fontsize=9)
axC.set_title('Liver AHSG ↔ Muscle CAV3', fontweight='bold', fontsize=10)
axC.text(0.05, 0.95, f'r = {r_ahsg:.3f}\n{sig_str_long(p_ahsg)}',
         transform=axC.transAxes, fontsize=7, va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
axC.legend(frameon=False, fontsize=7)

fig5.suptitle('Liver–Muscle Axis: CAV3 as a Key Signal Receiver',
              fontweight='bold', fontsize=12)
fig5.tight_layout()
fig5.savefig('figures/FigS9_cav3_expression_axis.png', dpi=300, facecolor=C_BG)
fig5.savefig('figures/FigS9_cav3_expression_axis.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig5)
print("  -> figures/FigS9_cav3_expression_axis.png|tiff")

# ============================================================
# Fig 6: Module-Trait Landscape Scatter (r_PD vs r_Urea)
# ============================================================
print("Generating Fig S2: Module-Trait Landscape...")

fig6, ax6 = plt.subplots(figsize=(8, 7))

plot_data = []
for tissue, mtc, mtp, gm in [
    ('Liver', liver_mtc, liver_mtp, liver_gm),
    ('Muscle', muscle_mtc, muscle_mtp, muscle_gm)
]:
    for mod in mtc.index:
        if mod == 'grey':
            continue
        r_pd   = mtc.loc[mod, 'PD']
        r_urea = mtc.loc[mod, 'Urea']
        p_pd   = mtp.loc[mod, 'PD']
        n      = (gm['Module'] == mod).sum()
        plot_data.append({
            'Tissue': tissue, 'Module': mod,
            'r_PD': r_pd, 'r_Urea': r_urea,
            'p_PD': p_pd, 'n_genes': n
        })
plot_df = pd.DataFrame(plot_data)

# Size proportional to sqrt(n_genes), color by r_PD direction and significance
plot_df['size'] = np.sqrt(plot_df['n_genes']) * 15
plot_df['color'] = np.where(
    (plot_df['r_PD'] > 0.3) & (plot_df['p_PD'] < 0.05), C_RED,
    np.where((plot_df['r_PD'] < -0.3) & (plot_df['p_PD'] < 0.05), C_BLUE,
    np.where(plot_df['r_PD'] > 0, '#F4A582', '#92C5DE')))

for tissue, marker, edge_style in [('Liver', 'D', 'solid'), ('Muscle', 'o', 'solid')]:
    sub = plot_df[plot_df['Tissue'] == tissue]
    ax6.scatter(sub['r_PD'], sub['r_Urea'], s=sub['size'], c=sub['color'],
                marker=marker, edgecolors='black', linewidth=0.4,
                alpha=0.85, label=tissue, zorder=3)

# Annotate key modules (|r_PD| > 0.4 or |r_Urea| > 0.3)
for _, row in plot_df.iterrows():
    if abs(row['r_PD']) > 0.4 or abs(row['r_Urea']) > 0.3:
        ax6.annotate(row['Module'], (row['r_PD'], row['r_Urea']),
                     fontsize=6.5, ha='center', va='bottom',
                     xytext=(0, 4), textcoords='offset points', fontweight='bold')

ax6.axhline(y=0, color='grey', linewidth=0.5, linestyle='--', alpha=0.5)
ax6.axvline(x=0, color='grey', linewidth=0.5, linestyle='--', alpha=0.5)

# Quadrant labels (mirrors quadrant annotation from your Step 2)
ax6.annotate('HIGH PD\nLow Urea\n(Anabolic)',
             xy=(0.65, -0.55), fontsize=7, ha='center', color='#1B7837',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#E5F5E0', alpha=0.7))
ax6.annotate('Low PD\nHigh Urea\n(Catabolic)',
             xy=(-0.55, 0.55), fontsize=7, ha='center', color='#762A83',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='#F4E6F7', alpha=0.7))

ax6.set_xlabel('Pearson r (Module ME ~ Protein Deposition)', fontsize=10)
ax6.set_ylabel('Pearson r (Module ME ~ Serum Urea)', fontsize=10)
ax6.set_title('Co-Expression Module Landscape\n(PD vs Urea Association)',
              fontweight='bold', fontsize=11)

# Legend
legend_elements = [
    plt.Line2D([0], [0], marker='D', color='w', markerfacecolor=C_GREY,
               markersize=8, label='Liver', markeredgecolor='black', markeredgewidth=0.4),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=C_GREY,
               markersize=8, label='Muscle', markeredgecolor='black', markeredgewidth=0.4),
    Patch(facecolor=C_RED, alpha=0.8, label='r_PD > 0.3, P < 0.05'),
    Patch(facecolor=C_BLUE, alpha=0.8, label='r_PD < -0.3, P < 0.05'),
    Patch(facecolor='#92C5DE', alpha=0.5, label='Non-significant'),
]
ax6.legend(handles=legend_elements, loc='lower left', fontsize=7, frameon=False)

fig6.tight_layout()
fig6.savefig('figures/FigS2_module_trait_landscape.png', dpi=300, facecolor=C_BG)
fig6.savefig('figures/FigS2_module_trait_landscape.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig6)
print("  -> figures/FigS2_module_trait_landscape.png|tiff")

# ============================================================
# Fig S8: Hub Gene Summary — ranked by GS_PD (direct PD association)
# ============================================================
print("Generating Fig S8: Hub Gene Summary...")

muscle_pos_mods = []
for mod in muscle_mtc.index:
    if mod == 'grey': continue
    if muscle_mtc.loc[mod, 'PD'] > 0.3 and muscle_mtp.loc[mod, 'PD'] < 0.05:
        muscle_pos_mods.append((mod, muscle_mtc.loc[mod, 'PD']))
muscle_pos_mods.sort(key=lambda x: x[1], reverse=True)
# Show top 3 PD-positive modules
top3 = [m[0] for m in muscle_pos_mods[:3]]

fig7, axes7 = plt.subplots(1, 3, figsize=(15, 5.5))

for idx, mod in enumerate(top3):
    ax = axes7[idx]
    mod_genes = muscle_gm[muscle_gm['Module'] == mod].copy()
    mod_genes = mod_genes.dropna(subset=['kME_module', 'GS_PD'])
    # Top 15 by GS_PD (not kME) — bars now represent PD association strength
    top15 = mod_genes.nlargest(15, 'GS_PD')

    # Bar color: red=PD+, intensity by |GS_PD|
    gs_vals = top15['GS_PD'].values
    kme_vals = top15['kME_module'].values

    bars = ax.barh(range(len(top15)), gs_vals,
                   color=[C_RED if g > 0 else C_BLUE for g in gs_vals],
                   alpha=0.85, edgecolor='white', height=0.7)
    # Annotate each bar with kME value
    for i, (gs, kme) in enumerate(zip(gs_vals, kme_vals)):
        ax.text(gs + 0.01, i, f'kME={kme:.2f}',
                fontsize=5.5, va='center', color='#333333')

    # Highlight SRPX
    for i, gene in enumerate(top15['Gene']):
        if gene == 'SRPX':
            bars[i].set_edgecolor('black')
            bars[i].set_linewidth(2.5)
            bars[i].set_hatch('////')
            bars[i].set_alpha(1.0)
            # Bold gene name
            ax.text(-0.02, i, gene, fontsize=8, fontweight='bold', color=C_RED,
                    va='center', ha='right')

    ax.set_yticks(range(len(top15)))
    ax.set_yticklabels([g if g != 'SRPX' else '' for g in top15['Gene']], fontsize=7)
    ax.set_xlabel('GS_PD (Gene Significance for PD)', fontsize=8, fontweight='bold')
    mod_r = muscle_mtc.loc[mod, 'PD']
    mod_p = muscle_mtp.loc[mod, 'PD']
    n_genes = len(mod_genes)
    # Count SRPX in module
    srpx_in_mod = 'SRPX' in mod_genes['Gene'].values
    ax.set_title(f'{mod} (r_PD={mod_r:+.2f}{sig_str(mod_p)})\nn={n_genes} genes',
                 fontweight='bold', fontsize=9)
    ax.invert_yaxis()

fig7.suptitle('Top PD-Associated Genes in Muscle PD-Positive Modules\n'
              '(Bars = GS_PD | Annotations = kME | SRPX highlighted with black border)',
              fontweight='bold', fontsize=11)
fig7.tight_layout()
fig7.savefig('figures/FigS8_hub_genes.png', dpi=300, facecolor=C_BG)
fig7.savefig('figures/FigS8_hub_genes.tiff', dpi=300, facecolor=C_BG,
             pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig7)
print("  -> figures/FigS8_hub_genes.png|tiff")

# ============================================================
# Fig S1: Soft Threshold Selection (already generated by R, redo in Python)
# ============================================================
print("Generating Fig S1: Soft Threshold Selection...")

figS1, axesS1 = plt.subplots(2, 2, figsize=(10, 8))

for t_idx, (tissue, prefix) in enumerate([
    ('Liver', 'liver'), ('Muscle', 'muscle')
]):
    expr = pd.read_csv(f'wgcna_output/{prefix}_expr.csv', index_col=0)
    powers = list(range(1, 21)) + list(range(22, 31, 2))

    # Compute SFT indices (simplified version)
    sft_r2 = []
    mean_k = []
    for pw in powers:
        try:
            # Simplified: use cor^power
            cor_mat = np.abs(np.corrcoef(expr.values.T))
            adj = cor_mat ** pw
            # Scale-free topology fit (simplified)
            k = adj.sum(axis=0)
            logk = np.log10(k[k > 0])
            if len(logk) > 10:
                hist, edges = np.histogram(logk, bins=20)
                centers = (edges[:-1] + edges[1:]) / 2
                hist_pos = hist[hist > 0]
                centers_pos = centers[hist > 0]
                if len(centers_pos) > 2:
                    z = np.polyfit(centers_pos, np.log10(hist_pos), 1)
                    r2 = np.corrcoef(centers_pos, np.log10(hist_pos))[0, 1] ** 2
                else:
                    r2 = 0
            else:
                r2 = 0
            sft_r2.append(r2)
            mean_k.append(np.mean(k))
        except:
            sft_r2.append(0)
            mean_k.append(100)

    ax_r2 = axesS1[t_idx, 0]
    ax_r2.scatter(powers, sft_r2, c='red', s=20, zorder=3)
    ax_r2.plot(powers, sft_r2, '-', color='red', alpha=0.4, linewidth=0.8)
    ax_r2.axhline(y=0.8, color='blue', linewidth=0.8, linestyle='--')
    ax_r2.set_xlabel('Soft Threshold (power)', fontsize=8)
    ax_r2.set_ylabel('Scale Free Topology Model Fit (R²)', fontsize=8)
    ax_r2.set_title(f'{tissue}: Scale Independence', fontsize=9, fontweight='bold')
    for pw, r2 in zip(powers, sft_r2):
        if r2 > 0.05:
            ax_r2.annotate(str(pw), (pw, r2), fontsize=5, ha='center', va='bottom')

    ax_mk = axesS1[t_idx, 1]
    ax_mk.scatter(powers, mean_k, c='#4575B4', s=20, zorder=3)
    ax_mk.plot(powers, mean_k, '-', color='#4575B4', alpha=0.4, linewidth=0.8)
    ax_mk.set_xlabel('Soft Threshold (power)', fontsize=8)
    ax_mk.set_ylabel('Mean Connectivity', fontsize=8)
    ax_mk.set_title(f'{tissue}: Mean Connectivity', fontsize=9, fontweight='bold')

figS1.suptitle('WGCNA Soft-Thresholding Power Selection', fontweight='bold', fontsize=12)
figS1.tight_layout()
figS1.savefig('figures/FigS1_soft_threshold.png', dpi=300, facecolor=C_BG)
plt.close(figS1)
print("  -> figures/FigS1_soft_threshold.png")

# ============================================================
# Fig S2: Liver-Muscle Bridging Candidate Summary
# ============================================================
print("Generating Fig S2: Bridging Candidates...")

bridge_df = pd.read_excel('wgcna_step2_integration.xlsx', sheet_name='Bridging_Candidates')
both_df = bridge_df[bridge_df['In_Both'] == True].head(15)

if len(both_df) > 0:
    figS2, axS2 = plt.subplots(figsize=(9, 5.5))

    genes   = both_df['Gene'].tolist()
    l_gs_pd = both_df['Liver_GS_PD'].fillna(0).tolist()
    m_gs_pd = both_df['Muscle_GS_PD'].fillna(0).tolist()

    x = np.arange(len(genes))
    w = 0.35

    axS2.barh(x - w/2, l_gs_pd, w, color=C_ORANGE, alpha=0.8,
              edgecolor='white', label='Liver GS_PD')
    axS2.barh(x + w/2, m_gs_pd, w, color=C_GREEN, alpha=0.8,
              edgecolor='white', label='Muscle GS_PD')

    axS2.set_yticks(x)
    axS2.set_yticklabels(genes, fontsize=8)
    axS2.axvline(x=0, color='black', linewidth=0.5)
    axS2.set_xlabel('Gene Significance (r with Protein Deposition)', fontsize=9)
    axS2.set_title('Liver–Muscle Axis Bridging Candidates\n'
                   '(Secreted / signaling factors present in BOTH tissues)',
                   fontweight='bold', fontsize=10)
    axS2.legend(loc='lower right', fontsize=8, frameon=False)
    axS2.invert_yaxis()

    figS2.tight_layout()
    figS2.savefig('figures/FigS11_bridging_candidates.png', dpi=300, facecolor=C_BG)
    plt.close(figS2)
    print("  -> figures/FigS11_bridging_candidates.png")

# ============================================================
# Fig S3: CAV3 vs Urea correlation (supporting cross-tissue logic)
# ============================================================
print("Generating Fig S10: CAV3 vs Serum Urea...")

figS3, axS3 = plt.subplots(figsize=(6, 5))

cav3_urea = []
for sample in common_samples:
    if 'CAV3' in muscle_expr.columns and sample in liver_traits.index:
        c = muscle_expr.loc[sample, 'CAV3']
        urea = liver_traits.loc[sample, 'Urea']
        if pd.notna(c) and pd.notna(urea):
            parts = sample.split('_')
            breed = parts[0]
            cav3_urea.append({'CAV3': c, 'Urea': urea, 'Breed': breed})
cav3_urea_df = pd.DataFrame(cav3_urea)

r_u, p_u = pearsonr(cav3_urea_df['CAV3'], cav3_urea_df['Urea'])

for breed, color, marker in [('DLY', C_RED, 'o'), ('TFB', C_BLUE, 's')]:
    sub = cav3_urea_df[cav3_urea_df['Breed'] == breed]
    axS3.scatter(sub['CAV3'], sub['Urea'], c=color, marker=marker, s=40,
                 edgecolors='black', linewidth=0.3, alpha=0.8, label=breed, zorder=3)

z3 = np.polyfit(cav3_urea_df['CAV3'], cav3_urea_df['Urea'], 1)
axS3.plot(np.sort(cav3_urea_df['CAV3']),
          np.poly1d(z3)(np.sort(cav3_urea_df['CAV3'])),
          '--', color='#333333', linewidth=1, zorder=2)

axS3.set_xlabel('Muscle CAV3 Expression (FPKM)', fontsize=10)
axS3.set_ylabel('Serum Urea (mmol/L)', fontsize=10)
axS3.set_title('Muscle CAV3 ↔ Serum Urea\n(Liver AA Catabolism Readout)',
               fontweight='bold', fontsize=10)
axS3.text(0.05, 0.95, f'r = {r_u:.3f}\n{sig_str_long(p_u)}',
          transform=axS3.transAxes, fontsize=8, va='top',
          bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
axS3.legend(frameon=False, fontsize=8)

figS3.tight_layout()
figS3.savefig('figures/FigS10_cav3_urea.png', dpi=300, facecolor=C_BG)
plt.close(figS3)
print("  -> figures/FigS10_cav3_urea.png")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("ALL FIGURES GENERATED")
print("=" * 60)
print("""
Main Figures:
  Fig 2: Module-Trait Correlation Heatmap (liver + muscle)
  Fig 3: Module-PD Association Barplot
  Fig 4: Cross-Tissue Module Eigengene Correlation Heatmap

Supplementary:
  Fig S1: Soft-Threshold Power Selection
  Fig S2: Module-Trait Landscape (r_PD vs r_Urea)
  Fig S8: Hub Gene Summary for Top PD-Positive Muscle Modules
  Fig S9: CAV3 Expression Pattern + AHSG-CAV3-PD Axis
  Fig S10: CAV3 vs Serum Urea Correlation
  Fig S11: Liver-Muscle Bridging Candidates
""")
print("All files in: figures/")
print("Done!")
