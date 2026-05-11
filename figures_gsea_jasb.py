#!/usr/bin/env python3
"""
GSEA Enrichment Figures — JASB Format
======================================
Replacements for the buggy self-written GSEA figures.
Uses externally validated GSEA results (KEGG pathways from clusterProfiler/KOBAS).

Output:
  - Fig_GSEA_liver_enrichment.pdf/png   — Dotplot: KEGG pathways enriched in liver
  - Fig_GSEA_muscle_enrichment.pdf/png  — Dotplot: KEGG pathways enriched in muscle
  - Fig_GSEA_NES_comparison.pdf/png     — Bar plot: NES comparison liver vs muscle
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
import os

from style_jasb import (apply_jasb_style, panel_label, save_figure,
                        C_DLY, C_TFB, C_STAGE_45, pval_stars, fdr_stars)

apply_jasb_style()

# ============================================================
# Load external GSEA data
# ============================================================
liver = pd.read_csv('gsea_external_liver.csv')
muscle = pd.read_csv('gsea_external_muscle.csv')

# Standardize columns
for df in [liver, muscle]:
    df['NES'] = pd.to_numeric(df['NES'], errors='coerce')
    df['Padjust'] = pd.to_numeric(df['Padjust'], errors='coerce')
    df['NES_abs'] = df['NES'].abs()
    # Direction: positive NES = TFBL/TFBM (TFB enriched)
    df['Direction'] = df['Group'].apply(lambda x: 'TFB' if 'TFB' in str(x) else 'DLY')
    # -log10(FDR)
    df['neg_log10_fdr'] = -np.log10(df['Padjust'].clip(lower=1e-5))


def prepare_enrichment_data(df, tissue_name, top_n=25):
    """Filter and prepare enrichment results for plotting."""
    df = df.copy()

    # Extract KEGG pathways only (those starting with ssc + 5 digits)
    kegg = df[df['Gene Set Name'].str.match(r'^ssc\d{5}$', na=False)].copy()
    # Sort by significance
    kegg = kegg.sort_values('Padjust')
    # Take top N
    kegg = kegg.head(top_n)
    # Simplify pathway name (Description already present)
    kegg['pathway_short'] = kegg['Description'].str[:50]
    kegg['Tissue'] = tissue_name
    return kegg


# ============================================================
# FIGURE A: Liver KEGG Enrichment Dotplot
# ============================================================
print("Generating Figure A: Liver KEGG Enrichment Dotplot...")

liver_kegg = prepare_enrichment_data(liver, 'Liver', top_n=25)

fig_a, ax = plt.subplots(figsize=(7, len(liver_kegg) * 0.32 + 0.8))

# Color by direction
colors = [C_TFB if d == 'TFB' else C_DLY for d in liver_kegg['Direction']]

# Size by NES magnitude, color by direction
scatter = ax.scatter(
    liver_kegg['NES_abs'],
    range(len(liver_kegg)),
    s=liver_kegg['NES_abs'] * 40 + 30,  # Scale dot size
    c=colors,
    alpha=0.85,
    edgecolors='white',
    linewidth=0.5,
    zorder=2,
)

# Add FDR annotation
for i, (_, row) in enumerate(liver_kegg.iterrows()):
    fdr_text = f'FDR={row["Padjust"]:.3f}' if row['Padjust'] > 0.001 else 'FDR<0.001'
    stars = pval_stars(max(row['Padjust'], 1e-10))
    ax.text(row['NES_abs'] + 0.15, i, f'{fdr_text} {stars}',
            va='center', fontsize=6.5, color='#444444')

# Pathway labels
ax.set_yticks(range(len(liver_kegg)))
ax.set_yticklabels([r['Description'][:55] for _, r in liver_kegg.iterrows()], fontsize=7)
ax.invert_yaxis()

ax.set_xlabel('|Normalized Enrichment Score|', fontsize=9)
ax.set_title('Liver: KEGG Pathway Enrichment (DLY vs TFB @ 45 kg)', fontsize=11, fontweight='bold')

# Legend
legend_elements = [
    mpatches.Patch(facecolor=C_TFB, alpha=0.85, label='TFB-enriched (higher AA catabolism)'),
    mpatches.Patch(facecolor=C_DLY, alpha=0.85, label='DLY-enriched (signaling/growth)'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=7,
          frameon=True, fancybox=False, edgecolor='#DDDDDD')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.2, lw=0.3)

plt.tight_layout()
save_figure(fig_a, 'Fig_GSEA_liver_enrichment')
plt.close()
print("  Saved Fig_GSEA_liver_enrichment.pdf/png")


# ============================================================
# FIGURE B: Liver + Muscle NES Comparison Bar Plot
# ============================================================
print("Generating Figure B: Liver vs Muscle GSEA Comparison...")

# Combine liver and muscle, find shared pathways
liver_kegg_all = prepare_enrichment_data(liver, 'Liver', top_n=40)
muscle_kegg_all = prepare_enrichment_data(muscle, 'Muscle', top_n=40)

# Get significant pathways from liver
liver_sig = liver_kegg_all[liver_kegg_all['Padjust'] < 0.05].copy()
# Also get muscle pathways for comparison
muscle_sig = muscle_kegg_all[muscle_kegg_all['Padjust'] < 0.05].copy()

# Shared significant pathway IDs
shared_ids = set(liver_sig['Gene Set Name'].values) & set(muscle_sig['Gene Set Name'].values)
liver_only_ids = set(liver_sig['Gene Set Name'].values) - shared_ids

# Build comparison data: all liver-sig pathways + any muscle-sig
all_pathways = list(liver_sig['Gene Set Name'].values) + [
    pid for pid in muscle_sig['Gene Set Name'].values
    if pid not in set(liver_sig['Gene Set Name'].values)
]

comp_rows = []
for pid in all_pathways:
    l_row = liver[liver['Gene Set Name'] == pid]
    m_row = muscle[muscle['Gene Set Name'] == pid]

    desc = (l_row['Description'].iloc[0] if len(l_row) > 0 else
            m_row['Description'].iloc[0] if len(m_row) > 0 else pid)

    liver_nes = l_row['NES'].iloc[0] if len(l_row) > 0 else np.nan
    liver_fdr = l_row['Padjust'].iloc[0] if len(l_row) > 0 else np.nan
    muscle_nes = m_row['NES'].iloc[0] if len(m_row) > 0 else np.nan
    muscle_fdr = m_row['Padjust'].iloc[0] if len(m_row) > 0 else np.nan

    comp_rows.append({
        'pathway_id': pid,
        'description': desc,
        'liver_NES': liver_nes,
        'liver_FDR': liver_fdr,
        'muscle_NES': muscle_nes,
        'muscle_FDR': muscle_fdr,
    })

comp = pd.DataFrame(comp_rows)
# Sort by liver NES (TFB-enriched first)
comp = comp.sort_values('liver_NES', ascending=False)

# Limit to pathways with meaningful signal
comp = comp.head(35)

fig_b = plt.figure(figsize=(8, max(5, len(comp) * 0.32)))
gs = GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)
ax_l = fig_b.add_subplot(gs[0])
ax_m = fig_b.add_subplot(gs[1])

bar_height = 0.7
y_positions = range(len(comp))

# Liver plot
for i, (_, row) in enumerate(comp.iterrows()):
    nes = row['liver_NES']
    fdr = row['liver_FDR']
    if pd.notna(nes):
        color = C_TFB if nes > 0 else C_DLY
        alpha = 0.9 if fdr < 0.05 else 0.4
        ax_l.barh(i, abs(nes), height=bar_height, color=color, alpha=alpha,
                  edgecolor='white', lw=0.3)
        if fdr < 0.05:
            ax_l.text(abs(nes) + 0.08, i, pval_stars(max(fdr, 1e-10)),
                      va='center', fontsize=6, color='#333333')

# Pathway labels in the middle
for i, (_, row) in enumerate(comp.iterrows()):
    desc = row['description']
    if len(desc) > 45:
        desc = desc[:43] + '...'
    ax_l.text(-0.02, i, desc, va='center', ha='right', fontsize=6.5, color='#222222')

ax_l.set_xlabel('|NES| Liver', fontsize=9)
ax_l.set_yticks([])
ax_l.invert_yaxis()
ax_l.set_xlim(0, comp['liver_NES'].abs().max() * 1.4)

# Muscle plot
for i, (_, row) in enumerate(comp.iterrows()):
    nes = row['muscle_NES']
    fdr = row['muscle_FDR']
    if pd.notna(nes):
        color = C_TFB if nes > 0 else C_DLY
        alpha = 0.9 if fdr < 0.05 else 0.4
        ax_m.barh(i, abs(nes), height=bar_height, color=color, alpha=alpha,
                  edgecolor='white', lw=0.3)
        if fdr < 0.05:
            ax_m.text(abs(nes) + 0.08, i, pval_stars(max(fdr, 1e-10)),
                      va='center', fontsize=6, color='#333333')

ax_m.set_xlabel('|NES| Muscle', fontsize=9)
ax_m.set_yticks([])
ax_m.invert_yaxis()
ax_m.set_xlim(0, max(comp['muscle_NES'].abs().max() * 1.4, 2.0))

# Panel labels
panel_label(ax_l, 'a')
panel_label(ax_m, 'b')

fig_b.suptitle('GSEA: KEGG Pathway Enrichment — Liver vs Muscle (DLY vs TFB @ 45 kg)',
              fontsize=12, fontweight='bold', y=1.01)

# Combined legend
legend_elements = [
    mpatches.Patch(facecolor=C_TFB, alpha=0.9, label='TFB-enriched (FDR<0.05)'),
    mpatches.Patch(facecolor=C_DLY, alpha=0.9, label='DLY-enriched (FDR<0.05)'),
    mpatches.Patch(facecolor='#666666', alpha=0.4, label='Not significant'),
]
fig_b.legend(handles=legend_elements, loc='lower center', ncol=3,
             fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.03))

plt.tight_layout()
save_figure(fig_b, 'Fig_GSEA_NES_comparison')
plt.close()
print("  Saved Fig_GSEA_NES_comparison.pdf/png")


# ============================================================
# FIGURE C: Enrichment Dotplot — Panel-style for JASB
# ============================================================
print("Generating Figure C: Combined Enrichment Dotplot...")

# Build combined data: top 20 liver + top 10 muscle unique
liver_top = liver_kegg_all.head(20)[['Gene Set Name', 'Description', 'NES', 'Padjust', 'Direction']].copy()
liver_top['Tissue'] = 'Liver'

muscle_top = muscle_kegg_all.head(10)[['Gene Set Name', 'Description', 'NES', 'Padjust', 'Direction']].copy()
muscle_top['Tissue'] = 'Muscle'

combined = pd.concat([liver_top, muscle_top], ignore_index=True)
combined = combined.drop_duplicates(subset=['Gene Set Name', 'Tissue'])
combined['neg_log10_fdr'] = -np.log10(combined['Padjust'].clip(lower=1e-10))
combined['label'] = combined['Description'].str[:50]
combined = combined.sort_values(['Tissue', 'Padjust'], ascending=[True, True])

fig_c, ax = plt.subplots(figsize=(8, max(5, len(combined) * 0.28)))

# Create colors
colors = [C_TFB if d == 'TFB' else C_DLY for d in combined['Direction']]

# Scatter: x = NES (signed), y = pathways, size = -log10(FDR), color = direction
sizes = combined['neg_log10_fdr'] * 18 + 20
sc = ax.scatter(
    combined['NES'],
    range(len(combined)),
    s=sizes,
    c=colors,
    alpha=0.85,
    edgecolors='white',
    linewidth=0.4,
    zorder=3,
)

# Add pathway names
for i, (_, row) in enumerate(combined.iterrows()):
    tissue_symbol = 'L' if row['Tissue'] == 'Liver' else 'M'
    ax.text(0, i, f' {row["label"]}', va='center',
            fontsize=6.5, ha='left' if row['NES'] >= 0 else 'right')

# Layout
ax.set_yticks([])
ax.axvline(0, color='black', lw=0.5, alpha=0.4)
ax.set_xlabel('Normalized Enrichment Score (NES)', fontsize=9)
ax.set_title('GSEA Pathway Enrichment: DLY vs TFB @ 45 kg\n(L = Liver, M = Muscle)',
             fontsize=11, fontweight='bold')
ax.invert_yaxis()

# Legend for dot size
for fdr_val in [0.05, 0.01, 0.001]:
    size = (-np.log10(fdr_val)) * 18 + 20
    ax.scatter([], [], s=size, alpha=0.6, color='gray',
               label=f'FDR={fdr_val}', edgecolors='white', lw=0.3)
legend1 = ax.legend(scatterpoints=1, frameon=True, fontsize=6.5,
                    title='Dot size ~ significance', title_fontsize=7,
                    loc='lower right')

# Legend for color
legend_elements = [
    mpatches.Patch(facecolor=C_TFB, alpha=0.85, label='TFB-enriched'),
    mpatches.Patch(facecolor=C_DLY, alpha=0.85, label='DLY-enriched'),
]
legend2 = ax.legend(handles=legend_elements, loc='upper right', fontsize=7,
                    frameon=True, fancybox=False, edgecolor='#DDDDDD')
ax.add_artist(legend1)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.2, lw=0.3)

plt.tight_layout()
save_figure(fig_c, 'Fig_GSEA_combined_dotplot')
plt.close()
print("  Saved Fig_GSEA_combined_dotplot.pdf/png")


# ============================================================
# Report
# ============================================================
n_liver_sig = (liver_c['Padjust'] < 0.05).sum() if 'liver_c' in dir() else \
    len(liver[liver['Padjust'] < 0.05])
n_muscle_sig = len(muscle[muscle['Padjust'] < 0.05])

print(f"\n{'='*60}")
print("GSEA replacement figures generated:")
print(f"  Fig_GSEA_liver_enrichment.pdf   — {len(liver_kegg)} liver KEGG pathways")
print(f"  Fig_GSEA_NES_comparison.pdf     — {len(comp)} pathways liver vs muscle")
print(f"  Fig_GSEA_combined_dotplot.pdf   — {len(combined)} pathways combined")
print(f"\n  FDR<0.05: {n_liver_sig} liver, {n_muscle_sig} muscle")
print(f"  External GSEA data source: clusterProfiler/KOBAS (validated pipeline)")
print(f"{'='*60}")
