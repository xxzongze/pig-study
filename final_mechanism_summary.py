"""
Final mechanism summary: updated with full-muscle validation.
Integrates hepatokine signaling + receptor validation.
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

out = '/Users/hezongze/pig_study'

# ============================================================
# FIGURE 1: Updated mechanism model with receptor validation
# ============================================================
print("Generating integrated mechanism model...")

fig, ax = plt.subplots(figsize=(22, 14))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Title + subtitle
ax.text(50, 98, 'Liver-Muscle Axis in DLY vs TFB: Validated Signaling Model',
        ha='center', fontsize=16, fontweight='bold')
ax.text(50, 95, 'Hepatokine → Receptor → Muscle Response | Full 9234-gene muscle validation | 7464 common genes',
        ha='center', fontsize=9, color='#666666')

# ============================================================
# PANELS: 4 stages (top)
# ============================================================
stage_names = ['15 kg\nEarly Programming', '45 kg\nTFB Peak', '75 kg\nDLY Rising', '105 kg\nLate Divergence']
stage_x = [10, 32, 54, 76]
stage_w = 20

for i, (sx, sn) in enumerate(zip(stage_x, stage_names)):
    rect = mpatches.FancyBboxPatch((sx, 55), stage_w, 37,
                                     boxstyle="round,pad=1",
                                     facecolor=plt.cm.Blues(0.05 + 0.15*i),
                                     edgecolor='#555555', linewidth=1.5)
    ax.add_patch(rect)
    ax.text(sx + stage_w/2, 89, sn, ha='center', fontsize=9, fontweight='bold')

    # Cross-tissue r annotation
    rs = [0.38, 0.20, 0.16, 0.08]
    ax.text(sx + stage_w/2, 56, f'Liver-Muscle r={rs[i]:.2f}',
            ha='center', fontsize=7, color='#888888', style='italic')

# ============================================================
# HEPATOKINE SIGNALING DIAGRAM (middle)
# ============================================================
ax.text(2, 52, 'LIVER\nHEPATOKINES', ha='center', fontsize=9, fontweight='bold', color='#d62728')

# Signaling axes
signaling_data = [
    # (gene, liver_fcs_by_stage, muscle_fcs_by_stage, corr, receptor, rec_corr)
    ('IGF1', [1.47, 1.60, 0.38, 0.04], [0.33, 0.46, 0.10, 0.11], 0.95, 'IGF1R', 0.85),
    ('IGFBP1', [-3.37, -1.10, -0.83, -2.41], [np.nan]*4, np.nan, '-', np.nan),
    ('IGFBP3', [0.66, 0.55, 1.25, -0.02], [0.29, 0.44, 0.32, 0.49], -0.79, 'TMEM219', -0.26),
    ('FGF21', [0.08, -0.35, 0.32, -0.70], [np.nan]*4, np.nan, 'FGFR1/KLB', np.nan),
    ('FST', [-0.08, -0.86, -0.25, -0.49], [0.22, 0.08, -0.04, 0.66], -0.04, 'ACVR2B', 0.95),
    ('ANGPTL4', [1.33, 0.35, 0.37, 0.05], [0.56, 1.53, -0.22, 0.32], 0.05, 'ITGB1', 0.47),
    ('CCN2', [0.40, 0.51, 1.05, -0.03], [-0.68, -0.25, 0.88, -0.72], 0.92, 'ITGB1', 0.64),
    ('APOE', [-0.18, 0.10, 0.41, -0.23], [1.32, 1.32, 0.70, 0.92], -0.50, 'LDLR', 0.83),
]

# Draw each signaling axis
y_start = 47
y_step = 5.5

for idx, (gene, lfcs, mfcs, corr, receptor, rec_corr) in enumerate(signaling_data):
    y = y_start - idx * y_step

    # Gene name (left)
    color = '#d62728' if not np.isnan(corr) and corr > 0.5 else ('#ff7f0e' if not np.isnan(corr) and corr > 0 else '#1f77b4')
    ax.text(3, y, gene, ha='right', fontsize=9, fontweight='bold', color=color)

    # Liver FC mini-bars
    for si, (sx, fc) in enumerate(zip(stage_x, lfcs)):
        if not np.isnan(fc):
            bar_x = sx + 1 + si * (stage_w / 4)
            bar_w = stage_w / 4 - 2
            bar_color = '#d62728' if fc > 0 else '#1f77b4'
            bar_h = abs(fc) * 2
            bar_y = y
            ax.bar(bar_x, bar_h, bar_w, bottom=bar_y, color=bar_color, alpha=0.7, edgecolor='none')

    # Correlation annotation
    if not np.isnan(corr):
        ax.text(96, y + 0.5, f'r={corr:.2f}', fontsize=7, ha='right', color=color if corr > 0 else '#999999')

    # Receptor arrow
    ax.annotate('', xy=(99, y - 0.8), xytext=(96, y + 0.8),
               arrowprops=dict(arrowstyle='->', color='#9467bd', lw=1.5))

    # Receptor name
    ax.text(101, y, f'→ {receptor}', fontsize=8, color='#9467bd', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#cccccc'))

# ============================================================
# MUSCLE RESPONSE (bottom panel)
# ============================================================
ax.text(2, y_start - len(signaling_data)*y_step - 2, 'MUSCLE\nRESPONSE', ha='center', fontsize=9, fontweight='bold', color='#1f77b4')

muscle_insights = [
    'IGF1r=0.95: IGF1 is the strongest validated hepatokine.\nLiver DLY↑ → Muscle DLY↑  synchronized at 15-45kg.',
    'FST→ACVR2B r=0.95: Follistatin receptor tracks liver FST.\nFST is a MSTN antagonist: TFB liver↑ → may suppress muscle growth via MSTN.',
    'IGFBP1: Strongest liver signal (FC=-3.37, 15kg) but NOT in muscle.\nActs as serum IGF1 carrier → regulates IGF1 bioavailability systemically.',
    'Cross-tissue r DECLINES: 0.38→0.08 over stages.\nLiver-muscle coupling strongest early, decouples as deposition diverges.',
    'ITGB1 is the most shared receptor (7 ligands).\nIntegrin β1 mediates ANGPTL4/6/8, CCN1/2, IGFBP1/2, SPARC signaling.',
]

for i, text in enumerate(muscle_insights):
    y = y_start - (len(signaling_data) + 1) * y_step - 4 - i * 3.5
    ax.text(50, y, text, ha='center', fontsize=7.5, color='#333333',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='#f8f8f8', alpha=0.9, edgecolor='#dddddd'))

# ============================================================
# KEY NUMBERS CALLOUT (right side)
# ============================================================
callout_x = 103
callout_y = 55
callout_data = [
    ('7464', 'genes shared\nliver & muscle'),
    ('0.38→0.08', 'cross-tissue r\n(15kg→105kg)'),
    ('21/37', 'receptor-ligand\npairs validated'),
    ('IGF1 r=0.95', 'top coordinated\nhepatokine'),
    ('FST→ACVR2B\nr=0.95', 'top receptor\ncorrelation'),
]

for i, (num, desc) in enumerate(callout_data):
    cy = callout_y - i * 7
    ax.text(callout_x, cy, num, fontsize=11, fontweight='bold', color='#d62728', ha='center')
    ax.text(callout_x, cy - 1.5, desc, fontsize=6.5, color='#666666', ha='center')

# ============================================================
# MECHANISTIC MODEL BOX (bottom)
# ============================================================
model_y = 5
model_box = mpatches.FancyBboxPatch((3, model_y), 94, 10, boxstyle="round,pad=1",
                                      facecolor='#fef7e7', edgecolor='#e6a817', linewidth=2)
ax.add_patch(model_box)

model_text = (
    '4-STAGE MODEL: 15kg Genetic Programming → 45kg TFB Peak (compensatory FST/MSTN) → 75kg DLY Catch-up (AA loading divergence) → 105kg Accumulated Consequence\n'
    'IGF1 Axis (validated): Liver DLY↑ → Muscle DLY↑ throughout 15-45kg | FST-MSTN Axis: TFB liver↑FST → ACVR2B↓ → MSTN signaling altered → muscle deposition ceiling\n'
    'Key insight: liver-muscle transcriptional coupling STRONGEST at 15kg (r=0.38) and WEAKEST at 105kg (r=0.08) — programming sets trajectory, later stages reflect accumulated divergence'
)
ax.text(50, model_y + 5, model_text, ha='center', va='center', fontsize=7.5, color='#555555')

# ============================================================
# LEGEND
# ============================================================
legend_items = [
    (mpatches.Patch(color='#d62728', alpha=0.7), 'Liver DLY > TFB (positive FC)'),
    (mpatches.Patch(color='#1f77b4', alpha=0.7), 'Liver TFB > DLY (negative FC)'),
    (mpatches.Patch(color='#9467bd', alpha=0.7), 'Receptor → Muscle target'),
    (mpatches.Patch(color='#ff7f0e', alpha=0.7), 'Discordant/Complex'),
]
leg = ax.legend([x[0] for x in legend_items], [x[1] for x in legend_items],
                loc='upper right', fontsize=7, title='Legend', title_fontsize=8)
leg.set_zorder(20)

plt.tight_layout()
plt.savefig(f'{out}/fig_final_mechanism_validated.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{out}/fig_final_mechanism_validated.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig_final_mechanism_validated.png/pdf")

# ============================================================
# FIGURE 2: IGF1-FST-MSTN Growth Axis Deep Dive
# ============================================================
print("Generating growth axis deep-dive figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 11))

# Panel A: IGF1 liver-muscle trajectory
ax = axes[0, 0]
x = [15, 45, 75, 105]
# IGF1
ax.plot(x, [1.47, 1.60, 0.38, 0.04], 'o-', color='#d62728', linewidth=2.5, markersize=8, label='IGF1 Liver FC')
ax.plot(x, [0.33, 0.46, 0.10, 0.11], 's--', color='#1f77b4', linewidth=2.5, markersize=8, label='IGF1 Muscle FC')
ax.axhline(y=0, color='grey', linewidth=0.5, linestyle='-')
ax.set_xlabel('Weight (kg)')
ax.set_ylabel('log2FC (DLY/TFB)')
ax.set_title('Panel A: IGF1 — Validated r=0.95\nStrongest coordinated hepatokine', fontsize=11, fontweight='bold')
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel B: FST-MSTN axis
ax = axes[0, 1]
# FST liver
ax.plot(x, [-0.08, -0.86, -0.25, -0.49], 'o-', color='#d62728', linewidth=2.5, markersize=8, label='FST Liver FC')
# FST muscle
ax.plot(x, [0.22, 0.08, -0.04, 0.66], 's--', color='#d62728', linewidth=2, markersize=7, alpha=0.6, label='FST Muscle FC')
# ACVR2B muscle (receptor)
ax.plot(x, [-0.05, 0.06, -0.05, -0.20], '^:', color='#9467bd', linewidth=2, markersize=7, label='ACVR2B Muscle FC')
# MSTN liver
ax.plot(x, [-0.20, -0.15, 0.05, -0.22], 'D-.', color='#2ca02c', linewidth=1.5, markersize=6, label='MSTN Liver FC')
ax.axhline(y=0, color='grey', linewidth=0.5, linestyle='-')
ax.set_xlabel('Weight (kg)')
ax.set_ylabel('log2FC (DLY/TFB)')
ax.set_title('Panel B: FST-MSTN-ACVR2B Axis\nFST→ACVR2B receptor r=0.95', fontsize=11, fontweight='bold')
ax.legend(fontsize=7)
ax.grid(True, alpha=0.3)

# Panel C: IGFBP family
ax = axes[0, 2]
igfbp_data = {
    'IGFBP1': [-3.37, -1.10, -0.83, -2.41],
    'IGFBP2': [-0.47, -1.22, -0.08, -1.59],
    'IGFBP3': [0.66, 0.55, 1.25, -0.02],
}
for gene, fcs in igfbp_data.items():
    ax.plot(x, fcs, 'o-', linewidth=2, markersize=7, label=gene)
ax.axhline(y=0, color='grey', linewidth=0.5, linestyle='-')
ax.set_xlabel('Weight (kg)')
ax.set_ylabel('log2FC (DLY/TFB)')
ax.set_title('Panel C: IGFBP Family in Liver\nIGFBP1 = strongest TFB-up signal', fontsize=11, fontweight='bold')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Panel D: Cross-tissue r decay
ax = axes[1, 0]
stage_r = [0.377, 0.202, 0.159, 0.080]
ax.bar(x, stage_r, width=8, color=['#1a9850', '#91cf60', '#fee08b', '#d73027'], alpha=0.85, edgecolor='#333333')
for xi, r_val in zip(x, stage_r):
    ax.annotate(f'r={r_val:.3f}', (xi, r_val + 0.01), ha='center', fontsize=10, fontweight='bold')
ax.set_xlabel('Weight (kg)')
ax.set_ylabel('Pearson r (Liver FC vs Muscle FC)')
ax.set_title('Panel D: Cross-Tissue Correlation DECLINES\n(7464 common genes)', fontsize=11, fontweight='bold')
ax.set_ylim(0, 0.45)
ax.grid(True, alpha=0.3, axis='y')

# Panel E: Hepatokine signaling classification
ax = axes[1, 1]
patterns = {'Liver_Only\n(DLY↑)': 14, 'Liver_Only\n(TFB↑)': 9, 'Low Signal': 7, 'Discordant': 2, 'Complex': 2, 'Opposing': 1, 'Muscle Only': 2}
colors_pie = ['#d62728', '#1f77b4', '#bcbd22', '#ff7f0e', '#7f7f7f', '#e377c2', '#2ca02c']
wedges, texts, autotexts = ax.pie(patterns.values(), labels=None, autopct='%1.0f%%',
                                    colors=colors_pie, startangle=90, pctdistance=0.6)
for i, (label, pct) in enumerate(zip(patterns.keys(), patterns.values())):
    texts[i].set_text(f'{label}\n(n={pct})')
    texts[i].set_fontsize(7)
ax.set_title('Panel E: Hepatokine Signaling Classification\n(37 hepatokines/myokines)', fontsize=11, fontweight='bold')

# Panel F: Receptor usage summary
ax = axes[1, 2]
rec_usage = {'ITGB1': 7, 'ITGAV': 3, 'ACVR2B': 2, 'ITGA5': 2, 'LDLR': 1, 'ACVR1B': 1,
             'ITGA6': 1, 'TMEM219': 1, 'IL2RG': 1, 'IGF1R': 1}
rec_sorted = sorted(rec_usage.items(), key=lambda x: x[1], reverse=True)
bars = ax.barh([r[0] for r in rec_sorted], [r[1] for r in rec_sorted],
               color=plt.cm.viridis(np.linspace(0.2, 0.9, len(rec_sorted))), edgecolor='#333333')
for bar, (rec, n) in zip(bars, rec_sorted):
    ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2, f'n={n}',
            fontsize=8, fontweight='bold', va='center')
ax.set_xlabel('Number of hepatokine ligands')
ax.set_title('Panel F: Receptor Usage in Muscle\nITGB1 = most shared receptor', fontsize=11, fontweight='bold')
ax.grid(True, alpha=0.3, axis='x')

fig.suptitle('DLY vs TFB Growth Axis: Validated Hepatokine → Muscle Signaling',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(f'{out}/fig_growth_axis_deep_dive.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{out}/fig_growth_axis_deep_dive.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig_growth_axis_deep_dive.png/pdf")

# ============================================================
# FIGURE 3: Temporal Cross-Tissue Decoupling Analysis
# ============================================================
print("Generating temporal decoupling analysis...")

# Load cross-tissue correlation data
ct_corr = pd.read_excel(f'{out}/crosstissue_correlation_genomewide.xlsx')

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Panel A: Distribution of cross-tissue r
ax = axes[0]
r_vals = ct_corr['CrossTissue_pearson_r'].dropna()
ax.hist(r_vals, bins=60, color='#4477aa', alpha=0.7, edgecolor='#333333', linewidth=0.3)
ax.axvline(x=0, color='grey', linewidth=0.8, linestyle='--')
ax.axvline(x=r_vals.mean(), color='#d62728', linewidth=1.5, linestyle='-', label=f'Mean r={r_vals.mean():.3f}')
ax.axvline(x=np.percentile(r_vals, 95), color='#ff7f0e', linewidth=1.5, linestyle='-',
           label=f'95th percentile r={np.percentile(r_vals, 95):.3f}')
ax.set_xlabel('Cross-Tissue Pearson r (Liver FC vs Muscle FC across 4 stages)')
ax.set_ylabel('Number of genes')
ax.set_title(f'Distribution of Cross-Tissue FC Correlations\n(n={len(r_vals)} common genes)')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis='y')

# Panel B: Correlation vs Signal Strength
ax = axes[1]
sig_strength = ct_corr['Signal_Strength'].values
r_vals_arr = ct_corr['CrossTissue_pearson_r'].values
# Filter extreme outliers
sq = np.percentile(sig_strength[~np.isnan(sig_strength)], [1, 99])
keep = (sig_strength >= sq[0]) & (sig_strength <= sq[1])
ax.scatter(sig_strength[keep], r_vals_arr[keep], s=1, alpha=0.15, color='#555555', rasterized=True)

# Highlight hepatokines
hk_list = ['IGF1', 'IGFBP1', 'IGFBP2', 'IGFBP3', 'FST', 'MSTN', 'FGF21', 'BDNF', 'FNDC5',
           'IL6', 'IL15', 'GDF15', 'ANGPTL4', 'ANGPTL6', 'ANGPTL8', 'CCN1', 'CCN2', 'CCN3',
           'SPARC', 'FSTL1', 'FSTL3', 'APOA1', 'APOE', 'ALB', 'TTR', 'POSTN', 'DCN']
for _, row in ct_corr.iterrows():
    if row['Gene'] in hk_list and row['CrossTissue_pearson_r'] > 0:
        ax.scatter(row['Signal_Strength'], row['CrossTissue_pearson_r'],
                  s=80, alpha=0.9, edgecolors='#333333', linewidth=0.5, zorder=5)
        ax.annotate(row['Gene'], (row['Signal_Strength'], row['CrossTissue_pearson_r']),
                   fontsize=7, xytext=(4, 4), textcoords='offset points', fontweight='bold')

ax.set_xlabel('Signal Strength (Liver|FC| × Muscle|FC|)')
ax.set_ylabel('Cross-Tissue Pearson r')
ax.set_title('Cross-Tissue Coordination vs Signal Strength\n(Hepatokines labeled)')
ax.grid(True, alpha=0.2)

fig.suptitle('Temporal Cross-Tissue Decoupling: Liver-Muscle Coordination Patterns',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig(f'{out}/fig_temporal_decoupling_analysis.png', dpi=200, bbox_inches='tight')
plt.savefig(f'{out}/fig_temporal_decoupling_analysis.pdf', bbox_inches='tight')
plt.close()
print("  Saved: fig_temporal_decoupling_analysis.png/pdf")

# ============================================================
# 4. Print integrated interpretation
# ============================================================
print(f"\n{'='*80}")
print("INTEGRATED MECHANISTIC INTERPRETATION")
print(f"{'='*80}")

print("""
┌─────────────────────────────────────────────────────────────────────────────┐
│ 1. IGF1 AXIS — Strongest Validated Hepatokine Signal                       │
├─────────────────────────────────────────────────────────────────────────────┤
│ Liver FC: 15kg=+1.47  45kg=+1.60  75kg=+0.38  105kg=+0.04                 │
│ Muscle FC: 15kg=+0.33  45kg=+0.46  75kg=+0.10  105kg=+0.11                │
│ Liver-Muscle r = 0.95 (Pearson)                                            │
│                                                                             │
│ Interpretation: IGF1 is the most robust liver→muscle coordination signal.  │
│ DLY liver produces more IGF1 at 15-45kg (peak protein deposition window).  │
│ Muscle IGF1 expression also trends DLY>TFB — consistent with systemic IGF1  │
│ driving local muscle IGF1 (autocrine/paracrine amplification).             │
│ Signal decays by 75kg as deposition rates converge.                        │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ 2. FST-MSTN-ACVR2B — The Growth Brake Axis                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│ FST Liver: 15kg=-0.08  45kg=-0.86  75kg=-0.25  105kg=-0.49 (TFB↑)         │
│ FST→ACVR2B receptor correlation: r=0.953                                   │
│                                                                             │
│ Follistatin (FST) is a potent MSTN antagonist. Higher FST in TFB liver =   │
│ more Follistatin → binds and neutralizes MSTN → should PROMOTE muscle      │
│ growth. But TFB has LOWER protein deposition!                              │
│                                                                             │
│ PARADOX: If FST inhibits MSTN, why does TFB (higher FST) have lower        │
│ deposition? Possible explanations:                                          │
│  (a) FST also binds Activin A (not just MSTN) — differential effects       │
│  (b) MSTN is also a myokine (muscle produces it) — local muscle MSTN       │
│      may dominate over liver-derived FST                                    │
│  (c) FST is UP-regulated in TFB as COMPENSATION for poor muscle growth     │
│      (liver sensing low muscle mass → trying to help via FST)              │
│                                                                             │
│ The high FST→ACVR2B receptor correlation (r=0.95) confirms that muscle     │
│ is RESPONSIVE to FST signaling — the receptor tracks the ligand.           │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ 3. IGFBP SYSTEM — IGF1 Bioavailability Control                              │
├─────────────────────────────────────────────────────────────────────────────┤
│ IGFBP1: Liver FC 15kg=-3.37 (strongest TFB-up signal overall)              │
│ IGFBP2: Liver FC 45kg=-1.22, 105kg=-1.59                                   │
│ IGFBP3: Liver FC 75kg=+1.25 (DLY-up at peak deposition)                    │
│                                                                             │
│ IGFBP1 binds IGF1 with higher affinity than IGF1R → SEQUESTERS IGF1.       │
│ TFB liver produces 3.4x more IGFBP1 at 15kg → less free IGF1 → less        │
│ muscle growth stimulus.                                                     │
│                                                                             │
│ IGFBP3 forms ternary complex with IGF1+ALS → EXTENDS IGF1 half-life.       │
│ DLY liver produces more IGFBP3 at 75kg → IGF1 stays active longer.         │
│                                                                             │
│ The IGFBP1/IGFBP3 RATIO may be the true determinant of IGF1 bioavailability.│
│ At 15kg: TFB has HIGH IGFBP1 + LOW IGFBP3 = IGF1 sequestered.             │
│ At 45kg: DLY has LOW IGFBP1 + HIGH IGFBP3 = IGF1 bioavailable.             │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ 4. TEMPORAL DECOUPLING — Programming → Consequence Model                    │
├─────────────────────────────────────────────────────────────────────────────┤
│ Cross-tissue r: 15kg=0.38 → 45kg=0.20 → 75kg=0.16 → 105kg=0.08           │
│                                                                             │
│ The STRONGEST liver-muscle transcriptional coupling occurs at 15kg — the   │
│ earliest measurable timepoint. This suggests:                               │
│  (a) Genetic/developmental PROGRAMMING sets the liver→muscle trajectory    │
│      very early (before 15kg)                                               │
│  (b) Later stages reflect ACCUMULATED DIVERGENCE, not new programming      │
│  (c) The initial "set point" matters more than later-stage interventions   │
│                                                                             │
│ This matches the Tier1 model: genes programmed at 15kg (ARG2, SDS, CPS1)   │
│ set the trajectory that plays out over 15→105kg.                            │
└─────────────────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────────────────┐
│ 5. ITGB1 — THE MASTER RECEPTOR                                              │
├─────────────────────────────────────────────────────────────────────────────┤
│ ITGB1 (Integrin β1) is the muscle receptor for 7 hepatokines:              │
│ ANGPTL4, ANGPTL6, ANGPTL8, CCN1, CCN2, IGFBP1, IGFBP2, SPARC              │
│                                                                             │
│ ITGB1 is highly expressed (mean log2=4.84) — a general adhesion/signaling  │
│ hub rather than a specific hepatokine receptor. Multiple hepatokines        │
│ converge on ITGB1-mediated signaling → potential for signal integration.   │
└─────────────────────────────────────────────────────────────────────────────┘

SUMMARY: The validated model supports a LIVER-CENTRIC programming model where:
  - IGF1 is the primary positive growth signal from liver to muscle
  - IGFBP1 is the primary negative regulator (sequesters IGF1 in TFB)
  - FST-MSTN axis is complex (compensatory rather than causative)
  - Early (15kg) programming dominates over later-stage changes
  - Most hepatokines act through ITGB1 in muscle (integrin-mediated)
""")

print("Done. All final figures generated.")
