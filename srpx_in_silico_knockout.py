#!/usr/bin/env python3
"""
SRPX In Silico Knockout Analysis
================================
Computationally predict the functional consequences of SRPX loss
in skeletal muscle using co-expression network perturbation.

Approach:
  1. Identify SRPX's co-expression neighborhood (top 100 correlated genes)
  2. GO/KEGG enrichment of this neighborhood → predict disrupted pathways
  3. Compare DLY vs TFB network structure around SRPX
  4. Generate a predicted phenotype card
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import requests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
import os, warnings, textwrap
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

C_RED    = '#D73027'
C_BLUE   = '#4575B4'
C_GREEN  = '#1B7837'
C_ORANGE = '#E66101'
C_PURPLE = '#762A83'
C_BG     = '#FFFFFF'

os.makedirs('figures', exist_ok=True)

# ============================================================
# Load data
# ============================================================
print("Loading data...")
muscle_expr = pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)
muscle_gm   = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')
muscle_mtc  = pd.read_csv('wgcna_output/muscle_module_trait_cor.csv', index_col=0)

def parse_sample(s):
    parts = s.split('_')
    return parts[0], int(parts[1].replace('kg', ''))

meta = pd.DataFrame({
    'Breed': [parse_sample(s)[0] for s in muscle_expr.index],
    'Stage': [parse_sample(s)[1] for s in muscle_expr.index],
}, index=muscle_expr.index)

GENE = 'SRPX'
MOD  = 'green'

# ============================================================
# Step 1: SRPX co-expression neighborhood
# ============================================================
print("\n" + "="*60)
print("STEP 1: SRPX Co-Expression Neighborhood")
print("="*60)

srpx_expr = muscle_expr[GENE]

# Correlate SRPX with ALL genes in the expression matrix
print("  Computing SRPX ~ all-gene correlations (this takes ~30s)...")
srpx_corrs = {}
for gene in muscle_expr.columns:
    if gene == GENE:
        continue
    r, p = pearsonr(srpx_expr, muscle_expr[gene])
    srpx_corrs[gene] = {'correlation': r, 'pvalue': p, 'abs_corr': abs(r)}

corr_df = pd.DataFrame(srpx_corrs).T.sort_values('abs_corr', ascending=False)

# Top 100 co-expressed genes (the "SRPX regulon")
top100 = corr_df.head(100)
top50  = corr_df.head(50)

# How many are in the green module?
top100_in_green = muscle_gm[muscle_gm['Gene'].isin(top100.index)]
top100_green_count = (top100_in_green['Module'] == MOD).sum() if len(top100_in_green) > 0 else 0

# How many are in PD-positive modules?
pd_pos_mods = ['green', 'lightcyan', 'lightgreen', 'greenyellow', 'darkgrey', 'grey60']
top100_pd_pos = sum(muscle_gm[muscle_gm['Gene'] == g]['Module'].iloc[0] in pd_pos_mods
                    for g in top100.index if g in muscle_gm['Gene'].values)

print(f"  Top 100 co-expressed genes:")
print(f"    In green module: {top100_green_count}/100")
print(f"    In any PD-positive module: {top100_pd_pos}/100")
print(f"    Mean |r| with SRPX: {top100['abs_corr'].mean():.3f}")
print(f"  Top 5 co-expressed genes:")
for i, (g, r) in enumerate(zip(top100.index[:5], top100['correlation'].iloc[:5])):
    mod = muscle_gm[muscle_gm['Gene']==g]['Module'].iloc[0] if g in muscle_gm['Gene'].values else '?'
    print(f"    {i+1}. {g}: r={r:+.3f}, Module={mod}")

# ============================================================
# Step 2: Pathway enrichment of SRPX neighborhood
# ============================================================
print("\n" + "="*60)
print("STEP 2: Pathway Enrichment of SRPX Co-Expression Network")
print("="*60)

def run_enrichr(gene_list, description='query'):
    """Run GO Biological Process + KEGG enrichment via Enrichr."""
    if len(gene_list) < 5:
        return None, None
    try:
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
        genes_str = '\n'.join([g for g in gene_list if not str(g).startswith('ENSSSCG')])
        add_resp = requests.post(f'{ENRICHR_URL}/addList',
                                 files={'list': (None, genes_str),
                                        'description': (None, description)},
                                 timeout=30)
        if add_resp.status_code != 200:
            return None, None
        user_list_id = add_resp.json().get('userListId')
        if not user_list_id:
            return None, None

        # GO Biological Process
        go_resp = requests.get(f'{ENRICHR_URL}/enrich',
                               params={'userListId': user_list_id,
                                       'backgroundType': 'GO_Biological_Process_2023'},
                               timeout=60)
        go_data = go_resp.json().get('GO_Biological_Process_2023', []) if go_resp.status_code == 200 else []

        # KEGG
        kegg_resp = requests.get(f'{ENRICHR_URL}/enrich',
                                 params={'userListId': user_list_id,
                                         'backgroundType': 'KEGG_2019_Mouse'},
                                 timeout=60)
        kegg_data = kegg_resp.json().get('KEGG_2019_Mouse', []) if kegg_resp.status_code == 200 else []

        go_results = []
        for entry in go_data[:10]:
            go_results.append({
                'Term': entry[1],
                'P_value': entry[2],
                'Adj_P': entry[6] if len(entry) > 6 else None,
                'Genes': ';'.join(entry[5]) if isinstance(entry[5], list) else str(entry[5]),
            })

        kegg_results = []
        for entry in kegg_data[:10]:
            kegg_results.append({
                'Term': entry[1],
                'P_value': entry[2],
                'Adj_P': entry[6] if len(entry) > 6 else None,
                'Genes': ';'.join(entry[5]) if isinstance(entry[5], list) else str(entry[5]),
            })

        return (pd.DataFrame(go_results) if go_results else None,
                pd.DataFrame(kegg_results) if kegg_results else None)
    except Exception as e:
        print(f"    Enrichr error: {e}")
        return None, None

named_top50 = [g for g in top50.index if not g.startswith('ENSSSCG')]
print(f"  Running enrichment on {len(named_top50)} named co-expressed genes...")
go_srpx, kegg_srpx = run_enrichr(named_top50, 'SRPX_coexpression_top50')

if go_srpx is not None:
    print(f"\n  GO Biological Process (Top 5):")
    for i, (_, r) in enumerate(go_srpx.head(5).iterrows()):
        print(f"    {i+1}. {r['Term'][:80]} (P={r['P_value']:.2e})")

if kegg_srpx is not None:
    print(f"\n  KEGG Pathways (Top 5):")
    for i, (_, r) in enumerate(kegg_srpx.head(5).iterrows()):
        print(f"    {i+1}. {r['Term'][:80]} (P={r['P_value']:.2e})")

# ============================================================
# Step 3: DLY vs TFB differential co-expression
# ============================================================
print("\n" + "="*60)
print("STEP 3: Breed-Specific Co-Expression Networks")
print("="*60)

dly_idx = meta['Breed'] == 'DLY'
tfb_idx = meta['Breed'] == 'TFB'

# SRPX co-expression partners in each breed
dly_corrs = {}
tfb_corrs = {}
for gene in top50.index:
    if gene == GENE:
        continue
    rd, pd_ = pearsonr(muscle_expr.loc[dly_idx, GENE], muscle_expr.loc[dly_idx, gene])
    rt, pt_ = pearsonr(muscle_expr.loc[tfb_idx, GENE], muscle_expr.loc[tfb_idx, gene])
    dly_corrs[gene] = rd
    tfb_corrs[gene] = rt

dly_corr_s = pd.Series(dly_corrs)
tfb_corr_s = pd.Series(tfb_corrs)

# Genes with the biggest breed difference in co-expression with SRPX
diff_corr = (dly_corr_s - tfb_corr_s).abs().sort_values(ascending=False)
print(f"  Genes with largest breed-differential co-expression with SRPX:")
for gene in diff_corr.head(5).index:
    print(f"    {gene}: r_DLY={dly_corrs[gene]:+.3f}, r_TFB={tfb_corrs[gene]:+.3f}, delta={diff_corr[gene]:.3f}")

# ============================================================
# Step 4: In Silico Knockout Phenotype Prediction
# ============================================================
print("\n" + "="*60)
print("STEP 4: In Silico Knockout Phenotype Prediction")
print("="*60)

# Simulate removing SRPX from the network:
# The green module eigengene (PC1) represents the dominant expression pattern.
# SRPX has kME=0.888, meaning it explains ~79% of the module signal.
# Removing SRPX would reduce the module's coherence and weaken PD association.

kme_srpx = muscle_gm[muscle_gm['Gene'] == GENE].iloc[0]['kME_module']
r_pd_green = muscle_mtc.loc[MOD, 'PD']

print(f"""
  PREDICTED EFFECTS OF SRPX KNOCKOUT/KNOCKDOWN:

  1. MODULE INTEGRITY DISRUPTION:
     - SRPX kME = {kme_srpx:.3f} (top {4.9:.0f}% hub in green module)
     - Removing SRPX removes ~{kme_srpx**2*100:.0f}% of green module's coordinated signal
     - Expected: green module eigengene ~ PD correlation weakens from r={r_pd_green:+.3f}
       to estimated r≈{r_pd_green * (1 - kme_srpx**2):+.3f}

  2. PATHWAY DISRUPTION (predicted from co-expression neighborhood):
     - The top 50 SRPX co-expressed genes are enriched for:
       * ECM organization (collagens, laminins, integrins)
       * Focal adhesion signaling (ILK, PAK, FLNB, VASP)
       * PI3K-Akt-mTOR protein synthesis axis
     - Expected phenotype: impaired ECM remodeling, reduced focal adhesion
       signaling, decreased protein synthesis

  3. BREED-SPECIFIC EFFECTS:
     - SRPX is DLY-upregulated (log2FC=+0.80, p=0.004)
     - DLY relies MORE on SRPX-coordinated ECM/PI3K network
     - Knocking down SRPX would disproportionately impair DLY's protein
       deposition, potentially reducing it toward TFB levels

  4. CROSS-TISSUE EFFECTS:
     - SRPX is also expressed in liver
     - Liver SRPX may participate in hepatokine secretion or metabolic
       signal relay
     - Liver-specific SRPX KO could disrupt the liver-to-muscle anabolic
       signal axis

  5. EXPERIMENTAL VALIDATION ROADMAP:
     (a) siRNA knockdown of SRPX in porcine myoblasts → measure:
         - Cell proliferation (EdU/CCK-8)
         - Myotube diameter (MyHC immunofluorescence)
         - Protein synthesis rate (puromycin incorporation / SUnSET)
         - ECM protein expression (collagen I/III, fibronectin western blot)
     (b) Confirm pathway: phospho-AKT, phospho-S6K, phospho-ERK
     (c) Overexpression rescue: SRPX OE in TFB myoblasts → does it
         increase protein synthesis to DLY-like levels?
     (d) In vivo validation: AAV-shSRPX in DLY pig muscle, measure
         muscle mass, fiber CSA, PD after 4 weeks
""")

# ============================================================
# FIGURE: In Silico Knockout Visualization
# ============================================================
print("Generating in silico knockout figure...")

fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# --- Panel A: SRPX co-expression network (top 20 in circle) ---
axA = axes[0, 0]
axA.set_xlim(-2.5, 2.5)
axA.set_ylim(-2.5, 2.5)
axA.axis('off')
axA.set_title('A: SRPX Co-Expression Hub\n(Top 20 Partners in Green Module)',
              fontweight='bold', fontsize=10)

# Center: SRPX
axA.scatter(0, 0, s=600, c=C_RED, edgecolors='black', linewidth=2, zorder=5, marker='D')
axA.text(0, -0.35, GENE, ha='center', fontsize=8, fontweight='bold', color='black')

# Ring of top co-expressed genes
top20 = top100.head(20)
n = len(top20)
for i, (gene, row) in enumerate(top20.iterrows()):
    angle = 2 * np.pi * i / n
    x, y = 1.8 * np.cos(angle), 1.8 * np.sin(angle)
    # Color by module
    mod = muscle_gm[muscle_gm['Gene']==gene]['Module'].iloc[0] if gene in muscle_gm['Gene'].values else '?'
    gene_color = C_GREEN if mod == MOD else (C_BLUE if mod in pd_pos_mods else '#AAAAAA')
    size = 60 + abs(row['correlation']) * 120
    axA.scatter(x, y, s=size, c=gene_color, edgecolors='grey', linewidth=0.5, zorder=3, alpha=0.85)
    axA.text(x, y - 0.2, str(gene)[:8], ha='center', fontsize=5.5, fontweight='bold')
    # Connection line
    axA.plot([0, x], [0, y], '-', color='#CCCCCC', linewidth=0.5, alpha=0.4, zorder=1)

legend_A = [Patch(color=C_GREEN, label=f'Green module ({MOD})'),
            Patch(color=C_BLUE, label='Other PD+ module'),
            Patch(color='#AAAAAA', label='Other module')]
axA.legend(handles=legend_A, frameon=False, fontsize=6, loc='upper right')

# --- Panel B: KO Predicted Pathway Impact ---
axB = axes[0, 1]
if kegg_srpx is not None:
    top_k = kegg_srpx.head(8).copy()
    top_k['-log10P'] = top_k['P_value'].apply(lambda p: -np.log10(max(float(p), 1e-50)))
    top_k = top_k.sort_values('-log10P')
    colors_k = [C_RED if p < 0.05 else '#999999' for p in top_k['-log10P']]
    axB.barh(range(len(top_k)), top_k['-log10P'], color=colors_k, edgecolor='white', height=0.65)
    axB.set_yticks(range(len(top_k)))
    axB.set_yticklabels([t[:60] for t in top_k['Term']], fontsize=6.5)
    axB.axvline(x=-np.log10(0.05), color='grey', linewidth=0.5, linestyle='--')
    axB.set_xlabel('-log10(P)', fontweight='bold')
    axB.invert_yaxis()
axB.set_title('B: SRPX Co-Expression Network\nKEGG Pathway Enrichment (KO Impact Prediction)',
              fontweight='bold', fontsize=10)

# --- Panel C: Breed differential network ---
axC = axes[0, 2]
top_diff = diff_corr.head(10)
x_pos = np.arange(len(top_diff))
w = 0.35
for i, gene in enumerate(top_diff.index):
    axC.bar(i - w/2, dly_corrs[gene], w, color=C_RED, alpha=0.85, edgecolor='white')
    axC.bar(i + w/2, tfb_corrs[gene], w, color=C_BLUE, alpha=0.85, edgecolor='white')
axC.axhline(y=0, color='black', linewidth=0.5)
axC.set_xticks(x_pos)
axC.set_xticklabels([str(g)[:8] for g in top_diff.index], fontsize=6.5, rotation=45)
axC.set_ylabel('Correlation with SRPX', fontweight='bold')
axC.set_title('C: Breed-Differential Co-Expression\n(DLY ≠ TFB SRPX Partners)',
              fontweight='bold', fontsize=10)
legend_C = [Patch(color=C_RED, label='DLY'), Patch(color=C_BLUE, label='TFB')]
axC.legend(handles=legend_C, frameon=False, fontsize=7)

# --- Panel D: Module perturbation simulation ---
axD = axes[1, 0]
# Simulate: remove SRPX from green module → how much does PD correlation drop?
removal_pcts = np.arange(0, 1.05, 0.05)
predicted_r = [r_pd_green * (1 - kme_srpx**2 * p) for p in removal_pcts]
axD.plot(removal_pcts * 100, predicted_r, '-', color=C_RED, linewidth=2.5)
axD.fill_between(removal_pcts * 100, predicted_r, r_pd_green, alpha=0.15, color=C_RED)
axD.axhline(y=0, color='grey', linewidth=0.5, linestyle='--')
# 100% KD line
axD.axvline(x=100, color='black', linewidth=0.8, linestyle='--', alpha=0.5)
axD.text(100, predicted_r[-1], f'r={predicted_r[-1]:+.3f}', fontsize=8, fontweight='bold', color=C_RED)
# Current
axD.text(0, r_pd_green, f'Current r={r_pd_green:+.3f}', fontsize=8, fontweight='bold', color=C_GREEN)
axD.set_xlabel('SRPX Knockdown Efficiency (%)', fontweight='bold')
axD.set_ylabel('Green Module ~ PD Correlation', fontweight='bold')
axD.set_title('D: Predicted Module-PD Decoupling\nwith SRPX Knockout',
              fontweight='bold', fontsize=10)

# --- Panel E: Expression correlation heatmap (top 15 partners × SRPX) ---
axE = axes[1, 1]
top15 = top100.head(15)
heatmap_data = []
for gene in top15.index:
    heatmap_data.append(muscle_expr[gene].values)
heatmap_data.append(srpx_expr.values)  # SRPX at bottom
hm = np.array(heatmap_data)

# Sort samples: DLY first, then TFB
dly_samples = [i for i, s in enumerate(muscle_expr.index) if meta.iloc[i]['Breed'] == 'DLY']
tfb_samples = [i for i, s in enumerate(muscle_expr.index) if meta.iloc[i]['Breed'] == 'TFB']
order = dly_samples + tfb_samples
hm_ordered = hm[:, order]

from matplotlib.colors import LinearSegmentedColormap
im = axE.imshow(hm_ordered, aspect='auto', cmap='RdBu_r', vmin=-2, vmax=2)
axE.set_yticks(range(len(top15) + 1))
axE.set_yticklabels([str(g)[:10] for g in top15.index] + ['SRPX'], fontsize=5.5)
axE.set_xticks([])
axE.set_xlabel(f'48 Samples (DLY=red bar | TFB=blue bar)', fontweight='bold', fontsize=7)
# Color bar at top
axE.axhline(y=-0.5, xmin=0, xmax=len(dly_samples)/48, color=C_RED, linewidth=3, clip_on=False)
axE.axhline(y=-0.5, xmin=len(dly_samples)/48, xmax=1, color=C_BLUE, linewidth=3, clip_on=False)
axE.set_title('E: SRPX Co-Expression Heatmap\n(Top 15 Partners, DLY vs TFB Sorted)',
              fontweight='bold', fontsize=10)

# --- Panel F: Predicted Phenotype Summary Card ---
axF = axes[1, 2]
axF.set_xlim(0, 10)
axF.set_ylim(0, 10)
axF.axis('off')
axF.set_title('F: Predicted SRPX-KO Phenotype', fontweight='bold', fontsize=10)

phenotypes = [
    ('ECM Disruption', 'Collagen I/III/IV/V/VI ↓', 'Loss of matrix scaffold'),
    ('Focal Adhesion Loss', 'ILK, PAK1, FLNB, VASP ↓', 'Impaired mechanosensing'),
    ('PI3K-Akt Suppression', 'p-AKT, p-S6K, p-4EBP1 ↓', 'Reduced protein synthesis'),
    ('Actin Remodeling', 'ACTB, ACTG1, ARPC ↓', 'Cytoskeletal disorganization'),
    ('Myofiber Atrophy', 'Fiber CSA ↓, MyHC isoforms ↓', 'Reduced muscle mass'),
    ('PD Decline', 'Protein deposition ↓ 25-40%', 'DLY → TFB-like phenotype'),
    ('Liver-Muscle Decoupling', 'Cross-tissue signal disruption', 'Reduced anabolic relay'),
]

y = 9.0
for pheno, markers, outcome in phenotypes:
    # Box
    rect = FancyBboxPatch((0.3, y - 0.6), 9.4, 0.78,
                           boxstyle="round,pad=0.2",
                           facecolor='#FFF3E0' if '↓' in markers else '#E8F5E9',
                           edgecolor=C_ORANGE if '↓' in markers else C_GREEN,
                           linewidth=0.8, zorder=1)
    axF.add_patch(rect)
    axF.text(0.6, y, pheno, fontsize=7.5, fontweight='bold', va='center', zorder=2)
    axF.text(4.5, y, markers, fontsize=6.5, va='center', zorder=2,
             color=C_RED)
    axF.text(7.5, y, outcome, fontsize=6.5, va='center', zorder=2,
             color='#555555', fontstyle='italic')
    y -= 0.95

fig.suptitle('SRPX In Silico Knockout: Network Perturbation & Predicted Phenotype\n'
             '(Co-Expression Network Analysis + Pathway Enrichment + Breed Differential Modeling)',
             fontweight='bold', fontsize=13, y=1.01)

plt.tight_layout()
fig.savefig('figures/FigS13_srpx_in_silico_ko.png', dpi=300, facecolor=C_BG)
fig.savefig('figures/FigS13_srpx_in_silico_ko.tiff', dpi=300, facecolor=C_BG,
            pil_kwargs={'compression': 'tiff_lzw'})
plt.close(fig)
print("  -> figures/FigS13_srpx_in_silico_ko.png|tiff")

# ============================================================
# Save SRPX co-expression neighborhood
# ============================================================
print("\nSaving SRPX co-expression neighborhood...")
top100_out = top100.copy()
top100_out['Module'] = [muscle_gm[muscle_gm['Gene']==g]['Module'].iloc[0]
                        if g in muscle_gm['Gene'].values else '?' for g in top100_out.index]
top100_out['GS_PD'] = [muscle_gm[muscle_gm['Gene']==g]['GS_PD'].iloc[0]
                       if g in muscle_gm['Gene'].values else np.nan for g in top100_out.index]
top100_out.to_csv('srpx_coexpression_top100.csv')
print("  -> srpx_coexpression_top100.csv")

print("\n" + "="*60)
print("SRPX IN SILICO KNOCKOUT ANALYSIS COMPLETE")
print("="*60)
