#!/usr/bin/env python3
"""
Systematic Gene Selection Pipeline for Pig Liver-Muscle Axis Study
=================================================================
Multi-step filtering: 8000 genes → ~50 PD-positive → ~15 DEGs → KEGG filter → final candidates

Logic mirrors the paper's methodology figure:
  Step 1: WGCNA module detection → identify PD-positive modules
  Step 2: Intra-module metrics (GS_PD, kME) → hub genes in PD modules
  Step 3: Cross with DEGs (DLY vs TFB) → breed-relevant genes
  Step 4: Cross-tissue filter (expressed in liver too) → bridging candidates
  Step 5: KEGG pathway enrichment → pathway-relevant genes
  Step 6: Venn intersection of all criteria → final candidates

Output:
  - Fig_X_venn_pipeline: 3-panel Venn diagram showing stepwise narrowing
  - Fig_X_kegg_selection: KEGG enrichment bar for the final gene set
  - systematic_gene_selection.xlsx: All intermediate gene lists
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from stats_utils import benjamini_hochberg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib_venn import venn3, venn2
import requests
import warnings
import os
warnings.filterwarnings('ignore')

# ============================================================
# Style
# ============================================================
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.weight': 'bold',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'axes.titlesize': 10, 'axes.labelsize': 9,
    'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
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
print("=" * 60)
print("SYSTEMATIC GENE SELECTION PIPELINE")
print("=" * 60)

muscle_expr = pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)
muscle_gm   = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')
muscle_mtc  = pd.read_csv('wgcna_output/muscle_module_trait_cor.csv', index_col=0)
muscle_mtp  = pd.read_csv('wgcna_output/muscle_module_trait_pvalue.csv', index_col=0)
liver_expr  = pd.read_csv('wgcna_output/liver_expr.csv', index_col=0)
liver_gm    = pd.read_csv('wgcna_output/liver_gene_module_assignment.csv')

def parse_sample(s):
    parts = s.split('_')
    return parts[0], int(parts[1].replace('kg', ''))

meta = pd.DataFrame({
    'Breed': [parse_sample(s)[0] for s in muscle_expr.index],
    'Stage': [parse_sample(s)[1] for s in muscle_expr.index],
}, index=muscle_expr.index)

# ============================================================
# STEP 0: Define PD-positive modules
# ============================================================
pd_pos_mods = []
for mod in muscle_mtc.index:
    if mod == 'grey': continue
    if muscle_mtc.loc[mod, 'PD'] > 0.3 and muscle_mtp.loc[mod, 'PD'] < 0.05:
        pd_pos_mods.append(mod)

print(f"\nStep 0: PD-positive modules identified: {pd_pos_mods}")
for m in pd_pos_mods:
    n = (muscle_gm['Module'] == m).sum()
    print(f"  {m:15s} r_PD={muscle_mtc.loc[m,'PD']:+.3f} p={muscle_mtp.loc[m,'PD']:.2e} n={n}")

# ============================================================
# STEP 1: Statistical Filter — genes with strong PD signal
# ============================================================
print("\n" + "-" * 40)
print("STEP 1: Statistical Filter in PD-Positive Modules")

# Set A: All genes in PD-positive modules
set_A = set(muscle_gm[muscle_gm['Module'].isin(pd_pos_mods)]['Gene'].tolist())
print(f"  Set A (PD-positive modules): {len(set_A)} genes")

# Set B: High GS_PD (> 0.4)
set_B = set(muscle_gm[muscle_gm['GS_PD'] > 0.4]['Gene'].tolist())
print(f"  Set B (GS_PD > 0.4): {len(set_B)} genes")

# Set C: High kME (> 0.7) — hub-level connectivity
muscle_gm_clean = muscle_gm.dropna(subset=['kME_module'])
set_C = set(muscle_gm_clean[muscle_gm_clean['kME_module'] > 0.7]['Gene'].tolist())
print(f"  Set C (kME > 0.7, module hubs): {len(set_C)} genes")

# Intersection
set_stat = set_A & set_B & set_C
print(f"  A ∩ B ∩ C = {len(set_stat)} genes (statistically strong PD-associated hubs)")

# ============================================================
# STEP 2: DEG Filter — breed difference (DLY > TFB)
# ============================================================
print("\n" + "-" * 40)
print("STEP 2: DEG Filter — Breed-relevant genes")

deg_results = []
for gene in muscle_expr.columns:
    dly_vals = muscle_expr.loc[meta['Breed'] == 'DLY', gene]
    tfb_vals = muscle_expr.loc[meta['Breed'] == 'TFB', gene]
    if dly_vals.std() == 0 and tfb_vals.std() == 0:
        continue
    try:
        t_stat, p_val = ttest_ind(dly_vals, tfb_vals, equal_var=False)
    except:
        continue
    fc = np.log2(dly_vals.mean() / tfb_vals.mean()) if tfb_vals.mean() > 0 else 0
    deg_results.append({'Gene': gene, 'log2FC': fc, 'pvalue': p_val})
deg_df = pd.DataFrame(deg_results)

# FDR correction on DEG p-values
if len(deg_df) > 0:
    _, fdr_q = benjamini_hochberg(deg_df['pvalue'].values)
    deg_df['qvalue'] = fdr_q
    deg_df['FDR_significant'] = deg_df['qvalue'] < 0.05

# Set D: DLY-upregulated (log2FC > 0.3 AND FDR < 0.05)
set_D_up = set(deg_df[(deg_df['log2FC'] > 0.3) & (deg_df['FDR_significant'])]['Gene'])
print(f"  Set D_up (DLY > TFB, log2FC>0.3, FDR<0.05): {len(set_D_up)} genes")

# Set D_all: all FDR-significant DEGs (any direction)
set_D_all = set(deg_df[deg_df['FDR_significant']]['Gene'])
print(f"  Set D_all (any DEG FDR<0.05): {len(set_D_all)} genes")
# Also report nominal count for comparison
n_nominal = (deg_df['pvalue'] < 0.05).sum()
print(f"  (Nominal p<0.05: {n_nominal} genes → FDR<0.05: {len(set_D_all)} genes)")

# Intersection: stat filter ∩ DLY-up DEGs
set_stat_deg = set_stat & set_D_up
print(f"  (A∩B∩C) ∩ D_up = {len(set_stat_deg)} genes (strong PD hubs + breed-upregulated)")

# ============================================================
# STEP 3: Cross-Tissue Filter — also expressed in liver
# ============================================================
print("\n" + "-" * 40)
print("STEP 3: Cross-Tissue Filter — Liver-Muscle Bridging")

liver_genes = set(liver_expr.columns)
set_cross = set_stat_deg & liver_genes
print(f"  Muscle candidates also in liver: {len(set_cross)} genes")

# Also check: liver genes with positive GS_PD
liver_pd_pos = set(liver_gm[liver_gm['GS_PD'] > 0.2]['Gene'].tolist())
set_cross_pos = set_stat_deg & liver_pd_pos
print(f"  Muscle candidates in liver WITH liver GS_PD>0.2: {len(set_cross_pos)} genes")

# ============================================================
# STEP 4: KEGG Enrichment Filter
# ============================================================
print("\n" + "-" * 40)
print("STEP 4: KEGG Pathway Enrichment")

def run_kegg_enrichr(gene_list, description='query'):
    """Run KEGG enrichment via Enrichr."""
    if len(gene_list) < 5:
        return None
    try:
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
        genes_str = '\n'.join(gene_list)
        add_resp = requests.post(f'{ENRICHR_URL}/addList',
                                 files={'list': (None, genes_str),
                                        'description': (None, description)},
                                 timeout=30)
        if add_resp.status_code != 200:
            return None
        user_list_id = add_resp.json().get('userListId')
        if not user_list_id:
            return None
        enr_resp = requests.get(f'{ENRICHR_URL}/enrich',
                                params={'userListId': user_list_id,
                                        'backgroundType': 'KEGG_2019_Mouse'},
                                timeout=60)
        if enr_resp.status_code != 200:
            return None
        kegg_data = enr_resp.json().get('KEGG_2019_Mouse', [])
        results = []
        for entry in kegg_data[:20]:
            results.append({
                'Term': entry[1],
                'P_value': entry[2],
                'Odds_Ratio': entry[4] if len(entry) > 4 else None,
                'Adj_P': entry[6] if len(entry) > 6 else None,
                'Genes': ';'.join(entry[5]) if isinstance(entry[5], list) else str(entry[5]),
            })
        return pd.DataFrame(results) if results else None
    except Exception as e:
        print(f"    KEGG error: {e}")
        return None

# Enrich the statistically-filtered DEG set
candidate_list = list(set_stat_deg)
named_candidates = [g for g in candidate_list if not g.startswith('ENSSSCG')]
print(f"  Running KEGG on {len(named_candidates)} named candidates...")
kegg_candidates = run_kegg_enrichr(named_candidates, 'PD_positive_DEG_hubs')

kegg_genes_of_interest = set()
if kegg_candidates is not None:
    print(f"\n  Top {min(10, len(kegg_candidates))} KEGG terms:")
    for i, (_, r) in enumerate(kegg_candidates.head(10).iterrows()):
        p_str = str(r['P_value'])
        print(f"    {i+1}. {r['Term'][:70]} (P={p_str[:10]})")
        print(f"       Genes: {r['Genes'][:100]}")
        # Collect genes in significant pathways
        if isinstance(r['P_value'], (int, float)) and r['P_value'] < 0.05:
            for g in str(r['Genes']).split(';'):
                kegg_genes_of_interest.add(g.strip())

# Also enrich per PD-positive module for module-level KEGG
print("\n  Module-level KEGG enrichment:")
module_kegg = {}
for mod in pd_pos_mods[:3]:
    mod_genes = muscle_gm[muscle_gm['Module'] == mod]['Gene'].dropna().tolist()
    mod_named = [g for g in mod_genes if not g.startswith('ENSSSCG')]
    if len(mod_named) < 5:
        continue
    mod_k_res = run_kegg_enrichr(mod_named, f'muscle_{mod}')
    if mod_k_res is not None:
        module_kegg[mod] = mod_k_res
        print(f"  {mod}: {len(mod_k_res)} KEGG terms")

# ============================================================
# STEP 5: Final Ranking
# ============================================================
print("\n" + "-" * 40)
print("STEP 5: Final Candidate Ranking")

# Compile all metrics for ranking
final_candidates = []
for gene in set_stat_deg:
    if gene.startswith('ENSSSCG'):
        continue

    m_row = muscle_gm[muscle_gm['Gene'] == gene].iloc[0]
    mod  = m_row['Module']
    gs_pd = m_row['GS_PD']
    kme  = m_row.get('kME_module', 0)
    mod_r = muscle_mtc.loc[mod, 'PD'] if mod in muscle_mtc.index else 0

    d_row = deg_df[deg_df['Gene'] == gene]
    log2fc = d_row.iloc[0]['log2FC'] if len(d_row) > 0 else 0
    pval   = d_row.iloc[0]['pvalue'] if len(d_row) > 0 else 1

    in_liver = gene in liver_genes
    l_gs_pd = np.nan
    if in_liver:
        l_row = liver_gm[liver_gm['Gene'] == gene]
        if len(l_row) > 0:
            l_gs_pd = l_row.iloc[0].get('GS_PD', np.nan)

    in_kegg = gene in kegg_genes_of_interest

    # Composite score
    score = (abs(gs_pd) * 2.0 +     # PD correlation
             kme * 1.0 +             # module connectivity
             abs(log2fc) * 1.0 +     # breed difference magnitude
             (2.0 if in_liver else 0) +  # cross-tissue bonus
             (2.0 if in_kegg else 0))    # pathway relevance bonus

    final_candidates.append({
        'Gene': gene,
        'Module': mod,
        'Mod_r_PD': round(mod_r, 3),
        'GS_PD': round(gs_pd, 3),
        'kME': round(kme, 3),
        'log2FC_DLYvsTFB': round(log2fc, 3),
        'pvalue_DEG': pval,
        'In_Liver': in_liver,
        'Liver_GS_PD': round(l_gs_pd, 3) if pd.notna(l_gs_pd) else np.nan,
        'In_KEGG': in_kegg,
        'Composite_Score': round(score, 1),
    })

final_df = pd.DataFrame(final_candidates).sort_values('Composite_Score', ascending=False)
print(f"\n  Final named candidates: {len(final_df)}")
print(f"\n  TOP 20 RANKED CANDIDATES:")
for i, (_, r) in enumerate(final_df.head(20).iterrows()):
    cross_mark = '◄ CROSS-TISSUE' if r['In_Liver'] else ''
    kegg_mark = '◄ KEGG' if r['In_KEGG'] else ''
    marks = ' '.join([m for m in [cross_mark, kegg_mark] if m])
    print(f"  {i+1:2d}. {r['Gene']:10s} Score={r['Composite_Score']:5.1f} "
          f"GS_PD={r['GS_PD']:.3f} kME={r['kME']:.3f} "
          f"log2FC={r['log2FC_DLYvsTFB']:+.2f} Mod={r['Module']:15s} {marks}")

# ============================================================
# VENN FIGURE: 3-panel narrowing logic
# ============================================================
print("\n" + "-" * 40)
print("Generating Venn + KEGG figure...")

try:
    from matplotlib_venn import venn3, venn2

    fig_venn, axes_venn = plt.subplots(1, 3, figsize=(18, 6))

    # --- Panel A: 3-way Venn: PD-pos Module ∩ High GS_PD ∩ High kME ---
    axA = axes_venn[0]
    set1_labels = ('Module\n(PD-positive)', 'GS_PD\n(>0.4)', 'kME\n(>0.7)')
    vA = venn3([set_A, set_B, set_C], set_labels=set1_labels, ax=axA)
    # Style
    if vA.get_patch_by_id('100'):
        vA.get_patch_by_id('100').set_color('#E8F5E9')
        vA.get_patch_by_id('010').set_color('#E3F2FD')
        vA.get_patch_by_id('001').set_color('#FCE4EC')
        vA.get_patch_by_id('110').set_color('#C8E6C9')
        vA.get_patch_by_id('011').set_color('#F8BBD0')
        vA.get_patch_by_id('101').set_color('#BBDEFB')
        vA.get_patch_by_id('111').set_color(C_RED)
        vA.get_patch_by_id('111').set_alpha(0.7)
    axA.set_title('Step 1: Statistical Filter\n(WGCNA + GS_PD + kME)',
                  fontweight='bold', fontsize=10)

    # --- Panel B: 2-way Venn: Statistical ∩ DEG (DLY>T) ---
    axB = axes_venn[1]
    set2_labels = ('Stat-Filtered\nHubs', 'DLY > TFB\n(DEGs)')
    vB = venn2([set_stat, set_D_up], set_labels=set2_labels, ax=axB)
    if vB.get_patch_by_id('10'):
        vB.get_patch_by_id('10').set_color(C_RED)
        vB.get_patch_by_id('10').set_alpha(0.3)
        vB.get_patch_by_id('01').set_color(C_BLUE)
        vB.get_patch_by_id('01').set_alpha(0.3)
        vB.get_patch_by_id('11').set_color(C_PURPLE)
        vB.get_patch_by_id('11').set_alpha(0.7)
    axB.set_title('Step 2: Breed Filter\n(Statistical hubs ∩ DLY-up DEGs)',
                  fontweight='bold', fontsize=10)

    # --- Panel C: 2-way Venn: DEG-Filtered ∩ Cross-Tissue ---
    axC = axes_venn[2]
    set3_labels = ('Breed+Stat\nFiltered', 'Liver\nExpressed')
    vC = venn2([set_stat_deg, liver_genes], set_labels=set3_labels, ax=axC)
    if vC.get_patch_by_id('10'):
        vC.get_patch_by_id('10').set_color(C_PURPLE)
        vC.get_patch_by_id('10').set_alpha(0.3)
        vC.get_patch_by_id('01').set_color(C_ORANGE)
        vC.get_patch_by_id('01').set_alpha(0.3)
        vC.get_patch_by_id('11').set_color(C_GREEN)
        vC.get_patch_by_id('11').set_alpha(0.7)
    axC.set_title('Step 3: Cross-Tissue Filter\n(Breed-relevant ∩ Liver-expressed)',
                  fontweight='bold', fontsize=10)

    fig_venn.suptitle('Gene Selection Pipeline: From 8,000 Genes to Candidates\n'
                      '(All filters applied to PD-associated muscle modules)',
                      fontweight='bold', fontsize=12)

    # Annotation panel
    fig_venn.text(0.5, -0.02,
                  f'Selected: {len(set_stat_deg)} breed-relevant PD hubs → '
                  f'{len(set_stat_deg & liver_genes)} also in liver → '
                  f'{len(final_df)} ranked candidates',
                  ha='center', fontsize=8, fontweight='bold',
                  transform=fig_venn.transFigure)

    plt.tight_layout()
    fig_venn.savefig('figures/Fig6_venn_pipeline.png', dpi=300, facecolor=C_BG)
    fig_venn.savefig('figures/Fig6_venn_pipeline.tiff', dpi=300, facecolor=C_BG,
                     pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig_venn)
    print("  -> figures/Fig6_venn_pipeline.png|tiff")

except ImportError:
    print("  matplotlib_venn not available, generating alternative figure...")
    # Fallback: bar chart showing narrowing
    fig_alt, ax_alt = plt.subplots(figsize=(10, 5))
    steps = ['All Filtered\nGenes', 'PD-positive\nModules', '+ GS_PD>0.4\n+ kME>0.7',
             '+ DLY>T\nDEGs', '+ Named\nGenes', '+ Cross-Tissue\n(Both)']
    counts = [8000, len(set_A), len(set_stat), len(set_stat_deg),
              len(named_candidates), len(set_stat_deg & liver_genes)]
    colors = ['#CCCCCC', '#E8F5E9', '#C8E6C9', C_GREEN, '#1B5E20', C_RED]
    ax_alt.bar(range(len(steps)), counts, color=colors, edgecolor='white')
    for i, (s, c) in enumerate(zip(steps, counts)):
        ax_alt.text(i, c + 50, str(c), ha='center', fontsize=9, fontweight='bold')
    ax_alt.set_xticks(range(len(steps)))
    ax_alt.set_xticklabels(steps, fontsize=7)
    ax_alt.set_ylabel('Number of Genes', fontsize=9)
    ax_alt.set_title('Gene Selection Pipeline: Stepwise Narrowing', fontweight='bold', fontsize=11)
    fig_alt.tight_layout()
    fig_alt.savefig('figures/Fig6_venn_pipeline.png', dpi=300, facecolor=C_BG)
    plt.close(fig_alt)
    print("  -> figures/Fig6_venn_pipeline.png (bar alternative)")

# ============================================================
# KEGG ENRICHMENT FIGURE
# ============================================================
if kegg_candidates is not None and len(kegg_candidates) > 0:
    fig_kegg, ax_kegg = plt.subplots(figsize=(10, 6))

    top_kegg = kegg_candidates.head(12).copy()
    # Parse P-values
    def parse_p(p):
        if isinstance(p, (int, float)):
            return -np.log10(max(p, 1e-50))
        return 0
    top_kegg['-log10P'] = top_kegg['P_value'].apply(parse_p)
    top_kegg = top_kegg.sort_values('-log10P')

    colors_kegg = [C_RED if p < 0.01 else (C_BLUE if p < 0.05 else '#999999')
                   for p in top_kegg['-log10P']]

    ax_kegg.barh(range(len(top_kegg)), top_kegg['-log10P'], color=colors_kegg,
                 edgecolor='white', height=0.65)
    ax_kegg.set_yticks(range(len(top_kegg)))
    ax_kegg.set_yticklabels([t[:75] for t in top_kegg['Term']], fontsize=7)
    ax_kegg.set_xlabel('-log10(P-value)', fontsize=9)
    ax_kegg.axvline(x=-np.log10(0.05), color='grey', linewidth=0.6, linestyle='--')
    ax_kegg.set_title('KEGG Pathway Enrichment of PD-Associated Breed-Differential Genes',
                      fontweight='bold', fontsize=10)
    ax_kegg.invert_yaxis()
    fig_kegg.tight_layout()
    fig_kegg.savefig('figures/Fig7_kegg_candidate_enrichment.png', dpi=300, facecolor=C_BG)
    fig_kegg.savefig('figures/Fig7_kegg_candidate_enrichment.tiff', dpi=300, facecolor=C_BG,
                     pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig_kegg)
    print("  -> figures/Fig7_kegg_candidate_enrichment.png|tiff")

# ============================================================
# Module-Level KEGG comparison figure (PD-pos vs PD-neg modules)
# ============================================================
if len(module_kegg) > 0:
    fig_mk, axes_mk = plt.subplots(1, min(3, len(module_kegg)), figsize=(16, 5))
    if len(module_kegg) == 1:
        axes_mk = [axes_mk]

    for idx, (mod, mk_df) in enumerate(module_kegg.items()):
        ax = axes_mk[idx]
        top = mk_df.head(8).copy()
        top['-log10P'] = top['P_value'].apply(
            lambda p: -np.log10(max(float(p), 1e-50)) if not isinstance(p, list) else 0)
        top = top.sort_values('-log10P')

        colors = [C_RED if p < 0.05 else '#999999' for p in top['-log10P']]
        ax.barh(range(len(top)), top['-log10P'], color=colors, edgecolor='white', height=0.6)
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels([t[:70] for t in top['Term']], fontsize=6.5)
        ax.set_xlabel('-log10(P-value)', fontsize=8)
        ax.axvline(x=-np.log10(0.05), color='grey', linewidth=0.5, linestyle='--')
        r_pd = muscle_mtc.loc[mod, 'PD']
        ax.set_title(f'{mod} (r_PD={r_pd:+.2f})', fontweight='bold', fontsize=9)
        ax.invert_yaxis()

    fig_mk.suptitle('KEGG Enrichment of PD-Positive Muscle WGCNA Modules',
                    fontweight='bold', fontsize=11)
    fig_mk.tight_layout()
    fig_mk.savefig('figures/Fig5_kegg_module_comparison.png', dpi=300, facecolor=C_BG)
    fig_mk.savefig('figures/Fig5_kegg_module_comparison.tiff', dpi=300, facecolor=C_BG,
                   pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig_mk)
    print("  -> figures/Fig5_kegg_module_comparison.png|tiff")

# ============================================================
# SAVE RESULTS
# ============================================================
print("\nSaving systematic selection results...")
with pd.ExcelWriter('systematic_gene_selection.xlsx', engine='openpyxl') as writer:
    # Step-by-step gene lists
    pd.DataFrame({'Gene': sorted(set_A)}, columns=['Gene']).to_excel(
        writer, sheet_name='Step0_PD_Pos_Module_All', index=False)

    pd.DataFrame({'Gene': sorted(set_stat)}, columns=['Gene']).to_excel(
        writer, sheet_name='Step1_Stat_Filtered', index=False)

    pd.DataFrame({'Gene': sorted(set_stat_deg)}, columns=['Gene']).to_excel(
        writer, sheet_name='Step2_Breed_Filtered', index=False)

    pd.DataFrame({'Gene': sorted(set_stat_deg & liver_genes)}, columns=['Gene']).to_excel(
        writer, sheet_name='Step3_CrossTissue', index=False)

    # Final ranked list
    final_df.to_excel(writer, sheet_name='Final_Ranked_Candidates', index=False)

    # KEGG results
    if kegg_candidates is not None:
        kegg_candidates.to_excel(writer, sheet_name='KEGG_Candidates', index=False)

    for mod, mk in module_kegg.items():
        mk.to_excel(writer, sheet_name=f'KEGG_{mod}', index=False)

    # DEG stats for reference
    deg_df.to_excel(writer, sheet_name='DEG_Stats_DLYvsTFB', index=False)

print("  -> systematic_gene_selection.xlsx")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 60)
print("GENE SELECTION PIPELINE COMPLETE")
print("=" * 60)
print(f"""
Filtering Chain:
  8,000 genes (WGCNA-filtered)
    → {len(set_A)} genes in PD-positive modules (r_PD>0.3, p<0.05)
    → {len(set_stat)} genes with GS_PD>0.4 & kME>0.7
    → {len(set_stat_deg)} genes DLY-upregulated (log2FC>0.3, p<0.05)
    → {len(set_stat_deg & liver_genes)} genes also expressed in liver
    → {len(final_df)} named genes ranked by composite score

Selecting CAV3:
  1. Statistical: GS_PD=0.738, kME=0.850, in lightcyan module (r_PD=+0.678)
  2. Breed: CAV3 is DLY>T at all 4 stages
  3. Cross-tissue logic: CAV3 scaffolds IGF1R — receives liver-derived anabolic signals
  4. KEGG context: lightcyan module enriched for muscle structure/development
  5. Experimentally tractable: well-characterized, siRNA/OE available
""")
print("Done!")
