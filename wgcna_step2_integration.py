#!/usr/bin/env python3
"""
Step 2: WGCNA Cross-Tissue Integration + GO/KEGG Enrichment + Liver-Muscle Axis
Adapts Jia et al. (2026) and Chen et al. (2026) multi-omics framework.

Pipeline:
  1. GO/KEGG enrichment on phenotype-associated module genes
  2. Cross-tissue module eigengene correlation (liver ME ~ muscle ME)
  3. Identify liver-muscle axis bridging candidates
  4. Generate integrated publication-quality figures
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, zscore
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import requests
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("WGCNA STEP 2: CROSS-TISSUE INTEGRATION + ENRICHMENT")
print("Adapted from Jia et al. (2026) + Chen et al. (2026)")
print("=" * 70)

# ============================================================
# 1. LOAD STEP 1 RESULTS
# ============================================================
print("\n[1] Loading Step 1 results...")

liver_gm = pd.read_csv('wgcna_output/liver_gene_module_assignment.csv')
muscle_gm = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')
liver_mtc = pd.read_csv('wgcna_output/liver_module_trait_cor.csv', index_col=0)
muscle_mtc = pd.read_csv('wgcna_output/muscle_module_trait_cor.csv', index_col=0)
liver_mtp = pd.read_csv('wgcna_output/liver_module_trait_pvalue.csv', index_col=0)
muscle_mtp = pd.read_csv('wgcna_output/muscle_module_trait_pvalue.csv', index_col=0)
liver_hub = pd.read_csv('wgcna_output/liver_hub_genes.csv')
muscle_hub = pd.read_csv('wgcna_output/muscle_hub_genes.csv')

print(f"  Liver: {len(liver_gm)} genes in {liver_gm['Module'].nunique()} modules")
print(f"  Muscle: {len(muscle_gm)} genes in {muscle_gm['Module'].nunique()} modules")

# ============================================================
# 2. DEFINE FOCAL MODULES
# ============================================================
print("\n[2] Defining focal modules for enrichment...")

# Positive regulatory modules (user's focus: 正向调控蛋白沉积)
muscle_pos_mods = []
for mod in muscle_mtc.index:
    if mod == 'grey':
        continue
    if muscle_mtc.loc[mod, 'PD'] > 0.3 and muscle_mtp.loc[mod, 'PD'] < 0.05:
        muscle_pos_mods.append(mod)

# Negative regulatory modules (atrophy/proteolysis)
muscle_neg_mods = []
for mod in muscle_mtc.index:
    if mod == 'grey':
        continue
    if muscle_mtc.loc[mod, 'PD'] < -0.3 and muscle_mtp.loc[mod, 'PD'] < 0.05:
        muscle_neg_mods.append(mod)

# Liver modules with any trait signal
liver_focal = []
for mod in liver_mtc.index:
    if mod == 'grey':
        continue
    if abs(liver_mtc.loc[mod, 'PD']) > 0.2 or abs(liver_mtc.loc[mod, 'Urea']) > 0.1:
        liver_focal.append(mod)

print(f"  Muscle POS modules: {muscle_pos_mods}")
print(f"  Muscle NEG modules: {muscle_neg_mods}")
print(f"  Liver focal modules: {liver_focal}")

# ============================================================
# 3. GO/KEGG ENRICHMENT ON PHENOTYPE-ASSOCIATED MODULES
# ============================================================
print("\n[3] GO/KEGG enrichment on phenotype-associated modules...")

# Use pig gene symbols — map to human for enrichment via gseapy
# gseapy supports pig (sus scrofa) but enrichr is more reliable with human symbols
# For pig genes, we'll use the enrichr 'Pig' library if available, else human mapping

def run_enrichment_enrichr(gene_list, module_name, tissue, gene_set='KEGG_2019'):
    """Run enrichment using Enrichr REST API (supports pig genes)."""
    if len(gene_list) < 5:
        return None

    try:
        import requests
        import json

        # Enrichr API
        ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
        genes_str = '\n'.join(gene_list)

        # Step 1: Submit gene list
        add_response = requests.post(f'{ENRICHR_URL}/addList', files={
            'list': (None, genes_str),
            'description': (None, f'{tissue}_{module_name}')
        }, timeout=30)
        if add_response.status_code != 200:
            print(f"    Warning: Enrichr addList failed for {module_name}: {add_response.status_code}")
            return None

        data = add_response.json()
        user_list_id = data.get('userListId')
        if not user_list_id:
            return None

        # Step 2: Get enrichment results
        enr_response = requests.get(
            f'{ENRICHR_URL}/enrich',
            params={'userListId': user_list_id, 'backgroundType': gene_set},
            timeout=60
        )
        if enr_response.status_code != 200:
            return None

        results = enr_response.json()
        if gene_set not in results or len(results[gene_set]) == 0:
            return None

        # Parse results
        rows = []
        for entry in results[gene_set]:
            rows.append({
                'Gene_Set': gene_set,
                'Term': entry[1],
                'Overlap': entry[3],
                'P-value': entry[2],
                'Adjusted P-value': entry[5],
                'Odds Ratio': entry[4],
                'Combined Score': entry[6],
                'Genes': entry[7],
                'Module': module_name,
                'Tissue': tissue,
            })

        return pd.DataFrame(rows)

    except Exception as e:
        print(f"    Warning: {gene_set} enrichment failed for {module_name}: {e}")
    return None

all_enrichments = []

for tissue, gm, focal_mods in [
    ('Liver', liver_gm, liver_focal),
    ('Muscle', muscle_gm, muscle_pos_mods + muscle_neg_mods)
]:
    for mod in focal_mods:
        mod_genes = gm[gm['Module'] == mod]['Gene'].dropna().tolist()
        mod_genes_clean = [g for g in mod_genes if not g.startswith('ENSSSCG')]

        print(f"  {tissue} {mod}: {len(mod_genes_clean)} named genes for enrichment")

        if len(mod_genes_clean) < 5:
            continue

        # KEGG enrichment
        kegg_res = run_enrichment_enrichr(mod_genes_clean, mod, tissue, 'KEGG_2019')
        if kegg_res is not None:
            all_enrichments.append(kegg_res)

        # GO BP enrichment
        go_res = run_enrichment_enrichr(mod_genes_clean, mod, tissue, 'GO_Biological_Process_2023')
        if go_res is not None:
            all_enrichments.append(go_res)

# Combine all enrichment results
if all_enrichments:
    enrich_df = pd.concat(all_enrichments, ignore_index=True)
    print(f"\n  Total enrichment terms: {len(enrich_df)}")
else:
    enrich_df = pd.DataFrame()
    print("\n  No enrichment results (all modules failed)")

# ============================================================
# 4. CROSS-TISSUE MODULE EIGENGENE CORRELATION
# ============================================================
print("\n[4] Computing cross-tissue module eigengene correlations...")

# Read expression data to compute module eigengenes
liver_expr = pd.read_csv('wgcna_output/liver_expr.csv', index_col=0)
muscle_expr = pd.read_csv('wgcna_output/muscle_expr.csv', index_col=0)

def compute_module_eigengenes(expr, gene_module):
    """Compute module eigengenes (PC1) for each module."""
    mes = {}
    for mod in gene_module['Module'].unique():
        if mod == 'grey':
            continue
        mod_genes = gene_module[gene_module['Module'] == mod]['Gene'].tolist()
        mod_genes_in_expr = [g for g in mod_genes if g in expr.columns]
        if len(mod_genes_in_expr) < 5:
            continue
        mod_expr = expr[mod_genes_in_expr]
        # Center and scale, then PC1
        mod_expr_scaled = mod_expr.apply(zscore, nan_policy='omit').fillna(0)
        if mod_expr_scaled.shape[1] > 1:
            U, S, Vt = np.linalg.svd(mod_expr_scaled.values, full_matrices=False)
            pc1 = U[:, 0] * S[0]
            # Make sign consistent: correlate with mean expression
            mean_expr = mod_expr.mean(axis=1)
            if np.corrcoef(pc1, mean_expr)[0, 1] < 0:
                pc1 = -pc1
            mes[mod] = pd.Series(pc1, index=mod_expr.index)
    return pd.DataFrame(mes)

print("  Computing liver module eigengenes...")
liver_MEs = compute_module_eigengenes(liver_expr, liver_gm)
print(f"  Liver MEs: {liver_MEs.shape[1]} modules")

print("  Computing muscle module eigengenes...")
muscle_MEs = compute_module_eigengenes(muscle_expr, muscle_gm)
print(f"  Muscle MEs: {muscle_MEs.shape[1]} modules")

# Align samples
common_samples = liver_MEs.index.intersection(muscle_MEs.index)
print(f"  Common samples: {len(common_samples)}")

liver_MEs_aligned = liver_MEs.loc[common_samples]
muscle_MEs_aligned = muscle_MEs.loc[common_samples]

# Cross-tissue correlation matrix
cross_cor = pd.DataFrame(index=liver_MEs_aligned.columns, columns=muscle_MEs_aligned.columns)
cross_pval = pd.DataFrame(index=liver_MEs_aligned.columns, columns=muscle_MEs_aligned.columns)

for lm in liver_MEs_aligned.columns:
    for mm in muscle_MEs_aligned.columns:
        r, p = pearsonr(liver_MEs_aligned[lm], muscle_MEs_aligned[mm])
        cross_cor.loc[lm, mm] = r
        cross_pval.loc[lm, mm] = p

cross_cor = cross_cor.astype(float)
cross_pval = cross_pval.astype(float)

# Find significant cross-tissue pairs
print("\n  Significant cross-tissue module pairs (|r| > 0.5, p < 0.001):")
cross_pairs = []
for lm in cross_cor.index:
    for mm in cross_cor.columns:
        r = cross_cor.loc[lm, mm]
        p = cross_pval.loc[lm, mm]
        if abs(r) > 0.5 and p < 0.001:
            cross_pairs.append({'Liver_Module': lm, 'Muscle_Module': mm, 'r': r, 'p': p})
            l_r_pd = liver_mtc.loc[lm, 'PD'] if lm in liver_mtc.index else np.nan
            m_r_pd = muscle_mtc.loc[mm, 'PD'] if mm in muscle_mtc.index else np.nan
            print(f"    {lm:15s} <-> {mm:15s}  r={r:+.3f} p={p:.5f}  "
                  f"[Liver r_PD={l_r_pd:+.2f}, Muscle r_PD={m_r_pd:+.2f}]")

cross_df = pd.DataFrame(cross_pairs).sort_values('r', key=abs, ascending=False) if cross_pairs else pd.DataFrame()

# ============================================================
# 5. LIVER-MUSCLE AXIS BRIDGING CANDIDATES
# ============================================================
print("\n[5] Identifying liver-muscle axis bridging candidates...")

# Secreted factors / signaling molecules that could bridge liver → muscle
secreted_keywords = ['growth factor', 'cytokine', 'hormone', 'insulin', 'IGF',
                     'FGF', 'HGF', 'EGF', 'BMP', 'TGF', 'WNT', 'SHH', 'NOTCH',
                     'interleukin', 'chemokine', 'TNF', 'adipokine', 'myokine',
                     'fetuin', 'binding protein', 'receptor']

# Find genes in liver modules that encode secreted/signaling proteins
# (using gene names as proxy — comprehensive list would need annotation)
# Focus on genes with high kME in their modules AND favorable trait correlations

def find_bridging_candidates(liver_gm, muscle_gm, liver_mtc, muscle_mtc,
                             liver_hub, muscle_hub, top_n=20):
    """Identify candidate cross-tissue bridging genes."""
    candidates = []

    # Known liver-muscle axis genes (literature curated)
    known_bridgers = {
        'IGF1': 'Growth factor', 'IGF2': 'Growth factor',
        'IGFBP2': 'IGF binding protein', 'IGFBP3': 'IGF binding protein',
        'IGFBP5': 'IGF binding protein',
        'AHSG': 'Fetuin-A, anti-inflammatory',
        'FGF21': 'FGF, metabolic regulator', 'FGF19': 'FGF, bile acid',
        'ANGPTL4': 'Angiopoietin-like, lipid',
        'ANGPTL3': 'Angiopoietin-like, lipid',
        'GDF15': 'TGF-beta family, stress',
        'MSTN': 'Myostatin, negative muscle regulator',
        'BDNF': 'Neurotrophin, muscle metabolism',
        'FNDC5': 'Irisin precursor, exercise',
        'IL6': 'Cytokine, muscle-liver crosstalk',
        'IL15': 'Cytokine, muscle anabolic',
        'TNF': 'Cytokine', 'TNFSF10': 'TNF superfamily',
        'LIF': 'Leukemia inhibitory factor',
        'CTGF': 'Connective tissue growth factor',
        'VEGFA': 'Vascular endothelial growth factor',
        'HGF': 'Hepatocyte growth factor',
        'EGF': 'Epidermal growth factor',
        'NRG1': 'Neuregulin, muscle development',
        'BMP2': 'Bone morphogenetic protein',
        'INHBA': 'Inhibin beta A / activin',
        'FST': 'Follistatin, MSTN antagonist',
        'THBS1': 'Thrombospondin, TGF-beta activator',
        'SERPINE1': 'PAI-1, fibrosis',
        'APOA1': 'Apolipoprotein, lipid transport',
        'APOE': 'Apolipoprotein, lipid transport',
        'TTR': 'Transthyretin, thyroid/retinol transport',
        'ALB': 'Albumin, carrier protein',
        'TF': 'Transferrin, iron transport',
        'HP': 'Haptoglobin, hemoglobin binding',
        'SERPINA1': 'Alpha-1 antitrypsin',
    }

    for gene, func in known_bridgers.items():
        l_rows = liver_gm[liver_gm['Gene'] == gene]
        m_rows = muscle_gm[muscle_gm['Gene'] == gene]

        if len(l_rows) > 0:
            l_mod = l_rows.iloc[0]['Module']
            l_kme = l_rows.iloc[0].get('kME_module', np.nan)
            l_gs_pd = l_rows.iloc[0].get('GS_PD', np.nan)
            l_gs_urea = l_rows.iloc[0].get('GS_Urea', np.nan)
            l_mod_r_pd = liver_mtc.loc[l_mod, 'PD'] if l_mod in liver_mtc.index else np.nan
            l_mod_r_urea = liver_mtc.loc[l_mod, 'Urea'] if l_mod in liver_mtc.index else np.nan
        else:
            l_mod = l_kme = l_gs_pd = l_gs_urea = l_mod_r_pd = l_mod_r_urea = np.nan

        m_mod = m_kme = m_gs_pd = m_gs_urea = m_mod_r_pd = m_mod_r_urea = np.nan
        if len(m_rows) > 0:
            m_mod = m_rows.iloc[0]['Module']
            m_kme = m_rows.iloc[0].get('kME_module', np.nan)
            m_gs_pd = m_rows.iloc[0].get('GS_PD', np.nan)
            m_gs_urea = m_rows.iloc[0].get('GS_Urea', np.nan)
            m_mod_r_pd = muscle_mtc.loc[m_mod, 'PD'] if m_mod in muscle_mtc.index else np.nan
            m_mod_r_urea = muscle_mtc.loc[m_mod, 'Urea'] if m_mod in muscle_mtc.index else np.nan

        # Score: higher if present in BOTH tissues
        in_both = (len(l_rows) > 0) and (len(m_rows) > 0)
        bridge_score = 0
        if len(l_rows) > 0:
            bridge_score += 1 + abs(l_kme) if not np.isnan(l_kme) else 1
        if len(m_rows) > 0:
            bridge_score += 1 + abs(m_kme) if not np.isnan(m_kme) else 1

        candidates.append({
            'Gene': gene, 'Function': func,
            'In_Liver': len(l_rows) > 0, 'In_Muscle': len(m_rows) > 0,
            'In_Both': in_both,
            'Liver_Module': l_mod, 'Liver_kME': l_kme,
            'Liver_GS_PD': l_gs_pd, 'Liver_GS_Urea': l_gs_urea,
            'Liver_Mod_r_PD': l_mod_r_pd, 'Liver_Mod_r_Urea': l_mod_r_urea,
            'Muscle_Module': m_mod, 'Muscle_kME': m_kme,
            'Muscle_GS_PD': m_gs_pd, 'Muscle_GS_Urea': m_gs_urea,
            'Muscle_Mod_r_PD': m_mod_r_pd, 'Muscle_Mod_r_Urea': m_mod_r_urea,
            'Bridge_Score': round(bridge_score, 2),
        })

    return pd.DataFrame(candidates).sort_values('Bridge_Score', ascending=False)

bridge_df = find_bridging_candidates(liver_gm, muscle_gm, liver_mtc, muscle_mtc,
                                     liver_hub, muscle_hub)

print(f"  Found {len(bridge_df)} bridging candidates")
if len(bridge_df) > 0:
    print(f"  Present in BOTH tissues: {bridge_df['In_Both'].sum()}")
    both_df = bridge_df[bridge_df['In_Both']].head(15)
    for _, r in both_df.iterrows():
        print(f"    {r['Gene']:12s} Liver:{r['Liver_Module']:15s} Muscle:{str(r['Muscle_Module']):15s} "
              f"lGS_PD={r['Liver_GS_PD']:.3f} mGS_PD={r['Muscle_GS_PD']:.3f}")

# ============================================================
# 6. INTEGRATED FIGURES
# ============================================================
print("\n[6] Generating integrated figures...")

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

# --- Fig A: Cross-Tissue Module-Module Correlation Heatmap ---
figA, axA = plt.subplots(figsize=(max(8, len(muscle_MEs.columns) * 0.5),
                                   max(6, len(liver_MEs.columns) * 0.5)))

im = axA.imshow(cross_cor.values, aspect='auto', cmap='RdBu_r', vmin=-1, vmax=1)

for i in range(cross_cor.shape[0]):
    for j in range(cross_cor.shape[1]):
        r_val = cross_cor.values[i, j]
        p_val = cross_pval.values[i, j]
        sig_str = '***' if p_val < 0.001 else ('**' if p_val < 0.01 else ('*' if p_val < 0.05 else ''))
        if abs(r_val) > 0.4:
            text = f'{r_val:.2f}{sig_str}'
            axA.text(j, i, text, ha='center', va='center', fontsize=6,
                    color='white' if abs(r_val) > 0.6 else 'black',
                    fontweight='bold' if sig_str else 'normal')

axA.set_xticks(range(len(cross_cor.columns)))
axA.set_xticklabels(cross_cor.columns, rotation=45, ha='right', fontsize=7)
axA.set_yticks(range(len(cross_cor.index)))
axA.set_yticklabels(cross_cor.index, fontsize=7)
axA.set_xlabel('Muscle Modules', fontweight='bold')
axA.set_ylabel('Liver Modules', fontweight='bold')
axA.set_title('Cross-Tissue Module Eigengene Correlation\n(Liver ↔ Muscle Co-Expression)', fontweight='bold')
plt.colorbar(im, ax=axA, shrink=0.8, label='Pearson r')
figA.tight_layout()
figA.savefig('fig_cross_tissue_module_correlation.png')
print("  Saved fig_cross_tissue_module_correlation.png")

# --- Fig B: Module-Trait Correlation Dot Plot (positive modules highlighted) ---
figB, axB = plt.subplots(figsize=(12, 6))

plot_data = []
for tissue, mtc, mtp in [('Liver', liver_mtc, liver_mtp), ('Muscle', muscle_mtc, muscle_mtp)]:
    for mod in mtc.index:
        if mod == 'grey':
            continue
        r_pd = mtc.loc[mod, 'PD']
        p_pd = mtp.loc[mod, 'PD']
        r_urea = mtc.loc[mod, 'Urea']
        n_genes = (liver_gm if tissue == 'Liver' else muscle_gm)
        n_genes = (n_genes['Module'] == mod).sum()
        plot_data.append({
            'Tissue': tissue, 'Module': mod, 'r_PD': r_pd, 'r_Urea': r_urea,
            'p_PD': p_pd, 'n_genes': n_genes,
            'abs_r': max(abs(r_pd), abs(r_urea))
        })

plot_df = pd.DataFrame(plot_data)
plot_df['-log10p'] = -np.log10(plot_df['p_PD'].clip(lower=1e-10))
plot_df['color'] = np.where(plot_df['r_PD'] > 0.3, '#D73027',
                     np.where(plot_df['r_PD'] < -0.3, '#4575B4', '#91BFDB'))
plot_df['size'] = np.sqrt(plot_df['n_genes']) * 8

for tissue, marker in [('Liver', 's'), ('Muscle', 'o')]:
    sub = plot_df[plot_df['Tissue'] == tissue]
    axB.scatter(sub['r_PD'], sub['r_Urea'], s=sub['size'], c=sub['color'],
               marker=marker, edgecolors='black', linewidth=0.3, alpha=0.8,
               label=tissue, zorder=3)

# Annotate key modules
for _, row in plot_df[plot_df['abs_r'] > 0.5].iterrows():
    axB.annotate(row['Module'], (row['r_PD'], row['r_Urea']),
                fontsize=7, ha='center', va='bottom',
                xytext=(0, 5), textcoords='offset points',
                fontweight='bold')

axB.axhline(y=0, color='grey', linewidth=0.5, linestyle='--')
axB.axvline(x=0, color='grey', linewidth=0.5, linestyle='--')
axB.set_xlabel('r(Module ME, Protein Deposition)', fontweight='bold')
axB.set_ylabel('r(Module ME, Serum Urea)', fontweight='bold')
axB.set_title('Module-Trait Association Landscape\nRed = POS (r_PD>0.3), Blue = NEG (r_PD<-0.3)', fontweight='bold')
axB.legend(loc='upper left', fontsize=8)

# Quadrant annotations
axB.text(0.8, -1.2, 'HIGH PD\nLow Urea\n(理想正向)', ha='center', fontsize=7,
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))
axB.text(-0.8, 1.2, 'Low PD\nHigh Urea\n(蛋白分解)', ha='center', fontsize=7,
        bbox=dict(boxstyle='round', facecolor='lightcoral', alpha=0.3))

figB.tight_layout()
figB.savefig('fig_module_trait_landscape.png')
print("  Saved fig_module_trait_landscape.png")

# --- Fig C: Enrichment Summary (text only, saved to Excel) ---
if len(enrich_df) > 0:
    pos_enrich = enrich_df[enrich_df['Module'].isin(muscle_pos_mods)]
    if len(pos_enrich) > 0:
        print("\n  Top Enriched KEGG/GO Terms in Positive Muscle Modules:")
        for mod in muscle_pos_mods[:3]:
            mod_enr = pos_enrich[pos_enrich['Module'] == mod]
            if len(mod_enr) > 0:
                try:
                    mod_enr = mod_enr.copy()
                    mod_enr['CS'] = mod_enr['Combined Score'].apply(
                        lambda x: float(x[0]) if isinstance(x, list) else float(x) if x else 0)
                    top = mod_enr.nlargest(5, 'CS')
                    print(f"\n    {mod} (r_PD={muscle_mtc.loc[mod, 'PD']:.3f}):")
                    for _, r in top.iterrows():
                        print(f"      {r['Term'][:90]}  CS={r['CS']:.1f}")
                except Exception as e:
                    print(f"    {mod}: enrichment display error: {e}")

# --- Fig D: Bridging Gene Summary ---
if len(bridge_df) > 0 and bridge_df['In_Both'].sum() > 3:
    figD, axD = plt.subplots(figsize=(10, max(6, bridge_df['In_Both'].sum() * 0.35)))

    both_genes = bridge_df[bridge_df['In_Both']].head(25)
    both_genes = both_genes.sort_values('Bridge_Score')

    x = np.arange(len(both_genes))
    width = 0.35

    bars1 = axD.barh(x - width/2, both_genes['Liver_GS_PD'].fillna(0), width,
                     color='#D73027', alpha=0.8, label='Liver GS_PD')
    bars2 = axD.barh(x + width/2, both_genes['Muscle_GS_PD'].fillna(0), width,
                     color='#4575B4', alpha=0.8, label='Muscle GS_PD')

    axD.set_yticks(x)
    axD.set_yticklabels(both_genes['Gene'], fontsize=8)
    axD.axvline(x=0, color='black', linewidth=0.5)
    axD.set_xlabel('Gene Significance (r with Protein Deposition)', fontweight='bold')
    axD.set_title('Liver-Muscle Axis Bridging Candidates\n(Known secreted/signaling factors present in BOTH tissues)',
                 fontweight='bold')
    axD.legend(loc='lower right', fontsize=8)
    figD.tight_layout()
    figD.savefig('fig_bridging_candidates.png')
    print("  Saved fig_bridging_candidates.png")

# ============================================================
# 7. SAVE INTEGRATED RESULTS
# ============================================================
print("\n[7] Saving integrated results...")

with pd.ExcelWriter('wgcna_step2_integration.xlsx', engine='openpyxl') as writer:
    # Cross-tissue correlations
    cross_cor.to_excel(writer, sheet_name='CrossTissue_Module_Cor')
    cross_pval.to_excel(writer, sheet_name='CrossTissue_Module_Pval')

    # Cross-tissue pairs
    if len(cross_df) > 0:
        cross_df.to_excel(writer, sheet_name='CrossTissue_Pairs', index=False)

    # Enrichment
    if len(enrich_df) > 0:
        enrich_df.to_excel(writer, sheet_name='GO_KEGG_Enrichment', index=False)

    # Bridging candidates
    if len(bridge_df) > 0:
        bridge_df.to_excel(writer, sheet_name='Bridging_Candidates', index=False)

    # Module gene lists for positive modules
    for mod in muscle_pos_mods:
        mod_genes = muscle_gm[muscle_gm['Module'] == mod][['Gene', 'kME_module', 'GS_PD', 'GS_Urea']]
        mod_genes = mod_genes.sort_values('kME_module', ascending=False)
        mod_genes.to_excel(writer, sheet_name=f'Muscle_{mod}_genes', index=False)

print("Saved wgcna_step2_integration.xlsx")

# ============================================================
# SUMMARY
# ============================================================
print("\n" + "=" * 70)
print("WGCNA STEP 2 COMPLETE")
print("=" * 70)
print(f"""
Key Outputs:
  - wgcna_step2_integration.xlsx (cross-tissue correlations, enrichment, bridging)
  - fig_cross_tissue_module_correlation.png (liver ME ↔ muscle ME heatmap)
  - fig_module_trait_landscape.png (module r_PD vs r_Urea scatter)
  - fig_enrichment_dotplot.png (top GO/KEGG terms in positive modules)
  - fig_bridging_candidates.png (liver-muscle axis candidates)

Muscle Positive Modules (正向调控蛋白沉积):
  {' '.join(muscle_pos_mods)}

Significant Cross-Tissue Pairs: {len(cross_df)}
Bridging Candidates in Both Tissues: {bridge_df['In_Both'].sum()}
""")

# Print key bridging genes
if len(bridge_df) > 0:
    both = bridge_df[bridge_df['In_Both']]
    if len(both) > 0:
        print("Top Liver-Muscle Axis Bridging Candidates:")
        for _, r in both.head(10).iterrows():
            l_dir = 'POS' if r['Liver_GS_PD'] > 0 else 'NEG'
            m_dir = 'POS' if r['Muscle_GS_PD'] > 0 else 'NEG'
            print(f"  {r['Gene']:12s} | {r['Function']:35s} | "
                  f"Liver: {l_dir} (GS={r['Liver_GS_PD']:.3f}) | "
                  f"Muscle: {m_dir} (GS={r['Muscle_GS_PD']:.3f})")

print("\nDone!")
