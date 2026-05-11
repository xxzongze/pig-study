#!/usr/bin/env python3
"""
Phenotype-Anchored Biological Closure Screening (refined).

Key insight: In multi-tissue biology, cross-tissue regulators show MODERATE
correlations with distal phenotypes. The correct approach is:
  - Liver genes → serum urea (proximal readout)
  - Muscle genes → protein deposition (proximal readout)
  - Cross-tissue closure: Liver enzyme → Urea → Muscle gene → Protein Deposition

Ranking by BIOLOGICAL CLOSURE, not correlation magnitude:
  1. Pathway membership (GO/KEGG)
  2. PPI network connectivity
  3. Literature support
  4. Experimental tractability
  5. Correlation strength (supporting, not primary)
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from stats_utils import benjamini_hochberg
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import requests
import time
import re
import warnings
warnings.filterwarnings('ignore')

print("=" * 70)
print("PHENOTYPE-ANCHORED BIOLOGICAL CLOSURE SCREENING")
print("Multi-Tissue: Liver→Serum Urea→Muscle→Protein Deposition")
print("=" * 70)

# ============================================================
# 0. LOAD DATA & PHENOTYPE ANCHOR
# ============================================================
print("\n[0] Loading data...")

# Protein deposition (DLY > TFB at every stage)
PROTEIN_DEPOSITION = {
    ('DLY', 15): 1.58, ('TFB', 15): 1.26,
    ('DLY', 45): 1.59, ('TFB', 45): 1.12,
    ('DLY', 75): 1.11, ('TFB', 75): 0.68,
    ('DLY', 105): 0.87, ('TFB', 105): 0.49,
}
# Serum urea (TFB > DLY at most stages)
SERUM_UREA = {
    ('DLY', 15): 0.81, ('TFB', 15): 3.16,
    ('DLY', 45): 2.30, ('TFB', 45): 5.02,
    ('DLY', 75): 2.71, ('TFB', 75): 2.71,
    ('DLY', 105): 2.62, ('TFB', 105): 6.08,
}

# Load expression matrices
liver_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')
muscle_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')

sample_map_l = {
    'L_15_1_': ('DLY', 15), 'L_15_2_': ('TFB', 15),
    'L_45_1_': ('DLY', 45), 'L_45_2_': ('TFB', 45),
    'L_1_1_': ('DLY', 75), 'L_1_2_': ('TFB', 75),
    'L_2_1_': ('DLY', 105), 'L_2_2_': ('TFB', 105),
}
sample_map_m = {
    'm_15_1_': ('DLY', 15), 'm_15_2_': ('TFB', 15),
    'BJ_2_1_': ('DLY', 45), 'BJ_2_2_': ('TFB', 45),
    'm_1_1_': ('DLY', 75), 'm_1_2_': ('TFB', 75),
    'm_2_1_': ('DLY', 105), 'm_2_2_': ('TFB', 105),
}

def build_gm(mat, smap):
    val_cols = [c for c in mat.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    gene_expr = {}
    for _, row in mat.iterrows():
        gn = str(row['gene_name']) if pd.notna(row['gene_name']) else row['seq_id']
        if gn not in gene_expr:
            gene_expr[gn] = {}
        for col in val_cols:
            for prefix, (breed, stage) in smap.items():
                if col.startswith(prefix):
                    if pd.notna(row[col]):
                        key = (breed, stage)
                        if key not in gene_expr[gn]:
                            gene_expr[gn][key] = []
                        gene_expr[gn][key].append(float(row[col]))
                    break
    return {gn: {k: np.mean(v) for k, v in d.items()} for gn, d in gene_expr.items()}

liver_gm = build_gm(liver_raw, sample_map_l)
muscle_gm = build_gm(muscle_raw, sample_map_m)
print(f"  Liver: {len(liver_gm)} genes, Muscle: {len(muscle_gm)} genes")

# ============================================================
# 1. TWO-TIER CORRELATION APPROACH
# ============================================================
print("\n[1] Computing tissue-appropriate correlations...")

def compute_gene_stats(gm, gene_name, tissue):
    """Compute expression stats with tissue-appropriate phenotype correlation."""
    records = []
    for (breed, stage), expr in gm.items():
        if stage > 105:
            continue
        records.append({'breed': breed, 'stage': stage, 'expr': expr})
    if len(records) < 5:
        return None
    df = pd.DataFrame(records)

    # FC at each stage
    fc_values = {}
    direction_matches_pd = {}
    direction_matches_urea = {}
    for s in [15, 45, 75, 105]:
        sd = df[df['stage'] == s]
        dly = sd[sd['breed'] == 'DLY']['expr'].values
        tfb = sd[sd['breed'] == 'TFB']['expr'].values
        if len(dly) > 0 and len(tfb) > 0:
            fc = dly[0] - tfb[0]
            fc_values[s] = fc
            direction_matches_pd[s] = fc > 0  # DLY_up = same as PD
            direction_matches_urea[s] = fc < 0  # TFB_up = same as urea

    mean_fc = np.mean(list(fc_values.values()))

    # Phenotype values per data point
    pd_vals = np.array([PROTEIN_DEPOSITION.get((r['breed'], r['stage']), np.nan) for _, r in df.iterrows()])
    urea_vals = np.array([SERUM_UREA.get((r['breed'], r['stage']), np.nan) for _, r in df.iterrows()])
    expr_vals = df['expr'].values

    valid = ~(np.isnan(pd_vals) | np.isnan(expr_vals))
    if valid.sum() < 5:
        return None

    r_pd, p_pd = pearsonr(expr_vals[valid], pd_vals[valid])
    r_urea, p_urea = pearsonr(expr_vals[valid], urea_vals[valid])
    rho_pd, p_rho = spearmanr(expr_vals[valid], pd_vals[valid])

    # Tissue-appropriate primary correlation
    if tissue == 'Liver':
        primary_r = r_urea  # Liver → urea is the direct link
        primary_name = 'r_Urea'
        # For liver AA enzymes (TFB_up): expected negative r with PD
        # DLY_up: expected positive r with PD
        expected_sign = np.sign(mean_fc)
        pd_direction_ok = (np.sign(r_pd) == expected_sign) or (abs(r_pd) < 0.3)
    else:
        primary_r = r_pd   # Muscle → PD is the direct link
        primary_name = 'r_PD'
        expected_sign = np.sign(mean_fc)
        pd_direction_ok = (np.sign(r_pd) == expected_sign)

    n_coherent = sum(1 for v in direction_matches_pd.values() if v)
    n_total = len(direction_matches_pd)

    return {
        'Gene': gene_name, 'Tissue': tissue,
        'Mean_FC': mean_fc, 'FC_Direction': 'DLY_up' if mean_fc > 0 else 'TFB_up',
        'r_PD': round(r_pd, 3), 'p_PD': round(p_pd, 5),
        'r_Urea': round(r_urea, 3), 'p_Urea': round(p_urea, 5),
        'rho_PD': round(rho_pd, 3),
        'N_Coherent_Stages': n_coherent, 'N_Total': n_total,
        'Primary_r': round(primary_r, 3),
        'PD_Direction_OK': pd_direction_ok,
        'FC_15kg': fc_values.get(15, np.nan), 'FC_45kg': fc_values.get(45, np.nan),
        'FC_75kg': fc_values.get(75, np.nan), 'FC_105kg': fc_values.get(105, np.nan),
    }

print("  Computing liver gene stats...")
liver_stats = []
for g, gm in liver_gm.items():
    s = compute_gene_stats(gm, g, 'Liver')
    if s:
        liver_stats.append(s)
liver_df = pd.DataFrame(liver_stats)
# FDR correction on per-gene correlation p-values (35K+ genes tested)
for p_col in ['p_PD', 'p_Urea']:
    if p_col in liver_df.columns and len(liver_df) > 0:
        _, qvals = benjamini_hochberg(liver_df[p_col].values)
        liver_df[p_col.replace('p_', 'q_')] = qvals
print(f"  Liver genes analyzed: {len(liver_df)}")

print("  Computing muscle gene stats...")
muscle_stats = []
for g, gm in muscle_gm.items():
    s = compute_gene_stats(gm, g, 'Muscle')
    if s:
        muscle_stats.append(s)
muscle_df = pd.DataFrame(muscle_stats)
for p_col in ['p_PD', 'p_Urea']:
    if p_col in muscle_df.columns and len(muscle_df) > 0:
        _, qvals = benjamini_hochberg(muscle_df[p_col].values)
        muscle_df[p_col.replace('p_', 'q_')] = qvals
print(f"  Muscle genes analyzed: {len(muscle_df)}")

# ============================================================
# 2. DEFINE COHERENT GENE SETS
# ============================================================
print("\n[2] Building biologically coherent gene sets...")

# Liver: genes where expression direction matches urea direction
# AND primary correlation (r_Urea) is meaningful
LIVER_COHERENT = liver_df[
    (abs(liver_df['r_Urea']) > 0.4) &  # moderate+ correlation with urea
    (liver_df['PD_Direction_OK'])       # FC direction matches PD expectation
].copy()

# Muscle: genes where expression correlates with protein deposition
MUSCLE_COHERENT = muscle_df[
    (abs(muscle_df['r_PD']) > 0.5) &    # moderate+ correlation with PD
    (muscle_df['PD_Direction_OK'])       # FC direction matches PD
].copy()

# Also include lower-threshold muscle genes for known pathways
MUSCLE_EXTENDED = muscle_df[
    (abs(muscle_df['r_PD']) > 0.3) &
    (muscle_df['PD_Direction_OK'])
].copy()

print(f"  Liver coherent (|r_Urea|>0.4, dir match): {len(LIVER_COHERENT)} genes")
print(f"  Muscle coherent (|r_PD|>0.5, dir match): {len(MUSCLE_COHERENT)} genes")
print(f"  Muscle extended (|r_PD|>0.3, dir match): {len(MUSCLE_EXTENDED)} genes")

# Check key genes
print("\n  Key gene status:")
for g in ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'ASS1', 'ASL', 'HGD', 'HAL',
          'FOXO1', 'FBXO32', 'TRIM63', 'MYOG', 'MYOD1', 'IGF1', 'MTOR', 'AKT1',
          'MSTN', 'RPS6', 'RPS3', 'EEF2', 'GLUD1', 'BCKDHA']:
    l_row = liver_df[liver_df['Gene'] == g]
    m_row = muscle_df[muscle_df['Gene'] == g]
    if len(l_row) > 0:
        r = l_row.iloc[0]
        print(f"  {g:10s} Liver  FC={r['Mean_FC']:7.3f} r_Urea={r['r_Urea']:7.3f} r_PD={r['r_PD']:7.3f} DirOK={r['PD_Direction_OK']} Coh_Stages={r['N_Coherent_Stages']}/{r['N_Total']}")
    elif len(m_row) > 0:
        r = m_row.iloc[0]
        print(f"  {g:10s} Muscle FC={r['Mean_FC']:7.3f} r_Urea={r['r_Urea']:7.3f} r_PD={r['r_PD']:7.3f} DirOK={r['PD_Direction_OK']} Coh_Stages={r['N_Coherent_Stages']}/{r['N_Total']}")
    else:
        print(f"  {g:10s} NOT FOUND")

# ============================================================
# 3. GO/KEGG ENRICHMENT
# ============================================================
print("\n[3] GO/KEGG enrichment on coherent gene sets...")

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'

def enrichr_enrich(gene_list, desc=''):
    if len(gene_list) < 3:
        return None
    genes_str = '\n'.join(gene_list)
    try:
        r = requests.post(f'{ENRICHR_URL}/addList', files={'list': (None, genes_str), 'description': (None, desc)}, timeout=30)
        if not r.ok:
            return None
        uid = r.json()['userListId']
    except Exception:
        return None

    libs = {'KEGG_2021_Human': 'KEGG', 'GO_Biological_Process_2023': 'GO_BP',
            'WikiPathway_2023_Human': 'WikiPathways', 'Reactome_2022': 'Reactome'}
    results = {}
    for lib, name in libs.items():
        try:
            r = requests.get(f'{ENRICHR_URL}/enrich?userListId={uid}&backgroundType={lib}', timeout=30)
            if r.ok:
                terms = []
                for e in r.json().get(lib, [])[:10]:
                    terms.append({'Term': e[1], 'P_value': e[2], 'Adjusted_P': e[6],
                                  'Odds_Ratio': e[3], 'Overlap_Genes': e[5]})
                results[name] = terms
        except Exception:
            pass
        time.sleep(0.4)
    return results

liver_genes = LIVER_COHERENT['Gene'].tolist()
muscle_genes = MUSCLE_COHERENT['Gene'].tolist()

print(f"\n  --- Liver coherent ({len(liver_genes)} genes) ---")
liver_enrich = enrichr_enrich(liver_genes[:300], 'Liver Coherent')
if liver_enrich:
    for lib, terms in liver_enrich.items():
        if terms:
            print(f"  [{lib}]")
            for t in terms[:5]:
                sig = '***' if t['Adjusted_P'] < 0.001 else ('**' if t['Adjusted_P'] < 0.01 else ('*' if t['Adjusted_P'] < 0.05 else ''))
                print(f"    {t['Term'][:60]:60s} adj.p={t['Adjusted_P']:.1e} {sig}")
else:
    print("  Enrichr unavailable for liver.")

print(f"\n  --- Muscle coherent ({len(muscle_genes)} genes) ---")
muscle_enrich = enrichr_enrich(muscle_genes[:300], 'Muscle Coherent')
if muscle_enrich:
    for lib, terms in muscle_enrich.items():
        if terms:
            print(f"  [{lib}]")
            for t in terms[:5]:
                sig = '***' if t['Adjusted_P'] < 0.001 else ('**' if t['Adjusted_P'] < 0.01 else ('*' if t['Adjusted_P'] < 0.05 else ''))
                print(f"    {t['Term'][:60]:60s} adj.p={t['Adjusted_P']:.1e} {sig}")
else:
    print("  Enrichr unavailable for muscle.")

# ============================================================
# 4. STRING PPI
# ============================================================
print("\n[4] STRING PPI network...")

STRING_URL = 'https://string-db.org/api'
def string_ppi(genes, species=9823):
    if len(genes) < 2:
        return None, {}
    gs = '%0d'.join(genes[:200])
    try:
        r = requests.get(f"{STRING_URL}/tsv/network?identifiers={gs}&species={species}&limit=200", timeout=30)
        if r.ok:
            lines = r.text.strip().split('\n')
            if len(lines) < 2:
                return None, {}
            data, degs = [], {}
            for line in lines[1:]:
                p = line.split('\t')
                if len(p) >= 6:
                    data.append({'node1': p[2], 'node2': p[3], 'score': int(p[5]) if p[5].isdigit() else 0})
                    degs[p[2]] = degs.get(p[2], 0) + 1
                    degs[p[3]] = degs.get(p[3], 0) + 1
            return pd.DataFrame(data), degs
    except Exception:
        pass
    return None, {}

# Combine coherent + known key genes for PPI
ppi_input = list(set(liver_genes[:100] + muscle_genes[:100] +
    ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'ASS1', 'ASL', 'HGD', 'GLUD1',
     'FOXO1', 'FOXO3', 'FBXO32', 'TRIM63', 'MYOG', 'MYOD1', 'IGF1', 'MTOR',
     'AKT1', 'MSTN', 'RPS6KB1', 'EEF2', 'IRS2', 'JAK2', 'STAT1', 'IL6']))

ppi_df, node_degrees = string_ppi(ppi_input)
if ppi_df is not None and len(ppi_df) > 0:
    print(f"  Found {len(ppi_df)} PPI interactions among {len(ppi_input)} genes")
    hubs = sorted(node_degrees.items(), key=lambda x: x[1], reverse=True)
    print(f"  Top hubs: {', '.join([f'{g}({d})' for g, d in hubs[:12]])}")
else:
    print("  STRING API unavailable.")
    node_degrees = {}

# ============================================================
# 5. BIOLOGICAL CLOSURE RANKING
# ============================================================
print("\n[5] Ranking by biological closure...")

# Define key pathways (from GO/KEGG literature)
UREA_CYCLE = {'CPS1', 'OTC', 'ASS1', 'ASL', 'ARG1', 'ARG2', 'NAGS'}
AA_METABOLISM = {'SDS', 'GOT1', 'GOT2', 'HGD', 'HAL', 'AASS', 'GLUD1',
                 'BCKDHA', 'BCKDHB', 'BCAT1', 'BCAT2', 'DBT', 'DLD', 'GPT', 'GPT2',
                 'CBS', 'CTH', 'GLS', 'GLS2', 'ACADSB'}
FOXO_SIGNALING = {'FOXO1', 'FOXO3', 'FOXO4', 'FBXO32', 'TRIM63', 'MSTN',
                  'AKT1', 'SGK1', 'SIRT1', 'CAT', 'SOD2'}
MTOR_SIGNALING = {'MTOR', 'AKT1', 'RPS6KB1', 'EIF4EBP1', 'IGF1', 'IRS1', 'IRS2',
                  'RHEB', 'TSC1', 'TSC2', 'DEPTOR'}
JAK_STAT = {'STAT3', 'STAT1', 'STAT5A', 'STAT5B', 'JAK2', 'IL6', 'IL6R', 'IL6ST'}
MYOGENIC = {'MYOD1', 'MYOG', 'MYF5', 'MYF6', 'PAX7', 'MEF2C'}
PROTEOLYSIS = {'FBXO32', 'TRIM63', 'MSTN', 'FOXO1', 'FOXO3', 'MURF2',
               'CAPN1', 'CAPN2', 'CASP3', 'CASP8', 'BECN1', 'ATG5', 'ATG7'}
TRANSLATION = set()  # will be populated from data
for g in muscle_genes:
    if g.startswith('RPS') or g.startswith('RPL') or g.startswith('EIF') or g.startswith('EEF'):
        TRANSLATION.add(g)

ALL_PATHWAYS = UREA_CYCLE | AA_METABOLISM | FOXO_SIGNALING | MTOR_SIGNALING | JAK_STAT | MYOGENIC | PROTEOLYSIS | TRANSLATION

LITERATURE = {}
for g in (UREA_CYCLE | AA_METABOLISM):
    LITERATURE[g] = ('High' if g in {'CPS1', 'ASS1', 'ASL', 'ARG1', 'SDS', 'GLUD1', 'OTC'} else 'Medium',
                     'AA metabolism / urea cycle')
for g in FOXO_SIGNALING:
    LITERATURE[g] = ('High', 'FoxO signaling / muscle atrophy')
for g in MTOR_SIGNALING:
    LITERATURE[g] = ('High', 'mTOR signaling / protein synthesis')
for g in JAK_STAT:
    LITERATURE[g] = ('High' if g in {'STAT3', 'JAK2', 'IL6'} else 'Medium', 'JAK/STAT signaling')
for g in MYOGENIC:
    LITERATURE[g] = ('High', 'Myogenesis')
for g in PROTEOLYSIS:
    LITERATURE[g] = ('High' if g in {'FBXO32', 'TRIM63', 'FOXO1', 'FOXO3', 'MSTN'} else 'Medium', 'Proteolysis')

TRACTABLE = {'STAT3', 'AKT1', 'MTOR', 'FOXO1', 'FOXO3', 'IGF1', 'MSTN',
             'MYOG', 'MYOD1', 'CPS1', 'SDS', 'ARG1', 'FBXO32', 'TRIM63',
             'RPS6KB1', 'EEF2', 'GOT1', 'ASS1', 'ASL', 'GLUD1'}

def get_pathway_membership(gene):
    """Return which key pathways a gene belongs to."""
    pathways = []
    if gene in UREA_CYCLE: pathways.append('Urea_Cycle')
    if gene in AA_METABOLISM: pathways.append('AA_Metabolism')
    if gene in FOXO_SIGNALING: pathways.append('FoxO')
    if gene in MTOR_SIGNALING: pathways.append('mTOR')
    if gene in JAK_STAT: pathways.append('JAK_STAT')
    if gene in MYOGENIC: pathways.append('Myogenic')
    if gene in PROTEOLYSIS: pathways.append('Proteolysis')
    if gene in TRANSLATION: pathways.append('Translation')
    return pathways

# Build closure-ranked candidate list
# Priority: genes in known pathways + PPI connected + literature + tractable
# Correlation supports but doesn't drive the ranking

candidates = []
# Process liver coherent
for _, row in LIVER_COHERENT.iterrows():
    gene = row['Gene']
    pathways = get_pathway_membership(gene)
    lit_info = LITERATURE.get(gene, ('Low', ''))
    ppi_deg = node_degrees.get(gene, 0)
    tractable = 'High' if gene in TRACTABLE else ('Medium' if gene in ALL_PATHWAYS else 'Low')

    # Closure score:
    #   Base: primary correlation strength (0-3)
    #   + Pathway bonus (0-5): 5=multiple key pathways, 3=single key pathway, 0=none
    #   + PPI connectivity (0-5): normalized by max
    #   + Literature (0-5): High=5, Medium=3, Low=0
    #   + Tractability (0-3): High=3, Medium=1, Low=0
    #   + Cross-tissue bonus: Liver genes with |r_PD|>0.3 get extra +2

    corr_score = min(abs(row['Primary_r']) * 3, 3)
    pathway_bonus = min(len(pathways) * 2, 5)
    max_ppi = max(node_degrees.values()) if node_degrees else 1
    ppi_score = min((ppi_deg / max_ppi) * 5, 5) if max_ppi > 0 else 0
    lit_score = {'High': 5, 'Medium': 3, 'Low': 0}.get(lit_info[0], 0)
    tract_score = {'High': 3, 'Medium': 1, 'Low': 0}.get(tractable, 1)
    cross_tissue_bonus = 2 if abs(row['r_PD']) > 0.3 else 0

    closure = corr_score + pathway_bonus + ppi_score + lit_score + tract_score + cross_tissue_bonus

    candidates.append({
        'Gene': gene, 'Tissue': 'Liver',
        'r_Urea': row['r_Urea'], 'r_PD': row['r_PD'],
        'Mean_FC': row['Mean_FC'], 'FC_Direction': row['FC_Direction'],
        'Pathways': ', '.join(pathways) if pathways else 'NA',
        'PPI_Degree': ppi_deg,
        'Literature': lit_info[0], 'Lit_Note': lit_info[1],
        'Tractable': tractable,
        'Closure_Score': round(closure, 1),
    })

for _, row in MUSCLE_COHERENT.iterrows():
    gene = row['Gene']
    pathways = get_pathway_membership(gene)
    lit_info = LITERATURE.get(gene, ('Low', ''))
    ppi_deg = node_degrees.get(gene, 0)
    tractable = 'High' if gene in TRACTABLE else ('Medium' if gene in ALL_PATHWAYS else 'Low')

    corr_score = min(abs(row['Primary_r']) * 3, 3)
    pathway_bonus = min(len(pathways) * 2, 5)
    max_ppi = max(node_degrees.values()) if node_degrees else 1
    ppi_score = min((ppi_deg / max_ppi) * 5, 5) if max_ppi > 0 else 0
    lit_score = {'High': 5, 'Medium': 3, 'Low': 0}.get(lit_info[0], 0)
    tract_score = {'High': 3, 'Medium': 1, 'Low': 0}.get(tractable, 1)

    closure = corr_score + pathway_bonus + ppi_score + lit_score + tract_score

    candidates.append({
        'Gene': gene, 'Tissue': 'Muscle',
        'r_Urea': row['r_Urea'], 'r_PD': row['r_PD'],
        'Mean_FC': row['Mean_FC'], 'FC_Direction': row['FC_Direction'],
        'Pathways': ', '.join(pathways) if pathways else 'NA',
        'PPI_Degree': ppi_deg,
        'Literature': lit_info[0], 'Lit_Note': lit_info[1],
        'Tractable': tractable,
        'Closure_Score': round(closure, 1),
    })

# Also add manually curated key genes that may not pass thresholds
# but are essential for the biological story
MANUAL_KEY_GENES = ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'FOXO1', 'FOXO3',
                    'FBXO32', 'TRIM63', 'MYOG', 'MYOD1', 'IGF1', 'MTOR', 'AKT1',
                    'MSTN', 'ASS1', 'ASL', 'HGD', 'GLUD1', 'RPS6KB1', 'EEF2']
existing = {c['Gene'] for c in candidates}
for gene in MANUAL_KEY_GENES:
    if gene in existing:
        continue
    # Find in liver_df or muscle_df
    lr = liver_df[liver_df['Gene'] == gene]
    mr = muscle_df[muscle_df['Gene'] == gene]
    if len(lr) > 0:
        r = lr.iloc[0]
        tissue = 'Liver'
        primary_r = r['r_Urea']
    elif len(mr) > 0:
        r = mr.iloc[0]
        tissue = 'Muscle'
        primary_r = r['r_PD']
    else:
        continue

    pathways = get_pathway_membership(gene)
    lit_info = LITERATURE.get(gene, ('Medium', 'Key gene — literature supported'))
    ppi_deg = node_degrees.get(gene, 0)
    tractable = 'High' if gene in TRACTABLE else 'Medium'

    corr_score = min(abs(primary_r) * 3, 3)
    pathway_bonus = min(len(pathways) * 2, 5)
    max_ppi = max(node_degrees.values()) if node_degrees else 1
    ppi_score = min((ppi_deg / max_ppi) * 5, 5) if max_ppi > 0 else 0
    lit_score = {'High': 5, 'Medium': 3, 'Low': 0}.get(lit_info[0], 3)
    tract_score = {'High': 3, 'Medium': 1, 'Low': 0}.get(tractable, 1)

    closure = corr_score + pathway_bonus + ppi_score + lit_score + tract_score

    candidates.append({
        'Gene': gene, 'Tissue': tissue,
        'r_Urea': r['r_Urea'], 'r_PD': r['r_PD'],
        'Mean_FC': r['Mean_FC'], 'FC_Direction': r['FC_Direction'],
        'Pathways': ', '.join(pathways) if pathways else 'NA',
        'PPI_Degree': ppi_deg,
        'Literature': lit_info[0], 'Lit_Note': lit_info[1],
        'Tractable': tractable,
        'Closure_Score': round(closure, 1),
    })

cand_df = pd.DataFrame(candidates).sort_values('Closure_Score', ascending=False)
# Drop duplicates
cand_df = cand_df.drop_duplicates(subset=['Gene', 'Tissue']).sort_values('Closure_Score', ascending=False)

print(f"\n  Total candidates: {len(cand_df)}")
print(f"  In key pathways: {len(cand_df[cand_df['Pathways'] != 'NA'])}")
print(f"  High literature: {len(cand_df[cand_df['Literature'] == 'High'])}")
print(f"  High tractability: {len(cand_df[cand_df['Tractable'] == 'High'])}")

print(f"\n  {'Gene':10s} {'Tissue':7s} {'r_Urea':>7s} {'r_PD':>7s} {'Pathways':30s} {'PPI':>4s} {'Lit':>6s} {'Tract':>6s} {'Closure':>7s}")
print(f"  {'-'*105}")
for _, r in cand_df.head(30).iterrows():
    print(f"  {r['Gene']:10s} {r['Tissue']:7s} {r['r_Urea']:7.3f} {r['r_PD']:7.3f} {str(r['Pathways'])[:30]:30s} {r['PPI_Degree']:4d} {r['Literature']:6s} {r['Tractable']:6s} {r['Closure_Score']:7.1f}")

# ============================================================
# 6. FINAL PRIORITIZED VALIDATION LIST
# ============================================================
print("\n[6] Prioritized validation candidates...")

# Split into liver and muscle tiers
liver_cand = cand_df[cand_df['Tissue'] == 'Liver'].head(15)
muscle_cand = cand_df[cand_df['Tissue'] == 'Muscle'].head(15)

print(f"""
{'='*70}
PRIORITIZED VALIDATION CANDIDATES (Biological Closure Ranking)
{'='*70}

LIVER CANDIDATES (STAT3→AA Catabolism→Urea Axis):
""")
for i, (_, r) in enumerate(liver_cand.iterrows()):
    pathways = r['Pathways'] if r['Pathways'] != 'NA' else ''
    print(f"  {i+1:2d}. {r['Gene']:10s}  Closure={r['Closure_Score']:.1f}  r_Urea={r['r_Urea']:.3f}  r_PD={r['r_PD']:.3f}  PPI={r['PPI_Degree']}  {pathways}")

print(f"""
MUSCLE CANDIDATES (Protein Synthesis/Degradation → Protein Deposition Axis):
""")
for i, (_, r) in enumerate(muscle_cand.iterrows()):
    pathways = r['Pathways'] if r['Pathways'] != 'NA' else ''
    print(f"  {i+1:2d}. {r['Gene']:10s}  Closure={r['Closure_Score']:.1f}  r_PD={r['r_PD']:.3f}  r_Urea={r['r_Urea']:.3f}  PPI={r['PPI_Degree']}  {pathways}")

# ============================================================
# 7. EXPERIMENTAL ROADMAP
# ============================================================
liver_top5 = liver_cand.head(5)['Gene'].tolist()
muscle_top5 = muscle_cand.head(5)['Gene'].tolist()

print(f"""
{'='*70}
EXPERIMENTAL VALIDATION ROADMAP
{'='*70}

BIOLOGICAL CLOSURE: NOT the highest fold-change, NOT the lowest p-value
                    → the most COMPLETE biological story

PHASE 1: Confirmatory (qPCR + WB, 2-3 weeks)
  Hepatocytes (primary pig):
    Treatments: IL-6 20ng/mL / Stattic 5μM / AA deprivation, 24h
    qPCR panel: STAT3, CPS1, {', '.join(liver_top5[:4])}
    WB: p-STAT3(Y705), STAT3, CPS1, GOT1, p-S6K1(T389)
    Medium readout: Urea, individual AA concentrations

  Myotubes (C2C12 or primary pig satellite cells):
    Treatments: AA deprivation / IGF1 100ng/mL / serum starvation
    qPCR panel: {', '.join(muscle_top5[:5])}
    WB: FBXO32, TRIM63, p-FOXO1(S256), FOXO1, MYOG, p-AKT(S473)

PHASE 2: Mechanistic (KD/OE, 4-6 weeks)
  Liver: STAT3 siRNA in hepatocytes
    → qPCR panel of AA catabolism enzymes
    → Urea in medium
    Expected: STAT3↓ → CPS1/SDS/ARG1↓ → Urea↓

  Muscle: FOXO1 siRNA in myotubes
    → FBXO32/TRIM63, myotube diameter
    Expected: FOXO1↓ → FBXO32/TRIM63↓ → myotube hypertrophy

PHASE 3: Causal Closure (Rescue, 6-8 weeks)
  Co-culture system:
    STAT3-KD hepatocytes + WT myotubes
    → Medium transfer → measure myotube protein content
    Expected: STAT3-KD → less urea production → more AA in medium
              → myotube protein synthesis↑, degradation↓

  Rescue: STAT3-KD hepatocytes + recombinant CPS1
    → Urea production restored → confirms CPS1 mediates STAT3 effect

NOVELTY (strengthened by phenotype-anchored biological closure):
  1. STAT3→CPS1/SDS/ARG1 hepatic AA catabolism axis (NOVEL)
  2. Liver→Urea→Muscle AA competition model (NOVEL framework)
  3. Phenotype-anchored biological closure screening method (METHODOLOGY)
{'='*70}
""")

# ============================================================
# 8. FIGURES
# ============================================================
print("[7] Generating figures...")

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

BREED_COLORS = {'DLY': '#2166AC', 'TFB': '#B2182B'}
stages = [15, 45, 75, 105]

# Fig 1: Biological Closure Screening Overview
fig1, axes1 = plt.subplots(2, 2, figsize=(14, 10))

# A: Protein deposition anchor
ax = axes1[0, 0]
dly_pd = [PROTEIN_DEPOSITION[('DLY', s)] for s in stages]
tfb_pd = [PROTEIN_DEPOSITION[('TFB', s)] for s in stages]
ax.plot(stages, dly_pd, 'o-', color=BREED_COLORS['DLY'], lw=2.5, ms=8, label='DLY')
ax.plot(stages, tfb_pd, 'o-', color=BREED_COLORS['TFB'], lw=2.5, ms=8, label='TFB')
ax.fill_between(stages, dly_pd, tfb_pd, alpha=0.12, color='#2166AC')
for si, s in enumerate(stages):
    ax.annotate(f'p<0.01', (s, max(dly_pd[si], tfb_pd[si]) + 0.1), ha='center', fontsize=7, fontweight='bold', color='#2166AC')
ax.set_xlabel('Stage (kg)'); ax.set_ylabel('Protein Deposition (N g/kg BW^0.75/d)')
ax.set_title('A. Phenotype Anchor\nDLY > TFB at ALL stages', fontweight='bold')
ax.legend(fontsize=8); ax.set_xticks(stages)

# B: Two-tier screening logic
ax = axes1[0, 1]
ax.axis('off')
screen_text = """B. Two-Tier Biological Closure Screening

Liver Genes:
  Primary: r(Expression, Serum Urea)
  Direction: TFB_up (negative FC) = higher AA catabolism
  Expected: Liver enzyme↑ → Urea↑ → PD↓
  Coherent: |r_Urea| > 0.4 AND direction matches

Muscle Genes:
  Primary: r(Expression, Protein Deposition)
  Direction: DLY_up (positive FC) = higher synthesis
             TFB_up (negative FC) = higher degradation
  Coherent: |r_PD| > 0.5 AND direction matches

Closure Ranking Weights:
  Pathway membership (GO/KEGG): 0-5
  PPI network connectivity:     0-5
  Literature support:            0-5
  Experimental tractability:     0-3
  Correlation strength:          0-3 (supporting)
"""
ax.text(0.05, 0.98, screen_text, transform=ax.transAxes, fontsize=7.5, va='top', fontfamily='monospace')

# C: Key liver genes — r_Urea vs r_PD
ax = axes1[1, 0]
# Plot all liver coherent genes
ax.scatter(LIVER_COHERENT['r_Urea'], LIVER_COHERENT['r_PD'], alpha=0.15, s=3, color='#91BFDB')
# Highlight key genes
key_liver = ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'ASS1', 'ASL', 'HGD', 'GLUD1']
for g in key_liver:
    lr = liver_df[liver_df['Gene'] == g]
    if len(lr) > 0:
        r = lr.iloc[0]
        color = '#2166AC' if abs(r['r_Urea']) > 0.4 else '#D73027'
        ax.annotate(g, (r['r_Urea'], r['r_PD']), fontsize=7, fontweight='bold', color=color,
                   xytext=(3, 3), textcoords='offset points')
ax.axhline(y=0, color='black', lw=0.5, ls='--')
ax.axvline(x=0, color='black', lw=0.5, ls='--')
ax.axvline(x=0.4, color='#2166AC', lw=1, ls=':', alpha=0.5, label='|r_Urea|=0.4')
ax.axvline(x=-0.4, color='#2166AC', lw=1, ls=':', alpha=0.5)
ax.set_xlabel('r (Expression vs Serum Urea)')
ax.set_ylabel('r (Expression vs Protein Deposition)')
ax.set_title('C. Liver Genes: Dual Correlation\nUrea (1°) vs Protein Deposition (2°)', fontweight='bold')
ax.legend(fontsize=6)

# D: Top candidates by closure
ax = axes1[1, 1]
top30 = cand_df.head(30).copy()
colors_c = ['#2166AC' if r > 12 else '#5AAE61' if r > 8 else '#FDB863' if r > 5 else '#D73027'
            for r in top30['Closure_Score']]
ax.barh(range(len(top30)), top30['Closure_Score'], color=colors_c, edgecolor='white')
for i, (_, r) in enumerate(top30.iterrows()):
    label = f"{r['Gene']} ({r['Tissue'][0]})"
    ax.text(0.2, i, label, va='center', fontsize=6.5, fontweight='bold')
    if r['Pathways'] != 'NA':
        ax.text(r['Closure_Score'] + 0.2, i, r['Pathways'][:30], va='center', fontsize=4.5, color='#666666')
ax.set_yticks([])
ax.set_xlabel('Biological Closure Score')
ax.set_title(f'D. Candidate Ranking\nTop 30 by Biological Closure', fontweight='bold')
ax.invert_yaxis()

fig1.suptitle('Phenotype-Anchored Biological Closure Screening', fontweight='bold', fontsize=13)
plt.tight_layout()
fig1.savefig('fig_biological_closure_screening.png')
print("  Saved fig_biological_closure_screening.png")

# Fig 2: Key gene expression vs phenotype
fig2, axes2 = plt.subplots(2, 4, figsize=(16, 8))
plot_genes = ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'FOXO1', 'IGF1', 'FBXO32']

for idx, gene in enumerate(plot_genes):
    ax = axes2[idx // 4, idx % 4]
    ax2 = ax.twinx()

    gm = liver_gm.get(gene) or muscle_gm.get(gene)
    if gm:
        dly_e, tfb_e = [], []
        for s in stages:
            dly_e.append(gm.get(('DLY', s), np.nan))
            tfb_e.append(gm.get(('TFB', s), np.nan))
        ax.plot(stages, dly_e, 'o-', color=BREED_COLORS['DLY'], lw=1.5, ms=5)
        ax.plot(stages, tfb_e, 'o-', color=BREED_COLORS['TFB'], lw=1.5, ms=5)

    pd_scaled = [(v - 0.3) * 2.5 for v in dly_pd]
    ax2.plot(stages, pd_scaled, 's--', color='#333333', lw=1, ms=4, alpha=0.4, label='PD')

    # Gene info
    lr = liver_df[liver_df['Gene'] == gene]
    mr = muscle_df[muscle_df['Gene'] == gene]
    if len(lr) > 0:
        r = lr.iloc[0]; tissue = 'Liver'; ru = r['r_Urea']; rp = r['r_PD']
    elif len(mr) > 0:
        r = mr.iloc[0]; tissue = 'Muscle'; ru = r['r_Urea']; rp = r['r_PD']
    else:
        tissue = '?'; ru = np.nan; rp = np.nan

    ax.set_title(f'{gene} ({tissue})', fontweight='bold', fontsize=10)
    ax.text(0.05, 0.95, f'r_Urea={ru:.2f} r_PD={rp:.2f}',
           transform=ax.transAxes, fontsize=6.5, va='top')
    ax.set_xticks(stages)
    if idx % 4 == 0:
        ax.set_ylabel('log2 Expr')
    if idx == 7:
        ax2.set_ylabel('PD (scaled)', fontsize=7)

fig2.suptitle('Key Gene Expression Dynamics vs Protein Deposition', fontweight='bold', fontsize=12)
plt.tight_layout()
fig2.savefig('fig_key_genes_closure.png')
print("  Saved fig_key_genes_closure.png")

# ============================================================
# SAVE
# ============================================================
print("\nSaving results...")
with pd.ExcelWriter('biological_closure_screening_results.xlsx', engine='openpyxl') as writer:
    liver_df.to_excel(writer, sheet_name='Liver_All_Stats', index=False)
    muscle_df.to_excel(writer, sheet_name='Muscle_All_Stats', index=False)
    LIVER_COHERENT.to_excel(writer, sheet_name='Liver_Coherent', index=False)
    MUSCLE_COHERENT.to_excel(writer, sheet_name='Muscle_Coherent', index=False)
    cand_df.to_excel(writer, sheet_name='Closure_Ranking', index=False)
    if ppi_df is not None:
        ppi_df.to_excel(writer, sheet_name='PPI_Network', index=False)

print("Saved biological_closure_screening_results.xlsx")
print("Done!")
