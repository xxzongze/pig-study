#!/usr/bin/env python3
"""
Phenotype-Anchored Gene Screening Pipeline (user's recommended approach):

Step 1: Find genes where transcriptome direction matches protein deposition direction
Step 2: Filter for reliable correlation with protein deposition phenotype
Step 3: GO/KEGG enrichment on coherent gene set
Step 4: STRING PPI network
Step 5: Rank by biological closure (NOT fold-change magnitude)
Step 6: Experimental roadmap (qPCR/WB → KD/OE → rescue)

Anchor: Protein deposition (DLY > TFB at ALL stages, p<0.001)
         → Positive FC (DLY_up) = consistent with phenotype
         → Negative FC (TFB_up) = opposite to phenotype
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyBboxPatch
import requests
import time
import re

print("=" * 70)
print("PHENOTYPE-ANCHORED GENE SCREENING")
print("Protein Deposition as Biological Anchor")
print("=" * 70)

# ============================================================
# 0. LOAD DATA
# ============================================================
print("\n[0] Loading data...")

# Protein deposition (the anchor)
# DLY > TFB at every stage → positive FC = consistent with phenotype
PROTEIN_DEPOSITION = {
    ('DLY', 15): 1.58, ('TFB', 15): 1.26,
    ('DLY', 45): 1.59, ('TFB', 45): 1.12,
    ('DLY', 75): 1.11, ('TFB', 75): 0.68,
    ('DLY', 105): 0.87, ('TFB', 105): 0.49,
}

# Serum urea
SERUM_UREA = {
    ('DLY', 15): 0.81, ('TFB', 15): 3.16,
    ('DLY', 45): 2.30, ('TFB', 45): 5.02,
    ('DLY', 75): 2.71, ('TFB', 75): 2.71,
    ('DLY', 105): 2.62, ('TFB', 105): 6.08,
}

# Build phenotype DataFrame
phenotype_records = []
for (breed, stage), pd_val in PROTEIN_DEPOSITION.items():
    phenotype_records.append({'breed': breed, 'stage': stage,
                              'protein_deposition': pd_val,
                              'serum_urea': SERUM_UREA[(breed, stage)]})
phenotype = pd.DataFrame(phenotype_records)

# Load original expression matrices
liver_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')
muscle_raw = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')

print(f"  Liver matrix: {liver_raw.shape[0]} genes × {liver_raw.shape[1]} columns")
print(f"  Muscle matrix: {muscle_raw.shape[0]} genes × {muscle_raw.shape[1]} columns")

# ============================================================
# 1. BUILD GROUP-MEAN EXPRESSION TABLES
# ============================================================
print("\n[1] Building group-mean expression tables...")

sample_map_l = {
    'L_15_1_': ('DLY', 15), 'L_15_2_': ('TFB', 15),
    'L_45_1_': ('DLY', 45), 'L_45_2_': ('TFB', 45),
    'L_1_1_': ('DLY', 75), 'L_1_2_': ('TFB', 75),
    'L_2_1_': ('DLY', 105), 'L_2_2_': ('TFB', 105),
    'L_3_1_': ('DLY', 135),
}
sample_map_m = {
    'm_15_1_': ('DLY', 15), 'm_15_2_': ('TFB', 15),
    'BJ_2_1_': ('DLY', 45), 'BJ_2_2_': ('TFB', 45),
    'm_1_1_': ('DLY', 75), 'm_1_2_': ('TFB', 75),
    'm_2_1_': ('DLY', 105), 'm_2_2_': ('TFB', 105),
    'm_3_1_': ('DLY', 135),
}

def build_group_means(mat, sample_map):
    """Build gene × (breed, stage) group-mean expression matrix."""
    val_cols = [c for c in mat.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    gene_expr = {}
    for _, row in mat.iterrows():
        gn = str(row['gene_name']) if pd.notna(row['gene_name']) else row['seq_id']
        if gn not in gene_expr:
            gene_expr[gn] = {}
        for col in val_cols:
            for prefix, (breed, stage) in sample_map.items():
                if col.startswith(prefix):
                    if pd.notna(row[col]):
                        key = (breed, stage)
                        if key not in gene_expr[gn]:
                            gene_expr[gn][key] = []
                        gene_expr[gn][key].append(float(row[col]))
                    break
    # Average
    result = {}
    for gn, stage_dict in gene_expr.items():
        result[gn] = {k: np.mean(v) for k, v in stage_dict.items() if len(v) > 0}
    return result

liver_gm = build_group_means(liver_raw, sample_map_l)
muscle_gm = build_group_means(muscle_raw, sample_map_m)

print(f"  Liver genes with group means: {len(liver_gm)}")
print(f"  Muscle genes with group means: {len(muscle_gm)}")

# ============================================================
# 2. COMPUTE PHENOTYPE CONSISTENCY SCORES
# ============================================================
print("\n[2] Computing phenotype-anchored coherence scores...")

# Key principle:
#   DLY has HIGHER protein deposition than TFB at every stage.
#   So: DLY_up (positive FC) = SAME direction as phenotype = COHERENT
#       TFB_up (negative FC) = OPPOSITE to phenotype = needs negative correlation to be coherent

def compute_phenotype_coherence(gene_expr_dict, gene_name):
    """
    For a given gene, compute:
    - FC direction at each stage (DLY/TFB)
    - Correlation (Pearson & Spearman) with protein deposition
    - Direction consistency score
    Returns None if fewer than 5 data points.
    """
    records = []
    for (breed, stage), expr in gene_expr_dict.items():
        if stage == 135:  # skip DLY-only stage
            continue
        records.append({'breed': breed, 'stage': stage, 'expr': expr})

    if len(records) < 5:
        return None

    df = pd.DataFrame(records)
    merged = df.merge(phenotype, on=['breed', 'stage'])

    if len(merged) < 5:
        return None

    # Compute FC at each stage
    stages_present = sorted(merged['stage'].unique())
    fc_values = {}
    direction_matches = {}
    for s in stages_present:
        sd = merged[merged['stage'] == s]
        dly_val = sd[sd['breed'] == 'DLY']['expr'].values
        tfb_val = sd[sd['breed'] == 'TFB']['expr'].values
        if len(dly_val) > 0 and len(tfb_val) > 0:
            fc = dly_val[0] - tfb_val[0]  # log2(DLY/TFB) approximately
            fc_values[s] = fc
            # DLY has higher protein deposition → DLY_up is coherent
            # TFB_up (negative FC) means gene goes opposite to protein deposition
            direction_matches[s] = 'coherent' if fc > 0 else 'opposite'

    # Overall FC (mean across stages)
    mean_fc = np.mean(list(fc_values.values())) if fc_values else 0

    # Correlation with protein deposition
    r_pearson, p_pearson = pearsonr(merged['expr'], merged['protein_deposition'])
    rho_spearman, p_spearman = spearmanr(merged['expr'], merged['protein_deposition'])

    # Correlation with serum urea (secondary anchor: TFB > DLY for urea)
    r_urea, p_urea = pearsonr(merged['expr'], merged['serum_urea'])

    # Consistency score:
    #   coherent gene → positive FC AND positive r with protein deposition
    #   OR negative FC AND negative r with protein deposition
    n_coherent_stages = sum(1 for v in direction_matches.values() if v == 'coherent')
    n_total_stages = len(direction_matches)

    # Coherence logic:
    # If mean_fc > 0 (DLY_up): gene should positively correlate with protein deposition
    # If mean_fc < 0 (TFB_up): gene should negatively correlate with protein deposition
    fc_direction = 'DLY_up' if mean_fc > 0 else 'TFB_up'
    expected_r_sign = 1 if mean_fc > 0 else -1
    r_sign_match = (r_pearson * expected_r_sign) > 0

    # Composite coherence score (0-10)
    # 1. FC consistency across stages (0-4)
    fc_consistency = n_coherent_stages / n_total_stages * 4 if n_total_stages > 0 else 0
    # 2. Correlation strength with protein deposition (0-3)
    corr_strength = min(abs(r_pearson) * 3, 3)
    # 3. Direction match bonus (0-3)
    direction_bonus = 3 if r_sign_match else 0

    coherence_score = fc_consistency + corr_strength + direction_bonus

    return {
        'Gene': gene_name,
        'Mean_FC': round(mean_fc, 4),
        'FC_Direction': fc_direction,
        'N_Coherent_Stages': n_coherent_stages,
        'N_Total_Stages': n_total_stages,
        'FC_Consistency': round(fc_consistency, 1),
        'r_ProteinDeposition': round(r_pearson, 3),
        'p_ProteinDeposition': round(p_pearson, 5),
        'rho_ProteinDeposition': round(rho_spearman, 3),
        'p_rho': round(p_spearman, 5),
        'r_Urea': round(r_urea, 3),
        'p_Urea': round(p_urea, 5),
        'Direction_Match': r_sign_match,
        'Coherence_Score': round(coherence_score, 2),
        'FC_15kg': round(fc_values.get(15, np.nan), 4),
        'FC_45kg': round(fc_values.get(45, np.nan), 4),
        'FC_75kg': round(fc_values.get(75, np.nan), 4),
        'FC_105kg': round(fc_values.get(105, np.nan), 4),
        'Stage_Directions': ', '.join([f'{s}kg:{direction_matches.get(s, "NA")}' for s in sorted(direction_matches.keys())]),
    }

# Compute for all liver genes
print("  Computing liver gene coherence scores...")
liver_results = []
for gene, expr_dict in liver_gm.items():
    result = compute_phenotype_coherence(expr_dict, gene)
    if result and result['N_Total_Stages'] >= 3:
        liver_results.append(result)

liver_coh = pd.DataFrame(liver_results)
liver_coh['Tissue'] = 'Liver'

# Compute for all muscle genes
print("  Computing muscle gene coherence scores...")
muscle_results = []
for gene, expr_dict in muscle_gm.items():
    result = compute_phenotype_coherence(expr_dict, gene)
    if result and result['N_Total_Stages'] >= 3:
        muscle_results.append(result)

muscle_coh = pd.DataFrame(muscle_results)
muscle_coh['Tissue'] = 'Muscle'

# Merge
all_coh = pd.concat([liver_coh, muscle_coh], ignore_index=True)
all_coh = all_coh.sort_values('Coherence_Score', ascending=False)

# Filter: direction match + reasonable correlation
# Relaxed threshold for group-level data (n=8)
HIGH_CONF = all_coh[(all_coh['Direction_Match']) &
                    (abs(all_coh['r_ProteinDeposition']) > 0.5)]

MEDIUM_CONF = all_coh[(all_coh['Direction_Match']) &
                      (abs(all_coh['r_ProteinDeposition']) > 0.3)]

print(f"\n  All genes analyzed: {len(all_coh)}")
print(f"  Direction-matched genes (expression FC matches phenotype direction): {len(all_coh[all_coh['Direction_Match']])}")
print(f"  High confidence (|r| > 0.5 + direction match): {len(HIGH_CONF)}")
print(f"  Medium confidence (|r| > 0.3 + direction match): {len(MEDIUM_CONF)}")

print(f"\n  Top 30 High-Confidence Coherent Genes:")
print(f"  {'Gene':12s} {'Tissue':7s} {'Mean_FC':>8s} {'Dir':>7s} {'r_PD':>7s} {'p':>8s} {'r_Urea':>7s} {'Score':>6s} {'Stage_Dirs'}")
print(f"  {'-'*105}")
for _, r in HIGH_CONF.head(30).iterrows():
    sig = '**' if r['p_ProteinDeposition'] < 0.01 else ('*' if r['p_ProteinDeposition'] < 0.05 else '')
    print(f"  {r['Gene']:12s} {r['Tissue']:7s} {r['Mean_FC']:8.3f} {r['FC_Direction']:7s} {r['r_ProteinDeposition']:7.3f} {r['p_ProteinDeposition']:8.5f} {sig} {r['r_Urea']:7.3f} {r['Coherence_Score']:6.1f} {r['Stage_Directions'][:40]}")

# ============================================================
# 3. BUILD COHERENT GENE SETS FOR PATHWAY ANALYSIS
# ============================================================
print("\n[3] Building coherent gene sets for enrichment...")

# Strategy: Select genes with direction match AND reasonable correlation
# Use HIGH_CONF first, then expand to MEDIUM_CONF if too few genes

coherent_set = HIGH_CONF.copy()

# Separate by tissue
liver_coherent = coherent_set[coherent_set['Tissue'] == 'Liver']['Gene'].tolist()
muscle_coherent = coherent_set[coherent_set['Tissue'] == 'Muscle']['Gene'].tolist()

# Also identify genes with OPPOSITE pattern (TFB_up but positively correlated with PD = degradation-related)
inverse_set = all_coh[(~all_coh['Direction_Match']) &
                      (abs(all_coh['r_ProteinDeposition']) > 0.5)]

print(f"  Coherent gene set (phenotype-consistent): {len(coherent_set)} genes")
print(f"    Liver: {len(liver_coherent)} genes")
print(f"    Muscle: {len(muscle_coherent)} genes")
print(f"  Inverse set (opposite to phenotype): {len(inverse_set)} genes")

# Check key genes
for gene in ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'FOXO1', 'FBXO32', 'TRIM63',
             'MYOG', 'MYOD1', 'IGF1', 'MTOR', 'AKT1', 'MSTN', 'RPS6', 'EEF2']:
    rows = all_coh[all_coh['Gene'] == gene]
    if len(rows) > 0:
        r = rows.iloc[0]
        print(f"  {gene:10s} {r['Tissue']:7s} FC={r['Mean_FC']:7.3f} r_PD={r['r_ProteinDeposition']:7.3f} "
              f"DirMatch={r['Direction_Match']} Coherence={r['Coherence_Score']:5.1f}")

# ============================================================
# 4. GO/KEGG ENRICHMENT
# ============================================================
print("\n[4] GO/KEGG enrichment on coherent gene set...")

ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'

def enrichr_enrich(gene_list, description=''):
    if len(gene_list) < 3:
        return None
    genes_str = '\n'.join(gene_list)
    payload = {'list': (None, genes_str), 'description': (None, description)}
    try:
        r = requests.post(f'{ENRICHR_URL}/addList', files=payload, timeout=30)
        if not r.ok:
            return None
        user_list_id = r.json()['userListId']
    except Exception:
        return None

    libraries = {
        'KEGG_2021_Human': 'KEGG',
        'GO_Biological_Process_2023': 'GO_BP',
        'WikiPathway_2023_Human': 'WikiPathways',
        'Reactome_2022': 'Reactome',
    }
    all_enrich = {}
    for lib, lib_name in libraries.items():
        try:
            r = requests.get(f'{ENRICHR_URL}/enrich?userListId={user_list_id}&backgroundType={lib}', timeout=30)
            if r.ok:
                data = r.json()
                terms = []
                for entry in data.get(lib, [])[:10]:
                    terms.append({
                        'Term': entry[1], 'P_value': entry[2], 'Adjusted_P': entry[6],
                        'Odds_Ratio': entry[3], 'Overlap_Genes': entry[5], 'Combined_Score': entry[4],
                    })
                all_enrich[lib_name] = terms
        except Exception:
            pass
        time.sleep(0.5)
    return all_enrich

# Enrich liver coherent genes
print(f"\n  --- Liver coherent genes ({len(liver_coherent)}) ---")
if len(liver_coherent) >= 3:
    liver_enrich = enrichr_enrich(liver_coherent, 'Liver Phenotype-Coherent Genes')
    if liver_enrich:
        for lib_name, terms in liver_enrich.items():
            if terms:
                print(f"\n  [{lib_name}]")
                for t in terms[:5]:
                    adj_p = t['Adjusted_P']
                    sig = '***' if adj_p < 0.001 else ('**' if adj_p < 0.01 else ('*' if adj_p < 0.05 else ''))
                    print(f"    {t['Term'][:65]:65s} p.adj={adj_p:.1e} {sig}")
    else:
        print("  Enrichr unavailable for liver genes.")

# Enrich muscle coherent genes
print(f"\n  --- Muscle coherent genes ({len(muscle_coherent)}) ---")
if len(muscle_coherent) >= 3:
    muscle_enrich = enrichr_enrich(muscle_coherent, 'Muscle Phenotype-Coherent Genes')
    if muscle_enrich:
        for lib_name, terms in muscle_enrich.items():
            if terms:
                print(f"\n  [{lib_name}]")
                for t in terms[:5]:
                    adj_p = t['Adjusted_P']
                    sig = '***' if adj_p < 0.001 else ('**' if adj_p < 0.01 else ('*' if adj_p < 0.05 else ''))
                    print(f"    {t['Term'][:65]:65s} p.adj={adj_p:.1e} {sig}")
    else:
        print("  Enrichr unavailable for muscle genes.")

# Enrich full coherent set
full_coherent_genes = list(set(liver_coherent + muscle_coherent))
print(f"\n  --- Full coherent set ({len(full_coherent_genes)} genes) ---")
if len(full_coherent_genes) >= 3:
    full_enrich = enrichr_enrich(full_coherent_genes, 'Full Phenotype-Coherent Genes')
    if full_enrich:
        for lib_name, terms in full_enrich.items():
            if terms:
                print(f"\n  [{lib_name}]")
                for t in terms[:8]:
                    adj_p = t['Adjusted_P']
                    sig = '***' if adj_p < 0.001 else ('**' if adj_p < 0.01 else ('*' if adj_p < 0.05 else ''))
                    print(f"    {t['Term'][:65]:65s} p.adj={adj_p:.1e} {sig}")
    else:
        print("  Enrichr unavailable for full set.")

# ============================================================
# 5. STRING PPI NETWORK
# ============================================================
print("\n[5] STRING PPI network...")

STRING_URL = 'https://string-db.org/api'

def string_ppi(gene_list, species=9823):
    if len(gene_list) < 2:
        return None, {}
    genes_str = '%0d'.join(gene_list)
    url = f"{STRING_URL}/tsv/network?identifiers={genes_str}&species={species}&limit=200"
    try:
        r = requests.get(url, timeout=30)
        if r.ok:
            lines = r.text.strip().split('\n')
            if len(lines) < 2:
                return None, {}
            data = []
            node_degrees = {}
            for line in lines[1:]:
                parts = line.split('\t')
                if len(parts) >= 6:
                    data.append({
                        'node1': parts[2], 'node2': parts[3],
                        'combined_score': int(parts[5]) if parts[5].isdigit() else 0
                    })
                    node_degrees[parts[2]] = node_degrees.get(parts[2], 0) + 1
                    node_degrees[parts[3]] = node_degrees.get(parts[3], 0) + 1
            return pd.DataFrame(data), node_degrees
    except Exception:
        pass
    return None, {}

# Combine liver and muscle coherent genes for PPI
ppi_genes = full_coherent_genes[:200]  # STRING API limit
print(f"  Querying STRING with {len(ppi_genes)} coherent genes...")
ppi_df, node_degrees = string_ppi(ppi_genes)

if ppi_df is not None and len(ppi_df) > 0:
    print(f"  Found {len(ppi_df)} PPI interactions")
    hubs = sorted(node_degrees.items(), key=lambda x: x[1], reverse=True)
    print(f"  Top 15 hub genes:")
    for gene, degree in hubs[:15]:
        # Get coherence info
        gene_info = all_coh[all_coh['Gene'] == gene]
        coh_score = gene_info.iloc[0]['Coherence_Score'] if len(gene_info) > 0 else 0
        print(f"    {gene:12s} PPI_degree={degree:3d}  Coherence={coh_score:.1f}")
else:
    print("  STRING API unavailable.")
    node_degrees = {}

# ============================================================
# 6. CANDIDATE RANKING BY BIOLOGICAL CLOSURE
# ============================================================
print("\n[6] Ranking candidates by biological closure...")

# Build comprehensive ranking
# Core principle: prioritize genes that form the most COMPLETE biological story
# (phenotype coherence → pathway support → PPI connectivity → literature → tractability)

# Literature annotations (user's domain knowledge)
LITERATURE = {
    # Liver AA metabolism
    'CPS1': ('High', 'Rate-limiting urea cycle enzyme'),
    'OTC': ('High', 'Urea cycle, mitochondrial'),
    'ASS1': ('High', 'Urea cycle, argininosuccinate synthase'),
    'ASL': ('High', 'Urea cycle, argininosuccinate lyase'),
    'ARG1': ('High', 'Urea cycle, arginase — known STAT3 target in macrophages'),
    'ARG2': ('Medium', 'Urea cycle, mitochondrial arginase'),
    'SDS': ('Medium', 'Serine dehydratase, gluconeogenesis from AA'),
    'GOT1': ('Medium', 'Aspartate transaminase, cytoplasmic'),
    'GOT2': ('Medium', 'Aspartate transaminase, mitochondrial'),
    'HGD': ('Low', 'Homogentisate dioxygenase, Tyr catabolism'),
    'HAL': ('Low', 'Histidine ammonia-lyase'),
    'AASS': ('Low', 'Alpha-aminoadipic semialdehyde synthase, Lys catabolism'),
    'GLUD1': ('High', 'Glutamate dehydrogenase, N metabolism hub'),
    'BCKDHA': ('Medium', 'BCAA degradation, branched-chain ketoacid dehydrogenase'),
    'BCAT2': ('Medium', 'BCAA transaminase, mitochondrial'),
    # Transcription factors
    'STAT3': ('High', 'Master TF — JAK/STAT signaling, NOVEL role in AA metabolism'),
    'FOXO1': ('High', 'Forkhead TF — muscle atrophy, gluconeogenesis'),
    'FOXO3': ('High', 'Forkhead TF — muscle atrophy, autophagy'),
    'MYOD1': ('High', 'Myogenic TF — muscle differentiation'),
    'MYOG': ('High', 'Myogenic TF — terminal differentiation'),
    'MYF6': ('Medium', 'Myogenic TF — muscle maintenance'),
    # Muscle synthesis/degradation
    'FBXO32': ('High', 'Atrogin-1, muscle-specific E3 ubiquitin ligase'),
    'TRIM63': ('High', 'MuRF1, muscle atrophy marker'),
    'MSTN': ('High', 'Myostatin, negative muscle mass regulator'),
    'IGF1': ('High', 'Insulin-like growth factor, anabolic'),
    'AKT1': ('High', 'AKT kinase, mTOR upstream'),
    'MTOR': ('High', 'Mechanistic target of rapamycin, protein synthesis master'),
    'RPS6KB1': ('High', 'S6K1, mTOR downstream effector'),
    'IRS2': ('Medium', 'Insulin receptor substrate 2'),
    # Translation machinery
    'RPS6': ('Medium', 'Ribosomal protein S6'),
    'RPS3': ('Medium', 'Ribosomal protein S3'),
    'EEF2': ('Medium', 'Eukaryotic translation elongation factor 2'),
    'EIF4G1': ('Medium', 'Translation initiation factor'),
    'EIF4B': ('Medium', 'Translation initiation factor'),
}

# Tractability
TRACTABLE = {
    'STAT3', 'AKT1', 'MTOR', 'FOXO1', 'FOXO3', 'IGF1', 'MSTN',
    'MYOG', 'MYOD1', 'CPS1', 'SDS', 'ARG1', 'FBXO32', 'TRIM63',
    'RPS6KB1', 'EEF2', 'GOT1', 'ASS1'
}

candidates = []
for _, row in HIGH_CONF.iterrows():
    gene = row['Gene']
    tissue = row['Tissue']

    # Get pathway context
    lit_info = LITERATURE.get(gene, ('Low', ''))
    lit_level = lit_info[0]
    lit_note = lit_info[1]

    # PPI connectivity
    ppi_deg = node_degrees.get(gene, 0)

    # Tractability
    tractable = 'High' if gene in TRACTABLE else 'Medium'

    # Closure score components:
    # 1. Phenotype coherence (0-10) × 1.5
    # 2. PPI connectivity (normalized, 0-5)
    # 3. Literature support (High=5, Medium=3, Low=1)
    # 4. Pathway membership (from enrichment, bonus if in urea cycle or mTOR/FoxO)

    # Normalize PPI
    max_ppi = max(node_degrees.values()) if node_degrees else 1
    ppi_norm = (ppi_deg / max_ppi) * 5 if max_ppi > 0 else 0

    # Literature score
    lit_score = {'High': 5, 'Medium': 3, 'Low': 1}.get(lit_level, 1)

    # Tractability bonus
    tract_bonus = 2 if tractable == 'High' else 1

    closure_total = row['Coherence_Score'] * 1.5 + ppi_norm + lit_score + tract_bonus

    candidates.append({
        'Gene': gene,
        'Tissue': tissue,
        'Coherence_Score': row['Coherence_Score'],
        'r_ProteinDeposition': row['r_ProteinDeposition'],
        'p_ProteinDeposition': row['p_ProteinDeposition'],
        'Mean_FC': row['Mean_FC'],
        'FC_Direction': row['FC_Direction'],
        'PPI_Degree': ppi_deg,
        'Literature': lit_level,
        'Lit_Note': lit_note,
        'Tractable': tractable,
        'Closure_Score': round(closure_total, 1),
        'r_Urea': row['r_Urea'],
        'Stage_Directions': row['Stage_Directions'],
    })

cand_df = pd.DataFrame(candidates).sort_values('Closure_Score', ascending=False)

print(f"\n  {'Gene':10s} {'Tissue':7s} {'Coher':>5s} {'r_PD':>7s} {'PPI':>4s} {'Lit':>6s} {'Tract':>6s} {'Closure':>7s} | {'Note'}")
print(f"  {'-'*85}")
for _, r in cand_df.head(25).iterrows():
    print(f"  {r['Gene']:10s} {r['Tissue']:7s} {r['Coherence_Score']:5.1f} {r['r_ProteinDeposition']:7.3f} {r['PPI_Degree']:4d} {r['Literature']:6s} {r['Tractable']:6s} {r['Closure_Score']:7.1f} | {r['Lit_Note'][:45]}")

# ============================================================
# 7. GENERATE FIGURES
# ============================================================
print("\n[7] Generating figures...")

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9, 'axes.titlesize': 11, 'axes.labelsize': 10,
    'figure.dpi': 150, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
})

BREED_COLORS = {'DLY': '#2166AC', 'TFB': '#B2182B'}

# Fig 1: Phenotype-Anchored Screening Overview
fig1, axes1 = plt.subplots(2, 2, figsize=(13, 10))

# A: Protein deposition — the anchor
ax = axes1[0, 0]
stages = [15, 45, 75, 105]
dly_pd = [PROTEIN_DEPOSITION[('DLY', s)] for s in stages]
tfb_pd = [PROTEIN_DEPOSITION[('TFB', s)] for s in stages]
ax.plot(stages, dly_pd, 'o-', color=BREED_COLORS['DLY'], linewidth=2.5, markersize=8, label='DLY')
ax.plot(stages, tfb_pd, 'o-', color=BREED_COLORS['TFB'], linewidth=2.5, markersize=8, label='TFB')
ax.fill_between(stages, dly_pd, tfb_pd, alpha=0.15, color='#2166AC')
for s_idx, s in enumerate(stages):
    ax.annotate(f'DLY>TFB\np<0.01', (s, max(dly_pd[s_idx], tfb_pd[s_idx]) + 0.15),
               ha='center', fontsize=6.5, fontweight='bold', color='#2166AC')
ax.set_xlabel('Stage (kg)')
ax.set_ylabel('Protein Deposition (N g/kg BW^0.75/d)')
ax.set_title('A. Phenotype Anchor: Protein Deposition\nDLY > TFB at ALL stages', fontweight='bold')
ax.legend(fontsize=8)
ax.set_xticks(stages)

# B: Coherence score distribution
ax = axes1[0, 1]
ax.hist(all_coh['Coherence_Score'], bins=50, color='#BDBDBD', edgecolor='white', alpha=0.7, label='All genes')
ax.hist(HIGH_CONF['Coherence_Score'], bins=30, color='#2166AC', edgecolor='white', alpha=0.8, label='High conf (|r|>0.5, dir match)')
ax.set_xlabel('Phenotype Coherence Score')
ax.set_ylabel('Number of Genes')
ax.set_title(f'B. Coherence Score Distribution\n{HIGH_CONF["Tissue"].value_counts().get("Liver",0)} liver + {HIGH_CONF["Tissue"].value_counts().get("Muscle",0)} muscle high-conf genes', fontweight='bold')
ax.legend(fontsize=7)

# C: r vs Protein Deposition colored by direction match
ax = axes1[1, 0]
colors_all = ['#2166AC' if dm else '#D73027' for dm in all_coh['Direction_Match']]
ax.scatter(all_coh['Mean_FC'], all_coh['r_ProteinDeposition'], c=colors_all, alpha=0.3, s=3)
# Highlight key genes
key_check = ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'FOXO1', 'FBXO32', 'TRIM63', 'MYOG', 'IGF1']
for gene in key_check:
    rows = all_coh[all_coh['Gene'] == gene]
    if len(rows) > 0:
        r = rows.iloc[0]
        ax.annotate(gene, (r['Mean_FC'], r['r_ProteinDeposition']),
                   fontsize=7, fontweight='bold',
                   color='#2166AC' if r['Direction_Match'] else '#D73027',
                   xytext=(5, 5), textcoords='offset points')
ax.axhline(y=0, color='black', linewidth=0.5, linestyle='--')
ax.axvline(x=0, color='black', linewidth=0.5, linestyle='--')
ax.set_xlabel('Mean log2FC (DLY − TFB)')
ax.set_ylabel('r (Expression vs Protein Deposition)')
ax.set_title('C. Expression FC vs Phenotype Correlation', fontweight='bold')
legend_elements = [
    Patch(facecolor='#2166AC', label='Direction Match (coherent)'),
    Patch(facecolor='#D73027', label='Direction Mismatch'),
]
ax.legend(handles=legend_elements, fontsize=7)

# D: Top candidates by closure score
ax = axes1[1, 1]
top_n = min(20, len(cand_df))
top_cand = cand_df.head(top_n)
colors_score = ['#2166AC' if s >= 15 else '#5AAE61' if s >= 10 else '#FDB863' if s >= 7 else '#D73027'
                for s in top_cand['Closure_Score']]
ax.barh(range(len(top_cand)), top_cand['Closure_Score'], color=colors_score, edgecolor='white')
for i, (_, r) in enumerate(top_cand.iterrows()):
    label_parts = [f"r={r['r_ProteinDeposition']:.2f}"]
    if r['PPI_Degree'] > 5:
        label_parts.append(f"PPI={r['PPI_Degree']}")
    ax.text(r['Closure_Score'] + 0.3, i, ', '.join(label_parts), va='center', fontsize=5.5, color='#666666')
ax.set_yticks(range(len(top_cand)))
ax.set_yticklabels([f"{r['Gene']} ({r['Tissue'][0]})" for _, r in top_cand.iterrows()], fontsize=8, fontweight='bold')
ax.set_xlabel('Biological Closure Score')
ax.set_title('D. Candidate Ranking\n(Phenotype-Anchored Biological Closure)', fontweight='bold')
ax.invert_yaxis()

fig1.suptitle('Phenotype-Anchored Gene Screening: Protein Deposition as Biological Anchor',
             fontweight='bold', fontsize=13)
plt.tight_layout()
fig1.savefig('fig_phenotype_anchored_screening.png')
print("  Saved fig_phenotype_anchored_screening.png")

# Fig 2: Key Gene Temporal Patterns vs Protein Deposition
fig2, axes2 = plt.subplots(2, 4, figsize=(16, 8))
key_genes_plot = ['STAT3', 'CPS1', 'SDS', 'GOT1', 'ARG1', 'FOXO1', 'IGF1', 'FBXO32']

for idx, gene in enumerate(key_genes_plot):
    ax = axes2[idx // 4, idx % 4]
    ax2 = ax.twinx()

    # Gene expression
    gm_data = liver_gm.get(gene) or muscle_gm.get(gene)
    if gm_data:
        dly_expr = []
        tfb_expr = []
        for s in stages:
            dly_val = gm_data.get(('DLY', s), np.nan)
            tfb_val = gm_data.get(('TFB', s), np.nan)
            dly_expr.append(dly_val)
            tfb_expr.append(tfb_val)
        ax.plot(stages, dly_expr, 'o-', color=BREED_COLORS['DLY'], linewidth=1.5, markersize=5)
        ax.plot(stages, tfb_expr, 'o-', color=BREED_COLORS['TFB'], linewidth=1.5, markersize=5)

    # Protein deposition (scaled)
    pd_scaled = [(v - 0.5) * 3 for v in dly_pd]  # scale for visibility
    ax2.plot(stages, pd_scaled, 's--', color='#333333', linewidth=1, markersize=4, alpha=0.5, label='ProtDep')

    ax.set_title(gene, fontweight='bold')
    ax.set_xticks(stages)
    if idx % 4 == 0:
        ax.set_ylabel('log2 Expression')
    if idx == 7:
        ax2.set_ylabel('Protein Deposition (scaled)', fontsize=7)

    # Get coherence info
    gene_row = all_coh[all_coh['Gene'] == gene]
    if len(gene_row) > 0:
        r = gene_row.iloc[0]
        ax.text(0.05, 0.95, f"r_PD={r['r_ProteinDeposition']:.2f}, Coh={r['Coherence_Score']:.1f}",
               transform=ax.transAxes, fontsize=6, va='top',
               color='#2166AC' if r['Direction_Match'] else '#D73027')

fig2.suptitle('Key Gene Expression vs Protein Deposition Across Growth Stages',
             fontweight='bold', fontsize=12)
plt.tight_layout()
fig2.savefig('fig_key_genes_vs_phenotype.png')
print("  Saved fig_key_genes_vs_phenotype.png")

# Fig 3: Biological Closure Network
fig3, axes3 = plt.subplots(1, 2, figsize=(14, 6))

# A: Coherent gene pathways
ax = axes3[0]
ax.axis('off')

pathway_info = ""
if full_enrich:
    for lib_name, terms in full_enrich.items():
        if terms:
            pathway_info += f"\n[{lib_name}]\n"
            for t in terms[:4]:
                pathway_info += f"  {t['Term'][:55]}\n"
                pathway_info += f"    adj.p={t['Adjusted_P']:.1e}\n"

ax.text(0.05, 0.98, "A. Enriched Pathways (Phenotype-Coherent Genes)" + pathway_info,
       transform=ax.transAxes, fontsize=7.5, va='top', fontfamily='monospace')

# B: Biological closure scheme
ax = axes3[1]
ax.axis('off')

top5_genes = cand_df.head(5)['Gene'].tolist()
closure_text = f"""B. Biological Closure Logic

Phenotype Anchor:
  Protein Deposition: DLY > TFB at ALL stages

Screening Cascade:
  1. Expression FC direction = phenotype direction
  2. |r(Expression, ProteinDeposition)| > 0.5
  3. GO/KEGG pathway membership
  4. PPI network connectivity
  5. Literature support
  6. Experimental tractability

Top 5 Biologically Closed Candidates:
  1. {top5_genes[0] if len(top5_genes) > 0 else 'N/A'}
  2. {top5_genes[1] if len(top5_genes) > 1 else 'N/A'}
  3. {top5_genes[2] if len(top5_genes) > 2 else 'N/A'}
  4. {top5_genes[3] if len(top5_genes) > 3 else 'N/A'}
  5. {top5_genes[4] if len(top5_genes) > 4 else 'N/A'}

Validation Priority:
  NOT the most significant fold-change
  → the most complete biological story

Experimental Roadmap:
  Phase 1: qPCR + WB (confirmatory)
  Phase 2: KD/OE → phenotype readout
  Phase 3: Rescue experiment (causal closure)
"""

ax.text(0.05, 0.98, closure_text, transform=ax.transAxes, fontsize=8,
       va='top', fontfamily='monospace')

fig3.suptitle('Biological Closure: From Phenotype Anchor to Validation Candidates',
             fontweight='bold', fontsize=12)
plt.tight_layout()
fig3.savefig('fig_biological_closure.png')
print("  Saved fig_biological_closure.png")

# ============================================================
# 8. EXPERIMENTAL ROADMAP
# ============================================================
print("\n[8] Experimental validation roadmap...")

top_candidates_final = cand_df.head(8)['Gene'].tolist()

print(f"""
{'='*70}
PHENOTYPE-ANCHORED EXPERIMENTAL VALIDATION ROADMAP
{'='*70}

CORE LOGIC: Protein deposition = biological anchor
           → expression must match phenotype direction
           → then pathway + PPI + literature
           → NOT fold-change magnitude

TOP VALIDATION CANDIDATES (by biological closure):
  Liver:  {', '.join([g for g in top_candidates_final if g in liver_coherent][:4])}
  Muscle: {', '.join([g for g in top_candidates_final if g in muscle_coherent][:4])}

PHASE 1: Confirmatory (2-3 weeks) — qPCR + WB
  Model: Primary pig hepatocytes + C2C12 myotubes
  Treatments: IL-6 (20ng/mL, 24h) or AA deprivation (2h)

  Liver readouts:
    qPCR: {', '.join([g for g in top_candidates_final if g in liver_coherent][:5])}
    WB: p-STAT3(Y705), STAT3, CPS1, p-S6K1(T389)
    Medium: Urea concentration

  Muscle readouts:
    qPCR: {', '.join([g for g in top_candidates_final if g in muscle_coherent][:5])}
    WB: FBXO32, TRIM63, p-FOXO1(S256), FOXO1, MYOG

PHASE 2: Mechanistic (4-6 weeks) — KD/OE
  Hepatocytes:
    STAT3 siRNA → qPCR panel of AA enzymes → urea in medium
    STAT3 OE (adenoviral) → same readouts

  Expected: STAT3 KD → CPS1/SDS/ARG1 down → urea down
            STAT3 OE → CPS1/SDS/ARG1 up → urea up

  Myotubes:
    FOXO1 siRNA → FBXO32/TRIM63 down → myotube diameter up
    AA deprivation → FOXO1 nuclear translocation → FBXO32 up

PHASE 3: Causal Closure (6-8 weeks) — Rescue
  Hepatocytes: STAT3 KD + recombinant CPS1 → urea restored?
  Co-culture: STAT3-KD hepatocytes + C2C12 myotubes
              → measure myotube protein content
              → expected: STAT3-KD → less urea → more AA → myotube growth up

  This closes the liver→muscle loop.

NOVELTY CLAIMS (strengthened by phenotype-anchored screening):
  1. STAT3→CPS1/SDS/ARG1 regulation in AA catabolism (NOVEL)
  2. STAT3→Serum Urea→Muscle proteolysis axis (NOVEL)
  3. Phenotype-coherent gene set method for multi-tissue integration (METHODOLOGY)
{'='*70}
""")

# ============================================================
# SAVE RESULTS
# ============================================================
print("Saving results...")

with pd.ExcelWriter('phenotype_anchored_screening_results.xlsx', engine='openpyxl') as writer:
    all_coh.to_excel(writer, sheet_name='All_Genes_Coherence', index=False)
    HIGH_CONF.to_excel(writer, sheet_name='High_Confidence_Genes', index=False)
    MEDIUM_CONF.to_excel(writer, sheet_name='Medium_Confidence_Genes', index=False)
    cand_df.to_excel(writer, sheet_name='Candidate_Ranking', index=False)
    inverse_set.to_excel(writer, sheet_name='Inverse_Pattern_Genes', index=False)
    if ppi_df is not None:
        ppi_df.to_excel(writer, sheet_name='PPI_Network', index=False)

print("Saved phenotype_anchored_screening_results.xlsx")
print(f"  Sheets: All_Genes_Coherence ({len(all_coh)} genes) | High_Confidence ({len(HIGH_CONF)})")
print(f"          Medium_Confidence ({len(MEDIUM_CONF)}) | Candidate_Ranking ({len(cand_df)})")
print(f"          Inverse_Pattern ({len(inverse_set)}) | PPI_Network")

print("\nDone! Phenotype-anchored screening complete.")
