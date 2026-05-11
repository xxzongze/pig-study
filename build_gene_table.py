#!/usr/bin/env python3
"""
Build comprehensive liver-muscle axis gene table.
Search key functional categories in both tissues, compute breed×stage log2FC.
"""
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import re

# ============================================================
# 0. Load data
# ============================================================
muscle = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')
liver = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')

sample_to_group_m = {
    'm_15_1_': 'DLYM15', 'm_15_2_': 'TFBM15',
    'BJ_2_1_': 'DLYM45', 'BJ_2_2_': 'TFBM45',
    'm_1_1_': 'DLYM75', 'm_1_2_': 'TFBM75',
    'm_2_1_': 'DLYM105', 'm_2_2_': 'TFBM105',
    'm_3_1_': 'DLYM135',
}
sample_to_group_l = {
    'L_15_1_': 'DLYL15', 'L_15_2_': 'TFBL15',
    'L_45_1_': 'DLYL45', 'L_45_2_': 'TFBL45',
    'L_1_1_': 'DLYL75', 'L_1_2_': 'TFBL75',
    'L_2_1_': 'DLYL105', 'L_2_2_': 'TFBL105',
    'L_3_1_': 'DLYL135',
}

def map_col_to_group(col, mapping):
    for prefix, group in mapping.items():
        if col.startswith(prefix):
            return group
    return None

def prepare_expr_long(df, sample_map):
    val_cols = [c for c in df.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    records = []
    for _, row in df.iterrows():
        gene_id = row['seq_id']
        gene_name = str(row['gene_name']) if pd.notna(row['gene_name']) else gene_id
        for col in val_cols:
            group = map_col_to_group(col, sample_map)
            if group is not None and pd.notna(row[col]):
                records.append({
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'group': group,
                    'expr': float(row[col])
                })
    return pd.DataFrame(records)

print("Preparing expression matrices...")
muscle_long = prepare_expr_long(muscle, sample_to_group_m)
liver_long = prepare_expr_long(liver, sample_to_group_l)

# Parse group info
def parse_group(g):
    breed = 'DLY' if 'DLY' in g else 'TFB'
    stage = int(re.search(r'(\d+)', g).group(1))
    tissue = 'Liver' if 'L' in g else 'Muscle'
    return breed, stage, tissue

for df_ in [muscle_long, liver_long]:
    parsed = df_['group'].apply(parse_group)
    df_['breed'] = [p[0] for p in parsed]
    df_['stage_kg'] = [p[1] for p in parsed]
    df_['tissue'] = [p[2] for p in parsed]

# Group means
liver_mean = liver_long.groupby(['gene_id', 'gene_name', 'group', 'breed', 'stage_kg'])['expr'].mean().reset_index()
muscle_mean = muscle_long.groupby(['gene_id', 'gene_name', 'group', 'breed', 'stage_kg'])['expr'].mean().reset_index()

# ============================================================
# 1. Define gene categories and search
# ============================================================
all_liver_genes = set(liver_long['gene_name'].unique())
all_liver_upper = {g.upper(): g for g in all_liver_genes}
all_muscle_genes = set(muscle_long['gene_name'].unique())
all_muscle_upper = {g.upper(): g for g in all_muscle_genes}

def find_gene(gene_list, gene_set, gene_upper_map):
    """Find genes, return dict of {query: actual_gene_name}."""
    found = {}
    for gene in gene_list:
        if gene in gene_set:
            found[gene] = gene
        elif gene.upper() in gene_upper_map:
            found[gene] = gene_upper_map[gene.upper()]
        elif gene.lower() in gene_set:
            found[gene] = gene.lower()
    return found

def find_by_pattern(patterns, gene_set):
    """Search genes by prefix. Returns {gene_name: category}."""
    found = {}
    for pattern, category in patterns.items():
        for gene in gene_set:
            if str(gene).upper().startswith(pattern.upper()):
                found[gene] = category
    return found

# ============================================================
# 2. LIVER gene categories
# ============================================================
print("\n=== Searching LIVER genes ===")

# 2a. AA catabolism enzymes (manually curated)
liver_aa_catabolism_genes = {
    # BCAA degradation
    'BCAT2': 'BCAA degradation',
    'BCKDHA': 'BCAA degradation',
    'BCKDHB': 'BCAA degradation',
    'DBT': 'BCAA degradation',
    'DLD': 'BCAA degradation',
    'ACADSB': 'BCAA/Ile/Val degradation',
    'BCAT1': 'BCAA degradation',
    'BCKDK': 'BCAA degradation (regulator)',
    'PPM1K': 'BCAA degradation (regulator)',
    'IVD': 'Leu degradation',
    'MCCC1': 'Leu degradation',
    'MCCC2': 'Leu degradation',
    'AUH': 'Leu degradation',
    'HIBCH': 'Val degradation',
    'HIBADH': 'Val degradation',
    'ALDH6A1': 'Val degradation',
    'ACADS': 'Ile/Val degradation',
    # Urea cycle
    'CPS1': 'Urea cycle',
    'OTC': 'Urea cycle',
    'ASS1': 'Urea cycle',
    'ASL': 'Urea cycle',
    'ARG1': 'Urea cycle',
    'ARG2': 'Urea cycle',
    'NAGS': 'Urea cycle',
    'SLC25A15': 'Urea cycle (ornithine transporter)',
    # Transaminases
    'GOT1': 'Transaminase',
    'GOT2': 'Transaminase',
    'GPT': 'Transaminase',
    'GPT2': 'Transaminase',
    'PSAT1': 'Transaminase',
    'TAT': 'Transaminase',
    # Specific AA catabolism
    'AASS': 'Lys catabolism',
    'HGD': 'Tyr catabolism',
    'HPD': 'Tyr catabolism',
    'FAH': 'Tyr catabolism',
    'TAT': 'Tyr catabolism',
    'PAH': 'Phe catabolism',
    'PCBD1': 'Phe catabolism',
    'SDS': 'Ser catabolism',
    'HAL': 'His catabolism',
    'GLUD1': 'Glu catabolism',
    'GLUD2': 'Glu catabolism',
    'GLUL': 'Gln synthesis',
    'GLS': 'Gln catabolism',
    'GLS2': 'Gln catabolism',
    # Sulfur AA
    'CBS': 'Cys/Met metabolism',
    'CTH': 'Cys/Met metabolism',
    'MAT1A': 'Met metabolism',
    'MAT2A': 'Met metabolism',
    'MAT2B': 'Met metabolism',
    'AHCY': 'Met metabolism',
    'BHMT': 'Met metabolism',
    'MTR': 'Met metabolism',
    'GNMT': 'Met metabolism',
    # Trp metabolism
    'TDO2': 'Trp catabolism',
    'IDO1': 'Trp catabolism',
    'IDO2': 'Trp catabolism',
    'KYNU': 'Trp catabolism',
    'KMO': 'Trp catabolism',
    'HAAO': 'Trp catabolism',
    'AFMID': 'Trp catabolism',
    # General N metabolism
    'GLDC': 'Gly/Ser/Thr metabolism',
    'AMT': 'Gly/Ser/Thr metabolism',
    'SHMT1': 'Gly/Ser/Thr metabolism',
    'SHMT2': 'Gly/Ser/Thr metabolism',
    'AGXT': 'Gly/Ser/Thr metabolism',
    'AGXT2': 'Gly/Ser/Thr metabolism',
}

liver_aa_found = find_gene(list(liver_aa_catabolism_genes.keys()), all_liver_genes, all_liver_upper)
print(f"AA catabolism genes found: {len(liver_aa_found)}/{len(liver_aa_catabolism_genes)}")

# 2b. IGF system
igf_genes = {
    'IGF1': 'IGF system',
    'IGF2': 'IGF system',
    'IGFBP1': 'IGF system',
    'IGFBP2': 'IGF system',
    'IGFBP3': 'IGF system',
    'IGFBP4': 'IGF system',
    'IGFBP5': 'IGF system',
    'IGFBP6': 'IGF system',
    'IGFBP7': 'IGF system',
    'IGF1R': 'IGF system',
    'IGF2R': 'IGF system',
    'INSR': 'IGF system',
    'IGFALS': 'IGF system',
    'PAPPA': 'IGF system',
    'PAPPA2': 'IGF system',
    'STC1': 'IGF system',
    'STC2': 'IGF system',
    'GH1': 'GH-IGF axis',
    'GHR': 'GH-IGF axis',
    'GHRHR': 'GH-IGF axis',
    'SST': 'GH-IGF axis',
    'SSTR2': 'GH-IGF axis',
}
igf_found = find_gene(list(igf_genes.keys()), all_liver_genes, all_liver_upper)
print(f"IGF system genes found: {len(igf_found)}/{len(igf_genes)}")

# 2c. PUFA/lipid synthesis & metabolism
pufa_genes = {
    'FADS1': 'PUFA synthesis',
    'FADS2': 'PUFA synthesis',
    'ELOVL2': 'PUFA synthesis',
    'ELOVL5': 'PUFA synthesis',
    'SCD': 'Lipid metabolism',
    'SCD5': 'Lipid metabolism',
    'FASN': 'Lipid synthesis',
    'ACACA': 'Lipid synthesis',
    'ACACB': 'Lipid synthesis',
    'DGAT1': 'Lipid synthesis',
    'DGAT2': 'Lipid synthesis',
    'GPAM': 'Lipid synthesis',
    'PPARA': 'Lipid regulation',
    'PPARG': 'Lipid regulation',
    'PPARGC1A': 'Lipid regulation',
    'SREBF1': 'Lipid regulation',
    'SREBF2': 'Lipid regulation',
    'MLXIPL': 'Lipid regulation',
    'INSIG1': 'Lipid regulation',
    'INSIG2': 'Lipid regulation',
    'CPT1A': 'FA oxidation',
    'CPT2': 'FA oxidation',
    'ACOX1': 'FA oxidation',
    'ACADM': 'FA oxidation',
    'ACADL': 'FA oxidation',
    'ACADVL': 'FA oxidation',
    'HADHA': 'FA oxidation',
    'HADHB': 'FA oxidation',
    'EHHADH': 'FA oxidation',
    'HMGCS2': 'Ketogenesis',
    'BDH1': 'Ketogenesis',
    'HMGCL': 'Ketogenesis',
    'ACAT1': 'Ketogenesis',
}
pufa_found = find_gene(list(pufa_genes.keys()), all_liver_genes, all_liver_upper)
print(f"PUFA/lipid genes found: {len(pufa_found)}/{len(pufa_genes)}")

# 2d. AA transporters
aat_genes = {
    'SLC1A1': 'AA transporter',
    'SLC1A2': 'AA transporter',
    'SLC1A3': 'AA transporter',
    'SLC1A4': 'AA transporter',
    'SLC1A5': 'AA transporter',
    'SLC3A1': 'AA transporter',
    'SLC3A2': 'AA transporter',
    'SLC7A1': 'AA transporter',
    'SLC7A2': 'AA transporter',
    'SLC7A3': 'AA transporter',
    'SLC7A5': 'AA transporter',
    'SLC7A6': 'AA transporter',
    'SLC7A7': 'AA transporter',
    'SLC7A8': 'AA transporter',
    'SLC7A9': 'AA transporter',
    'SLC7A10': 'AA transporter',
    'SLC7A11': 'AA transporter',
    'SLC16A10': 'AA transporter',
    'SLC25A15': 'AA transporter',
    'SLC36A1': 'AA transporter',
    'SLC36A4': 'AA transporter',
    'SLC38A1': 'AA transporter',
    'SLC38A2': 'AA transporter',
    'SLC38A3': 'AA transporter',
    'SLC38A4': 'AA transporter',
    'SLC38A5': 'AA transporter',
    'SLC38A7': 'AA transporter',
    'SLC38A9': 'AA transporter',
    'SLC43A1': 'AA transporter',
    'SLC43A2': 'AA transporter',
}
aat_found = find_gene(list(aat_genes.keys()), all_liver_genes, all_liver_upper)
print(f"AA transporter genes found: {len(aat_found)}/{len(aat_genes)}")

# 2e. mTOR & nutrient sensing (liver)
liver_signaling = {
    'MTOR': 'Nutrient signaling',
    'RPTOR': 'Nutrient signaling',
    'RICTOR': 'Nutrient signaling',
    'AKT1': 'Insulin/AKT signaling',
    'AKT2': 'Insulin/AKT signaling',
    'AKT3': 'Insulin/AKT signaling',
    'EIF4EBP1': 'Translation control',
    'RPS6KB1': 'Translation control',
    'RPS6KB2': 'Translation control',
    'EIF4E': 'Translation control',
    'EIF4G1': 'Translation control',
    'TSC1': 'mTOR regulation',
    'TSC2': 'mTOR regulation',
    'RHEB': 'mTOR regulation',
    'LAMTOR1': 'mTOR regulation',
    'SESN2': 'AA sensing',
    'SESN1': 'AA sensing',
    'GCN2': 'AA sensing',
    'ATF4': 'AA response',
    'DDIT3': 'AA response',
    'XBP1': 'ER stress/UPR',
    'ERN1': 'ER stress/UPR',
    'EIF2AK3': 'ER stress/UPR',
    'HSPA5': 'ER stress/UPR',
    'FOXO1': 'Transcription factor',
    'FOXO3': 'Transcription factor',
    'CREB1': 'Transcription factor',
    'CRTC2': 'Transcription factor',
    'NRF1': 'Transcription factor',
    'NRF2': 'Transcription factor',
    'TFAM': 'Mitochondrial',
    'PPARGC1A': 'Mitochondrial biogenesis',
}
liver_sig_found = find_gene(list(liver_signaling.keys()), all_liver_genes, all_liver_upper)
print(f"Liver signaling genes found: {len(liver_sig_found)}/{len(liver_signaling)}")

# Build liver annotations
liver_categories = {}
for query, actual in liver_aa_found.items():
    liver_categories[actual] = liver_aa_catabolism_genes.get(query, 'AA catabolism')
for query, actual in igf_found.items():
    liver_categories[actual] = igf_genes.get(query, 'IGF system')
for query, actual in pufa_found.items():
    liver_categories[actual] = pufa_genes.get(query, 'Lipid metabolism')
for query, actual in aat_found.items():
    liver_categories[actual] = aat_genes.get(query, 'AA transporter')
for query, actual in liver_sig_found.items():
    liver_categories[actual] = liver_signaling.get(query, 'Signaling')


# ============================================================
# 3. MUSCLE gene categories
# ============================================================
print("\n=== Searching MUSCLE genes ===")

# 3a. Ribosomal proteins
muscle_ribo_patterns = {
    'RPL': 'Ribosomal large subunit (RPL)',
    'RPS': 'Ribosomal small subunit (RPS)',
    'MRPL': 'Mito-ribosomal large (MRPL)',
    'MRPS': 'Mito-ribosomal small (MRPS)',
}
muscle_ribo = find_by_pattern(muscle_ribo_patterns, all_muscle_genes)
print(f"Ribosomal genes: {len(muscle_ribo)}")

# 3b. Translation factors
muscle_transl_patterns = {
    'EIF': 'Translation initiation (EIF)',
    'EEF': 'Translation elongation (EEF)',
    'ETF': 'Translation termination',
    'ERF': 'Translation regulation',
}
muscle_transl = find_by_pattern(muscle_transl_patterns, all_muscle_genes)
print(f"Translation factor genes: {len(muscle_transl)}")

# 3c. Myogenesis & muscle development
myogenic_genes = {
    'MYOD1': 'Myogenesis',
    'MYOG': 'Myogenesis',
    'MYF5': 'Myogenesis',
    'MYF6': 'Myogenesis',
    'PAX3': 'Myogenesis',
    'PAX7': 'Myogenesis',
    'MEF2A': 'Myogenesis',
    'MEF2C': 'Myogenesis',
    'MEF2D': 'Myogenesis',
    'MYMK': 'Myogenesis',
    'MYMX': 'Myogenesis',
    'MSTN': 'Myogenesis (negative)',
    'IGF1': 'Muscle growth',
    'IGF2': 'Muscle growth',
    'IGF1R': 'Muscle growth',
    'MSTN': 'Myostatin',
    'FST': 'Follistatin',
    'GDF11': 'TGF-beta family',
    'INHBA': 'TGF-beta family',
    'TGFB1': 'TGF-beta family',
    'ACVR2A': 'Activin receptor',
    'ACVR2B': 'Activin receptor',
    'SMAD2': 'TGF-beta signaling',
    'SMAD3': 'TGF-beta signaling',
}
myo_found = find_gene(list(myogenic_genes.keys()), all_muscle_genes, all_muscle_upper)
print(f"Myogenic genes found: {len(myo_found)}/{len(myogenic_genes)}")

# 3d. Protein degradation (Ubiquitin-proteasome + autophagy)
proteolysis_genes = {
    # Ubiquitin ligases
    'FBXO32': 'Ubiquitin-proteasome',
    'TRIM63': 'Ubiquitin-proteasome',
    'MURF1': 'Ubiquitin-proteasome',
    'FBXO30': 'Ubiquitin-proteasome',
    'FBXO40': 'Ubiquitin-proteasome',
    'TRIM54': 'Ubiquitin-proteasome',
    'TRIM55': 'Ubiquitin-proteasome',
    'TRIM72': 'Ubiquitin-proteasome',
    # Ubiquitin system
    'UBB': 'Ubiquitin-proteasome',
    'UBC': 'Ubiquitin-proteasome',
    'UBE2B': 'Ubiquitin-proteasome',
    'UBE2D1': 'Ubiquitin-proteasome',
    'UBE2D2': 'Ubiquitin-proteasome',
    'UBE2D3': 'Ubiquitin-proteasome',
    'PSMA1': 'Proteasome subunit',
    'PSMA2': 'Proteasome subunit',
    'PSMA3': 'Proteasome subunit',
    'PSMB1': 'Proteasome subunit',
    'PSMB2': 'Proteasome subunit',
    'PSMB5': 'Proteasome subunit',
    'PSMC1': 'Proteasome subunit',
    'PSMC2': 'Proteasome subunit',
    'PSMD1': 'Proteasome subunit',
    'PSMD2': 'Proteasome subunit',
    'PSMD4': 'Proteasome subunit',
    # Autophagy
    'ATG3': 'Autophagy',
    'ATG4B': 'Autophagy',
    'ATG5': 'Autophagy',
    'ATG7': 'Autophagy',
    'ATG12': 'Autophagy',
    'BECN1': 'Autophagy',
    'BECN2': 'Autophagy',
    'SQSTM1': 'Autophagy',
    'MAP1LC3A': 'Autophagy',
    'MAP1LC3B': 'Autophagy',
    'GABARAPL1': 'Autophagy',
    'GABARAPL2': 'Autophagy',
    'ULK1': 'Autophagy',
    'ULK2': 'Autophagy',
    'BNIP3': 'Autophagy',
    'BNIP3L': 'Autophagy',
    # Calpain system
    'CAPN1': 'Calpain system',
    'CAPN2': 'Calpain system',
    'CAPN3': 'Calpain system',
    'CAST': 'Calpain system',
}
prot_found = find_gene(list(proteolysis_genes.keys()), all_muscle_genes, all_muscle_upper)
print(f"Proteolysis genes found: {len(prot_found)}/{len(proteolysis_genes)}")

# 3e. Muscle structural & contractile
muscle_struct_genes = {
    'MYH1': 'Myosin heavy chain',
    'MYH2': 'Myosin heavy chain',
    'MYH3': 'Myosin heavy chain',
    'MYH4': 'Myosin heavy chain',
    'MYH7': 'Myosin heavy chain',
    'MYH8': 'Myosin heavy chain',
    'MYL1': 'Myosin light chain',
    'MYL2': 'Myosin light chain',
    'MYL3': 'Myosin light chain',
    'MYL4': 'Myosin light chain',
    'MYLPF': 'Myosin light chain',
    'ACTA1': 'Actin',
    'ACTC1': 'Actin',
    'ACTN2': 'Actinin',
    'ACTN3': 'Actinin',
    'TNNT1': 'Troponin',
    'TNNT2': 'Troponin',
    'TNNT3': 'Troponin',
    'TNNI1': 'Troponin',
    'TNNI2': 'Troponin',
    'TNNC1': 'Troponin',
    'TPM1': 'Tropomyosin',
    'TPM2': 'Tropomyosin',
    'TPM3': 'Tropomyosin',
    'TTN': 'Titin',
    'NEB': 'Nebulin',
    'DES': 'Desmin',
    'DMD': 'Dystrophin',
    'CSRP3': 'Muscle LIM protein',
    'MUSTN1': 'Musculoskeletal',
}
struct_found = find_gene(list(muscle_struct_genes.keys()), all_muscle_genes, all_muscle_upper)
print(f"Structural genes found: {len(struct_found)}/{len(muscle_struct_genes)}")

# 3f. Muscle mTOR/growth signaling
muscle_signaling = {
    'MTOR': 'mTOR signaling',
    'RPTOR': 'mTOR signaling',
    'RICTOR': 'mTOR signaling',
    'AKT1': 'AKT signaling',
    'AKT2': 'AKT signaling',
    'RPS6KB1': 'S6K/4E-BP',
    'RPS6KB2': 'S6K/4E-BP',
    'EIF4EBP1': 'S6K/4E-BP',
    'EIF4E': 'Cap-dependent translation',
    'EIF4G1': 'Cap-dependent translation',
    'EIF4A1': 'Cap-dependent translation',
    'EIF4A2': 'Cap-dependent translation',
    'MYC': 'Transcription factor',
    'MYCN': 'Transcription factor',
    'JUN': 'Transcription factor',
    'FOS': 'Transcription factor',
    'SRF': 'Transcription factor',
    'TEAD1': 'Transcription factor',
    'TEAD4': 'Transcription factor',
    'YAP1': 'Hippo pathway',
    'WWTR1': 'Hippo pathway',
    'LATS1': 'Hippo pathway',
    'LATS2': 'Hippo pathway',
    'MST1': 'Hippo pathway',
    'SAV1': 'Hippo pathway',
}
muscle_sig_found = find_gene(list(muscle_signaling.keys()), all_muscle_genes, all_muscle_upper)
print(f"Signaling genes found: {len(muscle_sig_found)}/{len(muscle_signaling)}")

# Build muscle annotations
muscle_categories = {}
for g, cat in muscle_ribo.items():
    muscle_categories[g] = cat
for g, cat in muscle_transl.items():
    muscle_categories[g] = cat
for query, actual in myo_found.items():
    muscle_categories[actual] = myogenic_genes.get(query, 'Myogenesis')
for query, actual in prot_found.items():
    muscle_categories[actual] = proteolysis_genes.get(query, 'Proteolysis')
for query, actual in struct_found.items():
    muscle_categories[actual] = muscle_struct_genes.get(query, 'Structure')
for query, actual in muscle_sig_found.items():
    muscle_categories[actual] = muscle_signaling.get(query, 'Signaling')


# ============================================================
# 4. Compute log2FC and build table
# ============================================================
def compute_gene_stats(df, gene_name, breed_col='breed', stage_col='stage_kg', expr_col='expr'):
    """For a gene, compute DLY and TFB means per stage, log2FC, and t-test p."""
    gene_df = df[df['gene_name'] == gene_name]
    results = {}
    for stage in sorted(gene_df[stage_col].unique()):
        dly = gene_df[(gene_df[breed_col] == 'DLY') & (gene_df[stage_col] == stage)][expr_col]
        tfb = gene_df[(gene_df[breed_col] == 'TFB') & (gene_df[stage_col] == stage)][expr_col]
        dly_mean = dly.mean() if len(dly) > 0 else np.nan
        tfb_mean = tfb.mean() if len(tfb) > 0 else np.nan
        if tfb_mean > 0 and not np.isnan(dly_mean):
            log2fc = np.log2(dly_mean / tfb_mean)
        else:
            log2fc = np.nan
        if len(dly) > 1 and len(tfb) > 1:
            _, p = ttest_ind(dly, tfb)
        else:
            p = np.nan
        results[stage] = {
            'DLY_mean': dly_mean,
            'TFB_mean': tfb_mean,
            'log2FC': log2fc,
            'p_value': p,
            'n_DLY': len(dly),
            'n_TFB': len(tfb),
        }
    return results

# ============================================================
# 5. Build the output table
# ============================================================
print("\n=== Building summary table ===")

rows = []
stages_of_interest = [15, 45, 75, 105]

# Liver genes
for gene_name, category in liver_categories.items():
    stats = compute_gene_stats(liver_mean, gene_name)
    row = {
        'Tissue': 'Liver',
        'Gene': gene_name,
        'Category': category,
    }
    for s in stages_of_interest:
        if s in stats:
            row[f'{s}kg_log2FC'] = round(stats[s]['log2FC'], 3) if not np.isnan(stats[s]['log2FC']) else ''
            row[f'{s}kg_DLY_mean'] = round(stats[s]['DLY_mean'], 2) if not np.isnan(stats[s]['DLY_mean']) else ''
            row[f'{s}kg_TFB_mean'] = round(stats[s]['TFB_mean'], 2) if not np.isnan(stats[s]['TFB_mean']) else ''
            row[f'{s}kg_p'] = f"{stats[s]['p_value']:.4f}" if not np.isnan(stats[s]['p_value']) else ''
        else:
            for suffix in ['_log2FC', '_DLY_mean', '_TFB_mean', '_p']:
                row[f'{s}kg{suffix}'] = ''
    # Mean abs log2FC across stages
    fcs = [abs(stats[s]['log2FC']) for s in stages_of_interest if s in stats and not np.isnan(stats[s]['log2FC'])]
    row['mean_abs_log2FC'] = round(np.mean(fcs), 3) if fcs else ''
    rows.append(row)

# Muscle genes
for gene_name, category in muscle_categories.items():
    stats = compute_gene_stats(muscle_mean, gene_name)
    row = {
        'Tissue': 'Muscle',
        'Gene': gene_name,
        'Category': category,
    }
    for s in stages_of_interest:
        if s in stats:
            row[f'{s}kg_log2FC'] = round(stats[s]['log2FC'], 3) if not np.isnan(stats[s]['log2FC']) else ''
            row[f'{s}kg_DLY_mean'] = round(stats[s]['DLY_mean'], 2) if not np.isnan(stats[s]['DLY_mean']) else ''
            row[f'{s}kg_TFB_mean'] = round(stats[s]['TFB_mean'], 2) if not np.isnan(stats[s]['TFB_mean']) else ''
            row[f'{s}kg_p'] = f"{stats[s]['p_value']:.4f}" if not np.isnan(stats[s]['p_value']) else ''
        else:
            for suffix in ['_log2FC', '_DLY_mean', '_TFB_mean', '_p']:
                row[f'{s}kg{suffix}'] = ''
    fcs = [abs(stats[s]['log2FC']) for s in stages_of_interest if s in stats and not np.isnan(stats[s]['log2FC'])]
    row['mean_abs_log2FC'] = round(np.mean(fcs), 3) if fcs else ''
    rows.append(row)

df_out = pd.DataFrame(rows)

# Sort by Tissue, then Category, then mean_abs_log2FC descending
cat_order = {
    # Liver categories in logical order
    'BCAA degradation': 1, 'BCAA/Ile/Val degradation': 2, 'Leu degradation': 3, 'Val degradation': 4,
    'Urea cycle': 5, 'Urea cycle (ornithine transporter)': 6,
    'Transaminase': 7,
    'Lys catabolism': 8, 'Tyr catabolism': 9, 'Phe catabolism': 10,
    'Ser catabolism': 11, 'His catabolism': 12,
    'Glu catabolism': 13, 'Gln synthesis': 14, 'Gln catabolism': 15,
    'Cys/Met metabolism': 16, 'Met metabolism': 17, 'Trp catabolism': 18,
    'Gly/Ser/Thr metabolism': 19,
    'IGF system': 20, 'GH-IGF axis': 21,
    'PUFA synthesis': 22, 'Lipid metabolism': 23, 'Lipid synthesis': 24,
    'FA oxidation': 25, 'Ketogenesis': 26, 'Lipid regulation': 27,
    'AA transporter': 28,
    'Nutrient signaling': 29, 'Insulin/AKT signaling': 30,
    'mTOR regulation': 31, 'AA sensing': 32, 'AA response': 33,
    'Translation control': 34, 'ER stress/UPR': 35,
    'Transcription factor': 36, 'Mitochondrial': 37, 'Mitochondrial biogenesis': 38,
}
# Muscle category order
muscle_cat_order = {
    'Ribosomal large subunit (RPL)': 100,
    'Ribosomal small subunit (RPS)': 101,
    'Mito-ribosomal large (MRPL)': 102,
    'Mito-ribosomal small (MRPS)': 103,
    'Translation initiation (EIF)': 104,
    'Translation elongation (EEF)': 105,
    'Translation termination': 106, 'Translation regulation': 107,
    'mTOR signaling': 108, 'AKT signaling': 109,
    'S6K/4E-BP': 110, 'Cap-dependent translation': 111,
    'Myogenesis': 112, 'Myogenesis (negative)': 113, 'Myostatin': 114,
    'Follistatin': 115, 'Muscle growth': 116,
    'TGF-beta family': 117, 'TGF-beta signaling': 118, 'Activin receptor': 119,
    'Ubiquitin-proteasome': 120, 'Proteasome subunit': 121,
    'Autophagy': 122, 'Calpain system': 123,
    'Myosin heavy chain': 124, 'Myosin light chain': 125,
    'Actin': 126, 'Actinin': 127,
    'Troponin': 128, 'Tropomyosin': 129,
    'Titin': 130, 'Nebulin': 131, 'Desmin': 132, 'Dystrophin': 133,
    'Muscle LIM protein': 134, 'Musculoskeletal': 135,
    'Transcription factor': 136, 'Hippo pathway': 137,
}
all_cat_order = {**cat_order, **muscle_cat_order}

df_out['cat_order'] = df_out['Category'].map(all_cat_order).fillna(999)
df_out['tissue_order'] = df_out['Tissue'].map({'Liver': 0, 'Muscle': 1})
df_out['abs_fc_num'] = pd.to_numeric(df_out['mean_abs_log2FC'], errors='coerce').fillna(0)

df_out = df_out.sort_values(['tissue_order', 'cat_order', 'abs_fc_num'],
                             ascending=[True, True, False])
df_out = df_out.drop(columns=['cat_order', 'tissue_order', 'abs_fc_num'])

# ============================================================
# 6. Save to Excel (with separate sheets for Liver and Muscle)
# ============================================================
print(f"\nTotal genes in table: {len(df_out)}")
print(f"  Liver: {len(df_out[df_out['Tissue']=='Liver'])}")
print(f"  Muscle: {len(df_out[df_out['Tissue']=='Muscle'])}")

# Simplify: for a clean summary table, keep only key columns
summary_cols = ['Tissue', 'Gene', 'Category']
for s in stages_of_interest:
    summary_cols.append(f'{s}kg_log2FC')
summary_cols.append('mean_abs_log2FC')

df_summary = df_out[summary_cols].copy()

with pd.ExcelWriter('liver_muscle_axis_genes.xlsx', engine='openpyxl') as writer:
    # Full table
    df_out.to_excel(writer, sheet_name='All_Genes_Detailed', index=False)

    # Summary (log2FC only)
    df_summary.to_excel(writer, sheet_name='Summary_log2FC', index=False)

    # Liver only
    liver_df = df_out[df_out['Tissue'] == 'Liver'].drop(columns=['Tissue'])
    liver_df.to_excel(writer, sheet_name='Liver_Genes', index=False)

    # Muscle only
    muscle_df = df_out[df_out['Tissue'] == 'Muscle'].drop(columns=['Tissue'])
    muscle_df.to_excel(writer, sheet_name='Muscle_Genes', index=False)

print(f"\nSaved: liver_muscle_axis_genes.xlsx")
print(f"  Sheet 1: All_Genes_Detailed — all genes with expression means and p-values")
print(f"  Sheet 2: Summary_log2FC — log2(DLY/TFB) for each stage + mean |FC|")
print(f"  Sheet 3: Liver_Genes — liver genes only")
print(f"  Sheet 4: Muscle_Genes — muscle genes only")

# Print summary stats
print("\n=== Gene counts by category ===")
cat_counts = df_out.groupby(['Tissue', 'Category']).size().reset_index(name='count')
for tissue in ['Liver', 'Muscle']:
    print(f"\n{tissue}:")
    for _, row in cat_counts[cat_counts['Tissue'] == tissue].iterrows():
        print(f"  {row['Category']}: {row['count']} genes")
