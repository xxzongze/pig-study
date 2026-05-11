#!/usr/bin/env python3
"""
Build liver-muscle cross-talk gene table.
Genes that mediate inter-tissue communication: myokines, hepatokines,
shared signaling, AA sensors, transcription factors, hormone axes.
"""
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import re

# ============================================================
# 0. Load data
# ============================================================
muscle = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix.xls', sep='\t')
liver = pd.read_csv('/Users/hezongze/Downloads/result_exp_matrix (4).xls', sep='\t')

sample_to_group_m = {
    'm_15_1_': ('DLY', 15, 'DLYM15'), 'm_15_2_': ('TFB', 15, 'TFBM15'),
    'BJ_2_1_': ('DLY', 45, 'DLYM45'), 'BJ_2_2_': ('TFB', 45, 'TFBM45'),
    'm_1_1_': ('DLY', 75, 'DLYM75'), 'm_1_2_': ('TFB', 75, 'TFBM75'),
    'm_2_1_': ('DLY', 105, 'DLYM105'), 'm_2_2_': ('TFB', 105, 'TFBM105'),
    'm_3_1_': ('DLY', 135, 'DLYM135'),
}
sample_to_group_l = {
    'L_15_1_': ('DLY', 15, 'DLYL15'), 'L_15_2_': ('TFB', 15, 'TFBL15'),
    'L_45_1_': ('DLY', 45, 'DLYL45'), 'L_45_2_': ('TFB', 45, 'TFBL45'),
    'L_1_1_': ('DLY', 75, 'DLYL75'), 'L_1_2_': ('TFB', 75, 'TFBL75'),
    'L_2_1_': ('DLY', 105, 'DLYL105'), 'L_2_2_': ('TFB', 105, 'TFBL105'),
    'L_3_1_': ('DLY', 135, 'DLYL135'),
}

def prepare_matrix(df, sample_map):
    val_cols = [c for c in df.columns if c not in ('seq_id', 'gene_name', 'length', 'description')]
    records = []
    for _, row in df.iterrows():
        gene_id = row['seq_id']
        gene_name = str(row['gene_name']) if pd.notna(row['gene_name']) else gene_id
        for col in val_cols:
            info = None
            for prefix, (breed, stage, group) in sample_map.items():
                if col.startswith(prefix):
                    info = (breed, stage, group)
                    break
            if info and pd.notna(row[col]):
                records.append({
                    'gene_name': gene_name,
                    'gene_id': gene_id,
                    'breed': info[0],
                    'stage_kg': info[1],
                    'group': info[2],
                    'expr': float(row[col])
                })
    return pd.DataFrame(records)

print("Preparing expression matrices...")
ml = prepare_matrix(muscle, sample_to_group_m)
ll = prepare_matrix(liver, sample_to_group_l)
ml['tissue'] = 'Muscle'
ll['tissue'] = 'Liver'

muscle_genes = set(ml['gene_name'].unique())
liver_genes = set(ll['gene_name'].unique())
muscle_upper = {g.upper(): g for g in muscle_genes}
liver_upper = {g.upper(): g for g in liver_genes}

# Group means
muscle_mean = ml.groupby(['gene_name', 'breed', 'stage_kg', 'tissue'])['expr'].mean().reset_index()
liver_mean = ll.groupby(['gene_name', 'breed', 'stage_kg', 'tissue'])['expr'].mean().reset_index()

# ============================================================
# 1. Define cross-talk gene categories
# ============================================================
def search_genes(gene_dict, muscle_set, liver_set, muscle_upper, liver_upper):
    """Search for genes in both tissues. Returns dict with tissue presence info."""
    results = {}
    for gene, info in gene_dict.items():
        in_muscle = gene in muscle_set or gene.upper() in muscle_upper
        in_liver = gene in liver_set or gene.upper() in liver_upper
        m_name = gene if gene in muscle_set else muscle_upper.get(gene.upper(), gene)
        l_name = gene if gene in liver_set else liver_upper.get(gene.upper(), gene)
        results[gene] = {
            'category': info.get('category', info) if isinstance(info, dict) else info,
            'subcategory': info.get('subcategory', '') if isinstance(info, dict) else '',
            'in_muscle': in_muscle,
            'in_liver': in_liver,
            'muscle_name': m_name,
            'liver_name': l_name,
        }
    return results

# ============================================================
# 1a. Myokines (muscle → liver)
# ============================================================
myokines = {
    'IL6': {'category': 'Myokine', 'subcategory': 'Cytokine'},
    'IL6R': {'category': 'Myokine', 'subcategory': 'Receptor'},
    'IL6ST': {'category': 'Myokine', 'subcategory': 'Receptor (gp130)'},
    'FNDC5': {'category': 'Myokine', 'subcategory': 'Irisin precursor'},
    'MSTN': {'category': 'Myokine', 'subcategory': 'Myostatin/GDF8'},
    'BDNF': {'category': 'Myokine', 'subcategory': 'Neurotrophin'},
    'SPARC': {'category': 'Myokine', 'subcategory': 'Osteonectin'},
    'DCN': {'category': 'Myokine', 'subcategory': 'Decorin'},
    'LIF': {'category': 'Myokine', 'subcategory': 'Cytokine'},
    'CTF1': {'category': 'Myokine', 'subcategory': 'Cytokine'},
    'VEGFA': {'category': 'Myokine', 'subcategory': 'Angiogenic'},
    'VEGFB': {'category': 'Myokine', 'subcategory': 'Angiogenic'},
    'NRG1': {'category': 'Myokine', 'subcategory': 'Neuregulin'},
    'NRG4': {'category': 'Myokine', 'subcategory': 'Neuregulin'},
    'MIF': {'category': 'Myokine', 'subcategory': 'Cytokine'},
    'FSTL1': {'category': 'Myokine', 'subcategory': 'Follistatin-like'},
    'METRNL': {'category': 'Myokine', 'subcategory': 'Meteorin-like'},
    'APLN': {'category': 'Myokine', 'subcategory': 'Apelin'},
    'CHI3L1': {'category': 'Myokine', 'subcategory': 'Chitinase-like'},
    'GDF15': {'category': 'Myokine', 'subcategory': 'TGF-beta family'},
    'IGF1': {'category': 'Myokine', 'subcategory': 'Growth factor (also hepatokine)'},
    'IGF2': {'category': 'Myokine', 'subcategory': 'Growth factor'},
    'FGF21': {'category': 'Myokine', 'subcategory': 'FGF family (also hepatokine)'},
}

# ============================================================
# 1b. Hepatokines (liver → muscle)
# ============================================================
hepatokines = {
    'ANGPTL3': {'category': 'Hepatokine', 'subcategory': 'Angiopoietin-like'},
    'ANGPTL4': {'category': 'Hepatokine', 'subcategory': 'Angiopoietin-like'},
    'ANGPTL6': {'category': 'Hepatokine', 'subcategory': 'Angiopoietin-like'},
    'ANGPTL8': {'category': 'Hepatokine', 'subcategory': 'Angiopoietin-like'},
    'FST': {'category': 'Hepatokine', 'subcategory': 'Follistatin'},
    'SHBG': {'category': 'Hepatokine', 'subcategory': 'Steroid binding'},
    'SELENOP': {'category': 'Hepatokine', 'subcategory': 'Selenium transport'},
    'AHSG': {'category': 'Hepatokine', 'subcategory': 'Fetuin-A'},
    'RBP4': {'category': 'Hepatokine', 'subcategory': 'Retinol binding'},
    'LECT2': {'category': 'Hepatokine', 'subcategory': 'Leukocyte chemotactic'},
    'HAMP': {'category': 'Hepatokine', 'subcategory': 'Hepcidin'},
    'TTR': {'category': 'Hepatokine', 'subcategory': 'Transthyretin'},
    'APOA1': {'category': 'Hepatokine', 'subcategory': 'Apolipoprotein'},
    'APOA2': {'category': 'Hepatokine', 'subcategory': 'Apolipoprotein'},
    'APOB': {'category': 'Hepatokine', 'subcategory': 'Apolipoprotein'},
    'APOC3': {'category': 'Hepatokine', 'subcategory': 'Apolipoprotein'},
    'APOE': {'category': 'Hepatokine', 'subcategory': 'Apolipoprotein'},
    'VTN': {'category': 'Hepatokine', 'subcategory': 'Vitronectin'},
    'ORM1': {'category': 'Hepatokine', 'subcategory': 'Acute phase'},
    'ORM2': {'category': 'Hepatokine', 'subcategory': 'Acute phase'},
    'HP': {'category': 'Hepatokine', 'subcategory': 'Haptoglobin'},
    'TTR': {'category': 'Hepatokine', 'subcategory': 'Transthyretin'},
    'C3': {'category': 'Hepatokine', 'subcategory': 'Complement'},
    'KNG1': {'category': 'Hepatokine', 'subcategory': 'Kininogen'},
    'SERPINA1': {'category': 'Hepatokine', 'subcategory': 'Serpin'},
    'SERPINC1': {'category': 'Hepatokine', 'subcategory': 'Serpin'},
    'F2': {'category': 'Hepatokine', 'subcategory': 'Coagulation'},
    'FGA': {'category': 'Hepatokine', 'subcategory': 'Fibrinogen'},
    'FGB': {'category': 'Hepatokine', 'subcategory': 'Fibrinogen'},
    'FGG': {'category': 'Hepatokine', 'subcategory': 'Fibrinogen'},
    'IGF1': {'category': 'Hepatokine', 'subcategory': 'IGF system'},
    'IGF2': {'category': 'Hepatokine', 'subcategory': 'IGF system'},
    'IGFBP1': {'category': 'Hepatokine', 'subcategory': 'IGFBP'},
    'IGFBP2': {'category': 'Hepatokine', 'subcategory': 'IGFBP'},
    'IGFBP3': {'category': 'Hepatokine', 'subcategory': 'IGFBP'},
    'IGFBP4': {'category': 'Hepatokine', 'subcategory': 'IGFBP'},
    'FGF21': {'category': 'Hepatokine', 'subcategory': 'FGF family'},
}

# ============================================================
# 1c. Shared AA / nutrient sensors
# ============================================================
aa_sensors = {
    'MTOR': {'category': 'AA/Nutrient sensing', 'subcategory': 'mTORC1'},
    'RPTOR': {'category': 'AA/Nutrient sensing', 'subcategory': 'mTORC1'},
    'RICTOR': {'category': 'AA/Nutrient sensing', 'subcategory': 'mTORC2'},
    'DEPTOR': {'category': 'AA/Nutrient sensing', 'subcategory': 'mTOR regulator'},
    'AKT1': {'category': 'AA/Nutrient sensing', 'subcategory': 'AKT pathway'},
    'AKT2': {'category': 'AA/Nutrient sensing', 'subcategory': 'AKT pathway'},
    'EIF2AK4': {'category': 'AA/Nutrient sensing', 'subcategory': 'GCN2/AA sensor'},
    'EIF2S1': {'category': 'AA/Nutrient sensing', 'subcategory': 'eIF2α'},
    'EIF2AK3': {'category': 'AA/Nutrient sensing', 'subcategory': 'PERK/ER stress'},
    'ATF4': {'category': 'AA/Nutrient sensing', 'subcategory': 'ISR transcription factor'},
    'DDIT3': {'category': 'AA/Nutrient sensing', 'subcategory': 'CHOP'},
    'SESN1': {'category': 'AA/Nutrient sensing', 'subcategory': 'Sestrin/Leu sensor'},
    'SESN2': {'category': 'AA/Nutrient sensing', 'subcategory': 'Sestrin/Leu sensor'},
    'SLC38A2': {'category': 'AA/Nutrient sensing', 'subcategory': 'SNAT2/AA transporter'},
    'SLC38A9': {'category': 'AA/Nutrient sensing', 'subcategory': 'Lysosomal Arg sensor'},
    'LAMTOR1': {'category': 'AA/Nutrient sensing', 'subcategory': 'Ragulator complex'},
    'LAMTOR2': {'category': 'AA/Nutrient sensing', 'subcategory': 'Ragulator complex'},
    'LAMTOR3': {'category': 'AA/Nutrient sensing', 'subcategory': 'Ragulator complex'},
    'LAMTOR4': {'category': 'AA/Nutrient sensing', 'subcategory': 'Ragulator complex'},
    'LAMTOR5': {'category': 'AA/Nutrient sensing', 'subcategory': 'Ragulator complex'},
    'RRAGA': {'category': 'AA/Nutrient sensing', 'subcategory': 'Rag GTPase'},
    'RRAGB': {'category': 'AA/Nutrient sensing', 'subcategory': 'Rag GTPase'},
    'RRAGC': {'category': 'AA/Nutrient sensing', 'subcategory': 'Rag GTPase'},
    'TSC1': {'category': 'AA/Nutrient sensing', 'subcategory': 'TSC complex'},
    'TSC2': {'category': 'AA/Nutrient sensing', 'subcategory': 'TSC complex'},
    'RHEB': {'category': 'AA/Nutrient sensing', 'subcategory': 'Rheb GTPase'},
    'RPS6KB1': {'category': 'AA/Nutrient sensing', 'subcategory': 'S6K1 effector'},
    'RPS6KB2': {'category': 'AA/Nutrient sensing', 'subcategory': 'S6K2 effector'},
    'EIF4EBP1': {'category': 'AA/Nutrient sensing', 'subcategory': '4E-BP1 effector'},
    'PRKAA1': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK α1'},
    'PRKAA2': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK α2'},
    'PRKAB1': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK β1'},
    'PRKAB2': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK β2'},
    'PRKAG1': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK γ1'},
    'PRKAG2': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK γ2'},
    'PRKAG3': {'category': 'AA/Nutrient sensing', 'subcategory': 'AMPK γ3'},
    'SIRT1': {'category': 'AA/Nutrient sensing', 'subcategory': 'Sirtuin/NAD+ sensor'},
    'SIRT3': {'category': 'AA/Nutrient sensing', 'subcategory': 'Sirtuin/mito'},
    'CAST': {'category': 'AA/Nutrient sensing', 'subcategory': 'mTOR-AA sensor complex'},
    'WDR59': {'category': 'AA/Nutrient sensing', 'subcategory': 'GATOR2 complex'},
    'WDR24': {'category': 'AA/Nutrient sensing', 'subcategory': 'GATOR2 complex'},
    'NPRL2': {'category': 'AA/Nutrient sensing', 'subcategory': 'GATOR1 complex'},
    'NPRL3': {'category': 'AA/Nutrient sensing', 'subcategory': 'GATOR1 complex'},
    'DEPDC5': {'category': 'AA/Nutrient sensing', 'subcategory': 'GATOR1 complex'},
}

# ============================================================
# 1d. Shared transcription factors / co-regulators
# ============================================================
shared_tfs = {
    'FOXO1': {'category': 'Shared TF/co-regulator', 'subcategory': 'FOXO family'},
    'FOXO3': {'category': 'Shared TF/co-regulator', 'subcategory': 'FOXO family'},
    'FOXO4': {'category': 'Shared TF/co-regulator', 'subcategory': 'FOXO family'},
    'PPARGC1A': {'category': 'Shared TF/co-regulator', 'subcategory': 'PGC1α/mito biogenesis'},
    'PPARA': {'category': 'Shared TF/co-regulator', 'subcategory': 'PPAR family'},
    'PPARD': {'category': 'Shared TF/co-regulator', 'subcategory': 'PPAR family'},
    'PPARG': {'category': 'Shared TF/co-regulator', 'subcategory': 'PPAR family'},
    'NFE2L2': {'category': 'Shared TF/co-regulator', 'subcategory': 'NRF2/oxidative stress'},
    'TFEB': {'category': 'Shared TF/co-regulator', 'subcategory': 'Lysosome/autophagy'},
    'TFE3': {'category': 'Shared TF/co-regulator', 'subcategory': 'Lysosome/autophagy'},
    'XBP1': {'category': 'Shared TF/co-regulator', 'subcategory': 'UPR/ER stress'},
    'ATF6': {'category': 'Shared TF/co-regulator', 'subcategory': 'UPR/ER stress'},
    'CREB1': {'category': 'Shared TF/co-regulator', 'subcategory': 'cAMP/PKA'},
    'CREBBP': {'category': 'Shared TF/co-regulator', 'subcategory': 'Co-activator'},
    'STAT3': {'category': 'Shared TF/co-regulator', 'subcategory': 'JAK/STAT'},
    'STAT5A': {'category': 'Shared TF/co-regulator', 'subcategory': 'JAK/STAT'},
    'STAT5B': {'category': 'Shared TF/co-regulator', 'subcategory': 'JAK/STAT'},
    'SMAD2': {'category': 'Shared TF/co-regulator', 'subcategory': 'TGF-β pathway'},
    'SMAD3': {'category': 'Shared TF/co-regulator', 'subcategory': 'TGF-β pathway'},
    'SMAD4': {'category': 'Shared TF/co-regulator', 'subcategory': 'TGF-β pathway'},
    'HIF1A': {'category': 'Shared TF/co-regulator', 'subcategory': 'Hypoxia'},
    'NRF1': {'category': 'Shared TF/co-regulator', 'subcategory': 'Mitochondrial biogenesis'},
    'TFAM': {'category': 'Shared TF/co-regulator', 'subcategory': 'mtDNA transcription'},
    'ESRRA': {'category': 'Shared TF/co-regulator', 'subcategory': 'ERRα/energy metabolism'},
    'RXRA': {'category': 'Shared TF/co-regulator', 'subcategory': 'Nuclear receptor'},
    'RXRG': {'category': 'Shared TF/co-regulator', 'subcategory': 'Nuclear receptor'},
    'KLF15': {'category': 'Shared TF/co-regulator', 'subcategory': 'AA/glucose metabolism'},
    'CEBPA': {'category': 'Shared TF/co-regulator', 'subcategory': 'Adipogenic/liver'},
    'CEBPB': {'category': 'Shared TF/co-regulator', 'subcategory': 'Adipogenic/liver'},
    'CEBPD': {'category': 'Shared TF/co-regulator', 'subcategory': 'Adipogenic/liver'},
    'MYC': {'category': 'Shared TF/co-regulator', 'subcategory': 'Growth/ribosome'},
    'MYCN': {'category': 'Shared TF/co-regulator', 'subcategory': 'Growth/ribosome'},
    'JUN': {'category': 'Shared TF/co-regulator', 'subcategory': 'AP-1'},
    'FOS': {'category': 'Shared TF/co-regulator', 'subcategory': 'AP-1'},
    'SRF': {'category': 'Shared TF/co-regulator', 'subcategory': 'Serum response'},
    'TEAD1': {'category': 'Shared TF/co-regulator', 'subcategory': 'Hippo/YAP'},
    'TEAD4': {'category': 'Shared TF/co-regulator', 'subcategory': 'Hippo/YAP'},
    'YAP1': {'category': 'Shared TF/co-regulator', 'subcategory': 'Hippo pathway'},
    'WWTR1': {'category': 'Shared TF/co-regulator', 'subcategory': 'Hippo/TAZ'},
}

# ============================================================
# 1e. Hormone receptors & endocrine axes
# ============================================================
hormone_axes = {
    'NR3C1': {'category': 'Hormone receptor/axis', 'subcategory': 'Glucocorticoid receptor'},
    'NR3C2': {'category': 'Hormone receptor/axis', 'subcategory': 'Mineralocorticoid receptor'},
    'HSD11B1': {'category': 'Hormone receptor/axis', 'subcategory': 'Cortisol activation'},
    'HSD11B2': {'category': 'Hormone receptor/axis', 'subcategory': 'Cortisol inactivation'},
    'THRA': {'category': 'Hormone receptor/axis', 'subcategory': 'Thyroid receptor α'},
    'THRB': {'category': 'Hormone receptor/axis', 'subcategory': 'Thyroid receptor β'},
    'INSR': {'category': 'Hormone receptor/axis', 'subcategory': 'Insulin receptor'},
    'IGF1R': {'category': 'Hormone receptor/axis', 'subcategory': 'IGF1 receptor'},
    'GHR': {'category': 'Hormone receptor/axis', 'subcategory': 'GH receptor'},
    'IRS1': {'category': 'Hormone receptor/axis', 'subcategory': 'Insulin signaling'},
    'IRS2': {'category': 'Hormone receptor/axis', 'subcategory': 'Insulin signaling'},
    'ADRB1': {'category': 'Hormone receptor/axis', 'subcategory': 'β-adrenergic'},
    'ADRB2': {'category': 'Hormone receptor/axis', 'subcategory': 'β-adrenergic'},
    'ADRA1A': {'category': 'Hormone receptor/axis', 'subcategory': 'α-adrenergic'},
    'ADRA1B': {'category': 'Hormone receptor/axis', 'subcategory': 'α-adrenergic'},
    'CNR1': {'category': 'Hormone receptor/axis', 'subcategory': 'Endocannabinoid'},
    'CNR2': {'category': 'Hormone receptor/axis', 'subcategory': 'Endocannabinoid'},
    'LEPR': {'category': 'Hormone receptor/axis', 'subcategory': 'Leptin receptor'},
    'ADIPOR1': {'category': 'Hormone receptor/axis', 'subcategory': 'Adiponectin receptor'},
    'ADIPOR2': {'category': 'Hormone receptor/axis', 'subcategory': 'Adiponectin receptor'},
    'SOCS1': {'category': 'Hormone receptor/axis', 'subcategory': 'Cytokine/JAK suppressor'},
    'SOCS2': {'category': 'Hormone receptor/axis', 'subcategory': 'Cytokine/JAK suppressor'},
    'SOCS3': {'category': 'Hormone receptor/axis', 'subcategory': 'Cytokine/JAK suppressor'},
    'CISH': {'category': 'Hormone receptor/axis', 'subcategory': 'Cytokine suppressor'},
}

# ============================================================
# 1f. Circulating AA transporters (multi-tissue)
# ============================================================
circulating_transporters = {
    'SLC1A1': {'category': 'AA transporter', 'subcategory': 'Glu/Asp (EAAT3)'},
    'SLC1A2': {'category': 'AA transporter', 'subcategory': 'Glu/Asp (EAAT2)'},
    'SLC1A3': {'category': 'AA transporter', 'subcategory': 'Glu/Asp (EAAT1)'},
    'SLC1A4': {'category': 'AA transporter', 'subcategory': 'Ala/Ser/Cys (ASCT1)'},
    'SLC1A5': {'category': 'AA transporter', 'subcategory': 'Gln (ASCT2)'},
    'SLC7A5': {'category': 'AA transporter', 'subcategory': 'LAT1/BCAA transporter'},
    'SLC3A2': {'category': 'AA transporter', 'subcategory': '4F2hc/LAT co-receptor'},
    'SLC7A1': {'category': 'AA transporter', 'subcategory': 'CAT1/Arg transporter'},
    'SLC7A2': {'category': 'AA transporter', 'subcategory': 'CAT2/Arg transporter'},
    'SLC7A6': {'category': 'AA transporter', 'subcategory': 'y+LAT2'},
    'SLC7A7': {'category': 'AA transporter', 'subcategory': 'y+LAT1'},
    'SLC7A8': {'category': 'AA transporter', 'subcategory': 'LAT2'},
    'SLC7A10': {'category': 'AA transporter', 'subcategory': 'Asc-1'},
    'SLC7A11': {'category': 'AA transporter', 'subcategory': 'xCT/Cys-Glu'},
    'SLC16A10': {'category': 'AA transporter', 'subcategory': 'TAT1/aromatic AA'},
    'SLC6A14': {'category': 'AA transporter', 'subcategory': 'ATB0+/broad specificity'},
    'SLC6A15': {'category': 'AA transporter', 'subcategory': 'B0AT2/BCAA'},
}

# ============================================================
# 1g. Proteases / ECM remodeling (affect inter-tissue signaling)
# ============================================================
ecm_signaling = {
    'MMP2': {'category': 'ECM/Protease', 'subcategory': 'Gelatinase'},
    'MMP9': {'category': 'ECM/Protease', 'subcategory': 'Gelatinase'},
    'MMP14': {'category': 'ECM/Protease', 'subcategory': 'MT1-MMP'},
    'TIMP1': {'category': 'ECM/Protease', 'subcategory': 'TIMP'},
    'TIMP2': {'category': 'ECM/Protease', 'subcategory': 'TIMP'},
    'TIMP3': {'category': 'ECM/Protease', 'subcategory': 'TIMP'},
    'PLAT': {'category': 'ECM/Protease', 'subcategory': 'tPA'},
    'PLAU': {'category': 'ECM/Protease', 'subcategory': 'uPA'},
    'SERPINE1': {'category': 'ECM/Protease', 'subcategory': 'PAI-1'},
    'CTGF': {'category': 'ECM/Protease', 'subcategory': 'CCN2/CTGF'},
    'CYR61': {'category': 'ECM/Protease', 'subcategory': 'CCN1'},
    'TGFB1': {'category': 'ECM/Protease', 'subcategory': 'TGF-β1'},
    'TGFBR1': {'category': 'ECM/Protease', 'subcategory': 'TGF-βR1'},
    'TGFBR2': {'category': 'ECM/Protease', 'subcategory': 'TGF-βR2'},
}

# ============================================================
# Combine all categories
# ============================================================
all_cross_talk = {}
for cat_dict in [myokines, hepatokines, aa_sensors, shared_tfs, hormone_axes,
                  circulating_transporters, ecm_signaling]:
    for gene, info in cat_dict.items():
        if gene not in all_cross_talk:
            all_cross_talk[gene] = info

found_genes = search_genes(all_cross_talk, muscle_genes, liver_genes, muscle_upper, liver_upper)

# ============================================================
# 2. Compute expression stats for each gene in each tissue
# ============================================================
stages = [15, 45, 75, 105]

def compute_stage_stats(df, gene_col='gene_name'):
    """For each gene, compute DLY/TFB mean, log2FC per stage."""
    results = {}
    for gene in df[gene_col].unique():
        gdf = df[df[gene_col] == gene]
        stats = {}
        for s in stages:
            dly = gdf[(gdf['breed']=='DLY') & (gdf['stage_kg']==s)]['expr']
            tfb = gdf[(gdf['breed']=='TFB') & (gdf['stage_kg']==s)]['expr']
            dly_m = dly.mean() if len(dly) > 0 else np.nan
            tfb_m = tfb.mean() if len(tfb) > 0 else np.nan
            if pd.notna(dly_m) and pd.notna(tfb_m) and tfb_m > 0:
                log2fc = np.log2(dly_m / tfb_m)
            else:
                log2fc = np.nan
            stats[s] = {'DLY': dly_m, 'TFB': tfb_m, 'log2FC': log2fc}
        results[gene] = stats
    return results

print("\nComputing expression stats...")
liver_stats = compute_stage_stats(liver_mean)
muscle_stats = compute_stage_stats(muscle_mean)

# ============================================================
# 3. Build output table
# ============================================================
print("Building cross-talk gene table...")
rows = []

for gene, info in found_genes.items():
    in_m = info['in_muscle']
    in_l = info['in_liver']

    if not in_m and not in_l:
        continue  # Gene not found in either tissue

    row = {
        'Gene': gene,
        'Category': info['category'],
        'Subcategory': info['subcategory'],
        'In_Muscle': 'Yes' if in_m else '—',
        'In_Liver': 'Yes' if in_l else '—',
    }

    # Muscle expression
    m_name = info['muscle_name']
    m_stat = muscle_stats.get(m_name, {})
    for s in stages:
        st = m_stat.get(s, {})
        row[f'M_{s}kg_log2FC'] = round(st.get('log2FC', np.nan), 3) if pd.notna(st.get('log2FC', np.nan)) else ''
        row[f'M_{s}kg_DLY'] = round(st.get('DLY', np.nan), 2) if pd.notna(st.get('DLY', np.nan)) else ''
        row[f'M_{s}kg_TFB'] = round(st.get('TFB', np.nan), 2) if pd.notna(st.get('TFB', np.nan)) else ''

    # Liver expression
    l_name = info['liver_name']
    l_stat = liver_stats.get(l_name, {})
    for s in stages:
        st = l_stat.get(s, {})
        row[f'L_{s}kg_log2FC'] = round(st.get('log2FC', np.nan), 3) if pd.notna(st.get('log2FC', np.nan)) else ''
        row[f'L_{s}kg_DLY'] = round(st.get('DLY', np.nan), 2) if pd.notna(st.get('DLY', np.nan)) else ''
        row[f'L_{s}kg_TFB'] = round(st.get('TFB', np.nan), 2) if pd.notna(st.get('TFB', np.nan)) else ''

    # Mean |log2FC| in each tissue
    m_fcs = []
    l_fcs = []
    for s in stages:
        if row[f'M_{s}kg_log2FC'] != '':
            m_fcs.append(abs(float(row[f'M_{s}kg_log2FC'])))
        if row[f'L_{s}kg_log2FC'] != '':
            l_fcs.append(abs(float(row[f'L_{s}kg_log2FC'])))
    row['Muscle_mean_abs_log2FC'] = round(np.mean(m_fcs), 3) if m_fcs else ''
    row['Liver_mean_abs_log2FC'] = round(np.mean(l_fcs), 3) if l_fcs else ''

    # Divergence score: |muscle_log2FC - liver_log2FC| at 105kg
    m105 = row.get('M_105kg_log2FC', '')
    l105 = row.get('L_105kg_log2FC', '')
    if m105 != '' and l105 != '':
        row['105kg_divergence'] = round(abs(float(m105) - float(l105)), 3)
    else:
        row['105kg_divergence'] = ''

    # Direction annotation
    if m105 != '' and l105 != '':
        m_sign = 'DLY↑' if float(m105) > 0.5 else ('TFB↑' if float(m105) < -0.5 else '≈')
        l_sign = 'DLY↑' if float(l105) > 0.5 else ('TFB↑' if float(l105) < -0.5 else '≈')
        if m_sign == l_sign:
            row['Concordance'] = f'Concordant ({m_sign} in both)'
        elif m_sign == '≈' or l_sign == '≈':
            row['Concordance'] = f'Partial ({m_sign}M/{l_sign}L)'
        else:
            row['Concordance'] = f'Discordant (M:{m_sign} L:{l_sign})'
    else:
        row['Concordance'] = ''

    rows.append(row)

df_out = pd.DataFrame(rows)

# Sort by category, then by 105kg_divergence
cat_order = {'Myokine': 1, 'Hepatokine': 2, 'AA/Nutrient sensing': 3,
             'Shared TF/co-regulator': 4, 'Hormone receptor/axis': 5,
             'AA transporter': 6, 'ECM/Protease': 7}
df_out['cat_sort'] = df_out['Category'].map(cat_order).fillna(99)
df_out['div_sort'] = pd.to_numeric(df_out['105kg_divergence'], errors='coerce').fillna(0)
df_out = df_out.sort_values(['cat_sort', 'div_sort'], ascending=[True, False]).drop(columns=['cat_sort', 'div_sort'])

# ============================================================
# 4. Also compute: genes found in BOTH tissues with highest |log2FC| difference
# ============================================================
# This identifies genes that behave very differently between tissues
both_tissue = df_out[(df_out['In_Muscle'] == 'Yes') & (df_out['In_Liver'] == 'Yes')].copy()
both_tissue['div_num'] = pd.to_numeric(both_tissue['105kg_divergence'], errors='coerce')
top_divergent = both_tissue.dropna(subset=['div_num']).nlargest(20, 'div_num')

print(f"\n=== Top 20 genes with most divergent liver vs muscle response at 105kg ===")
for _, r in top_divergent.iterrows():
    mfc = str(r['M_105kg_log2FC'])
    lfc = str(r['L_105kg_log2FC'])
    subcat = str(r['Subcategory'])[:30]
    print(f"  {str(r['Gene']):12s} {str(r['Category']):25s} {subcat:30s} "
          f"M_FC={mfc:>8s} L_FC={lfc:>8s} div={r['105kg_divergence']}")

# ============================================================
# 5. Save
# ============================================================
with pd.ExcelWriter('liver_muscle_crosstalk_genes.xlsx', engine='openpyxl') as writer:
    # Full detailed table
    df_out.to_excel(writer, sheet_name='All_Crosstalk_Genes', index=False)

    # Simplified: just log2FC columns
    slim_cols = ['Gene', 'Category', 'Subcategory', 'In_Muscle', 'In_Liver',
                 'Concordance', '105kg_divergence',
                 'M_15kg_log2FC', 'M_45kg_log2FC', 'M_75kg_log2FC', 'M_105kg_log2FC',
                 'L_15kg_log2FC', 'L_45kg_log2FC', 'L_75kg_log2FC', 'L_105kg_log2FC']
    slim = df_out[[c for c in slim_cols if c in df_out.columns]]
    slim.to_excel(writer, sheet_name='Log2FC_Summary', index=False)

    # Genes in BOTH tissues only
    both_cols = [c for c in slim_cols if c in df_out.columns]
    both_df = df_out[(df_out['In_Muscle']=='Yes') & (df_out['In_Liver']=='Yes')][both_cols]
    both_df.to_excel(writer, sheet_name='Both_Tissues', index=False)

    # Discriminant genes
    disc = df_out[df_out['Concordance'].str.contains('Discordant', na=False)]
    if len(disc) > 0:
        disc.to_excel(writer, sheet_name='Discordant_Patterns', index=False)

    # By category: one sheet per category
    for cat in sorted(df_out['Category'].unique()):
        cat_df = df_out[df_out['Category'] == cat]
        sheet_name = cat.replace('/', '-')[:31]
        cat_df.to_excel(writer, sheet_name=sheet_name, index=False)

# Summary stats
n_total = len(df_out)
n_both = len(df_out[(df_out['In_Muscle']=='Yes') & (df_out['In_Liver']=='Yes')])
n_muscle_only = len(df_out[(df_out['In_Muscle']=='Yes') & (df_out['In_Liver']=='—')])
n_liver_only = len(df_out[(df_out['In_Muscle']=='—') & (df_out['In_Liver']=='Yes')])

print(f"\n{'='*60}")
print(f"Total cross-talk genes: {n_total}")
print(f"  In BOTH tissues: {n_both}")
print(f"  Muscle only: {n_muscle_only}")
print(f"  Liver only: {n_liver_only}")
print(f"\nSaved: liver_muscle_crosstalk_genes.xlsx")
for s in ['All_Crosstalk_Genes', 'Log2FC_Summary', 'Both_Tissues']:
    print(f"  Sheet: {s}")
print(f"  Sheets by category: {', '.join(sorted(df_out['Category'].unique()))}")
