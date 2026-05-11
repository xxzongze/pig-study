#!/usr/bin/env python3
"""
整合外部 GSEA 结果与 pig-study 现有分析

输入:
  - GSEA分析统计表_Gsea_20260510_162059.csv  (肝脏, 外部网站)
  - GSEA分析统计表_Gsea_20260510_162140.csv  (肌肉, 外部网站)
  - key_results_summary.xlsx                (关键AA酶基因)
  - aa_candidate_genes_screened.csv         (AA候选基因, 肝-肌轴)
  - biological_closure_screening_results.xlsx

输出:
  - gsea_integration_summary.xlsx
  - 控制台打印整合报告
"""

import pandas as pd
import numpy as np
import os

# ============================================================
# 1. 加载外部 GSEA 结果
# ============================================================
GSEA_LIVER = 'gsea_external_liver.csv'
GSEA_MUSCLE = 'gsea_external_muscle.csv'

liver_raw = pd.read_csv(GSEA_LIVER)
muscle_raw = pd.read_csv(GSEA_MUSCLE)

liver_raw['Tissue'] = 'Liver'
muscle_raw['Tissue'] = 'Muscle'
gsea = pd.concat([liver_raw, muscle_raw], ignore_index=True)
gsea['Direction'] = gsea['Group'].apply(lambda x: 'TFB-up' if 'TFB' in str(x) else 'DLY-up')
gsea['NES_abs'] = gsea['NES'].abs()
gsea = gsea.sort_values('Padjust')

# ============================================================
# 2. 加载项目关键基因
# ============================================================
key_genes = pd.read_excel('key_results_summary.xlsx')
# 仅保留有基因名的行
key_genes_aa = key_genes[key_genes['Category'] == 'AA Catabolism Enzyme'].copy()

# ============================================================
# 3. KEGG pathway → 关键基因映射 (手工整理, 基于 KEGG 定义)
# ============================================================
PATHWAY_GENE_MAP = {
    'ssc00220': {  # Arginine biosynthesis
        'pathway_name': 'Arginine biosynthesis',
        'genes_in_pathway': ['ARG1', 'ARG2', 'ASL', 'ASS1', 'CPS1', 'OTC', 'NAGS', 'GOT1', 'GOT2', 'GLUD1'],
        'project_hit_genes': [],
    },
    'ssc00270': {  # Cysteine and methionine metabolism
        'pathway_name': 'Cysteine and methionine metabolism',
        'genes_in_pathway': ['SDS', 'CBS', 'CTH', 'CDO1', 'MAT1A', 'MAT2A', 'AHCY', 'GOT1', 'GOT2'],
        'project_hit_genes': [],
    },
    'ssc00250': {  # Alanine, aspartate and glutamate metabolism
        'pathway_name': 'Alanine, aspartate and glutamate metabolism',
        'genes_in_pathway': ['GOT1', 'GOT2', 'GPT', 'GPT2', 'GLUD1', 'GLUL', 'ASNS', 'ASL', 'ASS1'],
        'project_hit_genes': [],
    },
    'ssc00280': {  # Valine, leucine and isoleucine degradation
        'pathway_name': 'BCAA degradation',
        'genes_in_pathway': ['BCAT1', 'BCAT2', 'BCKDHA', 'BCKDHB', 'DBT', 'DLD', 'ACADSB', 'HIBCH', 'HMGCL'],
        'project_hit_genes': [],
    },
    'ssc00350': {  # Tyrosine metabolism
        'pathway_name': 'Tyrosine metabolism',
        'genes_in_pathway': ['HGD', 'TAT', 'FAH', 'PAH', 'HPD', 'GSTZ1', 'MAOA', 'MAOB'],
        'project_hit_genes': [],
    },
    'ssc00380': {  # Tryptophan metabolism
        'pathway_name': 'Tryptophan metabolism',
        'genes_in_pathway': ['TDO2', 'IDO1', 'IDO2', 'KYNU', 'HAAO', 'AFMID', 'AANAT', 'ASMT'],
        'project_hit_genes': [],
    },
    'ssc00400': {  # Phe, Tyr, Trp biosynthesis
        'pathway_name': 'Phe/Tyr/Trp biosynthesis',
        'genes_in_pathway': [],
        'project_hit_genes': [],
    },
    'ssc00020': {  # TCA cycle
        'pathway_name': 'TCA cycle',
        'genes_in_pathway': ['CS', 'ACO1', 'ACO2', 'IDH1', 'IDH2', 'IDH3A', 'OGDH', 'SUCLG1', 'SDHA', 'FH', 'MDH1'],
        'project_hit_genes': [],
    },
    'ssc00640': {  # Propanoate metabolism
        'pathway_name': 'Propanoate metabolism',
        'genes_in_pathway': ['ACSS1', 'ACSS2', 'ABAT', 'ALDH6A1', 'MCEE', 'MUT', 'PCCA', 'PCCB'],
        'project_hit_genes': [],
    },
}

# 匹配项目关键基因到通路
for k, v in PATHWAY_GENE_MAP.items():
    pathway_genes_upper = [g.upper() for g in v['genes_in_pathway']]
    hits = key_genes_aa[key_genes_aa['Gene'].str.upper().isin(pathway_genes_upper)]
    v['project_hit_genes'] = hits['Gene'].tolist() if len(hits) > 0 else []

# ============================================================
# 4. 构建整合表
# ============================================================
rows = []
for _, row in gsea.iterrows():
    gs_id = str(row['Gene Set Name'])
    pathway_info = PATHWAY_GENE_MAP.get(gs_id, {})
    project_hits = pathway_info.get('project_hit_genes', [])
    is_sig = row['Padjust'] < 0.05

    # 匹配 key_results_summary 中该通路相关基因的详细信息
    hit_details = ''
    if project_hits:
        for gene in project_hits:
            match = key_genes_aa[key_genes_aa['Gene'].str.upper() == gene.upper()]
            if len(match) > 0:
                m = match.iloc[0]
                hit_details += f"{gene}: log2FC@45={m['log2FC_45kg']:.1f}, r_SerumUrea={m['r_vs_SerumUrea']:.3f}, r_STAT3={m['r_vs_STAT3']:.3f}; "
        hit_details = hit_details.rstrip('; ')

    rows.append({
        'Tissue': row['Tissue'],
        'Pathway_ID': gs_id,
        'Description': row['Description'],
        'Direction': row['Direction'],
        'NES': row['NES'],
        'P_value': row['Pvalue'],
        'FDR': row['Padjust'],
        'FDR_significant': is_sig,
        'N_project_hit_genes': len(project_hits),
        'Project_hit_genes': ', '.join(project_hits) if project_hits else '',
        'Hit_gene_details': hit_details,
    })

integration = pd.DataFrame(rows)
integration['NES_abs'] = integration['NES'].abs()
integration = integration.sort_values(['FDR_significant', 'FDR', 'NES_abs'],
                                      ascending=[False, True, False])

# ============================================================
# 5. 保存
# ============================================================
OUTPUT = 'gsea_integration_summary.xlsx'
with pd.ExcelWriter(OUTPUT) as writer:
    # Sheet 1: 按 FDR 排序的完整整合表
    integration.to_excel(writer, sheet_name='Integrated_GSEA', index=False)

    # Sheet 2: 仅显著通路
    sig = integration[integration['FDR_significant'] == True]
    sig.to_excel(writer, sheet_name='FDR_significant_only', index=False)

    # Sheet 3: 关键 AA 酶基因汇总
    key_genes_aa_out = key_genes_aa[['Gene', 'Category', 'Direction', 'log2FC_45kg',
                                      'r_vs_SerumUrea', 'p_Urea', 'r_vs_STAT3', 'p_STAT3']].copy()
    key_genes_aa_out.to_excel(writer, sheet_name='Key_AA_enzymes', index=False)

print(f'输出: {OUTPUT}')
print(f'  总通路: {len(integration)}')
print(f'  FDR<0.05: {sig["Tissue"].count():.0f} ({sig[sig["Tissue"]=="Liver"].shape[0]} liver, {sig[sig["Tissue"]=="Muscle"].shape[0]} muscle)')
print(f'  有项目基因命中的通路: {(integration["N_project_hit_genes"] > 0).sum()}')

# ============================================================
# 6. 控制台报告
# ============================================================
print('\n' + '=' * 90)
print('整合报告: 外部 GSEA × pig-study 关键基因')
print('=' * 90)

print('\n--- 肝脏: TFB 上调的 AA 分解通路 (FDR<0.05) ---')
liver_aa = integration[(integration['Tissue'] == 'Liver') &
                        (integration['FDR_significant']) &
                        (integration['Direction'] == 'TFB-up')]
for _, r in liver_aa.iterrows():
    marker = ' [有项目基因]' if r['N_project_hit_genes'] > 0 else ''
    print(f"  {r['Description']:45s} NES={r['NES']:+6.3f}  FDR={r['FDR']:.4f}{marker}")
    if r['Hit_gene_details']:
        print(f"    → 项目已验证基因: {r['Hit_gene_details']}")

print('\n--- 肝脏: DLY 上调的信号/疾病通路 (FDR<0.05) ---')
liver_dly = integration[(integration['Tissue'] == 'Liver') &
                         (integration['FDR_significant']) &
                         (integration['Direction'] == 'DLY-up')]
for _, r in liver_dly.iterrows():
    print(f"  {r['Description']:45s} NES={r['NES']:+6.3f}  FDR={r['FDR']:.4f}")

print('\n--- 肌肉: 显著通路 (FDR<0.05) ---')
muscle_sig_df = integration[(integration['Tissue'] == 'Muscle') & (integration['FDR_significant'])]
for _, r in muscle_sig_df.iterrows():
    print(f"  {r['Description']:45s} NES={r['NES']:+6.3f}  FDR={r['FDR']:.4f}  {r['Direction']}")

print('\n--- 关键 AA 酶基因: 跨阶段表达 + 生化关联 ---')
for _, r in key_genes_aa.iterrows():
    urea_str = f'r_Urea={r["r_vs_SerumUrea"]:.3f}' if pd.notna(r['r_vs_SerumUrea']) else ''
    stat3_str = f'r_STAT3={r["r_vs_STAT3"]:.3f}' if pd.notna(r['r_vs_STAT3']) else ''
    print(f"  {r['Gene']:8s} {r['Direction'][:15] if pd.notna(r['Direction']) else '':15s} log2FC@45={float(r['log2FC_45kg']):+6.2f}  {urea_str:20s} {stat3_str}")

print('\n--- 氮平衡: DLY vs TFB @ 4 阶段 ---')
for _, r in key_genes[key_genes['Category'] == 'N Balance'].iterrows():
    print(f"  {r['Gene']}")

print('\n--- 生物学闭合性: 整合逻辑 ---')
print('''
  外部GSEA → 独立于项目内部筛选 → 作为外部验证证据
  TFB肝脏AA分解↑ (FDR<0.05)  ⇔  项目关键基因: SDS/GOT1/HGD/ARG1/ARG2/ASL 均为TFB↑
  BCAA降解 TFB↑ l(+肌)       ⇔  BCAT1 是唯一DLY↑的AA酶 (与整体AA分解方向相反)
  TCA cycle TFB↑              ⇔  AA碳骨架进入TCA的下游
  STAT3 (r_Urea=0.745)        ⇔  可能是AA分解的转录调控因子
  DLY 氮保留更高               ⇔  与TFB AA分解通路↑ 一致: 分解多 → 氮浪费多 → 尿素高
  DLY↑ 信号通路               ⇔  可能与更高效的蛋白质合成/细胞生长有关
''')
