#!/usr/bin/env python3
"""
花生四烯酸(AA)代谢——肝肌轴跨组织分析
Step 1: 肌肉AA受体/响应基因表达模式
Step 2: 肝脏AA酶 vs 肌肉AA受体 跨组织相关性
Step 3: 识别AA代谢介导肝肌crosstalk的候选信号轴
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from stats_utils import benjamini_hochberg
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 1. 数据加载
# ============================================================

# 肌肉基因矩阵
muscle = pd.read_csv('/Users/hezongze/pig_study/gene_expression/muscle_gene_matrix.xls',
                      sep=None, engine='python')

# 肝脏表达(已校正105kg)
liver = pd.read_csv('/Users/hezongze/pig_study/liver_expression_corrected.csv')

# 花生四烯酸通路基因(肝脏)
aa_liver = pd.read_csv('/Users/hezongze/Downloads/liver_arachidonic_pathway_genes.csv')

print("Muscle matrix shape:", muscle.shape)
print("Liver matrix shape:", liver.shape)
print("AA liver genes:", aa_liver.shape[0])

# ============================================================
# 2. 样本分组映射
# ============================================================

# 肌肉列名映射 (Sample Initial Name pattern)
muscle_groups = {}
for col in muscle.columns:
    if col == 'seq_id' or col == 'gene_name' or col == 'length' or col == 'description':
        continue
    # e.g., m_15_1_1 -> muscle, 15kg, breed1(DLY), rep1
    parts = col.split('_')
    if parts[0] == 'm':
        if parts[1] == '15':
            stage = '15'
            breed = 'DLY' if parts[2] == '1' else 'TFB'
        elif parts[1] == '1':
            stage = '75'
            breed = 'DLY' if parts[2] == '1' else 'TFB'
        elif parts[1] == '2':
            stage = '105'
            breed = 'DLY' if parts[2] == '1' else 'TFB'
        elif parts[1] == '3':
            stage = '135'
            breed = 'DLY' if parts[2] == '1' else 'TFB'
        else:
            continue
    elif parts[0] == 'BJ':
        if parts[1] == '2':
            stage = '45'
            breed = 'DLY' if parts[2] == '1' else 'TFB'
        else:
            continue
    else:
        continue
    muscle_groups[col] = f"{'DLY' if breed == 'DLY' else 'TFB'}M{stage}"

print(f"\nMuscle group mapping: {len(muscle_groups)} samples")
for g in sorted(set(muscle_groups.values())):
    n = sum(1 for v in muscle_groups.values() if v == g)
    print(f"  {g}: {n} replicates")

# 肝脏列名映射 (from corrected file)
liver_groups = {}
for col in liver.columns:
    if col == 'seq_id':
        continue
    parts = col.split('_')
    if parts[0] != 'L':
        continue
    stage_code = parts[1]
    breed_code = parts[2]

    if stage_code == '15':
        stage = '15'
    elif stage_code == '45':
        stage = '45'
    elif stage_code == '1':
        stage = '75'
    elif stage_code == '2':
        stage = '105'
    elif stage_code == '3':
        stage = '135'
    else:
        continue

    breed = 'DLY' if breed_code == '1' else 'TFB'
    liver_groups[col] = f"{breed}L{stage}"

print(f"\nLiver group mapping: {len(liver_groups)} samples")
for g in sorted(set(liver_groups.values())):
    n = sum(1 for v in liver_groups.values() if v == g)
    print(f"  {g}: {n} replicates")

# ============================================================
# 3. 定义基因集
# ============================================================

# 肌肉AA受体/响应基因（猪同源）
muscle_aa_receptors = {
    # 前列腺素受体
    'PTGER1': 'ENSSSCG00000016014',
    'PTGER2': 'ENSSSCG00000021862',
    'PTGER3': 'ENSSSCG00000003788',
    'PTGER4': 'ENSSSCG00000024439',
    'PTGIR': 'ENSSSCG00000026602',   # PGI2 receptor
    'PTGDR': 'ENSSSCG00000057054',   # PGD2 receptor
    'PTGDR2': 'ENSSSCG00000013106',
    'TBXA2R': 'ENSSSCG00000033759',  # Thromboxane receptor
    # 白三烯受体
    'CYSLTR1': 'ENSSSCG00000039203',
    'CYSLTR2': 'ENSSSCG00000009399',
    'LTB4R': 'ENSSSCG00000020941',   # LTB4R2 in file
    # PPAR家族（脂肪酸/类二十烷酸核受体）
    'PPARA': 'ENSSSCG00000010316',
    'PPARD': 'ENSSSCG00000015482',
    'PPARG': 'ENSSSCG00000005564',
    'PPARGC1A': 'ENSSSCG00000006575',
    # 脂肪酸转运/结合
    'CD36': 'ENSSSCG00000016347',
    'FABP3': 'ENSSSCG00000016497',   # muscle specific
    'FABP4': 'ENSSSCG00000011624',
    # AA信号相关的蛋白合成/降解调控
    'MTOR': 'ENSSSCG00000016169',
    'RPS6KB1': 'ENSSSCG00000012461',
    'FOXO1': 'ENSSSCG00000001779',
    'FOXO3': 'ENSSSCG00000009668',
    'FBXO32': 'ENSSSCG00000017676',  # Atrogin-1
    'TRIM63': 'ENSSSCG00000008312',  # MuRF1
    'MYOD1': 'ENSSSCG00000011878',
    'MYOG': 'ENSSSCG00000006763',
    # 已知的肝肌轴myokine/hepatokine受体
    'IL6R': 'ENSSSCG00000013829',
    'IGF1R': 'ENSSSCG00000029928',
    'INSR': 'ENSSSCG00000013683',
}

# 肝脏AA代谢酶（从已有文件中提取关键类别代表）
liver_aa_enzymes = {}
for _, row in aa_liver.iterrows():
    liver_aa_enzymes[row['Gene Name']] = {
        'gene_id': row['Gene ID'],
        'category': row['Category'],
        'pattern': row['Pattern'],
    }

# ============================================================
# 4. 计算每组的平均表达 (Z-score标准化)
# ============================================================

def calc_group_zscore(matrix, group_map, gene_ids_of_interest=None):
    """计算每个基因在每个组的平均表达，然后在基因内Z-score标准化"""
    # 提取表达值
    sample_cols = [c for c in matrix.columns if c in group_map]

    # 如果有gene_name列，用它建立索引
    if 'gene_name' in matrix.columns:
        expr = matrix.set_index('gene_name')[sample_cols]
    elif 'seq_id' in matrix.columns:
        expr = matrix.set_index('seq_id')[sample_cols]
    else:
        expr = matrix[sample_cols]

    # 按组取平均
    groups = sorted(set(group_map.values()))
    group_means = pd.DataFrame(index=expr.index)
    for g in groups:
        cols_in_g = [c for c in sample_cols if group_map[c] == g]
        if cols_in_g:
            group_means[g] = expr[cols_in_g].mean(axis=1)

    # 基因内Z-score
    zscore = group_means.subtract(group_means.mean(axis=1), axis=0).divide(
        group_means.std(axis=1), axis=0)

    return group_means, zscore

muscle_mean, muscle_z = calc_group_zscore(muscle, muscle_groups)
liver_mean, liver_z = calc_group_zscore(liver, liver_groups)

print(f"\nMuscle groups: {list(muscle_z.columns)}")
print(f"Liver groups: {list(liver_z.columns)}")

# ============================================================
# 5. 提取肌肉AA受体基因的表达
# ============================================================

# 肌肉基因名称可能在gene_name列，我们需要搜索
muscle_gene_col = 'gene_name' if 'gene_name' in muscle.columns else None

# 构建基因名->索引的映射
if muscle_gene_col:
    muscle_gene_to_idx = dict(zip(muscle[muscle_gene_col].str.strip().str.upper(),
                                   range(len(muscle))))
else:
    # 用seq_id
    muscle_gene_to_idx = dict(zip(muscle['seq_id'].str.strip().str.upper(),
                                   range(len(muscle))))

# 找到肌肉AA受体基因的实际Ensembl ID
# 先看看基因名称列的格式
if muscle_gene_col:
    print(f"\nMuscle gene_name examples: {muscle[muscle_gene_col].dropna().head(10).tolist()}")
print(f"Muscle seq_id examples: {muscle['seq_id'].head(5).tolist()}")

# 用seq_id搜索（Ensembl ID）
muscle_aa_found = {}
for gene_symbol, ensembl_id in muscle_aa_receptors.items():
    # 在muscle矩阵中搜索这个Ensembl ID
    mask = muscle['seq_id'].str.strip().str.upper() == ensembl_id.upper()
    if mask.any():
        idx = mask.values.argmax() if hasattr(mask, 'values') else mask.idxmax()
        muscle_aa_found[gene_symbol] = {
            'ensembl_id': ensembl_id,
            'row_idx': idx,
        }

print(f"\nFound {len(muscle_aa_found)}/{len(muscle_aa_receptors)} muscle AA receptor genes")

# ============================================================
# 6. 肌肉AA受体表达热图数据准备
# ============================================================

muscle_receptor_zscore = {}
muscle_receptor_mean = {}
for gene_symbol, info in muscle_aa_found.items():
    idx = info['row_idx']
    ensembl_id = info['ensembl_id']
    if ensembl_id in muscle_z.index:
        muscle_receptor_zscore[gene_symbol] = muscle_z.loc[ensembl_id]
        muscle_receptor_mean[gene_symbol] = muscle_mean.loc[ensembl_id]
    elif gene_symbol in muscle_z.index:
        muscle_receptor_zscore[gene_symbol] = muscle_z.loc[gene_symbol]
        muscle_receptor_mean[gene_symbol] = muscle_mean.loc[gene_symbol]

if muscle_receptor_zscore:
    mr_df = pd.DataFrame(muscle_receptor_zscore).T
    # 只保留15/45/75/105 kg
    target_stages = [c for c in mr_df.columns if '135' not in c and 'L' not in c]
    mr_df = mr_df[[c for c in mr_df.columns if c in target_stages or '135' not in c]]
    print(f"\nMuscle receptor Z-score matrix: {mr_df.shape}")
    print(mr_df.round(2).to_string())

# ============================================================
# 7. 肝脏AA酶表达数据提取
# ============================================================

# 从liver_z中提取AA酶基因
liver_aa_zscore = {}
liver_aa_mean = {}
for gene_symbol, info in liver_aa_enzymes.items():
    ensembl_id = info['gene_id']
    if ensembl_id in liver_z.index:
        liver_aa_zscore[gene_symbol] = liver_z.loc[ensembl_id]
        liver_aa_mean[gene_symbol] = liver_mean.loc[ensembl_id]

print(f"\nFound {len(liver_aa_zscore)} liver AA enzyme genes in expression matrix")

# ============================================================
# 8. 跨组织相关性分析：肝脏AA酶 vs 肌肉AA受体
# ============================================================

# 对齐组织样本：找到共同的时间点
# 肝脏: DLYL15, DLYL45, DLYL75, DLYL105, TFBL15, TFBL45, TFBL75, TFBL105
# 肌肉: DLYM15, DLYM45, DLYM75, DLYM105, TFBM15, TFBM45, TFBM75, TFBM105

def make_pair_name(tissue, breed, stage):
    return f"{breed}{tissue}{stage}"

correlation_results = []

for l_gene, l_info in liver_aa_enzymes.items():
    l_ensembl = l_info['gene_id']
    if l_ensembl not in liver_mean.index:
        continue

    for m_gene, m_ensembl in muscle_aa_receptors.items():
        if m_ensembl not in muscle_mean.index:
            # Try gene name
            continue

        # Collect paired expression values across all 8 breed-stage combos
        l_vals = []
        m_vals = []
        labels = []

        for breed in ['DLY', 'TFB']:
            for stage in ['15', '45', '75', '105']:
                l_key = f"{breed}L{stage}"
                m_key = f"{breed}M{stage}"
                if l_key in liver_mean.columns and m_key in muscle_mean.columns:
                    l_vals.append(liver_mean.loc[l_ensembl, l_key])
                    m_vals.append(muscle_mean.loc[m_ensembl, m_key])
                    labels.append(f"{breed}{stage}")

        if len(l_vals) >= 6:
            r, p = pearsonr(l_vals, m_vals)
            correlation_results.append({
                'Liver_Gene': l_gene,
                'Liver_Category': l_info['category'],
                'Liver_Pattern': l_info['pattern'],
                'Muscle_Gene': m_gene,
                'Pearson_r': r,
                'P_value': p,
                'N_pairs': len(l_vals),
                'Mean_Liver_Expr': np.mean(l_vals),
                'Mean_Muscle_Expr': np.mean(m_vals),
            })

corr_df = pd.DataFrame(correlation_results)
corr_df = corr_df.sort_values('P_value')

# FDR correction
if len(corr_df) > 0:
    _, fdr_q = benjamini_hochberg(corr_df['P_value'].values)
    corr_df['Q_value'] = fdr_q
    corr_df['FDR_significant'] = corr_df['Q_value'] < 0.05

n_nom = (corr_df['P_value'] < 0.05).sum()
n_fdr = (corr_df['FDR_significant']).sum()
print(f"\n=== 跨组织相关性分析: {len(corr_df)} pairs ===")
print(f"Nominal P<0.05: {n_nom} → FDR<0.05: {n_fdr}")

# FDR-significant
sig_corr = corr_df[corr_df['FDR_significant']]
print(f"\nFDR显著 (q<0.05): {len(sig_corr)} pairs")
if len(sig_corr) > 0:
    print(sig_corr[['Liver_Gene', 'Liver_Category', 'Muscle_Gene', 'Pearson_r', 'P_value', 'Q_value']].to_string())

# 强相关 (|r| > 0.7, FDR < 0.05)
strong_corr = sig_corr[sig_corr['Pearson_r'].abs() > 0.7]
print(f"\n强相关 (|r|>0.7, FDR<0.05): {len(strong_corr)} pairs")
if len(strong_corr) > 0:
    print(strong_corr[['Liver_Gene', 'Liver_Category', 'Muscle_Gene', 'Pearson_r', 'Q_value', 'N_pairs']].to_string())

# ============================================================
# 9. 品种特异性分析: DLY内 vs TFB内
# ============================================================

print("\n=== 品种内特异性相关 ===")
for breed in ['DLY', 'TFB']:
    breed_corr = []
    for l_gene, l_info in liver_aa_enzymes.items():
        l_ensembl = l_info['gene_id']
        if l_ensembl not in liver_mean.index:
            continue
        for m_gene, m_ensembl in muscle_aa_receptors.items():
            if m_ensembl not in muscle_mean.index:
                continue
            l_vals = []
            m_vals = []
            for stage in ['15', '45', '75', '105']:
                l_key = f"{breed}L{stage}"
                m_key = f"{breed}M{stage}"
                if l_key in liver_mean.columns and m_key in muscle_mean.columns:
                    l_vals.append(liver_mean.loc[l_ensembl, l_key])
                    m_vals.append(muscle_mean.loc[m_ensembl, m_key])
            if len(l_vals) >= 3:
                r, p = pearsonr(l_vals, m_vals)
                breed_corr.append({
                    'Liver_Gene': l_gene,
                    'Liver_Category': l_info['category'],
                    'Muscle_Gene': m_gene,
                    'Pearson_r': r,
                    'P_value': p,
                })

    bdf = pd.DataFrame(breed_corr).sort_values('P_value')
    if len(bdf) > 0:
        _, fdr_q = benjamini_hochberg(bdf['P_value'].values)
        bdf['Q_value'] = fdr_q
        bdf['FDR_significant'] = bdf['Q_value'] < 0.05
    bdf_sig = bdf[bdf['FDR_significant']]
    bdf_nom = (bdf['P_value'] < 0.05).sum()
    print(f"\n{breed}: nominal P<0.05={bdf_nom}, FDR<0.05={len(bdf_sig)} significant pairs")
    if len(bdf_sig) > 0:
        top = bdf_sig.head(10)
        print(top[['Liver_Gene', 'Muscle_Gene', 'Pearson_r', 'P_value', 'Q_value']].to_string())

# ============================================================
# 10. 保存结果
# ============================================================

corr_df.to_excel('/Users/hezongze/pig_study/aa_crosstissue_correlation.xlsx', index=False)
print("\nResults saved to: aa_crosstissue_correlation.xlsx")

# 输出显著相关对的详细数据用于绘图
if len(sig_corr) > 0:
    # 为每个显著对输出8点数据
    detail_rows = []
    for _, row in sig_corr.iterrows():
        l_gene = row['Liver_Gene']
        m_gene = row['Muscle_Gene']
        l_ensembl = liver_aa_enzymes[l_gene]['gene_id']
        m_ensembl = muscle_aa_receptors[m_gene]

        for breed in ['DLY', 'TFB']:
            for stage in ['15', '45', '75', '105']:
                l_key = f"{breed}L{stage}"
                m_key = f"{breed}M{stage}"
                if l_key in liver_mean.columns and m_key in muscle_mean.columns:
                    detail_rows.append({
                        'Liver_Gene': l_gene,
                        'Muscle_Gene': m_gene,
                        'Breed': breed,
                        'Stage': stage,
                        'Liver_Expr': liver_mean.loc[l_ensembl, l_key],
                        'Muscle_Expr': muscle_mean.loc[m_ensembl, m_key],
                    })

    detail_df = pd.DataFrame(detail_rows)
    detail_df.to_csv('/Users/hezongze/pig_study/aa_crosstalk_detail.csv', index=False)
    print("Detail data saved to: aa_crosstalk_detail.csv")

print("\nDone!")
