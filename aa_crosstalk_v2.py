#!/usr/bin/env python3
"""
花生四烯酸(AA)代谢——肝肌轴跨组织分析 (v2)
使用完整的39049基因TPM矩阵
"""

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from stats_utils import benjamini_hochberg
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# 1. 加载完整TPM矩阵
# ============================================================
tpm = pd.read_excel('/Users/hezongze/Downloads/gene.tpm.matrix.xlsx')
tpm = tpm.set_index('seq_id')
print(f"Full TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")

# ============================================================
# 2. 样本分组
# ============================================================

# 肝脏样本
liver_map = {}
# 肌肉样本
muscle_map = {}

for col in tpm.columns:
    parts = col.split('_')
    if parts[0] == 'L':
        # Liver: L_STAGE_BREED_REP
        stage_code = parts[1]
        breed_code = parts[2]
        if stage_code == '15': stage = '15'
        elif stage_code == '45': stage = '45'
        elif stage_code == '1': stage = '75'
        elif stage_code == '2': stage = '105'
        elif stage_code == '3': stage = '135'
        else: continue
        breed = 'DLY' if breed_code == '1' else 'TFB'
        liver_map[col] = f"{breed}_L_{stage}"
    elif parts[0] == 'm':
        # Muscle: m_STAGE_BREED_REP
        stage_code = parts[1]
        breed_code = parts[2]
        if stage_code == '15': stage = '15'
        elif stage_code == '1': stage = '75'
        elif stage_code == '2': stage = '105'
        elif stage_code == '3': stage = '135'
        else: continue
        breed = 'DLY' if breed_code == '1' else 'TFB'
        muscle_map[col] = f"{breed}_M_{stage}"
    elif parts[0] in ['DLYM', 'TFBM']:
        # DLYM_1-6 = DLY muscle 45kg, TFBM_1-6 = TFB muscle 45kg
        breed = 'DLY' if parts[0] == 'DLYM' else 'TFB'
        muscle_map[col] = f"{breed}_M_45"

print(f"\nLiver groups: {len(set(liver_map.values()))} | Muscle groups: {len(set(muscle_map.values()))}")
for g in sorted(set(liver_map.values())):
    print(f"  Liver {g}: {sum(1 for v in liver_map.values() if v==g)} reps")
for g in sorted(set(muscle_map.values())):
    print(f"  Muscle {g}: {sum(1 for v in muscle_map.values() if v==g)} reps")

# ============================================================
# 3. 计算组平均表达并Z-score标准化
# ============================================================
def group_zscore(tpm, group_map):
    groups = sorted(set(group_map.values()))
    means = pd.DataFrame(index=tpm.index)
    for g in groups:
        cols = [c for c in group_map if group_map[c] == g]
        if cols:
            means[g] = tpm[cols].mean(axis=1)
    z = means.subtract(means.mean(axis=1), axis=0).divide(means.std(axis=1).replace(0, np.nan), axis=0)
    return means, z

liver_mean, liver_z = group_zscore(tpm, liver_map)
muscle_mean, muscle_z = group_zscore(tpm, muscle_map)

# ============================================================
# 4. 加载肝脏AA通路基因列表
# ============================================================
aa_liver = pd.read_csv('/Users/hezongze/Downloads/liver_arachidonic_pathway_genes.csv')
print(f"\nAA pathway genes loaded: {len(aa_liver)}")

# ============================================================
# 5. 肌肉AA受体/响应基因（使用Ensembl ID直接搜索）
# ============================================================
muscle_aa_genes = {
    # 前列腺素受体 (细胞膜GPCR)
    'PTGER1': 'ENSSSCG00000016014',
    'PTGER2': 'ENSSSCG00000021862',
    'PTGER3': 'ENSSSCG00000003788',
    'PTGER4': 'ENSSSCG00000024439',
    'PTGIR': 'ENSSSCG00000026602',    # PGI2/IP receptor
    'PTGDR': 'ENSSSCG00000057054',    # PGD2/DP receptor
    'TBXA2R': 'ENSSSCG00000033759',   # Thromboxane/TP receptor
    # 白三烯受体
    'CYSLTR1': 'ENSSSCG00000039203',
    'CYSLTR2': 'ENSSSCG00000009399',
    'LTB4R': 'ENSSSCG00000020941',
    # 核受体 (AA/类二十烷酸/EET感应)
    'PPARA': 'ENSSSCG00000010316',
    'PPARD': 'ENSSSCG00000015482',
    'PPARG': 'ENSSSCG00000005564',
    'PPARGC1A': 'ENSSSCG00000006575',
    # 脂肪酸转运/膜受体
    'CD36': 'ENSSSCG00000016347',
    'FABP3': 'ENSSSCG00000016497',    # 肌肉型FABP
    # AA信号下游效应器 (蛋白合成/降解)
    'MTOR': 'ENSSSCG00000016169',
    'RPS6KB1': 'ENSSSCG00000012461',  # S6K1
    'FOXO1': 'ENSSSCG00000001779',
    'FOXO3': 'ENSSSCG00000009668',
    'FBXO32': 'ENSSSCG00000017676',   # Atrogin-1
    'TRIM63': 'ENSSSCG00000008312',   # MuRF1
    'MYOD1': 'ENSSSCG00000011878',
    'MYOG': 'ENSSSCG00000006763',
    # 肝肌轴相关受体
    'IL6R': 'ENSSSCG00000013829',
    'IGF1R': 'ENSSSCG00000029928',
}

# 验证哪些基因在矩阵中
valid_muscle = {}
for symbol, eid in muscle_aa_genes.items():
    if eid in tpm.index:
        valid_muscle[symbol] = eid
print(f"Muscle AA genes found in TPM: {len(valid_muscle)}/{len(muscle_aa_genes)}")

# ============================================================
# 6. 提取表达数据
# ============================================================

# 肝脏AA酶
liver_aa_data = {}
for _, row in aa_liver.iterrows():
    gid = row['Gene ID']
    if gid in liver_z.index:
        liver_aa_data[row['Gene Name']] = {
            'z': liver_z.loc[gid],
            'mean': liver_mean.loc[gid],
            'category': row['Category'],
            'pattern': row['Pattern'],
        }

# 肌肉AA受体
muscle_aa_data = {}
for symbol, eid in valid_muscle.items():
    if eid in muscle_z.index:
        muscle_aa_data[symbol] = {
            'z': muscle_z.loc[eid],
            'mean': muscle_mean.loc[eid],
        }

print(f"Liver AA enzymes with data: {len(liver_aa_data)}")
print(f"Muscle AA receptors with data: {len(muscle_aa_data)}")

# ============================================================
# 7. 肌肉AA受体表达热图
# ============================================================
print("\n" + "="*60)
print("肌肉AA受体/响应基因 Z-score 表达矩阵")
print("="*60)

mz_df = pd.DataFrame({k: v['z'] for k, v in muscle_aa_data.items()}).T
common_stages = ['15', '45', '75', '105']
# 排序: DLY各阶段 + TFB各阶段
col_order = []
for breed in ['DLY', 'TFB']:
    for s in common_stages:
        col = f"{breed}_M_{s}"
        if col in mz_df.columns:
            col_order.append(col)
mz_df = mz_df[col_order]
print(mz_df.round(2).to_string())

# ============================================================
# 8. 跨组织相关性分析
# ============================================================
print("\n" + "="*60)
print("肝脏AA酶 vs 肌肉AA受体 跨组织Pearson相关")
print("="*60)

corr_results = []
for l_gene, l_info in liver_aa_data.items():
    for m_gene, m_info in muscle_aa_data.items():
        l_vals, m_vals, labels = [], [], []
        for breed in ['DLY', 'TFB']:
            for s in common_stages:
                lk = f"{breed}_L_{s}"
                mk = f"{breed}_M_{s}"
                if lk in l_info['mean'].index and mk in m_info['mean'].index:
                    l_vals.append(l_info['mean'][lk])
                    m_vals.append(m_info['mean'][mk])
                    labels.append(f"{breed}{s}")

        if len(l_vals) >= 6:
            r, p = pearsonr(l_vals, m_vals)
            corr_results.append({
                'Liver_Gene': l_gene,
                'Liver_Category': l_info['category'],
                'Liver_Pattern': l_info['pattern'],
                'Muscle_Gene': m_gene,
                'Pearson_r': round(r, 4),
                'P_value': round(p, 4),
                'N_pairs': len(l_vals),
                'abs_r': abs(r),
            })

corr_df = pd.DataFrame(corr_results).sort_values('P_value')

# FDR correction
if len(corr_df) > 0:
    _, fdr_q = benjamini_hochberg(corr_df['P_value'].values)
    corr_df['Q_value'] = fdr_q
    corr_df['FDR_significant'] = corr_df['Q_value'] < 0.05

# FDR-filtered significant
sig = corr_df[corr_df['FDR_significant']].copy()
n_nom = (corr_df['P_value'] < 0.05).sum()
print(f"\n显著相关 (nominal P<0.05): {n_nom} → (FDR<0.05): {len(sig)}/{len(corr_df)} gene pairs")

# 强正相关
strong_pos = sig[sig['Pearson_r'] > 0.7].sort_values('Pearson_r', ascending=False)
print(f"\n强正相关 (r>0.7, FDR<0.05): {len(strong_pos)} pairs")
if len(strong_pos) > 0:
    print(strong_pos[['Liver_Gene', 'Liver_Category', 'Muscle_Gene', 'Pearson_r', 'P_value', 'Q_value']].to_string())

# 强负相关
strong_neg = sig[sig['Pearson_r'] < -0.7].sort_values('Pearson_r')
print(f"\n强负相关 (r<-0.7, FDR<0.05): {len(strong_neg)} pairs")
if len(strong_neg) > 0:
    print(strong_neg[['Liver_Gene', 'Liver_Category', 'Muscle_Gene', 'Pearson_r', 'P_value', 'Q_value']].to_string())

# ============================================================
# 9. 关键发现：分类分析
# ============================================================
print("\n" + "="*60)
print("按肝脏AA通路类别分组的跨组织相关性")
print("="*60)

for cat in aa_liver['Category'].unique():
    cat_genes = [g for g, info in liver_aa_data.items() if info['category'] == cat]
    if not cat_genes:
        continue
    cat_corr = sig[sig['Liver_Gene'].isin(cat_genes)]
    if len(cat_corr) > 0:
        mean_r = cat_corr['Pearson_r'].mean()
        print(f"\n{cat} ({len(cat_genes)} genes):")
        print(f"  Significant pairs: {len(cat_corr)}, Mean r = {mean_r:.3f}")
        top = cat_corr.sort_values('abs_r', ascending=False).head(5)
        for _, row in top.iterrows():
            print(f"    {row['Liver_Gene']} ↔ {row['Muscle_Gene']}: r={row['Pearson_r']:.3f}, P={row['P_value']:.4f}")

# ============================================================
# 10. 按Pattern分组(蛋白沉积窗口对齐)
# ============================================================
print("\n" + "="*60)
print("按蛋白沉积窗口对齐模式分组")
print("="*60)

for pat in aa_liver['Pattern'].unique():
    pat_genes = [g for g, info in liver_aa_data.items() if info['pattern'] == pat]
    if not pat_genes: continue
    pat_corr = sig[sig['Liver_Gene'].isin(pat_genes)]
    pos_r = len(pat_corr[pat_corr['Pearson_r'] > 0])
    neg_r = len(pat_corr[pat_corr['Pearson_r'] < 0])
    print(f"\n{pat} ({len(pat_genes)} genes):")
    print(f"  Significant pairs: {len(pat_corr)} (pos={pos_r}, neg={neg_r})")
    if len(pat_corr) > 0:
        top_abs = pat_corr.sort_values('abs_r', ascending=False).head(5)
        for _, row in top_abs.iterrows():
            print(f"    {row['Liver_Gene']} ↔ {row['Muscle_Gene']}: r={row['Pearson_r']:.3f}, P={row['P_value']:.4f}")

# ============================================================
# 11. 品种内特异性分析
# ============================================================
print("\n" + "="*60)
print("品种内跨组织相关性 (DLY-only vs TFB-only)")
print("="*60)

for breed in ['DLY', 'TFB']:
    breed_corr = []
    for l_gene, l_info in liver_aa_data.items():
        for m_gene, m_info in muscle_aa_data.items():
            l_vals, m_vals = [], []
            for s in common_stages:
                lk = f"{breed}_L_{s}"
                mk = f"{breed}_M_{s}"
                if lk in l_info['mean'].index and mk in m_info['mean'].index:
                    l_vals.append(l_info['mean'][lk])
                    m_vals.append(m_info['mean'][mk])
            if len(l_vals) >= 3:
                r, p = pearsonr(l_vals, m_vals)
                breed_corr.append({
                    'Liver_Gene': l_gene, 'Muscle_Gene': m_gene,
                    'Pearson_r': r, 'P_value': p
                })

    bdf = pd.DataFrame(breed_corr).sort_values('P_value')
    bsig = bdf[bdf['P_value'] < 0.05].sort_values('Pearson_r', ascending=False)
    print(f"\n{breed}: {len(bsig)} significant pairs / {len(bdf)} total")
    if len(bsig) > 0:
        print(bsig.head(10)[['Liver_Gene', 'Muscle_Gene', 'Pearson_r', 'P_value']].to_string())

# ============================================================
# 12. 核心发现：肝肌AA轴候选
# ============================================================
print("\n" + "="*60)
print("核心发现：花生四烯酸肝肌轴候选信号通路")
print("="*60)

# 找出同时满足以下条件的基因对:
# 1. P<0.05
# 2. |r|>0.6
# 3. 肝脏酶pattern为"TFB45与DLY75共同偏高"或"DLY75偏高"
key_pairs = sig[(sig['abs_r'] > 0.6) &
                (sig['Liver_Pattern'].isin(['TFB45与DLY75共同偏高', 'DLY75偏高，TFB45不突出']))]
key_pairs = key_pairs.sort_values('abs_r', ascending=False)

if len(key_pairs) > 0:
    print(f"\n高置信度肝肌AA轴候选 (|r|>0.6, PD窗口对齐): {len(key_pairs)} pairs")
    for _, row in key_pairs.iterrows():
        direction = "正" if row['Pearson_r'] > 0 else "负"
        print(f"  {row['Liver_Gene']}({row['Liver_Category']}) ↔ {row['Muscle_Gene']}: "
              f"r={row['Pearson_r']:.3f}, P={row['P_value']:.4f} [{direction}]")

# ============================================================
# 13. 保存结果
# ============================================================
corr_df.to_excel('/Users/hezongze/pig_study/aa_crosstissue_full_results.xlsx', index=False)
print(f"\nFull results saved: aa_crosstissue_full_results.xlsx ({len(corr_df)} pairs)")

# Detailed data for plotting
detail_rows = []
for _, row in sig[sig['abs_r'] > 0.5].iterrows():
    l_gene = row['Liver_Gene']
    m_gene = row['Muscle_Gene']
    l_info = liver_aa_data[l_gene]
    m_info = muscle_aa_data[m_gene]
    for breed in ['DLY', 'TFB']:
        for s in common_stages:
            lk = f"{breed}_L_{s}"
            mk = f"{breed}_M_{s}"
            if lk in l_info['mean'].index and mk in m_info['mean'].index:
                detail_rows.append({
                    'Liver_Gene': l_gene, 'Muscle_Gene': m_gene,
                    'Breed': breed, 'Stage': s,
                    'Liver_Mean_Expr': l_info['mean'][lk],
                    'Muscle_Mean_Expr': m_info['mean'][mk],
                    'Liver_Category': l_info['category'],
                    'Liver_Pattern': l_info['pattern'],
                    'Pearson_r': row['Pearson_r'], 'P_value': row['P_value'],
                })
pd.DataFrame(detail_rows).to_csv('/Users/hezongze/pig_study/aa_crosstalk_detail_v2.csv', index=False)
print("Detail data saved: aa_crosstalk_detail_v2.csv")

print("\nDone!")
