# Pig Study: 猪骨骼肌蛋白沉积/氮利用分析

DLY (杜长大) vs TFB (通城猪) 肝脏转录组 + 肌肉面板 + 血清代谢组分析。

## 新电脑上一键运行

```bash
# 1. 克隆
git clone <repo-url> && cd pig_study

# 2. 安装依赖
pip install pandas numpy scipy matplotlib seaborn gseapy openpyxl statsmodels

# 3. 确认依赖就绪
python3 -c "import gseapy; import statsmodels.api; from stats_utils import benjamini_hochberg; print('OK')"

# 4. 核心分析（需按顺序运行）
python3 gsea_multistage_pipeline.py    # 多阶段GSEA + ssGSEA (预计20-30分钟)
python3 hepatokine_crosstalk.py        # Hepatokine筛选 + 肝肌轴crosstalk (预计2-5分钟)
```

## 输出

| 脚本 | 输出文件 |
|------|---------|
| `gsea_multistage_pipeline.py` | `gsea_multistage_deg_results.xlsx` (15,823基因品种效应), `gsea_multistage_enrichment.xlsx` (442 FDR显著GSEA通路), `figures_final/fig_MS1-4.pdf` |
| `hepatokine_crosstalk.py` | `hepatokine_screening_results.xlsx` (11 FDR显著hepatokine), `hepatokine_muscle_crosstalk.xlsx`, `figures_final/fig_HK1-3.pdf` |

## 关键结果

| 分析 | 数字 |
|------|------|
| 多阶段线性模型 | n=48 (DLY/TFB × 15/45/75/105kg) |
| DEG FDR<0.05 | 534 genes (DLY-up 212 / TFB-up 322) |
| GSEA FDR<0.05 | 442 pathways (ALL TFB-enriched) |
| GSEA Top 1 | AA Metabolism (NES=-2.48, FDR=0.001) |
| Hepatokine FDR<0.05 | 11 (ARG2最强: log2FC=-2.82) |
| 核心发现 | TFB肝脏: AA分解↑ + 尿素循环↑ + 蛋白酶体↑ + 翻译↑ + IGFBP1/2↑ |

## 目录结构

```
pig_study/
├── gsea_multistage_pipeline.py    ← 多阶段GSEA (核心脚本1)
├── hepatokine_crosstalk.py        ← 肝肌轴crosstalk (核心脚本2)
├── gsea_45kg_pipeline.py          ← 单阶段GSEA (45kg only, 效力不足)
├── diagnostic_correlation.py      ← Group-mean vs individual相关诊断
├── stats_utils.py                 ← BH FDR + safe统计工具
├── master_analysis.py             ← 主分析
├── advanced_analysis.py           ← 高级分析
├── integrated_4stage_analysis.py  ← 4阶段DEG
├── gene_selection_pipeline.py     ← 基因筛选
├── biological_closure_screening.py← 生物闭合筛选
├── aa_crosstalk_v2.py             ← AA肝肌轴(TPM版)
├── aa_muscle_liver_crosstalk.py   ← AA肝肌轴(原版)
├── aa_mechanism_validation.py     ← AA机制验证
├── figures_*.py                   ← 各版本figure脚本
├── wgcna_*.py                     ← WGCNA (不推荐使用)
├── gene_expression/
│   ├── liver_gene_matrix.xls      ← 肝35,670基因×54样本
│   └── muscle_gene_matrix.xls     ← 肌305基因×58样本
├── serum_all_tidy.csv             ← 血清代谢组
└── figures_final/                 ← 输出图片目录(自动创建)
```

## 注意事项

1. **BH FDR修复**: `stats_utils.py`的q-value公式已修正。早期脚本如需正式qvalue列，请用修复后版本重跑
2. **运行顺序**: `hepatokine_crosstalk.py` 依赖 `gsea_multistage_pipeline.py` 的输出xlsx
3. **数据文件路径**: 所有路径均为repo内相对路径，无需修改
4. **gseapy首次运行**: 会自动下载MSigDB基因集缓存，需联网
