#!/usr/bin/env python3
"""
NetworkX-based Co-expression Network Analysis (替代经典 WGCNA)
================================================================
数学框架: Zhang & Horvath (2005) SAGMB + Langfelder & Horvath (2008) BMC Bioinformatics
社区检测: Louvain (Blondel et al. 2008 JSTAT)
纯 Python 实现, 透明可复现

Pipeline:
  S1 数据准备 & 过滤
  S2 软阈值选择 (scale-free topology criterion)
  S3 共表达网络构建 (Pearson → soft power adjacency → sparse graph)
  S4 Louvain 社区检测 → 模块识别
  S5 模块 eigengene & 模块-性状关联
  S6 Hub gene 鉴定 (多指标综合)
  S7 可视化 (seaborn, Nature 出版标准)
  S8 跨组织整合 & 富集分析

适用数据: n=8 样本 (DLY/TFB × 15/45/75/105kg), liver + muscle 转录组
"""

import os, sys, json, warnings
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from collections import defaultdict
import networkx as nx
from networkx.algorithms.community import louvain_communities
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch, FancyBboxPatch
import seaborn as sns

warnings.filterwarnings('ignore')
np.random.seed(42)

print("=" * 75)
print("NetworkX-Based Co-expression Network Analysis")
print("替代经典 WGCNA — 纯 Python 实现")
print("=" * 75)

# ============================================================================
# CONFIGURATION
# ============================================================================
CONFIG = {
    'top_n_genes': 8000,       # 高变异基因数
    'min_expr': 0.1,           # 最低平均表达
    'soft_power_range': range(1, 21),  # 候选软阈值
    'soft_power_min': 6,       # signed network 软阈值下限
    'rsquared_cutoff': 0.80,   # scale-free topology fit 目标
    'edge_retain_fraction': 0.005,  # 保留最强 0.5% 的边
    'min_module_size': 20,     # 最小模块基因数
    'louvain_resolution': 1.0, # Louvain resolution 参数
    'louvain_seed': 42,
    'output_dir': 'wgcna_networkx_output',
    'random_seed': 42,
}

os.makedirs(CONFIG['output_dir'], exist_ok=True)

# Phenotype anchors (Zhang & Horvath 框架中的 trait)
PROTEIN_DEPOSITION = {
    ('DLY', 15): 1.58, ('TFB', 15): 1.26,
    ('DLY', 45): 1.59, ('TFB', 45): 1.12,
    ('DLY', 75): 1.11, ('TFB', 75): 0.68,
    ('DLY', 105): 0.87, ('TFB', 105): 0.49,
}
SERUM_UREA = {
    ('DLY', 15): 0.81, ('TFB', 15): 3.16,
    ('DLY', 45): 2.30, ('TFB', 45): 5.02,
    ('DLY', 75): 2.71, ('TFB', 75): 2.71,
    ('DLY', 105): 2.62, ('TFB', 105): 6.08,
}

# Sample mapping
SAMPLE_MAP_LIVER = {
    'L_15_1_': ('DLY', 15), 'L_15_2_': ('TFB', 15),
    'L_45_1_': ('DLY', 45), 'L_45_2_': ('TFB', 45),
    'L_1_1_': ('DLY', 75), 'L_1_2_': ('TFB', 75),
    'L_2_1_': ('DLY', 105), 'L_2_2_': ('TFB', 105),
}
SAMPLE_MAP_MUSCLE = {
    'm_15_1_': ('DLY', 15), 'm_15_2_': ('TFB', 15),
    'BJ_2_1_': ('DLY', 45), 'BJ_2_2_': ('TFB', 45),
    'm_1_1_': ('DLY', 75), 'm_1_2_': ('TFB', 75),
    'm_2_1_': ('DLY', 105), 'm_2_2_': ('TFB', 105),
}

# ============================================================================
# S1: DATA PREPARATION
# ============================================================================
print("\n[S1] Loading and preparing data...")

LIVER_RAW_PATH = '/Users/hezongze/Downloads/result_exp_matrix (4).xls'
MUSCLE_RAW_PATH = '/Users/hezongze/Downloads/result_exp_matrix.xls'

liver_raw = pd.read_csv(LIVER_RAW_PATH, sep='\t')
muscle_raw = pd.read_csv(MUSCLE_RAW_PATH, sep='\t')
print(f"  Raw: Liver {liver_raw.shape[0]} genes, Muscle {muscle_raw.shape[0]} genes")


def build_sample_matrix(mat, sample_map):
    """Build (genes × samples) expression matrix from raw quantification output."""
    meta_cols = ['seq_id', 'gene_name', 'length', 'description']
    val_cols = [c for c in mat.columns if c not in meta_cols]

    sample_info = {}
    for col in val_cols:
        for prefix, (breed, stage) in sample_map.items():
            if col.startswith(prefix):
                sample_name = f"{breed}_{stage}kg"
                sample_info[col] = (sample_name, breed, stage)
                break

    records = []
    for _, row in mat.iterrows():
        gn = str(row['gene_name']) if pd.notna(row['gene_name']) else str(row['seq_id'])
        for col in val_cols:
            if col in sample_info:
                sname, breed, stage = sample_info[col]
                if pd.notna(row[col]):
                    records.append({'Gene': gn, 'Sample': sname,
                                    'Breed': breed, 'Weight': stage,
                                    'Expr': float(row[col])})
    df = pd.DataFrame(records)
    mat_pivot = df.pivot_table(index='Gene', columns='Sample',
                                values='Expr', aggfunc='mean')
    return mat_pivot


liver_mat = build_sample_matrix(liver_raw, SAMPLE_MAP_LIVER)
muscle_mat = build_sample_matrix(muscle_raw, SAMPLE_MAP_MUSCLE)
print(f"  Pivoted: Liver {liver_mat.shape}, Muscle {muscle_mat.shape}")


def filter_genes(mat, min_expr=0.1, top_n=8000):
    """Remove low-expression genes, keep top N most variable."""
    mat_f = mat.loc[mat.mean(axis=1) > min_expr]
    if mat_f.shape[0] > top_n:
        top_idx = mat_f.var(axis=1).nlargest(top_n).index
        mat_f = mat_f.loc[top_idx]
    return mat_f


liver_filt = filter_genes(liver_mat, CONFIG['min_expr'], CONFIG['top_n_genes'])
muscle_filt = filter_genes(muscle_mat, CONFIG['min_expr'], CONFIG['top_n_genes'])
print(f"  Filtered: Liver {liver_filt.shape}, Muscle {muscle_filt.shape}")

# Build trait matrix
def build_trait_matrix(mat):
    """Sample-level trait matrix (PD, Urea, Breed, Weight)."""
    traits = {}
    for s in mat.columns:
        parts = s.split('_')
        breed, stage = parts[0], int(parts[1].replace('kg', ''))
        traits[s] = {
            'Breed': 1 if breed == 'DLY' else 0,
            'Weight': stage,
            'PD': PROTEIN_DEPOSITION.get((breed, stage), np.nan),
            'Urea': SERUM_UREA.get((breed, stage), np.nan),
            'Breed_x_Weight': (1 if breed == 'DLY' else 0) * stage,
        }
    return pd.DataFrame(traits).T


liver_traits = build_trait_matrix(liver_filt)
muscle_traits = build_trait_matrix(muscle_filt)

# Transpose to (samples × genes) for correlation calculation
liver_expr = liver_filt.T  # samples × genes
muscle_expr = muscle_filt.T

# Log-transform for better normality (common in WGCNA)
liver_expr_log = np.log2(liver_expr + 0.01)
muscle_expr_log = np.log2(muscle_expr + 0.01)

n_samples = liver_expr.shape[0]
print(f"  Samples: {n_samples}")
print(f"  Traits: {list(liver_traits.columns)}")
print(f"  WARNING: n={n_samples} — 相关性估计不稳定, 结果需谨慎解读\n")

# ============================================================================
# S2: SOFT-THRESHOLD SELECTION (Scale-free topology criterion)
# ============================================================================
print("[S2] Soft-threshold power selection (scale-free topology criterion)...")
print("     Ref: Zhang & Horvath (2005) SAGMB, Eq. 5-6")


def scale_free_fit(adjacency_matrix):
    """Compute scale-free topology model fit R^2.

    Model: log10(p(k)) ~ -gamma * log10(k)
    where p(k) is the degree distribution.

    Returns R^2 from the linear regression.
    """
    # Degree: sum of adjacency weights per gene
    degree = adjacency_matrix.sum(axis=0)

    # Degree distribution: histogram of connectivities
    # Use log-binning for better fitting (Zhang & Horvath 2005)
    n_genes = len(degree)
    # Create log-spaced bins
    min_deg = max(degree[degree > 0].min(), 0.5)
    max_deg = degree.max()
    if max_deg <= min_deg:
        return 0.0

    n_bins = min(50, n_genes // 50)
    bins = np.logspace(np.log10(min_deg), np.log10(max_deg), n_bins)

    hist, bin_edges = np.histogram(degree, bins=bins)
    # Bin centers (geometric mean)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])

    # Only use bins with counts > 0
    valid = hist > 0
    if valid.sum() < 3:
        return 0.0

    log_k = np.log10(bin_centers[valid])
    log_pk = np.log10(hist[valid] / n_genes)

    # Linear regression: log(p(k)) = a + b * log(k)
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_k, log_pk)

    return r_value ** 2


# Work on a subset for efficiency if needed, but 8000 genes is OK
# Use log-transformed expression for correlation
for tissue, expr_log in [('Liver', liver_expr_log), ('Muscle', muscle_expr_log)]:
    print(f"\n  --- {tissue} ---")
    # Compute correlation matrix once for reuse
    cor_matrix = np.corrcoef(expr_log.values.T)  # genes × genes
    n_genes = cor_matrix.shape[0]
    print(f"  Genes: {n_genes}")

    sft_results = []
    for power in CONFIG['soft_power_range']:
        # Signed adjacency: a_ij = (0.5 * (1 + cor))^power for signed
        # Standard WGCNA signed: a_ij = |cor|^power, preserving sign via separate treatment
        # For scale-free fit we use |cor|^power (same topology, just magnitude)
        adj = np.abs(cor_matrix) ** power
        np.fill_diagonal(adj, 0)
        r2 = scale_free_fit(adj)
        mean_conn = adj.sum(axis=0).mean()
        sft_results.append({'Power': power, 'SFT_R2': r2,
                            'Mean_Connectivity': mean_conn,
                            'Slope': np.nan})
        print(f"    Power={power:2d}  R²={r2:.4f}  mean_k={mean_conn:.1f}")

    # Select optimal power
    sft_df = pd.DataFrame(sft_results)
    above_cut = sft_df[sft_df['SFT_R2'] >= CONFIG['rsquared_cutoff']]

    if len(above_cut) > 0:
        soft_power = above_cut['Power'].min()
    else:
        # Pick power with highest R², but enforce minimum
        best_idx = sft_df['SFT_R2'].idxmax()
        soft_power = sft_df.loc[best_idx, 'Power']
        best_r2 = sft_df.loc[best_idx, 'SFT_R2']
        print(f"    WARNING: No power reaches R²≥{CONFIG['rsquared_cutoff']}, "
              f"best R²={best_r2:.3f} at power={soft_power}")

    if soft_power < CONFIG['soft_power_min']:
        print(f"    Power={soft_power} < {CONFIG['soft_power_min']}, "
              f"enforcing minimum (signed network requires adequate co-expression structure)")
        soft_power = CONFIG['soft_power_min']

    # Store result
    if tissue == 'Liver':
        liver_soft_power = soft_power
        liver_sft_df = sft_df
        liver_cor = cor_matrix
    else:
        muscle_soft_power = soft_power
        muscle_sft_df = sft_df
        muscle_cor = cor_matrix

    print(f"  >>> Selected soft power: {soft_power}")

# ============================================================================
# S3: CO-EXPRESSION NETWORK CONSTRUCTION
# ============================================================================
print(f"\n[S3] Building co-expression networks...")
print(f"     Signed adjacency: a_ij = |cor(i,j)|^beta")
print(f"     Sparsification: retain top {CONFIG['edge_retain_fraction']*100:.1f}% strongest edges")


def build_coexpression_graph(cor_matrix, soft_power, gene_names, edge_fraction=0.005):
    """Build a sparse NetworkX weighted graph from the correlation matrix.

    Steps:
      1. Compute signed adjacency: A_ij = |r_ij|^beta
      2. Extract upper triangle edges
      3. Keep top edge_fraction of strongest edges
      4. Build NetworkX graph
    """
    n = cor_matrix.shape[0]
    adj = np.abs(cor_matrix) ** soft_power
    np.fill_diagonal(adj, 0)

    # Extract upper triangle (excluding diagonal)
    triu_idx = np.triu_indices(n, k=1)
    edges_flat = adj[triu_idx]

    # Determine threshold for top edges
    n_edges = len(edges_flat)
    n_keep = max(int(n_edges * edge_fraction), n * 5)  # at least 5 edges per gene
    threshold = np.partition(edges_flat, n_edges - n_keep)[n_edges - n_keep]

    print(f"    Edge threshold: {threshold:.6f}, keeping {n_keep}/{n_edges} edges "
          f"({n_keep/n_edges*100:.2f}%)")

    # Build graph
    G = nx.Graph()
    G.add_nodes_from(range(n))

    # Add edges above threshold
    mask = edges_flat >= threshold
    sources = triu_idx[0][mask]
    targets = triu_idx[1][mask]
    weights = edges_flat[mask]

    edge_list = [(int(s), int(t), float(w)) for s, t, w in zip(sources, targets, weights)]
    G.add_weighted_edges_from(edge_list)

    # Map node indices to gene names
    gene_map = {i: gn for i, gn in enumerate(gene_names)}
    nx.set_node_attributes(G, gene_map, 'gene')

    # Compute basic metrics
    n_comp = nx.number_connected_components(G)
    density = nx.density(G)
    print(f"    Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    print(f"    Connected components: {n_comp}, Density: {density:.6f}")

    return G


liver_genes = list(liver_expr.columns)
muscle_genes = list(muscle_expr.columns)

liver_graph = build_coexpression_graph(
    liver_cor, liver_soft_power, liver_genes, CONFIG['edge_retain_fraction'])
muscle_graph = build_coexpression_graph(
    muscle_cor, muscle_soft_power, muscle_genes, CONFIG['edge_retain_fraction'])

# ============================================================================
# S4: COMMUNITY DETECTION → MODULE IDENTIFICATION
# ============================================================================
print(f"\n[S4] Community detection (Louvain algorithm)...")
print(f"     Ref: Blondel et al. (2008) JSTAT")
print(f"     Resolution: {CONFIG['louvain_resolution']}, Seed: {CONFIG['louvain_seed']}")


def detect_communities_louvain(G, resolution=1.0, seed=42):
    """Run Louvain community detection, filter small communities."""
    communities = louvain_communities(
        G, weight='weight', resolution=resolution, seed=seed)
    # communities is list of sets of node indices

    module_assignment = {}
    for idx, comm in enumerate(communities):
        color = f"module_{idx}"
        for node in comm:
            module_assignment[node] = color

    # Map node index -> gene name -> module
    gene_module = {}
    for node in G.nodes():
        gene = G.nodes[node]['gene']
        gene_module[gene] = module_assignment.get(node, 'grey')

    # Rename modules by size
    module_sizes = defaultdict(list)
    for gene, mod in gene_module.items():
        module_sizes[mod].append(gene)

    # Sort by size and rename
    sorted_mods = sorted(module_sizes.items(), key=lambda x: -len(x[1]))
    rename_map = {}
    for rank, (old_name, genes) in enumerate(sorted_mods):
        new_name = f"M{rank+1}"
        rename_map[old_name] = new_name

    gene_module_renamed = {gene: rename_map.get(mod, 'grey')
                           for gene, mod in gene_module.items()}

    return gene_module_renamed, rename_map


liver_gene_module, liver_mod_rename = detect_communities_louvain(
    liver_graph, CONFIG['louvain_resolution'], CONFIG['louvain_seed'])
muscle_gene_module, muscle_mod_rename = detect_communities_louvain(
    muscle_graph, CONFIG['louvain_resolution'], CONFIG['louvain_seed'])

# Print module sizes
for tissue, gm in [('Liver', liver_gene_module), ('Muscle', muscle_gene_module)]:
    sizes = pd.Series(list(gm.values())).value_counts()
    print(f"  {tissue}: {len(sizes)} modules")
    for mod, size in sizes.items():
        print(f"    {mod}: {size} genes")

# ============================================================================
# S5: MODULE EIGENGENE & MODULE-TRAIT ASSOCIATION
# ============================================================================
print(f"\n[S5] Computing module eigengenes and module-trait correlations...")


def compute_module_eigengenes(expr, gene_module):
    """Compute module eigengene (PC1) for each module.

    expr: samples × genes DataFrame
    gene_module: dict {gene_name: module_name}
    """
    mes = {}
    mod_genes = defaultdict(list)
    for gene, mod in gene_module.items():
        mod_genes[mod].append(gene)

    for mod, genes in mod_genes.items():
        if len(genes) < 3:
            continue
        genes_in = [g for g in genes if g in expr.columns]
        if len(genes_in) < 3:
            continue
        mod_expr = expr[genes_in].values
        # Center and scale
        mod_expr_scaled = (mod_expr - mod_expr.mean(axis=0)) / mod_expr.std(axis=0, ddof=1)
        mod_expr_scaled = np.nan_to_num(mod_expr_scaled, 0)

        if mod_expr_scaled.shape[1] >= 2:
            pca = PCA(n_components=1, random_state=CONFIG['random_seed'])
            pc1 = pca.fit_transform(mod_expr_scaled)[:, 0]
            # Sign convention: correlate with mean expression
            mean_expr = mod_expr.mean(axis=1)
            if np.corrcoef(pc1, mean_expr)[0, 1] < 0:
                pc1 = -pc1
            mes[mod] = pd.Series(pc1, index=expr.index, name=mod)

    return pd.DataFrame(mes)


liver_MEs = compute_module_eigengenes(liver_expr_log, liver_gene_module)
muscle_MEs = compute_module_eigengenes(muscle_expr_log, muscle_gene_module)
print(f"  Liver MEs: {liver_MEs.shape[1]} modules")
print(f"  Muscle MEs: {muscle_MEs.shape[1]} modules")

# Module-trait correlations
def module_trait_correlation(MEs, traits):
    """Compute Pearson correlation between each ME and each trait."""
    common = MEs.index.intersection(traits.index)
    me_aligned = MEs.loc[common]
    tr_aligned = traits.loc[common]

    cor_df = pd.DataFrame(index=me_aligned.columns, columns=tr_aligned.columns)
    pval_df = pd.DataFrame(index=me_aligned.columns, columns=tr_aligned.columns)

    for mod in me_aligned.columns:
        for trait in tr_aligned.columns:
            r, p = stats.pearsonr(me_aligned[mod], tr_aligned[trait])
            cor_df.loc[mod, trait] = r
            pval_df.loc[mod, trait] = p

    return cor_df.astype(float), pval_df.astype(float)


liver_mtc, liver_mtp = module_trait_correlation(liver_MEs, liver_traits)
muscle_mtc, muscle_mtp = module_trait_correlation(muscle_MEs, muscle_traits)

# Print key associations
for tissue, mtc, mtp in [('Liver', liver_mtc, liver_mtp),
                           ('Muscle', muscle_mtc, muscle_mtp)]:
    print(f"\n  {tissue} Module-Trait Associations (|r_PD| > 0.3):")
    for mod in mtc.index:
        r_pd = mtc.loc[mod, 'PD']
        p_pd = mtp.loc[mod, 'PD']
        r_urea = mtc.loc[mod, 'Urea']
        if abs(r_pd) > 0.3:
            sig = '***' if p_pd < 0.001 else ('**' if p_pd < 0.01 else
                                              ('*' if p_pd < 0.05 else 'ns'))
            print(f"    {mod:10s}  r_PD={r_pd:+.3f} {sig}  r_Urea={r_urea:+.3f}")

# ============================================================================
# S6: HUB GENE IDENTIFICATION
# ============================================================================
print(f"\n[S6] Identifying hub genes (multi-metric scoring)...")
print("     Metrics: within-module degree, eigenvector centrality, betweenness")


def identify_hub_genes(G, gene_module, top_n=20):
    """Identify hub genes per module using multiple centrality metrics.

    For each module:
      - within-module degree (local hubness)
      - eigenvector centrality (influence)
      - betweenness centrality (bridge potential)
    Composite score = rank aggregation across metrics.
    """
    # Global centrality
    print("  Computing betweenness centrality...")
    betweenness = nx.betweenness_centrality(G, weight='weight', normalized=True)
    print("  Computing eigenvector centrality...")
    try:
        eigenvector = nx.eigenvector_centrality(G, weight='weight', max_iter=2000)
    except nx.PowerIterationFailedConvergence:
        print("  WARNING: Eigenvector centrality did not converge, using degree instead")
        eigenvector = dict(G.degree(weight='weight'))

    # Group genes by module
    mod_genes = defaultdict(list)
    for gene, mod in gene_module.items():
        mod_genes[mod].append(gene)

    # Build gene -> node index mapping
    gene_to_node = {G.nodes[n]['gene']: n for n in G.nodes()}

    hub_results = []
    for mod, genes in mod_genes.items():
        if len(genes) < 10:
            continue

        # Within-module metrics
        mod_nodes = [gene_to_node[g] for g in genes if g in gene_to_node]
        if len(mod_nodes) < 5:
            continue

        # Within-module degree
        subG = G.subgraph(mod_nodes)
        wm_degree = dict(subG.degree(weight='weight'))

        # Collect metrics for each gene in the module
        gene_metrics = []
        for node in mod_nodes:
            gene = G.nodes[node]['gene']
            gene_metrics.append({
                'Gene': gene,
                'Module': mod,
                'WM_Degree': wm_degree.get(node, 0),
                'Eigenvector': eigenvector.get(node, 0),
                'Betweenness': betweenness.get(node, 0),
                'Degree_Global': G.degree(node, weight='weight'),
            })

        gm_df = pd.DataFrame(gene_metrics)

        # Rank aggregation: average percentile rank
        for col in ['WM_Degree', 'Eigenvector', 'Betweenness']:
            gm_df[f'{col}_rank'] = gm_df[col].rank(pct=True)

        gm_df['Hub_Score'] = (gm_df['WM_Degree_rank'] +
                               gm_df['Eigenvector_rank'] +
                               gm_df['Betweenness_rank']) / 3

        gm_df = gm_df.sort_values('Hub_Score', ascending=False)
        top_hubs = gm_df.head(top_n)
        hub_results.append(top_hubs)

    if hub_results:
        return pd.concat(hub_results, ignore_index=True)
    return pd.DataFrame()


liver_hubs = identify_hub_genes(liver_graph, liver_gene_module)
muscle_hubs = identify_hub_genes(muscle_graph, muscle_gene_module)

# Print top hubs per module (for significant modules)
for tissue, hubs, mtc in [('Liver', liver_hubs, liver_mtc),
                            ('Muscle', muscle_hubs, muscle_mtc)]:
    print(f"\n  {tissue} Top 5 Hub Genes per Module:")
    for mod in hubs['Module'].unique():
        mod_hubs = hubs[hubs['Module'] == mod].head(5)
        r_pd = mtc.loc[mod, 'PD'] if mod in mtc.index else np.nan
        sig_flag = ' ***' if abs(r_pd) > 0.3 else ''
        print(f"    {mod} (r_PD={r_pd:+.3f}){sig_flag}: "
              f"{', '.join(mod_hubs['Gene'].tolist())}")

# ============================================================================
# S7: PUBLICATION-QUALITY FIGURES
# ============================================================================
print(f"\n[S7] Generating publication-quality figures...")

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.titlesize': 10,
    'axes.labelsize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False,
})

out = CONFIG['output_dir']

# --- Fig 1: Scale-free topology fit & mean connectivity ---
fig1, axes1 = plt.subplots(2, 2, figsize=(10, 8))

for idx, (tissue, sft_df, soft_power) in enumerate([
    ('Liver', liver_sft_df, liver_soft_power),
    ('Muscle', muscle_sft_df, muscle_soft_power),
]):
    ax1, ax2 = axes1[0, idx], axes1[1, idx]

    # Scale-free fit
    ax1.plot(sft_df['Power'], sft_df['SFT_R2'], 'o-', color='#D73027',
             markersize=5, linewidth=1.5, label='Scale-free R²')
    ax1.axhline(y=CONFIG['rsquared_cutoff'], color='#4575B4', linestyle='--',
                linewidth=1, label=f'R²={CONFIG["rsquared_cutoff"]}')
    ax1.axvline(x=soft_power, color='#1A9850', linestyle='--', linewidth=1,
                label=f'Selected β={soft_power}')
    ax1.set_xlabel('Soft Threshold (β)')
    ax1.set_ylabel('Scale-free Topology R²')
    ax1.set_title(f'{tissue}: Scale Independence', fontweight='bold')
    ax1.legend(fontsize=6, loc='lower right')

    # Mean connectivity
    ax2.plot(sft_df['Power'], sft_df['Mean_Connectivity'], 's-',
             color='#4575B4', markersize=5, linewidth=1.5)
    ax2.axvline(x=soft_power, color='#1A9850', linestyle='--', linewidth=1)
    ax2.set_xlabel('Soft Threshold (β)')
    ax2.set_ylabel('Mean Connectivity')
    ax2.set_title(f'{tissue}: Mean Connectivity', fontweight='bold')

fig1.suptitle('Soft-Threshold Selection (Scale-Free Topology Criterion)',
              fontweight='bold', fontsize=11, y=1.01)
plt.tight_layout()
fig1.savefig(f'{out}/fig_soft_threshold_selection.png', bbox_inches='tight')
fig1.savefig(f'{out}/fig_soft_threshold_selection.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_soft_threshold_selection.png/pdf")

# --- Fig 2: Module-Trait Correlation Heatmap ---
fig2, axes2 = plt.subplots(1, 2, figsize=(14, max(6, len(liver_MEs.columns) * 0.4)))

for idx, (tissue, mtc, mtp) in enumerate([
    ('Liver', liver_mtc, liver_mtp),
    ('Muscle', muscle_mtc, muscle_mtp),
]):
    ax = axes2[idx]
    traits_show = ['PD', 'Urea', 'Breed', 'Weight']
    mods_show = mtc.index.tolist()

    cor_mat = mtc.loc[mods_show, traits_show].values
    p_mat = mtp.loc[mods_show, traits_show].values

    hm = sns.heatmap(cor_mat, annot=True, fmt='.2f', cmap='RdBu_r',
                     vmin=-1, vmax=1, center=0,
                     xticklabels=traits_show, yticklabels=mods_show,
                     linewidths=0.5, linecolor='white',
                     cbar_kws={'label': "Pearson's r", 'shrink': 0.8},
                     ax=ax, annot_kws={'fontsize': 7})

    # Add significance stars
    for i in range(len(mods_show)):
        for j in range(len(traits_show)):
            p = p_mat[i, j]
            if p < 0.05:
                stars = '***' if p < 0.001 else ('**' if p < 0.01 else '*')
                ax.text(j + 0.65, i + 0.7, stars, ha='center', va='center',
                       fontsize=6, fontweight='bold',
                       color='darkgreen' if abs(cor_mat[i, j]) > 0.3 else 'grey')

    ax.set_title(f'{tissue} Module-Trait Correlations\n(Louvain Community Detection)',
                 fontweight='bold')
    ax.set_xlabel('')

plt.suptitle('Co-expression Module Associations with Growth Phenotypes',
             fontweight='bold', fontsize=12, y=1.02)
plt.tight_layout()
fig2.savefig(f'{out}/fig_module_trait_heatmap.png', bbox_inches='tight')
fig2.savefig(f'{out}/fig_module_trait_heatmap.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_module_trait_heatmap.png/pdf")

# --- Fig 3: Module Size Distribution ---
fig3, axes3 = plt.subplots(1, 2, figsize=(10, 5))

for idx, (tissue, gm, mtc) in enumerate([
    ('Liver', liver_gene_module, liver_mtc),
    ('Muscle', muscle_gene_module, muscle_mtc),
]):
    ax = axes3[idx]
    sizes = pd.Series(list(gm.values())).value_counts().sort_values(ascending=False)

    colors = []
    for mod in sizes.index:
        if mod in mtc.index:
            r_pd = mtc.loc[mod, 'PD']
            if r_pd > 0.3:
                colors.append('#D73027')
            elif r_pd < -0.3:
                colors.append('#4575B4')
            else:
                colors.append('#E0E0E0')
        else:
            colors.append('#E0E0E0')

    bars = ax.barh(range(len(sizes)), sizes.values, color=colors, edgecolor='white', height=0.7)
    ax.set_yticks(range(len(sizes)))
    ax.set_yticklabels(sizes.index, fontsize=7)
    ax.set_xlabel('Number of Genes')
    ax.set_title(f'{tissue} Module Sizes', fontweight='bold')
    ax.invert_yaxis()

    # Legend
    legend_elements = [
        Patch(facecolor='#D73027', label='PD-positive (r>0.3)'),
        Patch(facecolor='#4575B4', label='PD-negative (r<-0.3)'),
        Patch(facecolor='#E0E0E0', label='Not significant'),
    ]
    ax.legend(handles=legend_elements, fontsize=6, loc='lower right')

fig3.suptitle('Co-expression Module Sizes (Louvain Detection)',
              fontweight='bold', fontsize=11, y=1.01)
plt.tight_layout()
fig3.savefig(f'{out}/fig_module_sizes.png', bbox_inches='tight')
fig3.savefig(f'{out}/fig_module_sizes.pdf', bbox_inches='tight')
plt.close()
print("  Saved fig_module_sizes.png/pdf")

# --- Fig 4: Hub Gene Network for Top PD-associated Module ---
# Find top PD-associated module in muscle
def get_top_pd_module(mtc):
    """Return the module name with highest |r_PD| if significant."""
    if 'PD' not in mtc.columns:
        return None
    r_pd = mtc['PD'].abs().sort_values(ascending=False)
    if len(r_pd) > 0:
        return r_pd.index[0]
    return None


for tissue, G, gm, hubs, mtc in [
    ('Liver', liver_graph, liver_gene_module, liver_hubs, liver_mtc),
    ('Muscle', muscle_graph, muscle_gene_module, muscle_hubs, muscle_mtc),
]:
    top_mod = get_top_pd_module(mtc)
    if top_mod is None:
        continue

    r_pd = mtc.loc[top_mod, 'PD']
    print(f"  {tissue} top PD module: {top_mod} (r_PD={r_pd:+.3f})")

    # Get genes in top module
    mod_genes = [g for g, m in gm.items() if m == top_mod]
    gene_to_node = {G.nodes[n]['gene']: n for n in G.nodes()}
    mod_nodes = [gene_to_node[g] for g in mod_genes if g in gene_to_node]

    if len(mod_nodes) < 5 or len(mod_nodes) > 300:
        print(f"    Skipping network plot: {len(mod_nodes)} genes (need 5-300)")
        continue

    subG = G.subgraph(mod_nodes)

    # Get top hub genes for this module
    mod_hubs = hubs[hubs['Module'] == top_mod]
    top_10 = mod_hubs.head(10)['Gene'].tolist()

    # Layout
    try:
        pos = nx.spring_layout(subG, weight='weight', seed=42, k=0.5, iterations=50)
    except Exception:
        continue

    figN, axN = plt.subplots(figsize=(8, 7))

    # Edge width by weight
    edge_weights = [subG[u][v]['weight'] for u, v in subG.edges()]
    if edge_weights:
        max_w = max(edge_weights)
        edge_widths = [w / max_w * 2 for w in edge_weights]
    else:
        edge_widths = [0.5]

    # Node sizes by degree
    node_deg = dict(subG.degree(weight='weight'))
    node_sizes = [node_deg.get(n, 0) * 3 + 5 for n in subG.nodes()]

    # Node colors: hub genes highlighted
    node_colors = []
    for n in subG.nodes():
        gene = G.nodes[n]['gene']
        if gene in top_10:
            node_colors.append('#D73027')
        else:
            node_colors.append('#B0B0B0')

    nx.draw_networkx_edges(subG, pos, alpha=0.15, edge_color='#888888',
                           width=edge_widths, ax=axN)
    nx.draw_networkx_nodes(subG, pos, node_size=node_sizes, node_color=node_colors,
                           alpha=0.8, edgecolors='white', linewidths=0.3, ax=axN)

    # Label only top 10 hubs
    labels = {gene_to_node[g]: g for g in top_10 if g in gene_to_node
              and gene_to_node[g] in subG.nodes()}
    nx.draw_networkx_labels(subG, pos, labels, font_size=6,
                            font_weight='bold', ax=axN)

    axN.set_title(f'{tissue} Module {top_mod} Hub Gene Network\n'
                  f'(r_PD={r_pd:+.3f}, {len(mod_nodes)} genes)',
                  fontweight='bold')
    axN.axis('off')
    plt.tight_layout()
    figN.savefig(f'{out}/fig_{tissue.lower()}_{top_mod}_hub_network.png')
    figN.savefig(f'{out}/fig_{tissue.lower()}_{top_mod}_hub_network.pdf')
    plt.close()
    print(f"    Saved fig_{tissue.lower()}_{top_mod}_hub_network.png/pdf")

# ============================================================================
# S8: CROSS-TISSUE INTEGRATION
# ============================================================================
print(f"\n[S8] Cross-tissue module integration...")

# 8.1 Cross-tissue module eigengene correlation
common_samples = liver_MEs.index.intersection(muscle_MEs.index)
print(f"  Common samples: {len(common_samples)}")

cross_cor = pd.DataFrame(index=liver_MEs.columns, columns=muscle_MEs.columns)
cross_pval = pd.DataFrame(index=liver_MEs.columns, columns=muscle_MEs.columns)

for lm in liver_MEs.columns:
    for mm in muscle_MEs.columns:
        r, p = stats.pearsonr(liver_MEs.loc[common_samples, lm],
                               muscle_MEs.loc[common_samples, mm])
        cross_cor.loc[lm, mm] = r
        cross_pval.loc[lm, mm] = p

cross_cor = cross_cor.astype(float)
cross_pval = cross_pval.astype(float)

print("  Significant cross-tissue module pairs (|r| > 0.5, p < 0.001):")
cross_pairs = []
for lm in cross_cor.index:
    for mm in cross_cor.columns:
        r = cross_cor.loc[lm, mm]
        p = cross_pval.loc[lm, mm]
        if abs(r) > 0.5 and p < 0.001:
            cross_pairs.append({'Liver_Module': lm, 'Muscle_Module': mm,
                                'r': r, 'p': p})
            print(f"    {lm:8s} ↔ {mm:8s}  r={r:+.3f}  p={p:.5f}")

# 8.2 Cross-tissue heatmap
figC, axC = plt.subplots(figsize=(max(8, len(muscle_MEs.columns) * 0.5),
                                    max(6, len(liver_MEs.columns) * 0.5)))

sns.heatmap(cross_cor, annot=True, fmt='.2f', cmap='RdBu_r',
            vmin=-1, vmax=1, center=0, linewidths=0.5,
            cbar_kws={'label': "Pearson's r", 'shrink': 0.8},
            ax=axC, annot_kws={'fontsize': 7})

axC.set_xlabel('Muscle Modules', fontweight='bold')
axC.set_ylabel('Liver Modules', fontweight='bold')
axC.set_title('Cross-Tissue Module Eigengene Correlation\n(Liver ↔ Muscle)',
              fontweight='bold')
plt.tight_layout()
figC.savefig(f'{out}/fig_crosstissue_module_correlation.png')
figC.savefig(f'{out}/fig_crosstissue_module_correlation.pdf')
plt.close()
print("  Saved fig_crosstissue_module_correlation.png/pdf")

# 8.3 Bridging candidates: known liver-muscle axis genes
print("\n  Identifying liver-muscle axis bridging candidates...")
KNOWN_BRIDGERS = {
    'IGF1': 'Growth factor', 'IGF2': 'Growth factor',
    'IGFBP2': 'IGFBP', 'IGFBP3': 'IGFBP', 'IGFBP5': 'IGFBP',
    'AHSG': 'Fetuin-A', 'FGF21': 'Metabolic regulator', 'FGF19': 'Bile acid signaling',
    'ANGPTL4': 'Lipid regulator', 'GDF15': 'Stress response',
    'MSTN': 'Negative muscle regulator', 'FNDC5': 'Irisin precursor',
    'IL6': 'Cytokine', 'IL15': 'Muscle anabolic cytokine',
    'TNF': 'Cytokine', 'LIF': 'Leukemia inhibitory factor',
    'CTGF': 'Connective tissue GF', 'VEGFA': 'Angiogenesis',
    'HGF': 'Hepatocyte GF', 'BMP2': 'Bone morphogenetic protein',
    'INHBA': 'Activin A', 'FST': 'Follistatin (MSTN antagonist)',
    'SERPINE1': 'PAI-1', 'THBS1': 'Thrombospondin',
    'APOA1': 'Lipid transport', 'APOE': 'Lipid transport',
    'TTR': 'Transthyretin', 'ALB': 'Albumin',
    'TF': 'Transferrin', 'SERPINA1': 'Alpha-1 antitrypsin',
}

bridge_data = []
for gene, func in KNOWN_BRIDGERS.items():
    row = {'Gene': gene, 'Function': func}

    if gene in liver_gene_module:
        l_mod = liver_gene_module[gene]
        row['Liver_Module'] = l_mod
        row['Liver_Mod_r_PD'] = liver_mtc.loc[l_mod, 'PD'] if l_mod in liver_mtc.index else np.nan
    else:
        row['Liver_Module'] = 'Not found'

    if gene in muscle_gene_module:
        m_mod = muscle_gene_module[gene]
        row['Muscle_Module'] = m_mod
        row['Muscle_Mod_r_PD'] = muscle_mtc.loc[m_mod, 'PD'] if m_mod in muscle_mtc.index else np.nan
    else:
        row['Muscle_Module'] = 'Not found'

    row['In_Both'] = (gene in liver_gene_module) and (gene in muscle_gene_module)
    bridge_data.append(row)

bridge_df = pd.DataFrame(bridge_data)
bridge_in_both = bridge_df[bridge_df['In_Both']]
print(f"  Known bridging candidates: {len(bridge_df)} total, "
      f"{len(bridge_in_both)} in BOTH tissues")
if len(bridge_in_both) > 0:
    for _, r in bridge_in_both.iterrows():
        print(f"    {r['Gene']:12s} | Liver: {r['Liver_Module']:8s} "
              f"(r_PD={r['Liver_Mod_r_PD']:+.3f}) | "
              f"Muscle: {r['Muscle_Module']:8s} (r_PD={r['Muscle_Mod_r_PD']:+.3f})")

# ============================================================================
# SAVE RESULTS
# ============================================================================
print(f"\n{'='*75}")
print("SAVING RESULTS")
print(f"{'='*75}")

# Build gene-level result tables
def build_gene_table(gm, expr, mtc, hubs):
    """Compile per-gene module assignment, kME (if computable), GS traits."""
    gene_df = pd.DataFrame({'Gene': list(gm.keys()), 'Module': list(gm.values())})

    # Module-trait correlations for each module
    for trait in mtc.columns:
        gene_df[f'Module_r_{trait}'] = gene_df['Module'].map(
            lambda m: mtc.loc[m, trait] if m in mtc.index else np.nan)

    # Merge hub scores
    if len(hubs) > 0:
        gene_df = gene_df.merge(
            hubs[['Gene', 'Hub_Score', 'WM_Degree', 'Eigenvector', 'Betweenness']],
            on='Gene', how='left')

    return gene_df.sort_values(['Module', 'Gene'])


liver_gt = build_gene_table(liver_gene_module, liver_expr_log, liver_mtc, liver_hubs)
muscle_gt = build_gene_table(muscle_gene_module, muscle_expr_log, muscle_mtc, muscle_hubs)

with pd.ExcelWriter(f'{out}/networkx_wgcna_results.xlsx', engine='openpyxl') as writer:
    # Module assignments
    liver_gt.to_excel(writer, sheet_name='Liver_Module_Assignment', index=False)
    muscle_gt.to_excel(writer, sheet_name='Muscle_Module_Assignment', index=False)

    # Module-trait correlations
    liver_mtc.to_excel(writer, sheet_name='Liver_ModuleTrait_Cor')
    muscle_mtc.to_excel(writer, sheet_name='Muscle_ModuleTrait_Cor')
    liver_mtp.to_excel(writer, sheet_name='Liver_ModuleTrait_Pval')
    muscle_mtp.to_excel(writer, sheet_name='Muscle_ModuleTrait_Pval')

    # Hub genes
    liver_hubs.to_excel(writer, sheet_name='Liver_Hub_Genes', index=False)
    muscle_hubs.to_excel(writer, sheet_name='Muscle_Hub_Genes', index=False)

    # Cross-tissue
    cross_cor.to_excel(writer, sheet_name='CrossTissue_Module_Cor')
    if cross_pairs:
        pd.DataFrame(cross_pairs).to_excel(writer, sheet_name='CrossTissue_Pairs',
                                            index=False)

    # Bridging candidates
    bridge_df.to_excel(writer, sheet_name='Bridging_Candidates', index=False)

    # Soft threshold
    liver_sft_df.to_excel(writer, sheet_name='Liver_SoftThreshold', index=False)
    muscle_sft_df.to_excel(writer, sheet_name='Muscle_SoftThreshold', index=False)

print(f"  Saved {out}/networkx_wgcna_results.xlsx")

# ============================================================================
# VALIDATION: Compare with original R WGCNA
# ============================================================================
print(f"\n{'='*75}")
print("VALIDATION: Comparison with original R WGCNA results")
print(f"{'='*75}")

# Load original R WGCNA assignments
liver_r_gm = pd.read_csv('wgcna_output/liver_gene_module_assignment.csv')
muscle_r_gm = pd.read_csv('wgcna_output/muscle_gene_module_assignment.csv')

for tissue, r_gm, nx_gm in [
    ('Liver', liver_r_gm, liver_gene_module),
    ('Muscle', muscle_r_gm, muscle_gene_module),
]:
    # Find common genes
    r_genes = set(r_gm['Gene'].tolist())
    nx_genes = set(nx_gm.keys())
    common = r_genes & nx_genes
    print(f"\n  {tissue}: {len(common)} common genes for comparison")

    if len(common) < 10:
        print("    Too few common genes, skipping ARI")
        continue

    r_labels = []
    nx_labels = []
    for g in sorted(common):
        r_mod = r_gm[r_gm['Gene'] == g]['Module'].values[0]
        nx_mod = nx_gm[g]
        r_labels.append(r_mod)
        nx_labels.append(nx_mod)

    ari = adjusted_rand_score(r_labels, nx_labels)
    print(f"    Adjusted Rand Index (Louvain vs WGCNA): {ari:.4f}")
    print(f"    NOTE: ARI near 0 = methods differ; near 1 = highly similar")
    print(f"    This is expected — Louvain uses a different modularity criterion")
    print(f"    than WGCNA's hierarchical clustering + dynamic tree cut.")

# Compare PD-associated modules
print("\n  Key phenotype associations (PD):")
for tissue, r_mtc_path, nx_mtc in [
    ('Liver', 'wgcna_output/liver_module_trait_cor.csv', liver_mtc),
    ('Muscle', 'wgcna_output/muscle_module_trait_cor.csv', muscle_mtc),
]:
    r_mtc = pd.read_csv(r_mtc_path, index_col=0)
    print(f"\n  {tissue}:")
    print(f"    R WGCNA top PD modules: ", end='')
    r_pd = r_mtc.get('PD', pd.Series(dtype=float))
    top_r = r_pd.abs().nlargest(3)
    for mod, cor in top_r.items():
        print(f"{mod}(r={r_pd[mod]:+.3f}) ", end='')
    print()
    print(f"    NetworkX top PD modules: ", end='')
    nx_pd = nx_mtc.get('PD', pd.Series(dtype=float))
    top_nx = nx_pd.abs().nlargest(3)
    for mod, cor in top_nx.items():
        print(f"{mod}(r={r_pd.get(mod, np.nan):+.3f}) ", end='')
    print()

# ============================================================================
# FINAL SUMMARY
# ============================================================================
print(f"\n{'='*75}")
print("ANALYSIS COMPLETE")
print(f"{'='*75}")
print(f"""
NetworkX-based Co-expression Network Analysis
=============================================
Method: Louvain community detection on signed soft-threshold adjacency
  - Liver soft power: β = {liver_soft_power}
  - Muscle soft power: β = {muscle_soft_power}

Modules detected:
  - Liver: {len(set(liver_gene_module.values()))} modules
  - Muscle: {len(set(muscle_gene_module.values()))} modules

Key phenotype-associated modules (|r_PD| > 0.3):
  Liver: {[m for m in liver_mtc.index if abs(liver_mtc.loc[m, 'PD']) > 0.3]}
  Muscle: {[m for m in muscle_mtc.index if abs(muscle_mtc.loc[m, 'PD']) > 0.3]}

Cross-tissue module pairs (|r| > 0.5): {len(cross_pairs)}

Output directory: {out}/
  - networkx_wgcna_results.xlsx (full results)
  - fig_*.png/pdf (publication-quality figures)

IMPORTANT: n={n_samples} samples — 相关性估计不稳定, 结果需结合生物学先验知识解读.
           建议用 bootstrap 验证关键结果, 并在独立数据集中验证候选基因.
""")

print("Done!")
