#!/usr/bin/env python3
"""
Shared statistical utilities for the pig muscle protein deposition study.
Provides: FDR correction (Benjamini-Hochberg), robust correlation helpers,
and sample-size-aware statistical wrappers.
"""
import numpy as np
from scipy.stats import pearsonr, spearmanr, ttest_ind


def benjamini_hochberg(p_values, alpha=0.05):
    """
    Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values : array-like
        Raw p-values.
    alpha : float
        FDR threshold (default 0.05).

    Returns
    -------
    rejected : np.ndarray (bool)
        Whether each hypothesis is rejected at FDR < alpha.
    p_corrected : np.ndarray
        FDR-adjusted p-values (q-values).
    """
    p = np.asarray(p_values, dtype=float)
    mask = ~np.isnan(p)
    n = mask.sum()
    if n == 0:
        return np.zeros_like(p, dtype=bool), np.full_like(p, np.nan)

    # Sort p-values
    sorted_idx = np.argsort(p[mask])
    sorted_p = p[mask][sorted_idx]

    # BH procedure
    ranks = np.arange(1, n + 1)
    bh_thresholds = alpha * ranks / n

    # Find largest k where p(k) <= alpha * k / n
    below = sorted_p <= bh_thresholds
    if below.any():
        max_k = np.max(np.where(below)[0])
        threshold = sorted_p[max_k]
    else:
        threshold = 0.0

    rejected = np.zeros(len(p), dtype=bool)
    p_corrected = np.full_like(p, np.nan)

    # Compute q-values: q_(i) = min_{j >= i} p_(j) * n / j
    # Take cumulative min from right (largest rank) to left (smallest rank)
    q_raw = sorted_p * n / ranks
    q_values = np.minimum.accumulate(q_raw[::-1])[::-1]

    p_corrected[mask] = q_values[np.argsort(sorted_idx)]
    rejected[mask] = p[mask] <= threshold
    rejected = np.asarray(rejected)

    return rejected, p_corrected


def safe_pearsonr(x, y):
    """Pearson correlation with safety against constant arrays and NaN."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 3:
        return np.nan, np.nan
    x, y = x[mask], y[mask]
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan
    return pearsonr(x, y)


def safe_spearmanr(x, y):
    """Spearman correlation with safety against constant arrays and NaN."""
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    if mask.sum() < 3:
        return np.nan, np.nan
    x, y = x[mask], y[mask]
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan
    rho, p = spearmanr(x, y)
    return rho, p


def safe_ttest(a, b):
    """Welch's t-test with safety against small samples and zero variance."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    if np.std(a) == 0 and np.std(b) == 0:
        return np.nan, np.nan
    return ttest_ind(a, b, equal_var=False)


def apply_fdr_to_dataframe(df, p_col='p_value', alpha=0.05):
    """
    Add FDR correction columns to a DataFrame containing raw p-values.

    Adds columns: 'FDR_significant' (bool), 'q_value' (FDR-adjusted p).
    Returns the modified DataFrame.
    """
    if p_col not in df.columns:
        raise ValueError(f"Column '{p_col}' not found in DataFrame")
    p_vals = df[p_col].values
    rejected, q_vals = benjamini_hochberg(p_vals, alpha=alpha)
    df = df.copy()
    df['FDR_significant'] = rejected
    df['q_value'] = q_vals
    return df
