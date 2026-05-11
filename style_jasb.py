#!/usr/bin/env python3
"""
Unified JASB (Journal of Animal Science and Biotechnology) figure style.
Single import, consistent formatting across all figures.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import matplotlib.ticker as ticker
import seaborn as sns
import os
import numpy as np

# ============================================================
# Color Palette — colorblind-friendly + JASB compatible
# ============================================================
C_DLY = '#2166AC'        # DLY breed (blue)
C_TFB = '#B2182B'        # TFB breed (red)
C_DLY_LIGHT = '#92C5DE'
C_TFB_LIGHT = '#F4A582'
C_NS = '#BDBDBD'         # not significant
C_UP = '#2166AC'         # upregulated
C_DN = '#B2182B'         # downregulated
C_STAGE_15  = '#FEE08B'  # stage colors (light to dark)
C_STAGE_45  = '#FDAE61'
C_STAGE_75  = '#F46D43'
C_STAGE_105 = '#A50026'

# Categorical palette for modules/pathways
C_MODULE = ['#1F78B4', '#33A02C', '#FB9A99', '#E31A1C', '#FF7F00',
            '#6A3D9A', '#B15928', '#A6CEE3', '#B2DF8A', '#FDBF6F']

# Diverging colormaps
CMAP_DIVERGING = sns.diverging_palette(240, 10, as_cmap=True)  # blue-red
CMAP_CORRELATION = sns.diverging_palette(250, 15, s=75, l=40, as_cmap=True)

# ============================================================
# Figure Dimensions
# ============================================================
# JASB: single-column = 8.5 cm, double-column = 17 cm
FIG_W_SINGLE = 3.35    # inches (~8.5 cm)
FIG_W_DOUBLE = 6.69    # inches (~17 cm)
FIG_W_1_5 = 5.0        # 1.5 column

# ============================================================
# Global RC Params
# ============================================================
def apply_jasb_style():
    """Apply JASB-compatible matplotlib rcParams globally."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 8,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 7.5,
        'ytick.labelsize': 7.5,
        'legend.fontsize': 7,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.05,
        'axes.linewidth': 0.6,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'grid.alpha': 0.3,
        'grid.linewidth': 0.3,
    })


def panel_label(ax, label, x=-0.08, y=1.05, fontsize=11, fontweight='bold', **kwargs):
    """Add a journal-style panel label (a), (b), etc. to an axes."""
    ax.text(x, y, f'({label})', transform=ax.transAxes,
            fontsize=fontsize, fontweight=fontweight, va='bottom', ha='left',
            **kwargs)


def styled_legend(ax, handles=None, labels=None, loc='best', frameon=True,
                  fancybox=False, edgecolor='white', **kwargs):
    """Create a clean, minimal legend."""
    if handles is None:
        legend = ax.legend(loc=loc, frameon=frameon, fancybox=fancybox,
                           edgecolor=edgecolor, **kwargs)
    else:
        legend = ax.legend(handles, labels, loc=loc, frameon=frameon,
                           fancybox=fancybox, edgecolor=edgecolor, **kwargs)
    if legend:
        legend.get_frame().set_linewidth(0.3)
    return legend


def save_figure(fig, name, outdir='figures_final', formats=('pdf', 'png')):
    """Save figure in multiple formats with consistent naming."""
    os.makedirs(outdir, exist_ok=True)
    for fmt in formats:
        path = os.path.join(outdir, f'{name}.{fmt}')
        fig.savefig(path, dpi=300)


# ============================================================
# Colorbar helpers
# ============================================================
def add_colorbar(fig, im, label='', shrink=0.8, aspect=20, pad=0.02):
    """Add a clean colorbar."""
    cbar = fig.colorbar(im, shrink=shrink, aspect=aspect, pad=pad)
    cbar.ax.tick_params(labelsize=6.5, width=0.4)
    cbar.outline.set_linewidth(0.4)
    if label:
        cbar.set_label(label, fontsize=7)
    return cbar


# ============================================================
# Significance annotation
# ============================================================
def pval_stars(p):
    """Return significance stars."""
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'ns'


def fdr_stars(fdr):
    """Return FDR significance stars with value."""
    if fdr < 0.001:
        return f'FDR<0.001***'
    elif fdr < 0.01:
        return f'FDR={fdr:.3f}**'
    elif fdr < 0.05:
        return f'FDR={fdr:.3f}*'
    else:
        return f'FDR={fdr:.3f}'
