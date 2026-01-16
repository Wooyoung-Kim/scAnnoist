"""QC Visualization Tools

Comprehensive quality control visualization for scRNA-seq data.
Provides pre/post QC plots and comparison visualizations.
"""

from typing import Annotated, Optional, Dict, List
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import scanpy as sc
from langchain_core.tools import tool


# Global state for storing plots
_qc_plots = {}


def _setup_plot_style():
    """Setup matplotlib style for publication-quality plots."""
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']


@tool
def plot_qc_metrics_before(
    adata_path: Annotated[str, "Path to h5ad file before QC"],
    output_dir: Annotated[str, "Output directory for plots"] = "./qc_plots",
) -> Annotated[str, "Summary of QC plots created"]:
    """
    Generate comprehensive QC plots for raw data BEFORE filtering.

    Creates:
    - Violin plots for n_genes, n_counts, pct_counts_mt
    - Scatter plots showing relationships
    - Histogram distributions
    - Summary statistics

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = sc.read_h5ad(adata_path)

    # Calculate QC metrics if not already present
    if 'n_genes' not in adata.obs.columns:
        adata.var_names_make_unique()
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
    if 'n_counts' not in adata.obs.columns:
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)

    # Calculate MT percentage if not present
    if 'pct_counts_mt' not in adata.obs.columns:
        adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    n_cells_before = adata.n_obs

    # Create comprehensive QC plots
    fig = plt.figure(figsize=(20, 12))

    # 1. Violin plots
    ax1 = plt.subplot(3, 4, 1)
    sc.pl.violin(adata, ['n_genes'], jitter=0.4, ax=ax1, show=False)
    ax1.set_title('Genes per Cell (Before QC)', fontsize=12, fontweight='bold')

    ax2 = plt.subplot(3, 4, 2)
    sc.pl.violin(adata, ['n_counts'], jitter=0.4, ax=ax2, show=False)
    ax2.set_title('UMI Counts per Cell (Before QC)', fontsize=12, fontweight='bold')

    ax3 = plt.subplot(3, 4, 3)
    sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, ax=ax3, show=False)
    ax3.set_title('MT% per Cell (Before QC)', fontsize=12, fontweight='bold')

    # 2. Histograms with statistics
    ax4 = plt.subplot(3, 4, 4)
    n_genes = adata.obs['n_genes'].values
    ax4.hist(n_genes, bins=50, edgecolor='black', alpha=0.7)
    ax4.axvline(np.median(n_genes), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(n_genes):.0f}')
    ax4.axvline(np.mean(n_genes), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(n_genes):.0f}')
    ax4.set_xlabel('Genes per Cell')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution: Genes per Cell', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    ax5 = plt.subplot(3, 4, 5)
    n_counts = adata.obs['n_counts'].values
    ax5.hist(n_counts, bins=50, edgecolor='black', alpha=0.7, color='orange')
    ax5.axvline(np.median(n_counts), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(n_counts):.0f}')
    ax5.axvline(np.mean(n_counts), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(n_counts):.0f}')
    ax5.set_xlabel('UMI Counts per Cell')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Distribution: UMI Counts', fontsize=12, fontweight='bold')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    ax6 = plt.subplot(3, 4, 6)
    mt_pct = adata.obs['pct_counts_mt'].values
    ax6.hist(mt_pct, bins=50, edgecolor='black', alpha=0.7, color='green')
    ax6.axvline(np.median(mt_pct), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(mt_pct):.1f}%')
    ax6.axvline(np.mean(mt_pct), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(mt_pct):.1f}%')
    ax6.set_xlabel('Mitochondrial %')
    ax6.set_ylabel('Frequency')
    ax6.set_title('Distribution: MT%', fontsize=12, fontweight='bold')
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    # 3. Scatter plots showing relationships
    ax7 = plt.subplot(3, 4, 7)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', ax=ax7, show=False)
    ax7.set_title('Counts vs Genes (Before QC)', fontsize=12, fontweight='bold')

    ax8 = plt.subplot(3, 4, 8)
    sc.pl.scatter(adata, x='n_counts', y='pct_counts_mt', ax=ax8, show=False)
    ax8.set_title('Counts vs MT% (Before QC)', fontsize=12, fontweight='bold')

    ax9 = plt.subplot(3, 4, 9)
    sc.pl.scatter(adata, x='n_genes', y='pct_counts_mt', ax=ax9, show=False)
    ax9.set_title('Genes vs MT% (Before QC)', fontsize=12, fontweight='bold')

    # 4. Box plots for batch comparison (if batch key exists)
    if 'batch' in adata.obs.columns or 'sample' in adata.obs.columns:
        batch_key = 'batch' if 'batch' in adata.obs.columns else 'sample'

        ax10 = plt.subplot(3, 4, 10)
        adata.obs.boxplot(column='n_genes', by=batch_key, ax=ax10)
        ax10.set_title(f'Genes per Cell by {batch_key}', fontsize=12, fontweight='bold')
        ax10.set_xlabel(batch_key)
        plt.suptitle('')  # Remove default title

        ax11 = plt.subplot(3, 4, 11)
        adata.obs.boxplot(column='n_counts', by=batch_key, ax=ax11)
        ax11.set_title(f'UMI Counts by {batch_key}', fontsize=12, fontweight='bold')
        ax11.set_xlabel(batch_key)
        plt.suptitle('')

        ax12 = plt.subplot(3, 4, 12)
        adata.obs.boxplot(column='pct_counts_mt', by=batch_key, ax=ax12)
        ax12.set_title(f'MT% by {batch_key}', fontsize=12, fontweight='bold')
        ax12.set_xlabel(batch_key)
        plt.suptitle('')

    plt.tight_layout()

    # Save figure
    output_file = output_path / "qc_before_filtering.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _qc_plots['qc_before'] = str(output_file)

    # Summary statistics
    summary = f"""
QC Metrics (Before Filtering):
- Total cells: {n_cells_before:,}
- Genes per cell: {np.median(n_genes):.0f} (median), {np.mean(n_genes):.0f} (mean)
- UMI counts per cell: {np.median(n_counts):.0f} (median), {np.mean(n_counts):.0f} (mean)
- MT% per cell: {np.median(mt_pct):.2f}% (median), {np.mean(mt_pct):.2f}% (mean)

Plot saved: {output_file}
"""

    return summary


@tool
def plot_qc_metrics_after(
    adata_path: Annotated[str, "Path to h5ad file after QC"],
    output_dir: Annotated[str, "Output directory for plots"] = "./qc_plots",
) -> Annotated[str, "Summary of QC plots created"]:
    """
    Generate comprehensive QC plots for data AFTER filtering.

    Creates the same plots as plot_qc_metrics_before for comparison.

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = sc.read_h5ad(adata_path)

    n_cells_after = adata.n_obs

    # Ensure QC metrics are present
    if 'n_genes' not in adata.obs.columns:
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
    if 'n_counts' not in adata.obs.columns:
        adata.obs['n_counts'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)
    if 'pct_counts_mt' not in adata.obs.columns:
        adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('mt-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Create comprehensive QC plots (same structure as before)
    fig = plt.figure(figsize=(20, 12))

    # 1. Violin plots
    ax1 = plt.subplot(3, 4, 1)
    sc.pl.violin(adata, ['n_genes'], jitter=0.4, ax=ax1, show=False)
    ax1.set_title('Genes per Cell (After QC)', fontsize=12, fontweight='bold')

    ax2 = plt.subplot(3, 4, 2)
    sc.pl.violin(adata, ['n_counts'], jitter=0.4, ax=ax2, show=False)
    ax2.set_title('UMI Counts per Cell (After QC)', fontsize=12, fontweight='bold')

    ax3 = plt.subplot(3, 4, 3)
    sc.pl.violin(adata, ['pct_counts_mt'], jitter=0.4, ax=ax3, show=False)
    ax3.set_title('MT% per Cell (After QC)', fontsize=12, fontweight='bold')

    # 2. Histograms
    ax4 = plt.subplot(3, 4, 4)
    n_genes = adata.obs['n_genes'].values
    ax4.hist(n_genes, bins=50, edgecolor='black', alpha=0.7)
    ax4.axvline(np.median(n_genes), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(n_genes):.0f}')
    ax4.axvline(np.mean(n_genes), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(n_genes):.0f}')
    ax4.set_xlabel('Genes per Cell')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Distribution: Genes per Cell', fontsize=12, fontweight='bold')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    ax5 = plt.subplot(3, 4, 5)
    n_counts = adata.obs['n_counts'].values
    ax5.hist(n_counts, bins=50, edgecolor='black', alpha=0.7, color='orange')
    ax5.axvline(np.median(n_counts), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(n_counts):.0f}')
    ax5.axvline(np.mean(n_counts), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(n_counts):.0f}')
    ax5.set_xlabel('UMI Counts per Cell')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Distribution: UMI Counts', fontsize=12, fontweight='bold')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    ax6 = plt.subplot(3, 4, 6)
    mt_pct = adata.obs['pct_counts_mt'].values
    ax6.hist(mt_pct, bins=50, edgecolor='black', alpha=0.7, color='green')
    ax6.axvline(np.median(mt_pct), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(mt_pct):.1f}%')
    ax6.axvline(np.mean(mt_pct), color='blue', linestyle='--', linewidth=2, label=f'Mean: {np.mean(mt_pct):.1f}%')
    ax6.set_xlabel('Mitochondrial %')
    ax6.set_ylabel('Frequency')
    ax6.set_title('Distribution: MT%', fontsize=12, fontweight='bold')
    ax6.legend()
    ax6.grid(True, alpha=0.3)

    # 3. Scatter plots
    ax7 = plt.subplot(3, 4, 7)
    sc.pl.scatter(adata, x='n_counts', y='n_genes', ax=ax7, show=False)
    ax7.set_title('Counts vs Genes (After QC)', fontsize=12, fontweight='bold')

    ax8 = plt.subplot(3, 4, 8)
    sc.pl.scatter(adata, x='n_counts', y='pct_counts_mt', ax=ax8, show=False)
    ax8.set_title('Counts vs MT% (After QC)', fontsize=12, fontweight='bold')

    ax9 = plt.subplot(3, 4, 9)
    sc.pl.scatter(adata, x='n_genes', y='pct_counts_mt', ax=ax9, show=False)
    ax9.set_title('Genes vs MT% (After QC)', fontsize=12, fontweight='bold')

    # 4. Box plots for batch comparison
    if 'batch' in adata.obs.columns or 'sample' in adata.obs.columns:
        batch_key = 'batch' if 'batch' in adata.obs.columns else 'sample'

        ax10 = plt.subplot(3, 4, 10)
        adata.obs.boxplot(column='n_genes', by=batch_key, ax=ax10)
        ax10.set_title(f'Genes per Cell by {batch_key}', fontsize=12, fontweight='bold')
        ax10.set_xlabel(batch_key)
        plt.suptitle('')

        ax11 = plt.subplot(3, 4, 11)
        adata.obs.boxplot(column='n_counts', by=batch_key, ax=ax11)
        ax11.set_title(f'UMI Counts by {batch_key}', fontsize=12, fontweight='bold')
        ax11.set_xlabel(batch_key)
        plt.suptitle('')

        ax12 = plt.subplot(3, 4, 12)
        adata.obs.boxplot(column='pct_counts_mt', by=batch_key, ax=ax12)
        ax12.set_title(f'MT% by {batch_key}', fontsize=12, fontweight='bold')
        ax12.set_xlabel(batch_key)
        plt.suptitle('')

    plt.tight_layout()

    # Save figure
    output_file = output_path / "qc_after_filtering.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _qc_plots['qc_after'] = str(output_file)

    # Summary statistics
    summary = f"""
QC Metrics (After Filtering):
- Total cells: {n_cells_after:,}
- Genes per cell: {np.median(n_genes):.0f} (median), {np.mean(n_genes):.0f} (mean)
- UMI counts per cell: {np.median(n_counts):.0f} (median), {np.mean(n_counts):.0f} (mean)
- MT% per cell: {np.median(mt_pct):.2f}% (median), {np.mean(mt_pct):.2f}% (mean)

Plot saved: {output_file}
"""

    return summary


@tool
def plot_qc_comparison(
    adata_before_path: Annotated[str, "Path to h5ad file before QC"],
    adata_after_path: Annotated[str, "Path to h5ad file after QC"],
    output_dir: Annotated[str, "Output directory for plots"] = "./qc_plots",
) -> Annotated[str, "Summary of comparison plots created"]:
    """
    Create side-by-side comparison plots of Before vs After QC.

    Shows the effect of QC filtering on data quality metrics.

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata_before = sc.read_h5ad(adata_before_path)
    adata_after = sc.read_h5ad(adata_after_path)

    # Ensure metrics are calculated
    for adata in [adata_before, adata_after]:
        if 'n_genes' not in adata.obs.columns:
            adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
        if 'n_counts' not in adata.obs.columns:
            adata.obs['n_counts'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)
        if 'pct_counts_mt' not in adata.obs.columns:
            adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('Mt-') | adata.var_names.str.startswith('mt-')
            sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    n_cells_before = adata_before.n_obs
    n_cells_after = adata_after.n_obs
    pct_retained = (n_cells_after / n_cells_before) * 100

    # Create comparison figure
    fig = plt.figure(figsize=(18, 10))

    # 1. Cell count comparison
    ax1 = plt.subplot(2, 4, 1)
    categories = ['Before QC', 'After QC']
    cell_counts = [n_cells_before, n_cells_after]
    colors = ['lightcoral', 'lightgreen']
    bars = ax1.bar(categories, cell_counts, color=colors, edgecolor='black', linewidth=1.5)
    ax1.set_ylabel('Number of Cells')
    ax1.set_title('Cell Count Comparison', fontsize=12, fontweight='bold')
    for bar, count in zip(bars, cell_counts):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}',
                ha='center', va='bottom', fontweight='bold')
    ax1.text(0.5, max(cell_counts)*0.5, f'{pct_retained:.1f}% retained',
            ha='center', fontsize=14, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))
    ax1.grid(True, alpha=0.3, axis='y')

    # 2-4. Violin plots comparison
    metrics = ['n_genes', 'n_counts', 'pct_counts_mt']
    titles = ['Genes per Cell', 'UMI Counts per Cell', 'MT% per Cell']

    for idx, (metric, title) in enumerate(zip(metrics, titles), start=2):
        ax = plt.subplot(2, 4, idx)

        # Combine data for violin plot
        before_data = adata_before.obs[metric].values
        after_data = adata_after.obs[metric].values

        data_combined = [before_data, after_data]
        positions = [1, 2]

        parts = ax.violinplot(data_combined, positions=positions, widths=0.7,
                             showmeans=True, showextrema=True, showmedians=True)

        for pc in parts['bodies']:
            pc.set_facecolor('lightblue')
            pc.set_alpha(0.7)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(['Before', 'After'])
        ax.set_ylabel(metric.replace('_', ' ').title())
        ax.set_title(f'{title} Comparison', fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, axis='y')

        # Add median values as text
        median_before = np.median(before_data)
        median_after = np.median(after_data)
        ax.text(1, median_before, f'{median_before:.1f}', ha='right', va='center',
               fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        ax.text(2, median_after, f'{median_after:.1f}', ha='left', va='center',
               fontweight='bold', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # 5-7. Histogram overlays
    for idx, (metric, title) in enumerate(zip(metrics, titles), start=5):
        ax = plt.subplot(2, 4, idx)

        before_data = adata_before.obs[metric].values
        after_data = adata_after.obs[metric].values

        ax.hist(before_data, bins=50, alpha=0.5, label='Before QC', color='red', edgecolor='black')
        ax.hist(after_data, bins=50, alpha=0.5, label='After QC', color='green', edgecolor='black')
        ax.set_xlabel(metric.replace('_', ' ').title())
        ax.set_ylabel('Frequency')
        ax.set_title(f'{title} Distribution', fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # 8. Summary statistics table
    ax8 = plt.subplot(2, 4, 8)
    ax8.axis('tight')
    ax8.axis('off')

    table_data = [
        ['Metric', 'Before QC', 'After QC', 'Change'],
        ['Cells', f'{n_cells_before:,}', f'{n_cells_after:,}', f'-{n_cells_before-n_cells_after:,}'],
        ['Median Genes', f'{np.median(adata_before.obs["n_genes"]):.0f}',
         f'{np.median(adata_after.obs["n_genes"]):.0f}', ''],
        ['Median Counts', f'{np.median(adata_before.obs["n_counts"]):.0f}',
         f'{np.median(adata_after.obs["n_counts"]):.0f}', ''],
        ['Median MT%', f'{np.median(adata_before.obs["pct_counts_mt"]):.2f}%',
         f'{np.median(adata_after.obs["pct_counts_mt"]):.2f}%', ''],
    ]

    table = ax8.table(cellText=table_data, cellLoc='center', loc='center',
                     colWidths=[0.3, 0.25, 0.25, 0.2])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)

    # Style header row
    for i in range(4):
        table[(0, i)].set_facecolor('#4CAF50')
        table[(0, i)].set_text_props(weight='bold', color='white')

    ax8.set_title('QC Summary Statistics', fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()

    # Save figure
    output_file = output_path / "qc_before_after_comparison.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _qc_plots['qc_comparison'] = str(output_file)

    summary = f"""
QC Comparison Summary:
- Cells before QC: {n_cells_before:,}
- Cells after QC: {n_cells_after:,}
- Cells removed: {n_cells_before - n_cells_after:,} ({100-pct_retained:.1f}%)
- Cells retained: {pct_retained:.1f}%

Plot saved: {output_file}

All QC plots generated:
{chr(10).join(f'- {k}: {v}' for k, v in _qc_plots.items())}
"""

    return summary


@tool
def get_qc_plot_paths() -> Annotated[Dict[str, str], "Dictionary of plot names and file paths"]:
    """
    Get paths to all QC plots that have been generated.

    Returns:
        Dictionary mapping plot names to file paths.
    """
    return _qc_plots.copy()
