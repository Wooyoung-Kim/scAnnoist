"""Integration Visualization Tools

Comprehensive batch integration visualization for scRNA-seq data.
Provides pre/post integration plots, batch mixing analysis, and condition-wise cell type distribution.
"""

from typing import Annotated, Optional, Dict, List
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
from langchain_core.tools import tool


# Global state for storing plots
_integration_plots = {}


def _setup_plot_style():
    """Setup matplotlib style for publication-quality plots."""
    plt.style.use('seaborn-v0_8-whitegrid')
    sns.set_context("paper", font_scale=1.2)
    plt.rcParams['figure.dpi'] = 150
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']


@tool
def plot_integration_before(
    adata_path: Annotated[str, "Path to h5ad file before integration"],
    batch_key: Annotated[str, "Column name for batch information"] = "batch",
    output_dir: Annotated[str, "Output directory for plots"] = "./integration_plots",
) -> Annotated[str, "Summary of integration plots created"]:
    """
    Generate comprehensive plots for data BEFORE batch integration.

    Creates:
    - UMAP colored by batch (shows batch effect)
    - UMAP colored by clusters (if available)
    - Batch composition bar plots
    - PCA plot colored by batch

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = sc.read_h5ad(adata_path)

    if batch_key not in adata.obs.columns:
        return f"Error: batch_key '{batch_key}' not found in adata.obs. Available keys: {list(adata.obs.columns)}"

    # Ensure we have PCA and UMAP
    if 'X_pca' not in adata.obsm.keys():
        sc.tl.pca(adata, n_comps=50)
    if 'X_umap' not in adata.obsm.keys():
        sc.pp.neighbors(adata, use_rep='X_pca')
        sc.tl.umap(adata)

    # Create comprehensive figure
    fig = plt.figure(figsize=(20, 12))

    # 1. UMAP colored by batch
    ax1 = plt.subplot(2, 4, 1)
    sc.pl.umap(adata, color=batch_key, ax=ax1, show=False, frameon=False, title='UMAP by Batch (Before Integration)')
    ax1.set_title('UMAP by Batch (Before Integration)', fontsize=12, fontweight='bold')

    # 2. UMAP colored by leiden clusters (if available)
    if 'leiden' in adata.obs.columns:
        ax2 = plt.subplot(2, 4, 2)
        sc.pl.umap(adata, color='leiden', ax=ax2, show=False, frameon=False, title='UMAP by Cluster')
        ax2.set_title('UMAP by Cluster (Before Integration)', fontsize=12, fontweight='bold')
    else:
        # Run leiden clustering if not present
        sc.tl.leiden(adata, resolution=0.5)
        ax2 = plt.subplot(2, 4, 2)
        sc.pl.umap(adata, color='leiden', ax=ax2, show=False, frameon=False, title='UMAP by Cluster')
        ax2.set_title('UMAP by Cluster (Before Integration)', fontsize=12, fontweight='bold')

    # 3. PCA colored by batch
    ax3 = plt.subplot(2, 4, 3)
    sc.pl.pca(adata, color=batch_key, ax=ax3, show=False, frameon=False, title='PCA by Batch')
    ax3.set_title('PCA by Batch (Before Integration)', fontsize=12, fontweight='bold')

    # 4. Batch composition
    ax4 = plt.subplot(2, 4, 4)
    batch_counts = adata.obs[batch_key].value_counts().sort_index()
    colors_batch = plt.cm.tab20(np.linspace(0, 1, len(batch_counts)))
    bars = ax4.bar(range(len(batch_counts)), batch_counts.values, color=colors_batch, edgecolor='black')
    ax4.set_xticks(range(len(batch_counts)))
    ax4.set_xticklabels(batch_counts.index, rotation=45, ha='right')
    ax4.set_ylabel('Number of Cells')
    ax4.set_title('Cell Count by Batch', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    # Add count labels on bars
    for bar, count in zip(bars, batch_counts.values):
        height = bar.get_height()
        ax4.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}', ha='center', va='bottom', fontsize=9)

    # 5. Batch by cluster heatmap (if leiden available)
    if 'leiden' in adata.obs.columns:
        ax5 = plt.subplot(2, 4, 5)
        cross_tab = pd.crosstab(adata.obs['leiden'], adata.obs[batch_key], normalize='index')
        sns.heatmap(cross_tab, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax5, cbar_kws={'label': 'Proportion'})
        ax5.set_title('Batch Proportion per Cluster', fontsize=12, fontweight='bold')
        ax5.set_xlabel('Batch')
        ax5.set_ylabel('Cluster')

    # 6. UMAP density per batch (separate subplots)
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)

    if n_batches <= 6:
        for idx, batch in enumerate(batches[:6], start=6):
            ax = plt.subplot(2, 4, idx)
            mask = adata.obs[batch_key] == batch
            ax.scatter(adata.obsm['X_umap'][~mask, 0], adata.obsm['X_umap'][~mask, 1],
                      c='lightgray', s=1, alpha=0.3, rasterized=True)
            ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                      c='red', s=5, alpha=0.6, rasterized=True)
            ax.set_title(f'Batch: {batch}', fontsize=10, fontweight='bold')
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.axis('off')

    plt.tight_layout()

    # Save figure
    output_file = output_path / "integration_before.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _integration_plots['integration_before'] = str(output_file)

    summary = f"""
Integration Analysis (Before):
- Total cells: {adata.n_obs:,}
- Number of batches: {n_batches}
- Batch sizes: {dict(batch_counts)}

Plot saved: {output_file}

⚠️  Expect to see batch clustering in UMAP before integration.
"""

    return summary


@tool
def plot_integration_after(
    adata_path: Annotated[str, "Path to h5ad file after integration"],
    batch_key: Annotated[str, "Column name for batch information"] = "batch",
    celltype_key: Annotated[str, "Column name for cell type annotation"] = "cell_type",
    embedding_key: Annotated[str, "Embedding to use (X_scVI, X_harmony, etc.)"] = "X_scVI",
    output_dir: Annotated[str, "Output directory for plots"] = "./integration_plots",
) -> Annotated[str, "Summary of integration plots created"]:
    """
    Generate comprehensive plots for data AFTER batch integration.

    Creates:
    - UMAP colored by batch (should show good mixing)
    - UMAP colored by cell type
    - Batch composition by cell type
    - Integration quality metrics visualization

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = sc.read_h5ad(adata_path)

    if batch_key not in adata.obs.columns:
        return f"Error: batch_key '{batch_key}' not found in adata.obs"

    # Check if embedding exists, if not use X_pca
    if embedding_key not in adata.obsm.keys():
        print(f"Warning: {embedding_key} not found, available: {list(adata.obsm.keys())}")
        if 'X_pca' in adata.obsm.keys():
            embedding_key = 'X_pca'
        else:
            return f"Error: No suitable embedding found in adata.obsm"

    # Compute UMAP from integrated embedding if not already present
    if 'X_umap' not in adata.obsm.keys() or f'X_umap_{embedding_key}' not in adata.obsm.keys():
        sc.pp.neighbors(adata, use_rep=embedding_key)
        sc.tl.umap(adata)

    # Create comprehensive figure
    fig = plt.figure(figsize=(20, 14))

    # 1. UMAP colored by batch
    ax1 = plt.subplot(3, 4, 1)
    sc.pl.umap(adata, color=batch_key, ax=ax1, show=False, frameon=False, title='UMAP by Batch (After Integration)')
    ax1.set_title('UMAP by Batch (After Integration)', fontsize=12, fontweight='bold')

    # 2. UMAP colored by cell type
    if celltype_key in adata.obs.columns:
        ax2 = plt.subplot(3, 4, 2)
        sc.pl.umap(adata, color=celltype_key, ax=ax2, show=False, frameon=False, title='UMAP by Cell Type')
        ax2.set_title('UMAP by Cell Type (After Integration)', fontsize=12, fontweight='bold')
    else:
        print(f"Warning: celltype_key '{celltype_key}' not found. Available: {list(adata.obs.columns)}")
        # Use leiden instead
        if 'leiden' not in adata.obs.columns:
            sc.tl.leiden(adata, resolution=0.5)
        ax2 = plt.subplot(3, 4, 2)
        sc.pl.umap(adata, color='leiden', ax=ax2, show=False, frameon=False, title='UMAP by Cluster')
        ax2.set_title('UMAP by Cluster (After Integration)', fontsize=12, fontweight='bold')
        celltype_key = 'leiden'  # Use leiden as fallback

    # 3. Batch composition
    ax3 = plt.subplot(3, 4, 3)
    batch_counts = adata.obs[batch_key].value_counts().sort_index()
    colors_batch = plt.cm.tab20(np.linspace(0, 1, len(batch_counts)))
    bars = ax3.bar(range(len(batch_counts)), batch_counts.values, color=colors_batch, edgecolor='black')
    ax3.set_xticks(range(len(batch_counts)))
    ax3.set_xticklabels(batch_counts.index, rotation=45, ha='right')
    ax3.set_ylabel('Number of Cells')
    ax3.set_title('Cell Count by Batch', fontsize=12, fontweight='bold')
    ax3.grid(True, alpha=0.3, axis='y')
    for bar, count in zip(bars, batch_counts.values):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}', ha='center', va='bottom', fontsize=9)

    # 4. Cell type composition
    ax4 = plt.subplot(3, 4, 4)
    celltype_counts = adata.obs[celltype_key].value_counts().sort_index()
    colors_celltype = plt.cm.tab20b(np.linspace(0, 1, len(celltype_counts)))
    bars = ax4.barh(range(len(celltype_counts)), celltype_counts.values, color=colors_celltype, edgecolor='black')
    ax4.set_yticks(range(len(celltype_counts)))
    ax4.set_yticklabels(celltype_counts.index, fontsize=9)
    ax4.set_xlabel('Number of Cells')
    ax4.set_title('Cell Count by Cell Type', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='x')
    for bar, count in zip(bars, celltype_counts.values):
        width = bar.get_width()
        ax4.text(width, bar.get_y() + bar.get_height()/2.,
                f'{count:,}', ha='left', va='center', fontsize=8)

    # 5. Batch by cell type heatmap (stacked bar)
    ax5 = plt.subplot(3, 4, 5)
    cross_tab = pd.crosstab(adata.obs[celltype_key], adata.obs[batch_key])
    cross_tab_pct = cross_tab.div(cross_tab.sum(axis=1), axis=0) * 100  # Percentage per cell type

    # Stacked bar plot
    cross_tab_pct.T.plot(kind='bar', stacked=True, ax=ax5, legend=False, colormap='tab20')
    ax5.set_xlabel('Batch')
    ax5.set_ylabel('Percentage (%)')
    ax5.set_title('Cell Type Composition per Batch', fontsize=12, fontweight='bold')
    ax5.set_xticklabels(ax5.get_xticklabels(), rotation=45, ha='right')
    ax5.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax5.grid(True, alpha=0.3, axis='y')

    # 6. Batch by cell type heatmap (absolute counts)
    ax6 = plt.subplot(3, 4, 6)
    sns.heatmap(cross_tab, annot=True, fmt='d', cmap='Blues', ax=ax6, cbar_kws={'label': 'Cell Count'})
    ax6.set_title('Cell Count: Cell Type × Batch', fontsize=12, fontweight='bold')
    ax6.set_xlabel('Batch')
    ax6.set_ylabel('Cell Type')

    # 7. Batch proportion per cell type (normalized by cell type)
    ax7 = plt.subplot(3, 4, 7)
    cross_tab_norm = cross_tab.div(cross_tab.sum(axis=1), axis=0)
    sns.heatmap(cross_tab_norm, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax7, cbar_kws={'label': 'Proportion'})
    ax7.set_title('Batch Proportion per Cell Type', fontsize=12, fontweight='bold')
    ax7.set_xlabel('Batch')
    ax7.set_ylabel('Cell Type')

    # 8. UMAP density per batch (overlay)
    batches = adata.obs[batch_key].unique()
    n_batches = len(batches)

    if n_batches <= 5:
        for idx, batch in enumerate(batches[:5], start=8):
            ax = plt.subplot(3, 4, idx)
            mask = adata.obs[batch_key] == batch
            ax.scatter(adata.obsm['X_umap'][~mask, 0], adata.obsm['X_umap'][~mask, 1],
                      c='lightgray', s=1, alpha=0.2, rasterized=True)
            ax.scatter(adata.obsm['X_umap'][mask, 0], adata.obsm['X_umap'][mask, 1],
                      c='blue', s=5, alpha=0.7, rasterized=True)
            ax.set_title(f'Batch: {batch}', fontsize=10, fontweight='bold')
            ax.set_xlabel('UMAP 1')
            ax.set_ylabel('UMAP 2')
            ax.axis('off')

    plt.tight_layout()

    # Save figure
    output_file = output_path / "integration_after.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _integration_plots['integration_after'] = str(output_file)

    summary = f"""
Integration Analysis (After):
- Total cells: {adata.n_obs:,}
- Number of batches: {n_batches}
- Number of cell types: {len(celltype_counts)}
- Embedding used: {embedding_key}

Batch distribution:
{batch_counts.to_dict()}

Cell type distribution:
{celltype_counts.to_dict()}

Plot saved: {output_file}

✓ Batches should be well-mixed in UMAP after integration.
"""

    return summary


@tool
def plot_condition_celltype_distribution(
    adata_path: Annotated[str, "Path to h5ad file with annotations"],
    condition_key: Annotated[str, "Column name for experimental condition (e.g., 'condition', 'treatment', 'age')"],
    celltype_key: Annotated[str, "Column name for cell type annotation"] = "cell_type",
    output_dir: Annotated[str, "Output directory for plots"] = "./integration_plots",
) -> Annotated[str, "Summary of condition-wise cell type distribution plots"]:
    """
    Create detailed condition-wise cell type distribution plots.

    Useful for analyzing how cell type composition changes across conditions
    (e.g., young vs aged, treated vs control, etc.).

    Creates:
    - Stacked bar plots showing cell type proportions per condition
    - Heatmaps of cell counts and proportions
    - Statistical comparisons between conditions

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata = sc.read_h5ad(adata_path)

    if condition_key not in adata.obs.columns:
        return f"Error: condition_key '{condition_key}' not found in adata.obs. Available: {list(adata.obs.columns)}"
    if celltype_key not in adata.obs.columns:
        return f"Error: celltype_key '{celltype_key}' not found in adata.obs. Available: {list(adata.obs.columns)}"

    # Create cross-tabulation
    cross_tab = pd.crosstab(adata.obs[celltype_key], adata.obs[condition_key])
    cross_tab_pct = cross_tab.div(cross_tab.sum(axis=0), axis=1) * 100  # Percentage per condition

    # Create comprehensive figure
    fig = plt.figure(figsize=(18, 12))

    # 1. Stacked bar plot (percentage per condition)
    ax1 = plt.subplot(2, 3, 1)
    cross_tab_pct.T.plot(kind='bar', stacked=True, ax=ax1, colormap='tab20', edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Percentage (%)', fontsize=12, fontweight='bold')
    ax1.set_title('Cell Type Proportion per Condition', fontsize=14, fontweight='bold')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45, ha='right')
    ax1.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax1.grid(True, alpha=0.3, axis='y')

    # 2. Grouped bar plot (absolute counts)
    ax2 = plt.subplot(2, 3, 2)
    cross_tab.T.plot(kind='bar', ax=ax2, colormap='tab20', edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Cell Count', fontsize=12, fontweight='bold')
    ax2.set_title('Cell Count per Condition', fontsize=14, fontweight='bold')
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=45, ha='right')
    ax2.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')

    # 3. Heatmap (absolute counts)
    ax3 = plt.subplot(2, 3, 3)
    sns.heatmap(cross_tab, annot=True, fmt='d', cmap='Blues', ax=ax3,
                cbar_kws={'label': 'Cell Count'}, linewidths=0.5, linecolor='gray')
    ax3.set_title('Cell Count Heatmap', fontsize=14, fontweight='bold')
    ax3.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Cell Type', fontsize=12, fontweight='bold')

    # 4. Heatmap (percentage per condition)
    ax4 = plt.subplot(2, 3, 4)
    sns.heatmap(cross_tab_pct, annot=True, fmt='.1f', cmap='YlOrRd', ax=ax4,
                cbar_kws={'label': 'Percentage (%)'}, linewidths=0.5, linecolor='gray')
    ax4.set_title('Cell Type Proportion Heatmap (%)', fontsize=14, fontweight='bold')
    ax4.set_xlabel('Condition', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Cell Type', fontsize=12, fontweight='bold')

    # 5. Fold change heatmap (if 2 conditions, compare to first condition as baseline)
    conditions = adata.obs[condition_key].unique()
    if len(conditions) == 2:
        ax5 = plt.subplot(2, 3, 5)
        baseline = conditions[0]
        comparison = conditions[1]

        fold_change = cross_tab_pct[comparison] / (cross_tab_pct[baseline] + 1e-6)  # Avoid division by zero
        fold_change_log2 = np.log2(fold_change)

        # Create dataframe for heatmap
        fc_df = pd.DataFrame({
            f'{comparison} vs {baseline} (log2 FC)': fold_change_log2
        })

        sns.heatmap(fc_df, annot=True, fmt='.2f', cmap='RdBu_r', center=0, ax=ax5,
                   cbar_kws={'label': 'log2(Fold Change)'}, linewidths=0.5, linecolor='gray')
        ax5.set_title(f'Cell Type Enrichment\n{comparison} vs {baseline}', fontsize=14, fontweight='bold')
        ax5.set_ylabel('Cell Type', fontsize=12, fontweight='bold')
        ax5.set_xlabel('')

    elif len(conditions) > 2:
        # For multiple conditions, show difference from mean
        ax5 = plt.subplot(2, 3, 5)
        mean_pct = cross_tab_pct.mean(axis=1)
        deviation = cross_tab_pct.sub(mean_pct, axis=0)

        sns.heatmap(deviation, annot=True, fmt='.1f', cmap='RdBu_r', center=0, ax=ax5,
                   cbar_kws={'label': 'Deviation from Mean (%)'}, linewidths=0.5, linecolor='gray')
        ax5.set_title('Deviation from Mean Proportion', fontsize=14, fontweight='bold')
        ax5.set_xlabel('Condition', fontsize=12, fontweight='bold')
        ax5.set_ylabel('Cell Type', fontsize=12, fontweight='bold')

    # 6. Total cell count per condition
    ax6 = plt.subplot(2, 3, 6)
    condition_counts = adata.obs[condition_key].value_counts().sort_index()
    colors_cond = plt.cm.Set2(np.linspace(0, 1, len(condition_counts)))
    bars = ax6.bar(range(len(condition_counts)), condition_counts.values,
                   color=colors_cond, edgecolor='black', linewidth=1.5)
    ax6.set_xticks(range(len(condition_counts)))
    ax6.set_xticklabels(condition_counts.index, rotation=45, ha='right')
    ax6.set_ylabel('Total Cell Count', fontsize=12, fontweight='bold')
    ax6.set_title('Total Cells per Condition', fontsize=14, fontweight='bold')
    ax6.grid(True, alpha=0.3, axis='y')

    for bar, count in zip(bars, condition_counts.values):
        height = bar.get_height()
        ax6.text(bar.get_x() + bar.get_width()/2., height,
                f'{count:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')

    plt.tight_layout()

    # Save figure
    output_file = output_path / f"condition_{condition_key}_celltype_distribution.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _integration_plots[f'condition_{condition_key}_distribution'] = str(output_file)

    # Generate summary statistics
    summary_text = f"""
Condition-wise Cell Type Distribution Analysis:
- Condition variable: {condition_key}
- Cell type variable: {celltype_key}
- Number of conditions: {len(conditions)}
- Number of cell types: {len(cross_tab)}

Condition distribution:
{condition_counts.to_dict()}

Cell type counts per condition:
{cross_tab.to_string()}

Cell type proportions (%) per condition:
{cross_tab_pct.round(2).to_string()}

Plot saved: {output_file}
"""

    return summary_text


@tool
def plot_integration_comparison(
    adata_before_path: Annotated[str, "Path to h5ad file before integration"],
    adata_after_path: Annotated[str, "Path to h5ad file after integration"],
    batch_key: Annotated[str, "Column name for batch information"] = "batch",
    output_dir: Annotated[str, "Output directory for plots"] = "./integration_plots",
) -> Annotated[str, "Summary of before/after integration comparison"]:
    """
    Create side-by-side comparison of before and after integration.

    Shows the effect of batch correction on data structure and batch mixing.

    Returns paths to saved figures.
    """
    import scanpy as sc

    _setup_plot_style()
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load data
    adata_before = sc.read_h5ad(adata_before_path)
    adata_after = sc.read_h5ad(adata_after_path)

    if batch_key not in adata_before.obs.columns or batch_key not in adata_after.obs.columns:
        return f"Error: batch_key '{batch_key}' not found in one of the datasets"

    # Ensure both have UMAP
    for adata, label in [(adata_before, 'before'), (adata_after, 'after')]:
        if 'X_umap' not in adata.obsm.keys():
            if 'X_pca' not in adata.obsm.keys():
                sc.tl.pca(adata)
            sc.pp.neighbors(adata)
            sc.tl.umap(adata)

    # Create comparison figure
    fig = plt.figure(figsize=(18, 8))

    # Row 1: UMAP by batch
    ax1 = plt.subplot(2, 3, 1)
    sc.pl.umap(adata_before, color=batch_key, ax=ax1, show=False, frameon=False,
              title='Before Integration - Batch')
    ax1.set_title('Before Integration - Batch', fontsize=14, fontweight='bold')

    ax2 = plt.subplot(2, 3, 2)
    sc.pl.umap(adata_after, color=batch_key, ax=ax2, show=False, frameon=False,
              title='After Integration - Batch')
    ax2.set_title('After Integration - Batch', fontsize=14, fontweight='bold')

    # Row 2: UMAP by cluster (if available)
    if 'leiden' in adata_before.obs.columns and 'leiden' in adata_after.obs.columns:
        ax3 = plt.subplot(2, 3, 4)
        sc.pl.umap(adata_before, color='leiden', ax=ax3, show=False, frameon=False,
                  title='Before Integration - Cluster')
        ax3.set_title('Before Integration - Cluster', fontsize=14, fontweight='bold')

        ax4 = plt.subplot(2, 3, 5)
        sc.pl.umap(adata_after, color='leiden', ax=ax4, show=False, frameon=False,
                  title='After Integration - Cluster')
        ax4.set_title('After Integration - Cluster', fontsize=14, fontweight='bold')

    # Summary comparison (text)
    ax5 = plt.subplot(2, 3, 3)
    ax5.axis('off')

    batch_counts_before = adata_before.obs[batch_key].value_counts().sort_index()
    batch_counts_after = adata_after.obs[batch_key].value_counts().sort_index()

    summary_lines = [
        "Integration Summary",
        "=" * 40,
        f"Cells before: {adata_before.n_obs:,}",
        f"Cells after: {adata_after.n_obs:,}",
        "",
        "Batch distribution (before):",
    ]
    for batch, count in batch_counts_before.items():
        summary_lines.append(f"  {batch}: {count:,}")

    summary_lines.append("")
    summary_lines.append("Batch distribution (after):")
    for batch, count in batch_counts_after.items():
        summary_lines.append(f"  {batch}: {count:,}")

    ax5.text(0.1, 0.9, '\n'.join(summary_lines), transform=ax5.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

    # Batch mixing score comparison (if multiple batches)
    batches = adata_before.obs[batch_key].unique()
    if len(batches) > 1:
        ax6 = plt.subplot(2, 3, 6)

        # Simple mixing score: average minimum distance to other batch in UMAP space
        # Higher score = better mixing
        from scipy.spatial.distance import cdist

        def compute_mixing_score(adata, batch_key):
            """Compute simple batch mixing score in UMAP space."""
            scores = []
            for batch in adata.obs[batch_key].unique():
                batch_mask = adata.obs[batch_key] == batch
                batch_cells = adata.obsm['X_umap'][batch_mask]
                other_cells = adata.obsm['X_umap'][~batch_mask]

                if len(batch_cells) > 0 and len(other_cells) > 0:
                    # Sample if too many cells (for speed)
                    if len(batch_cells) > 1000:
                        idx = np.random.choice(len(batch_cells), 1000, replace=False)
                        batch_cells = batch_cells[idx]
                    if len(other_cells) > 1000:
                        idx = np.random.choice(len(other_cells), 1000, replace=False)
                        other_cells = other_cells[idx]

                    # Compute minimum distance to other batch
                    dists = cdist(batch_cells, other_cells)
                    min_dists = dists.min(axis=1)
                    scores.append(np.mean(min_dists))

            return np.mean(scores) if scores else 0

        mixing_before = compute_mixing_score(adata_before, batch_key)
        mixing_after = compute_mixing_score(adata_after, batch_key)

        categories = ['Before\nIntegration', 'After\nIntegration']
        mixing_scores = [mixing_before, mixing_after]
        colors = ['lightcoral', 'lightgreen']

        bars = ax6.bar(categories, mixing_scores, color=colors, edgecolor='black', linewidth=2)
        ax6.set_ylabel('Avg Min Distance to Other Batch', fontsize=11, fontweight='bold')
        ax6.set_title('Batch Mixing Score\n(Lower = Better Mixing)', fontsize=12, fontweight='bold')
        ax6.grid(True, alpha=0.3, axis='y')

        for bar, score in zip(bars, mixing_scores):
            height = bar.get_height()
            ax6.text(bar.get_x() + bar.get_width()/2., height,
                    f'{score:.3f}', ha='center', va='bottom', fontsize=11, fontweight='bold')

        # Add interpretation
        improvement = ((mixing_before - mixing_after) / mixing_before) * 100 if mixing_before > 0 else 0
        ax6.text(0.5, 0.5, f'Improvement: {improvement:.1f}%',
                transform=ax6.transAxes, ha='center', fontsize=11,
                bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.5))

    plt.tight_layout()

    # Save figure
    output_file = output_path / "integration_before_after_comparison.png"
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    plt.close()

    _integration_plots['integration_comparison'] = str(output_file)

    summary = f"""
Integration Before/After Comparison:
- Cells: {adata_before.n_obs:,} → {adata_after.n_obs:,}
- Batches: {len(batches)}

Plot saved: {output_file}

Visual inspection:
- Before integration: Expect batch clustering
- After integration: Expect batch mixing with preserved biological structure
"""

    return summary


@tool
def get_integration_plot_paths() -> Annotated[Dict[str, str], "Dictionary of plot names and file paths"]:
    """
    Get paths to all integration plots that have been generated.

    Returns:
        Dictionary mapping plot names to file paths.
    """
    return _integration_plots.copy()
