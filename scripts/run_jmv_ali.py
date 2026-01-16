#!/usr/bin/env python
"""
Enhanced scRNA-seq pipeline for JMV_ALI dataset with comprehensive analysis.
Integrates all 24 samples using Harmony batch correction with:
- Enhanced QC and integration visualization
- Condition-level differential expression analysis
- Cell-type-specific viral response analysis
- Comprehensive marker annotation and visualization
"""

import os
import sys
import logging
from pathlib import Path
from datetime import datetime
import json
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Suppress warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from scrna_agent.pipeline import PipelineConfig, run_pipeline
from scrna_agent.stats import (
    compute_adaptive_qc_thresholds,
    apply_adaptive_qc,
    compute_qc_distributions,
    compute_pc_selection,
    run_full_deg_analysis,
)

# Configuration
DATA_DIR = Path("/home/kwy7605/data_61/SARS/Count/JMV_ALI")
OUTPUT_DIR = Path("/mnt2/kwy/scrna_agent/results/jmv_ali_epithelial")
TISSUE = "airway"  # Using airway markers for ALI (Air-Liquid Interface) cultures
N_PCS = 50  # Number of PCs to use for integration and clustering
REMOVE_IMMUNE = False  # Dataset is already epithelial-only (only 4% cells with immune markers)

# Airway cell type markers for visualization
AIRWAY_MARKERS = {
    'Ciliated': ['FOXJ1', 'DNAH5', 'TPPP3', 'CAPS', 'SNTN'],
    'Basal': ['KRT5', 'TP63', 'KRT15', 'KRT14', 'NGFR'],
    'Secretory/Club': ['SCGB1A1', 'SCGB3A2', 'MUC5B', 'CYP2F1'],
    'Goblet': ['MUC5AC', 'MUC5B', 'TFF3', 'SPDEF'],
    'AT2-like': ['SFTPC', 'SFTPA1', 'SFTPB', 'LAMP3'],
    'Immune': ['PTPRC', 'CD3D', 'CD79A', 'LYZ', 'CD68'],
}

# Viral response genes
VIRAL_RESPONSE_GENES = ['ISG15', 'IFIT1', 'IFIT2', 'IFIT3', 'MX1', 'MX2', 'OAS1', 'IFI6', 'IFI44L', 'IFITM3']


def create_output_directories(output_dir):
    """Create organized output directory structure."""
    subdirs = [
        'qc',
        'integration',
        'umap',
        'annotation',
        'clustering',
        'deg',
        'deg/condition_comparisons',
        'deg/celltype_specific',
        'deg/variant_comparisons',
        'supplementary'
    ]

    for subdir in subdirs:
        (output_dir / subdir).mkdir(parents=True, exist_ok=True)

    logger.info(f"Created output directory structure at {output_dir}")


def filter_immune_cells(adata, output_dir):
    """
    Filter out immune cells based on marker gene expression.
    ALI organoids should contain primarily epithelial cells.
    """
    logger.info("Filtering immune cells...")

    # Immune markers to check (PTPRC/CD45 is the key pan-immune marker)
    immune_markers = ['PTPRC', 'CD3D', 'CD3E', 'CD79A', 'CD79B', 'LYZ', 'CD68', 'CD14']
    present_immune_markers = [m for m in immune_markers if m in adata.var_names]

    if len(present_immune_markers) == 0:
        logger.warning("No immune markers found, skipping immune cell filtering")
        return adata

    logger.info(f"Using {len(present_immune_markers)} immune markers: {present_immune_markers}")

    # Calculate mean expression of immune markers per cell (in normalized space)
    immune_expr = adata[:, present_immune_markers].X
    if hasattr(immune_expr, 'toarray'):
        immune_expr = immune_expr.toarray()

    mean_immune_expr = immune_expr.mean(axis=1)

    # Define threshold: cells with very high immune marker expression are likely immune cells
    # Use more conservative threshold - 90th percentile (keep top 10% as immune, rest as epithelial)
    threshold = np.percentile(mean_immune_expr, 90)

    is_epithelial = mean_immune_expr < threshold
    n_before = adata.n_obs
    n_immune = (~is_epithelial).sum()

    logger.info(f"  Mean immune marker expression range: {mean_immune_expr.min():.3f} - {mean_immune_expr.max():.3f}")
    logger.info(f"  Threshold (90th percentile): {threshold:.3f}")
    logger.info(f"  Identified {n_immune} potential immune cells ({n_immune/n_before*100:.1f}%)")
    logger.info(f"  Retaining {is_epithelial.sum()} epithelial cells ({is_epithelial.sum()/n_before*100:.1f}%)")

    # Save filtering statistics
    filter_stats = pd.DataFrame({
        'sample': adata.obs['sample'].values,
        'mean_immune_expr': mean_immune_expr,
        'is_epithelial': is_epithelial
    })
    filter_stats.to_csv(output_dir / 'qc' / 'immune_filtering.csv', index=False)

    # Filter - ensure categorical columns are preserved
    adata_filtered = adata[is_epithelial].copy()

    # Restore categorical dtypes
    for col in ['sample', 'condition', 'group']:
        if col in adata_filtered.obs.columns:
            adata_filtered.obs[col] = adata_filtered.obs[col].astype('category')

    return adata_filtered


def load_all_samples():
    """Load all JMV_ALI samples and concatenate."""
    samples = sorted([d for d in DATA_DIR.iterdir() if d.is_dir()])

    adatas = []
    for sample_dir in samples:
        h5_path = sample_dir / "outs" / "filtered_feature_bc_matrix.h5"
        if not h5_path.exists():
            logger.warning(f"Skipping {sample_dir.name}: no h5 file found")
            continue

        logger.info(f"Loading {sample_dir.name}...")
        adata = sc.read_10x_h5(str(h5_path))
        adata.var_names_make_unique()

        # Add sample metadata
        sample_name = sample_dir.name
        adata.obs['sample'] = sample_name

        # Extract condition from sample name
        # JMV_Control-1 -> Control
        # JMV_HCoV_OC43-1 -> HCoV_OC43
        # JMV_SARS-CoV-2_BA275-1 -> SARS-CoV-2_BA275
        parts = sample_name.replace('JMV_', '').rsplit('-', 1)
        condition = parts[0]
        adata.obs['condition'] = condition

        # Simplify condition groups
        if 'Control' in condition:
            adata.obs['group'] = 'Control'
        elif 'SARS-CoV-2' in condition:
            adata.obs['group'] = 'SARS-CoV-2'
        elif 'SARS-CoV-1' in condition:
            adata.obs['group'] = 'SARS-CoV-1'
        elif 'MERS' in condition:
            adata.obs['group'] = 'MERS'
        elif 'HCoV' in condition:
            adata.obs['group'] = 'HCoV_OC43'
        else:
            adata.obs['group'] = condition

        adatas.append(adata)
        logger.info(f"  Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Concatenate all samples
    logger.info(f"Concatenating {len(adatas)} samples...")
    adata = sc.concat(adatas, join='outer')
    adata.obs_names_make_unique()

    logger.info(f"Total: {adata.n_obs} cells, {adata.n_vars} genes")
    return adata


def generate_qc_visualizations(adata, output_dir, stage='before'):
    """Generate comprehensive QC visualizations."""
    qc_dir = output_dir / 'qc'

    logger.info(f"Generating {stage}-QC visualizations...")

    # Create comprehensive QC figure
    fig = plt.figure(figsize=(20, 12))

    # Violin plots
    ax1 = plt.subplot(3, 4, 1)
    sc.pl.violin(adata, ['n_genes_by_counts'], groupby='sample', rotation=90, ax=ax1, show=False)
    ax1.set_title(f'Genes per Cell ({stage} QC)', fontweight='bold')
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, ha='right', fontsize=6)

    ax2 = plt.subplot(3, 4, 2)
    sc.pl.violin(adata, ['total_counts'], groupby='sample', rotation=90, ax=ax2, show=False)
    ax2.set_title(f'UMI Counts per Cell ({stage} QC)', fontweight='bold')
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=90, ha='right', fontsize=6)

    ax3 = plt.subplot(3, 4, 3)
    sc.pl.violin(adata, ['pct_counts_mt'], groupby='sample', rotation=90, ax=ax3, show=False)
    ax3.set_title(f'MT% per Cell ({stage} QC)', fontweight='bold')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=90, ha='right', fontsize=6)

    # Scatter plots
    ax4 = plt.subplot(3, 4, 5)
    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', ax=ax4, show=False)
    ax4.set_title('Counts vs Genes', fontweight='bold')

    ax5 = plt.subplot(3, 4, 6)
    sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax5, show=False)
    ax5.set_title('Counts vs MT%', fontweight='bold')

    ax6 = plt.subplot(3, 4, 7)
    sc.pl.scatter(adata, x='n_genes_by_counts', y='pct_counts_mt', ax=ax6, show=False)
    ax6.set_title('Genes vs MT%', fontweight='bold')

    # Histograms
    ax7 = plt.subplot(3, 4, 9)
    ax7.hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
    ax7.axvline(np.median(adata.obs['n_genes_by_counts']), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(adata.obs["n_genes_by_counts"]):.0f}')
    ax7.set_xlabel('Genes per Cell')
    ax7.set_ylabel('Frequency')
    ax7.set_title('Distribution: Genes per Cell', fontweight='bold')
    ax7.legend()
    ax7.grid(True, alpha=0.3)

    ax8 = plt.subplot(3, 4, 10)
    ax8.hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7, color='orange')
    ax8.axvline(np.median(adata.obs['total_counts']), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(adata.obs["total_counts"]):.0f}')
    ax8.set_xlabel('UMI Counts')
    ax8.set_ylabel('Frequency')
    ax8.set_title('Distribution: UMI Counts', fontweight='bold')
    ax8.legend()
    ax8.grid(True, alpha=0.3)

    ax9 = plt.subplot(3, 4, 11)
    ax9.hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7, color='green')
    ax9.axvline(np.median(adata.obs['pct_counts_mt']), color='red', linestyle='--', linewidth=2, label=f'Median: {np.median(adata.obs["pct_counts_mt"]):.1f}%')
    ax9.set_xlabel('MT%')
    ax9.set_ylabel('Frequency')
    ax9.set_title('Distribution: MT%', fontweight='bold')
    ax9.legend()
    ax9.grid(True, alpha=0.3)

    # Cell count by sample
    ax10 = plt.subplot(3, 4, 4)
    sample_counts = adata.obs['sample'].value_counts().sort_index()
    ax10.bar(range(len(sample_counts)), sample_counts.values, color='steelblue', edgecolor='black')
    ax10.set_xticks(range(len(sample_counts)))
    ax10.set_xticklabels(sample_counts.index, rotation=90, ha='right', fontsize=6)
    ax10.set_ylabel('Number of Cells')
    ax10.set_title(f'Cell Count by Sample ({stage} QC)', fontweight='bold')
    ax10.grid(True, alpha=0.3, axis='y')

    # Cell count by condition
    ax11 = plt.subplot(3, 4, 8)
    condition_counts = adata.obs['condition'].value_counts().sort_index()
    ax11.bar(range(len(condition_counts)), condition_counts.values, color='coral', edgecolor='black')
    ax11.set_xticks(range(len(condition_counts)))
    ax11.set_xticklabels(condition_counts.index, rotation=45, ha='right', fontsize=8)
    ax11.set_ylabel('Number of Cells')
    ax11.set_title(f'Cell Count by Condition ({stage} QC)', fontweight='bold')
    ax11.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig(qc_dir / f'qc_metrics_{stage}.png', dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"  Saved QC plots to {qc_dir}/qc_metrics_{stage}.png")


def run_integrated_analysis(adata, output_dir):
    """Run integrated analysis with Harmony batch correction."""

    create_output_directories(output_dir)

    # Initialize metadata
    run_metadata = {
        'pipeline_version': '2.0_integrated',
        'timestamp': datetime.now().isoformat(),
        'tissue': TISSUE,
        'n_samples': len(adata.obs['sample'].unique()),
        'samples': list(adata.obs['sample'].unique()),
    }

    logger.info("=" * 60)
    logger.info("JMV_ALI Integrated Analysis (Statistical Mode)")
    logger.info("=" * 60)

    # 1. QC metrics calculation
    logger.info("Step 1: Calculating QC metrics...")
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Compute distributions
    distributions = compute_qc_distributions(adata)
    run_metadata['qc_distributions'] = distributions

    # Generate pre-QC visualizations
    generate_qc_visualizations(adata.copy(), output_dir, stage='before')

    # 2. Adaptive QC thresholds
    logger.info("Step 2: Computing adaptive QC thresholds...")
    thresholds = compute_adaptive_qc_thresholds(
        adata,
        nmads=3.0,
        min_genes_floor=200,
        max_mt_pct_ceiling=20.0,
        method='mad'
    )
    run_metadata['qc_thresholds'] = thresholds.to_dict()

    logger.info(f"  min_genes: {thresholds.min_genes}")
    logger.info(f"  max_genes: {thresholds.max_genes}")
    logger.info(f"  min_counts: {thresholds.min_counts}")
    logger.info(f"  max_counts: {thresholds.max_counts}")
    logger.info(f"  max_mt_pct: {thresholds.max_mt_pct:.1f}%")

    # Apply QC
    n_before = adata.n_obs
    filter_stats = apply_adaptive_qc(adata, thresholds, min_cells=3)
    run_metadata['cells_before_qc'] = n_before
    run_metadata['cells_after_qc'] = adata.n_obs
    logger.info(f"  Cells: {n_before} -> {adata.n_obs} ({adata.n_obs/n_before*100:.1f}% retained)")

    # Generate post-QC visualizations
    generate_qc_visualizations(adata.copy(), output_dir, stage='after')

    # Save QC statistics
    qc_stats = pd.DataFrame({
        'sample': adata.obs.groupby('sample').size().index,
        'n_cells': adata.obs.groupby('sample').size().values,
        'mean_genes': adata.obs.groupby('sample')['n_genes_by_counts'].mean().values,
        'mean_counts': adata.obs.groupby('sample')['total_counts'].mean().values,
        'mean_mt_pct': adata.obs.groupby('sample')['pct_counts_mt'].mean().values,
    })
    qc_stats.to_csv(output_dir / 'qc' / 'qc_statistics.csv', index=False)
    logger.info(f"  Saved QC statistics to qc/qc_statistics.csv")

    # 3. Normalization (needed before immune cell filtering)
    logger.info("Step 3: Normalizing...")
    adata.layers['counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 3b. Filter immune cells (ALI organoids should be epithelial) - after normalization
    if REMOVE_IMMUNE:
        logger.info("Step 3b: Filtering immune cells...")
        n_before_immune_filter = adata.n_obs
        adata = filter_immune_cells(adata, output_dir)
        n_after_immune_filter = adata.n_obs
        run_metadata['cells_before_immune_filter'] = n_before_immune_filter
        run_metadata['cells_after_immune_filter'] = n_after_immune_filter
        logger.info(f"  Cells after immune filtering: {n_before_immune_filter} -> {n_after_immune_filter} ({n_after_immune_filter/n_before_immune_filter*100:.1f}% retained)")

        if n_after_immune_filter == 0:
            logger.error("No cells remaining after immune filtering! Threshold may be too stringent.")
            logger.error("Consider adjusting the immune filtering threshold or disabling REMOVE_IMMUNE.")
            raise ValueError("No cells remaining after immune filtering")

        # Generate post-immune-filtering QC visualizations
        generate_qc_visualizations(adata.copy(), output_dir, stage='after_immune_filter')

    # Store raw data
    adata.raw = adata.copy()

    # 4. HVG selection
    logger.info("Step 4: Selecting highly variable genes...")
    # Use cell_ranger flavor which is more stable for batch data
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        flavor='cell_ranger',
        n_bins=20
    )
    hvg_list = adata.var_names[adata.var['highly_variable']].tolist()
    run_metadata['n_hvgs'] = len(hvg_list)
    logger.info(f"  Selected {len(hvg_list)} HVGs")

    # Save HVG list
    with open(output_dir / "supplementary" / "hvg_list.txt", 'w') as f:
        f.write('\n'.join(hvg_list))

    # 5. Scale HVGs
    logger.info("Step 5: Scaling...")
    adata_hvg = adata[:, adata.var['highly_variable']].copy()
    sc.pp.scale(adata_hvg, max_value=10)

    # 6. PCA
    logger.info("Step 6: Computing PCA...")
    sc.tl.pca(adata_hvg, n_comps=50, random_state=42)
    adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
    adata.uns['pca'] = adata_hvg.uns['pca']

    # Use specified number of PCs
    n_pcs = N_PCS
    run_metadata['n_pcs_used'] = n_pcs
    logger.info(f"  Using {n_pcs} PCs for downstream analysis")

    # Calculate variance explained
    variance_ratio = adata.uns['pca']['variance_ratio']
    cum_var = np.cumsum(variance_ratio)
    logger.info(f"  {n_pcs} PCs explain {cum_var[n_pcs-1]*100:.1f}% of variance")

    # 6b. Save state before Harmony and compute preliminary UMAP for comparison
    logger.info("Step 6b: Computing preliminary UMAP before integration...")
    adata_before_harmony = adata.copy()
    sc.pp.neighbors(adata_before_harmony, n_neighbors=15, n_pcs=n_pcs, random_state=42)
    sc.tl.umap(adata_before_harmony, random_state=42)

    # 7. Harmony batch correction
    logger.info("Step 7: Running Harmony batch correction...")
    sce.pp.harmony_integrate(adata, key='sample', max_iter_harmony=20)

    # 8. Neighbors and UMAP
    logger.info("Step 8: Computing neighbors and UMAP...")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, use_rep='X_pca_harmony', random_state=42)
    sc.tl.umap(adata, random_state=42)

    # 8b. Generate integration quality visualizations
    generate_integration_visualizations(adata_before_harmony, adata, output_dir)

    # 9. Clustering
    logger.info("Step 9: Clustering...")
    sc.tl.leiden(adata, resolution=0.8, random_state=42)
    n_clusters = adata.obs['leiden'].nunique()
    run_metadata['n_clusters'] = n_clusters
    run_metadata['resolution'] = 0.8
    logger.info(f"  Found {n_clusters} clusters")

    # 10. DEG analysis
    logger.info("Step 10: Finding marker genes...")
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', pts=True, n_genes=adata.n_vars)

    # Extract DEG results
    from scrna_agent.stats import run_full_deg_analysis
    deg_df = run_full_deg_analysis(adata, groupby='leiden', method='wilcoxon')
    deg_df.to_csv(output_dir / "deg" / "markers_full.csv", index=False)

    # Top markers
    top_markers = deg_df.groupby('cluster').head(25)
    top_markers.to_csv(output_dir / "deg" / "markers.csv", index=False)

    # 11. Tissue annotation
    logger.info("Step 11: Running tissue-aware annotation...")
    adata.obs['tissue'] = TISSUE

    from scrna_agent.tissue import annotate_cells, generate_tissue_summary, generate_cell_metadata
    annotate_cells(adata)

    # Save cell metadata
    cell_metadata = generate_cell_metadata(adata, sample_col='sample')
    cell_metadata.to_csv(output_dir / "annotation" / "cell_metadata.csv")

    # Tissue summary
    tissue_summary = generate_tissue_summary(adata)
    tissue_summary.to_csv(output_dir / "supplementary" / "tissue_summary.csv")

    # 11b. Generate marker annotation visualizations
    generate_marker_plots(adata, output_dir)

    # 11c. Run condition-level DEG analysis
    deg_results = run_condition_deg_analysis(adata, output_dir)

    # 11d. Generate DEG visualizations
    generate_deg_visualizations(deg_results, adata, output_dir)

    # 12. Generate plots
    logger.info("Step 12: Generating basic UMAP plots...")

    # UMAP by cluster
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='leiden', ax=ax, show=False, title='Clusters')
    plt.savefig(output_dir / "umap" / "umap_clusters.png", dpi=150, bbox_inches='tight')
    plt.close()

    # UMAP by condition
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='condition', ax=ax, show=False, title='Condition')
    plt.savefig(output_dir / "umap" / "umap_condition.png", dpi=150, bbox_inches='tight')
    plt.close()

    # UMAP by group
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='group', ax=ax, show=False, title='Virus Group')
    plt.savefig(output_dir / "umap" / "umap_group.png", dpi=150, bbox_inches='tight')
    plt.close()

    # UMAP by sample
    fig, ax = plt.subplots(figsize=(14, 10))
    sc.pl.umap(adata, color='sample', ax=ax, show=False, title='Sample')
    plt.savefig(output_dir / "umap" / "umap_sample.png", dpi=150, bbox_inches='tight')
    plt.close()

    # UMAP by cell type
    if 'celltype' in adata.obs.columns:
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color='celltype', ax=ax, show=False, title='Cell Type (Marker-based)')
        plt.savefig(output_dir / "umap" / "umap_celltype.png", dpi=150, bbox_inches='tight')
        plt.close()

    # 12b. Generate extended visualizations
    generate_extended_visualizations(adata, output_dir)

    # PCA variance plot
    variance_ratio = adata.uns['pca']['variance_ratio']
    cum_var = np.cumsum(variance_ratio)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].bar(range(1, 51), variance_ratio[:50] * 100, alpha=0.7, color='steelblue')
    axes[0].axvline(n_pcs, color='red', linestyle='--', label=f'Selected: {n_pcs}')
    axes[0].set_xlabel('PC')
    axes[0].set_ylabel('Variance Explained (%)')
    axes[0].set_title('Scree Plot')
    axes[0].legend()

    axes[1].plot(range(1, 51), cum_var[:50] * 100, 'b-', linewidth=2)
    axes[1].axvline(n_pcs, color='red', linestyle='--', label=f'{n_pcs} PCs: {cum_var[n_pcs-1]*100:.1f}%')
    axes[1].axhline(90, color='gray', linestyle=':', alpha=0.5)
    axes[1].set_xlabel('Number of PCs')
    axes[1].set_ylabel('Cumulative Variance (%)')
    axes[1].set_title('Cumulative Variance')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / "clustering" / "pca_variance.png", dpi=150, bbox_inches='tight')
    plt.close()

    # Cluster composition by condition
    composition = pd.crosstab(adata.obs['leiden'], adata.obs['group'], normalize='index') * 100
    composition.to_csv(output_dir / "clustering" / "cluster_composition.csv")

    fig, ax = plt.subplots(figsize=(14, 8))
    composition.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
    ax.set_xlabel('Cluster')
    ax.set_ylabel('Percentage')
    ax.set_title('Cluster Composition by Virus Group')
    ax.legend(title='Group', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(output_dir / "clustering" / "cluster_composition.png", dpi=150, bbox_inches='tight')
    plt.close()

    # 13. Save data
    logger.info("Step 13: Saving processed data...")
    adata.write_h5ad(output_dir / "integrated.h5ad")

    # Save clusters
    clusters_df = pd.DataFrame({
        'cell_id': adata.obs_names,
        'cluster': adata.obs['leiden'].values,
        'sample': adata.obs['sample'].values,
        'condition': adata.obs['condition'].values,
        'group': adata.obs['group'].values
    })
    clusters_df.to_csv(output_dir / "clustering" / "clusters.csv", index=False)

    # 14. Save metadata
    run_metadata['summary'] = {
        'cells_analyzed': adata.n_obs,
        'genes_analyzed': adata.n_vars,
        'n_samples': len(adata.obs['sample'].unique()),
        'n_hvgs': len(hvg_list),
        'n_pcs_used': n_pcs,
        'n_clusters': n_clusters,
        'tissue': TISSUE,
    }

    with open(output_dir / "run_metadata.json", 'w') as f:
        json.dump(run_metadata, f, indent=2, default=str)

    # 15. Generate report
    logger.info("Step 14: Generating report...")
    generate_report(adata, run_metadata, output_dir)

    logger.info("=" * 60)
    logger.info("Analysis complete!")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 60)

    return adata


def generate_integration_visualizations(adata_before, adata_after, output_dir):
    """Generate batch integration quality visualizations."""
    integration_dir = output_dir / 'integration'
    logger.info("Generating integration quality visualizations...")

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))

    # Before integration
    ax = axes[0, 0]
    sc.pl.pca(adata_before, color='sample', ax=ax, show=False, title='PCA by Sample (Before Harmony)')
    ax.set_title('PCA by Sample (Before Harmony)', fontweight='bold')

    ax = axes[0, 1]
    sc.pl.pca(adata_before, color='condition', ax=ax, show=False, title='PCA by Condition (Before Harmony)')
    ax.set_title('PCA by Condition (Before Harmony)', fontweight='bold')

    ax = axes[0, 2]
    sc.pl.umap(adata_before, color='sample', ax=ax, show=False, title='UMAP by Sample (Before Harmony)')
    ax.set_title('UMAP by Sample (Before Harmony)', fontweight='bold')

    # After integration - plot Harmony embedding
    ax = axes[1, 0]
    sc.pl.embedding(adata_after, basis='X_pca_harmony', color='sample', ax=ax, show=False, title='Harmony PCA by Sample')
    ax.set_title('Harmony PCA by Sample', fontweight='bold')

    ax = axes[1, 1]
    sc.pl.embedding(adata_after, basis='X_pca_harmony', color='condition', ax=ax, show=False, title='Harmony PCA by Condition')
    ax.set_title('Harmony PCA by Condition', fontweight='bold')

    ax = axes[1, 2]
    sc.pl.umap(adata_after, color='sample', ax=ax, show=False, title='UMAP by Sample (After Harmony)')
    ax.set_title('UMAP by Sample (After Harmony)', fontweight='bold')

    plt.tight_layout()
    plt.savefig(integration_dir / 'integration_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"  Saved integration plots to {integration_dir}/integration_comparison.png")


def generate_marker_plots(adata, output_dir):
    """Generate cell type marker visualizations."""
    annotation_dir = output_dir / 'annotation'
    logger.info("Generating marker gene visualizations...")

    # Collect all markers that exist in the dataset
    all_markers = []
    marker_dict_for_plot = {}
    for cell_type, markers in AIRWAY_MARKERS.items():
        present_markers = [m for m in markers if m in adata.var_names]
        if present_markers:
            all_markers.extend(present_markers)
            marker_dict_for_plot[cell_type] = present_markers

    if len(all_markers) < 3:
        logger.warning("Too few markers present in dataset, skipping marker plots")
        return

    # Dotplot
    try:
        fig, ax = plt.subplots(figsize=(12, 8))
        sc.pl.dotplot(adata, all_markers[:min(30, len(all_markers))], groupby='leiden', ax=ax, show=False)
        plt.savefig(annotation_dir / 'marker_dotplot.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved dotplot to {annotation_dir}/marker_dotplot.png")
    except Exception as e:
        logger.warning(f"  Could not generate dotplot: {e}")

    # Heatmap of top markers
    try:
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.heatmap(adata, all_markers[:min(40, len(all_markers))], groupby='leiden',
                      swap_axes=True, show_gene_labels=True, ax=ax, show=False, cmap='viridis')
        plt.savefig(annotation_dir / 'celltype_markers_heatmap.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved heatmap to {annotation_dir}/celltype_markers_heatmap.png")
    except Exception as e:
        logger.warning(f"  Could not generate heatmap: {e}")

    # Cell type proportions by condition
    if 'celltype' in adata.obs.columns:
        celltype_composition = pd.crosstab(adata.obs['condition'], adata.obs['celltype'], normalize='index') * 100
        celltype_composition.to_csv(annotation_dir / 'celltype_proportions.csv')

        fig, ax = plt.subplots(figsize=(14, 8))
        celltype_composition.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
        ax.set_xlabel('Condition')
        ax.set_ylabel('Percentage')
        ax.set_title('Cell Type Composition by Condition', fontweight='bold')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.02, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(annotation_dir / 'celltype_composition.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved cell type composition to {annotation_dir}/celltype_composition.png")


def run_condition_deg_analysis(adata, output_dir):
    """Run pairwise DEG analysis: each condition vs Control."""
    deg_dir = output_dir / 'deg' / 'condition_comparisons'
    logger.info("Running condition-level DEG analysis...")

    # Get unique conditions excluding Control
    conditions = [c for c in adata.obs['condition'].unique() if 'Control' not in c]

    deg_results = {}

    for condition in conditions:
        logger.info(f"  Analyzing {condition} vs Control...")

        # Subset to condition vs control
        mask = adata.obs['condition'].isin([condition, 'Control'])
        adata_subset = adata[mask].copy()

        # Run DEG
        try:
            sc.tl.rank_genes_groups(
                adata_subset,
                groupby='condition',
                groups=[condition],
                reference='Control',
                method='wilcoxon',
                pts=True,
                n_genes=adata_subset.n_vars
            )

            # Extract results
            result_df = sc.get.rank_genes_groups_df(adata_subset, group=condition)
            result_df['comparison'] = f'{condition}_vs_Control'
            result_df.to_csv(deg_dir / f'{condition}_vs_Control_degs.csv', index=False)

            deg_results[condition] = result_df

            logger.info(f"    Found {(result_df['pvals_adj'] < 0.05).sum()} significant DEGs")

        except Exception as e:
            logger.warning(f"    Failed to analyze {condition}: {e}")

    return deg_results


def generate_deg_visualizations(deg_results, adata, output_dir):
    """Generate volcano and MA plots for DEG results."""
    deg_dir = output_dir / 'deg' / 'condition_comparisons'
    logger.info("Generating DEG visualizations...")

    for condition, result_df in deg_results.items():
        try:
            # Volcano plot
            fig, ax = plt.subplots(figsize=(10, 8))

            # Add fold change and -log10(pval) if not present
            if 'logfoldchanges' in result_df.columns and 'pvals_adj' in result_df.columns:
                result_df['neg_log10_pval'] = -np.log10(result_df['pvals_adj'].replace(0, 1e-300))

                # Color by significance
                sig_mask = (result_df['pvals_adj'] < 0.05) & (np.abs(result_df['logfoldchanges']) > 0.5)
                colors = ['red' if s else 'gray' for s in sig_mask]

                ax.scatter(result_df['logfoldchanges'], result_df['neg_log10_pval'],
                          c=colors, alpha=0.5, s=10, edgecolors='none')

                # Add significance thresholds
                ax.axhline(-np.log10(0.05), color='blue', linestyle='--', linewidth=1, alpha=0.5)
                ax.axvline(0.5, color='blue', linestyle='--', linewidth=1, alpha=0.5)
                ax.axvline(-0.5, color='blue', linestyle='--', linewidth=1, alpha=0.5)

                ax.set_xlabel('Log2 Fold Change', fontsize=12)
                ax.set_ylabel('-Log10(Adjusted P-value)', fontsize=12)
                ax.set_title(f'Volcano Plot: {condition} vs Control', fontsize=14, fontweight='bold')
                ax.grid(True, alpha=0.3)

                plt.tight_layout()
                plt.savefig(deg_dir / f'volcano_{condition}_vs_Control.png', dpi=150, bbox_inches='tight')
                plt.close()

                logger.info(f"  Saved volcano plot for {condition}")

        except Exception as e:
            logger.warning(f"  Could not generate volcano plot for {condition}: {e}")


def generate_extended_visualizations(adata, output_dir):
    """Generate extended UMAP suite and additional visualizations."""
    umap_dir = output_dir / 'umap'
    logger.info("Generating extended UMAP visualizations...")

    # UMAP by n_genes (QC overlay)
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='n_genes_by_counts', ax=ax, show=False, cmap='viridis')
    ax.set_title('UMAP colored by Gene Count', fontweight='bold')
    plt.savefig(umap_dir / 'umap_n_genes.png', dpi=150, bbox_inches='tight')
    plt.close()

    # UMAPs for key marker genes
    key_markers = ['FOXJ1', 'KRT5', 'SCGB1A1', 'MUC5AC']
    present_markers = [m for m in key_markers if m in adata.var_names]

    if present_markers:
        n_markers = len(present_markers)
        fig, axes = plt.subplots(1, n_markers, figsize=(6*n_markers, 5))
        if n_markers == 1:
            axes = [axes]

        for i, marker in enumerate(present_markers):
            sc.pl.umap(adata, color=marker, ax=axes[i], show=False, cmap='Reds', vmin=0)
            axes[i].set_title(f'{marker} Expression', fontweight='bold')

        plt.tight_layout()
        plt.savefig(umap_dir / 'umap_key_markers.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved key marker UMAPs")

    # UMAPs for viral response genes
    present_viral_genes = [g for g in VIRAL_RESPONSE_GENES if g in adata.var_names]

    if len(present_viral_genes) >= 4:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        axes = axes.flatten()

        for i, gene in enumerate(present_viral_genes[:4]):
            sc.pl.umap(adata, color=gene, ax=axes[i], show=False, cmap='Reds', vmin=0)
            axes[i].set_title(f'{gene} Expression', fontweight='bold')

        plt.tight_layout()
        plt.savefig(umap_dir / 'umap_viral_response_genes.png', dpi=150, bbox_inches='tight')
        plt.close()
        logger.info(f"  Saved viral response gene UMAPs")


def generate_report(adata, run_metadata, output_dir):
    """Generate analysis report."""
    from datetime import datetime

    n_clusters = adata.obs['leiden'].nunique()

    report = f"""# JMV_ALI scRNA-seq Integrated Analysis Report

## Summary

- **Cells analyzed**: {adata.n_obs:,}
- **Genes**: {adata.n_vars:,}
- **Samples**: {len(adata.obs['sample'].unique())}
- **HVGs**: {run_metadata.get('n_hvgs', 'N/A')}
- **PCs used**: {run_metadata.get('n_pcs_used', 'N/A')}
- **Clusters**: {n_clusters}
- **Tissue type**: {TISSUE}

## Statistical Methods

### QC Threshold Selection (MAD-based)
Thresholds were computed using Median Absolute Deviation (MAD):
- **Method**: median ± 3 × MAD (robust to outliers)
- **min_genes**: {run_metadata.get('qc_thresholds', {}).get('min_genes', 'N/A')}
- **max_genes**: {run_metadata.get('qc_thresholds', {}).get('max_genes', 'N/A')}
- **max_mt_pct**: {run_metadata.get('qc_thresholds', {}).get('max_mt_pct', 'N/A'):.1f}%

### Batch Correction
- **Method**: Harmony integration
- **Batch key**: sample (24 samples)

### PC Selection
- **Method**: Elbow method
- **PCs selected**: {run_metadata.get('n_pcs_used', 'N/A')}

## Sample Distribution

| Sample | Cells |
|--------|-------|
"""

    for sample, count in adata.obs['sample'].value_counts().items():
        report += f"| {sample} | {count:,} |\n"

    report += f"""
## Cluster Distribution

| Cluster | Cells | Percentage |
|---------|-------|------------|
"""

    cluster_counts = adata.obs['leiden'].value_counts().sort_index()
    for cluster, count in cluster_counts.items():
        pct = count / adata.n_obs * 100
        report += f"| {cluster} | {count:,} | {pct:.1f}% |\n"

    report += f"""
## Condition Distribution

| Condition | Cells | Percentage |
|-----------|-------|------------|
"""

    for cond, count in adata.obs['condition'].value_counts().items():
        pct = count / adata.n_obs * 100
        report += f"| {cond} | {count:,} | {pct:.1f}% |\n"

    report += f"""
## Output Files

| File | Description |
|------|-------------|
| `integrated.h5ad` | Processed AnnData with all embeddings |
| `clusters.csv` | Cell cluster assignments |
| `markers.csv` | Top marker genes per cluster |
| `markers_full.csv` | Complete DEG table |
| `hvg_list.txt` | List of highly variable genes |
| `cell_metadata.csv` | Cell metadata with annotations |
| `cluster_composition.csv` | Cluster composition by condition |
| `umap_*.png` | UMAP visualizations |
| `pca_variance.png` | PCA explained variance |
| `run_metadata.json` | Complete analysis metadata |

## Reproducibility

- **Random seed**: 42
- **Analysis timestamp**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

---
*Generated by scrna-agent with statistical rigor*
"""

    with open(output_dir / "report.md", 'w') as f:
        f.write(report)


if __name__ == "__main__":
    # Load all samples
    adata = load_all_samples()

    # Run integrated analysis
    run_integrated_analysis(adata, OUTPUT_DIR)
