"""
MVP scRNA-seq Analysis Pipeline
Standalone pipeline that works without LangChain/LLM dependencies

This pipeline implements statistically rigorous methods:
- Adaptive QC thresholds using robust statistics (MAD)
- Variance-aware normalization with size-factor correction
- Explained variance-based PC selection
- Full differential expression with multiple testing correction
"""

import os
import json
import logging
from pathlib import Path
from typing import Optional, Dict, Any, Union, List, Tuple
from dataclasses import dataclass, field
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Configuration for the scRNA-seq analysis pipeline."""
    # QC parameters - adaptive by default
    adaptive_qc: bool = True  # Use adaptive MAD-based thresholds
    qc_method: str = "mad"  # "mad" or "quantile"
    qc_nmads: float = 3.0  # Number of MADs for outlier detection
    min_genes: int = 200  # Fallback/floor for min_genes
    min_cells: int = 3  # Minimum cells per gene
    max_mt_pct: float = 20.0  # Ceiling for max MT%

    # Normalization
    target_sum: float = 1e4  # Size factor target
    n_top_genes: int = 2000  # Number of HVGs
    hvg_flavor: str = "seurat_v3"  # HVG selection method

    # Dimensionality reduction
    n_pcs: int = 50  # Maximum PCs to compute
    n_neighbors: int = 15  # kNN neighbors
    pc_selection_method: str = "elbow"  # "elbow" or "variance_ratio"
    variance_threshold: float = 0.90  # For variance_ratio method

    # Clustering
    resolution: float = 0.5
    clustering_method: str = "leiden"  # "leiden" or "louvain"

    # DEG detection
    n_genes_per_cluster: int = 25  # Top genes to report per cluster
    deg_method: str = "wilcoxon"  # "wilcoxon", "t-test", "logreg"
    deg_correction: str = "benjamini-hochberg"  # Multiple testing correction
    min_fold_change: float = 0.25  # Minimum log2FC for significance

    # Tissue and annotation
    tissue: Optional[str] = None  # Single tissue for all cells
    tissue_map_path: Optional[str] = None  # Path to tissue-map CSV
    sample_col: str = "sample"  # Column for sample IDs
    allow_missing_tissue: bool = False  # Allow samples without tissue
    markers_path: Optional[str] = None  # Custom markers YAML
    run_annotation: bool = True  # Run marker-based annotation

    # Biological validation (NEW - from JMV_ALI session learnings)
    sample_type: str = "primary_tissue"  # "organoid", "primary_tissue", "cell_line", "pbmc"
    use_multi_method_annotation: bool = False  # Use CellTypist + marker-based consensus
    celltypist_model: Optional[str] = None  # CellTypist model (auto-selected if None)

    # Reproducibility
    seed: int = 42

    # Output format
    report_format: str = "md"  # "md" or "html"

    @classmethod
    def from_yaml(cls, path: str) -> "PipelineConfig":
        """Load configuration from YAML file."""
        import yaml
        with open(path, 'r') as f:
            config_dict = yaml.safe_load(f)
        # Filter to only valid fields
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in (config_dict or {}).items() if k in valid_fields}
        return cls(**filtered) if filtered else cls()

    def to_yaml(self, path: str):
        """Save configuration to YAML file."""
        import yaml
        with open(path, 'w') as f:
            yaml.dump(self.__dict__, f, default_flow_style=False)


@dataclass
class PipelineResults:
    """Results from the pipeline run."""
    adata: Any = None
    qc_summary: pd.DataFrame = None
    clusters: pd.DataFrame = None
    markers: pd.DataFrame = None
    markers_full: pd.DataFrame = None  # Full DEG table
    cell_metadata: pd.DataFrame = None
    tissue_summary: pd.DataFrame = None
    hvg_list: List[str] = field(default_factory=list)
    run_metadata: Dict[str, Any] = field(default_factory=dict)
    output_dir: str = None
    files_created: list = field(default_factory=list)


def load_data(
    input_path: str,
    input_format: str = "auto"
) -> sc.AnnData:
    """
    Load single-cell data from various formats.

    Args:
        input_path: Path to input file or directory
        input_format: "h5ad", "10x", "10x_h5", or "auto" (auto-detect)

    Returns:
        AnnData object
    """
    input_path = Path(input_path)

    if input_format == "auto":
        if input_path.suffix == ".h5ad":
            input_format = "h5ad"
        elif input_path.suffix == ".h5":
            input_format = "10x_h5"
        elif input_path.is_dir():
            # Check for 10x files
            mtx_files = list(input_path.glob("*matrix*.mtx*"))
            if mtx_files:
                input_format = "10x"
            else:
                raise ValueError(f"Cannot auto-detect format for directory: {input_path}")
        else:
            raise ValueError(f"Cannot auto-detect format for: {input_path}")

    logger.info(f"Loading data from {input_path} (format: {input_format})")

    if input_format == "h5ad":
        adata = sc.read_h5ad(input_path)
    elif input_format == "10x_h5":
        # 10x CellRanger h5 format
        adata = sc.read_10x_h5(str(input_path))
    elif input_format == "10x":
        # Handle 10x CellRanger output (mtx format)
        if input_path.is_dir():
            # Check for filtered_feature_bc_matrix subdirectory
            filtered_dir = input_path / "filtered_feature_bc_matrix"
            if filtered_dir.exists():
                input_path = filtered_dir
            # Also check outs/filtered_feature_bc_matrix
            outs_filtered = input_path / "outs" / "filtered_feature_bc_matrix"
            if outs_filtered.exists():
                input_path = outs_filtered
        adata = sc.read_10x_mtx(str(input_path), var_names='gene_symbols', cache=True)
    else:
        raise ValueError(f"Unsupported format: {input_format}")

    # Ensure var_names are unique
    adata.var_names_make_unique()

    logger.info(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def run_qc(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
    """
    Run quality control filtering with statistically justified thresholds.

    Uses adaptive MAD-based thresholds when config.adaptive_qc=True:
    - Thresholds are computed as median ± nmads * MAD
    - MAD (Median Absolute Deviation) is robust to outliers
    - This approach adapts to dataset-specific distributions

    Args:
        adata: AnnData object (modified in place)
        config: Pipeline configuration
        run_metadata: Optional dict to store QC metadata

    Returns:
        DataFrame with QC summary including statistical rationale
    """
    from scrna_agent.stats import (
        compute_adaptive_qc_thresholds,
        apply_adaptive_qc,
        compute_qc_distributions,
        generate_qc_rationale
    )

    logger.info("Running QC with statistically justified thresholds...")

    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )

    # Compute distribution statistics for reporting
    distributions = compute_qc_distributions(adata)

    # Store pre-filter stats
    qc_data = {
        'metric': [],
        'value': [],
        'description': []
    }

    qc_data['metric'].append('n_cells_raw')
    qc_data['value'].append(n_cells_before)
    qc_data['description'].append('Total cells before QC')

    qc_data['metric'].append('n_genes_raw')
    qc_data['value'].append(n_genes_before)
    qc_data['description'].append('Total genes before QC')

    # Add distribution statistics
    for metric_name, dist in distributions.items():
        short_name = metric_name.replace('_by_counts', '').replace('pct_counts_', '')
        qc_data['metric'].append(f'{short_name}_median')
        qc_data['value'].append(dist['median'])
        qc_data['description'].append(f'Median {short_name} (before QC)')

        qc_data['metric'].append(f'{short_name}_mad')
        qc_data['value'].append(dist['mad'])
        qc_data['description'].append(f'MAD of {short_name} (robust spread)')

        qc_data['metric'].append(f'{short_name}_q05')
        qc_data['value'].append(dist['q05'])
        qc_data['description'].append(f'5th percentile {short_name}')

        qc_data['metric'].append(f'{short_name}_q95')
        qc_data['value'].append(dist['q95'])
        qc_data['description'].append(f'95th percentile {short_name}')

    if config.adaptive_qc:
        # Compute adaptive thresholds using robust statistics
        thresholds = compute_adaptive_qc_thresholds(
            adata,
            nmads=config.qc_nmads,
            min_genes_floor=config.min_genes,
            max_mt_pct_ceiling=config.max_mt_pct,
            method=config.qc_method
        )

        # Apply adaptive filtering
        filter_stats = apply_adaptive_qc(adata, thresholds, config.min_cells)

        # Record thresholds
        qc_data['metric'].append('qc_method')
        qc_data['value'].append(config.qc_method)
        qc_data['description'].append('QC threshold selection method')

        qc_data['metric'].append('qc_nmads')
        qc_data['value'].append(config.qc_nmads)
        qc_data['description'].append('Number of MADs for outlier detection')

        qc_data['metric'].append('min_genes_threshold')
        qc_data['value'].append(thresholds.min_genes)
        qc_data['description'].append('Adaptive min genes threshold')

        qc_data['metric'].append('max_genes_threshold')
        qc_data['value'].append(thresholds.max_genes)
        qc_data['description'].append('Adaptive max genes threshold (doublet filter)')

        qc_data['metric'].append('min_counts_threshold')
        qc_data['value'].append(thresholds.min_counts)
        qc_data['description'].append('Adaptive min UMI counts threshold')

        qc_data['metric'].append('max_counts_threshold')
        qc_data['value'].append(thresholds.max_counts)
        qc_data['description'].append('Adaptive max UMI counts threshold')

        qc_data['metric'].append('max_mt_pct_threshold')
        qc_data['value'].append(thresholds.max_mt_pct)
        qc_data['description'].append('Adaptive max MT% threshold')

        # Store in run_metadata
        if run_metadata is not None:
            run_metadata['qc_thresholds'] = thresholds.to_dict()
            run_metadata['qc_distributions'] = distributions
            run_metadata['qc_rationale'] = generate_qc_rationale(thresholds, distributions)

    else:
        # Use fixed thresholds (legacy behavior)
        sc.pp.filter_cells(adata, min_genes=config.min_genes)
        sc.pp.filter_genes(adata, min_cells=config.min_cells)

        mt_mask = adata.obs['pct_counts_mt'] < config.max_mt_pct
        adata._inplace_subset_obs(mt_mask)

        qc_data['metric'].append('min_genes_threshold')
        qc_data['value'].append(config.min_genes)
        qc_data['description'].append('Fixed min genes threshold')

        qc_data['metric'].append('max_mt_pct_threshold')
        qc_data['value'].append(config.max_mt_pct)
        qc_data['description'].append('Fixed max MT% threshold')

    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars

    qc_data['metric'].append('n_cells_filtered')
    qc_data['value'].append(n_cells_after)
    qc_data['description'].append('Cells after QC')

    qc_data['metric'].append('n_genes_filtered')
    qc_data['value'].append(n_genes_after)
    qc_data['description'].append('Genes after QC')

    qc_data['metric'].append('cells_removed_total')
    qc_data['value'].append(n_cells_before - n_cells_after)
    qc_data['description'].append('Total cells removed by QC')

    qc_data['metric'].append('retention_rate')
    qc_data['value'].append(n_cells_after / n_cells_before * 100)
    qc_data['description'].append('Cell retention rate (%)')

    qc_summary = pd.DataFrame(qc_data)

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['cells_before_qc'] = n_cells_before
        run_metadata['cells_after_qc'] = n_cells_after
        run_metadata['genes_before_qc'] = n_genes_before
        run_metadata['genes_after_qc'] = n_genes_after

    logger.info(f"QC complete: {n_cells_before} -> {n_cells_after} cells "
                f"({n_cells_after/n_cells_before*100:.1f}% retained), "
                f"{n_genes_before} -> {n_genes_after} genes")

    return qc_summary


def normalize_and_scale(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> List[str]:
    """
    Normalize, log-transform, and identify highly variable genes.

    Normalization approach:
    1. Size-factor normalization: counts / total_counts * target_sum
       - This corrects for sequencing depth differences between cells
       - target_sum is typically 10,000 (counts per 10k)

    2. Log1p transformation: log(x + 1)
       - Stabilizes variance (variance no longer proportional to mean)
       - Makes data more normally distributed

    3. HVG selection using Seurat v3 method:
       - Fits mean-variance relationship
       - Identifies genes with high variance relative to their mean
       - These genes capture biological variation

    Args:
        adata: AnnData object (modified in place)
        config: Pipeline configuration
        run_metadata: Optional dict to store normalization metadata

    Returns:
        List of highly variable gene names
    """
    logger.info("Normalizing with size-factor correction and log1p transformation...")

    # Store raw counts
    adata.layers['counts'] = adata.X.copy()

    # Size-factor normalization
    # Mathematical basis: X_norm = X / sum(X) * target_sum
    # This normalizes for sequencing depth
    sc.pp.normalize_total(adata, target_sum=config.target_sum)

    # Log1p transformation
    # Mathematical basis: X_log = log(X + 1)
    # Variance stabilization: Var(log(X)) ≈ constant for Poisson data
    sc.pp.log1p(adata)

    # Store normalized data for later use (e.g., DEG analysis)
    adata.raw = adata.copy()

    # Find highly variable genes
    # Seurat v3 method: variance stabilizing transformation approach
    # Fits a LOESS regression to mean-variance relationship
    # Selects genes with high standardized variance
    logger.info(f"Selecting HVGs using {config.hvg_flavor} method...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config.n_top_genes,
        flavor=config.hvg_flavor,
        layer='counts'
    )

    # Get HVG list
    hvg_list = adata.var_names[adata.var['highly_variable']].tolist()
    n_hvg = len(hvg_list)

    # Subset to HVGs for downstream analysis
    adata_hvg = adata[:, adata.var['highly_variable']].copy()

    # Scale (z-score normalization per gene)
    # Mathematical basis: X_scaled = (X - mean(X)) / std(X)
    # max_value clips extreme values to prevent outlier influence
    sc.pp.scale(adata_hvg, max_value=10)

    # Copy scaled data back
    adata.X = adata.raw.X.copy()
    adata.obsm['X_hvg_scaled'] = adata_hvg.X

    # Store HVG statistics in var
    if 'highly_variable_rank' not in adata.var.columns:
        # Create ranking based on variance
        hvg_var = adata.var.loc[adata.var['highly_variable'], 'highly_variable_nbatches']
        if 'highly_variable_nbatches' in adata.var.columns:
            adata.var.loc[adata.var['highly_variable'], 'hvg_rank'] = range(1, n_hvg + 1)

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['normalization_method'] = 'size_factor_log1p'
        run_metadata['target_sum'] = config.target_sum
        run_metadata['hvg_method'] = config.hvg_flavor
        run_metadata['n_hvgs'] = n_hvg
        run_metadata['hvg_list'] = hvg_list

    logger.info(f"Identified {n_hvg} highly variable genes using {config.hvg_flavor} method")

    return hvg_list


def run_dimensionality_reduction(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> int:
    """
    Run PCA and UMAP with model-driven PC selection.

    PCA approach:
    1. Fit PCA on HVGs only (these capture biological variation)
    2. Select optimal number of PCs using elbow or variance threshold
    3. Build kNN graph in PC space

    UMAP approach:
    1. Non-linear dimensionality reduction preserving local structure
    2. Uses kNN graph from PCA space
    3. Produces 2D visualization

    Mathematical basis for PC selection:
    - Elbow method: Find point of maximum curvature in variance explained curve
    - Variance threshold: Select PCs explaining X% of total variance

    Args:
        adata: AnnData object (modified in place)
        config: Pipeline configuration
        run_metadata: Optional dict to store dimensionality reduction metadata

    Returns:
        Number of PCs selected for downstream analysis
    """
    from scrna_agent.stats import compute_pc_selection, generate_pc_selection_rationale

    logger.info("Running dimensionality reduction with model-driven PC selection...")

    # Set random seed for reproducibility
    np.random.seed(config.seed)

    # PCA on HVGs
    # Mathematical basis: Find orthogonal directions of maximum variance
    # Using HVGs ensures we capture biological variation, not technical noise
    logger.info(f"Computing PCA on highly variable genes (max {config.n_pcs} components)...")
    sc.tl.pca(adata, n_comps=config.n_pcs, use_highly_variable=True, random_state=config.seed)

    # Model-driven PC selection
    pc_result = compute_pc_selection(
        adata,
        method=config.pc_selection_method,
        variance_threshold=config.variance_threshold,
        min_pcs=10,
        max_pcs=config.n_pcs
    )

    n_pcs_use = pc_result.n_pcs_selected

    logger.info(f"PCA complete: {adata.obsm['X_pca'].shape}")
    logger.info(f"Selected {n_pcs_use} PCs using {config.pc_selection_method} method "
                f"(explains {pc_result.cumulative_variance[n_pcs_use-1]*100:.1f}% variance)")

    # Build kNN graph from PCA space
    # Mathematical basis: Approximate high-dimensional distances using k nearest neighbors
    # This captures local structure while being computationally efficient
    logger.info(f"Building kNN graph (k={config.n_neighbors}) from {n_pcs_use} PCs...")
    sc.pp.neighbors(
        adata,
        n_neighbors=config.n_neighbors,
        n_pcs=n_pcs_use,
        random_state=config.seed
    )

    # UMAP
    # Mathematical basis: Preserves local structure by modeling fuzzy simplicial sets
    # Optimizes low-dimensional embedding to match high-dimensional neighborhoods
    logger.info("Computing UMAP embedding...")
    sc.tl.umap(adata, random_state=config.seed)

    # Store PC selection results in adata.uns
    adata.uns['pc_selection'] = {
        'method': config.pc_selection_method,
        'n_pcs_selected': n_pcs_use,
        'variance_explained': pc_result.variance_explained[:n_pcs_use],
        'cumulative_variance': pc_result.cumulative_variance[:n_pcs_use],
        'elbow_point': pc_result.elbow_point
    }

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['pc_selection'] = pc_result.to_dict()
        run_metadata['n_pcs_used'] = n_pcs_use
        run_metadata['n_neighbors'] = config.n_neighbors
        run_metadata['pc_rationale'] = generate_pc_selection_rationale(pc_result)

    logger.info(f"UMAP complete: {adata.obsm['X_umap'].shape}")

    return n_pcs_use


def run_clustering(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> pd.DataFrame:
    """
    Run graph-based clustering on kNN graph.

    Uses community detection algorithms on the kNN graph built from PCA:
    - Leiden: Improved modularity optimization with guaranteed connectivity
    - Louvain: Classic modularity-based community detection

    Mathematical basis:
    - Modularity Q = (1/2m) * sum[(A_ij - k_i*k_j/2m) * delta(c_i, c_j)]
    - Higher resolution -> more clusters (finer granularity)
    - Lower resolution -> fewer clusters (coarser grouping)

    Args:
        adata: AnnData with kNN graph computed
        config: Pipeline configuration
        run_metadata: Optional dict to store clustering metadata

    Returns:
        DataFrame with cluster assignments
    """
    logger.info(f"Running {config.clustering_method} clustering "
                f"(resolution={config.resolution})...")

    np.random.seed(config.seed)

    if config.clustering_method == "leiden":
        sc.tl.leiden(adata, resolution=config.resolution, random_state=config.seed)
        cluster_key = 'leiden'
    else:
        sc.tl.louvain(adata, resolution=config.resolution, random_state=config.seed)
        cluster_key = 'louvain'

    # Create clusters dataframe
    clusters_df = pd.DataFrame({
        'cell_id': adata.obs_names,
        'cluster': adata.obs[cluster_key].values
    })

    n_clusters = adata.obs[cluster_key].nunique()

    # Calculate cluster statistics
    cluster_sizes = adata.obs[cluster_key].value_counts()

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['clustering_method'] = config.clustering_method
        run_metadata['resolution'] = config.resolution
        run_metadata['n_clusters'] = n_clusters
        run_metadata['cluster_sizes'] = cluster_sizes.to_dict()

    logger.info(f"Found {n_clusters} clusters (sizes range: "
                f"{cluster_sizes.min()} - {cluster_sizes.max()})")

    return clusters_df


def find_markers(
    adata: sc.AnnData,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Find marker genes for each cluster using statistical DEG analysis.

    Statistical approach:
    1. Wilcoxon rank-sum test: Non-parametric, robust to non-normality
       - Tests if expression differs between cluster vs rest
       - Appropriate for scRNA-seq's zero-inflated distributions

    2. Multiple testing correction: Benjamini-Hochberg procedure
       - Controls False Discovery Rate (FDR)
       - Adjusts p-values for multiple comparisons

    3. Effect size: Log2 fold change
       - log2(mean_in / mean_out)
       - Interpretable magnitude of differential expression

    Args:
        adata: AnnData with normalized expression
        config: Pipeline configuration
        run_metadata: Optional dict to store DEG metadata

    Returns:
        Tuple of (top_markers_df, full_deg_df)
        - top_markers_df: Top n genes per cluster for reporting
        - full_deg_df: Complete DEG table for all genes
    """
    from scrna_agent.stats import run_full_deg_analysis

    logger.info(f"Finding marker genes using {config.deg_method} test...")

    cluster_key = 'leiden' if 'leiden' in adata.obs else 'louvain'

    # Run full DEG analysis with proper statistics
    full_deg_df = run_full_deg_analysis(
        adata,
        groupby=cluster_key,
        method=config.deg_method,
        corr_method=config.deg_correction,
        min_fold_change=config.min_fold_change
    )

    # Create top markers table (for report and quick reference)
    top_markers_list = []
    for cluster in full_deg_df['cluster'].unique():
        cluster_df = full_deg_df[full_deg_df['cluster'] == cluster]
        # Filter to significant genes and take top n
        sig_genes = cluster_df[cluster_df['significant']].head(config.n_genes_per_cluster)

        if len(sig_genes) < config.n_genes_per_cluster:
            # If not enough significant, take top by p-value
            sig_genes = cluster_df.head(config.n_genes_per_cluster)

        for rank, (_, row) in enumerate(sig_genes.iterrows(), 1):
            top_markers_list.append({
                'cluster': cluster,
                'gene': row['gene'],
                'log2fc': row['log2fc'],
                'pval': row['pval'],
                'pval_adj': row['pval_adj'],
                'score': row['score'],
                'pct_in': row['pct_in'],
                'pct_out': row['pct_out'],
                'significant': row['significant'],
                'rank': rank
            })

    top_markers_df = pd.DataFrame(top_markers_list)

    # Summary statistics
    n_sig = full_deg_df['significant'].sum()
    n_genes_tested = len(full_deg_df['gene'].unique())

    # Store in run_metadata
    if run_metadata is not None:
        run_metadata['deg_method'] = config.deg_method
        run_metadata['deg_correction'] = config.deg_correction
        run_metadata['min_fold_change'] = config.min_fold_change
        run_metadata['n_significant_degs'] = int(n_sig)
        run_metadata['n_genes_tested'] = n_genes_tested

    logger.info(f"DEG analysis complete: {n_sig} significant genes out of "
                f"{n_genes_tested} tested (FDR < 0.05, |log2FC| >= {config.min_fold_change})")

    return top_markers_df, full_deg_df


def generate_plots(
    adata: sc.AnnData,
    output_dir: str,
    config: PipelineConfig
) -> list:
    """Generate UMAP, QC, and PC selection plots."""
    logger.info("Generating plots...")

    plots_created = []
    cluster_key = 'leiden' if 'leiden' in adata.obs else 'louvain'

    # UMAP colored by clusters
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=cluster_key, ax=ax, show=False, title='Clusters')
    umap_path = os.path.join(output_dir, 'umap.png')
    plt.savefig(umap_path, dpi=150, bbox_inches='tight')
    plt.close()
    plots_created.append(umap_path)

    # QC distribution plots
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    ax = axes[0]
    ax.hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Genes per cell')
    ax.set_ylabel('Count')
    ax.set_title('Genes per Cell (post-QC)')
    ax.axvline(adata.obs['n_genes_by_counts'].median(), color='red', linestyle='--',
               label=f"Median: {adata.obs['n_genes_by_counts'].median():.0f}")
    ax.legend()

    ax = axes[1]
    ax.hist(adata.obs['total_counts'], bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('Total counts')
    ax.set_ylabel('Count')
    ax.set_title('Total Counts per Cell (post-QC)')
    ax.axvline(adata.obs['total_counts'].median(), color='red', linestyle='--',
               label=f"Median: {adata.obs['total_counts'].median():.0f}")
    ax.legend()

    ax = axes[2]
    ax.hist(adata.obs['pct_counts_mt'], bins=50, edgecolor='black', alpha=0.7)
    ax.set_xlabel('MT %')
    ax.set_ylabel('Count')
    ax.set_title('Mitochondrial % (post-QC)')
    ax.axvline(adata.obs['pct_counts_mt'].median(), color='red', linestyle='--',
               label=f"Median: {adata.obs['pct_counts_mt'].median():.1f}%")
    ax.legend()

    plt.tight_layout()
    qc_path = os.path.join(output_dir, 'qc_plots.png')
    plt.savefig(qc_path, dpi=150, bbox_inches='tight')
    plt.close()
    plots_created.append(qc_path)

    # PC variance explained plot (elbow plot)
    if 'pca' in adata.uns and 'variance_ratio' in adata.uns['pca']:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        variance_ratio = adata.uns['pca']['variance_ratio']
        n_pcs = len(variance_ratio)
        cum_var = np.cumsum(variance_ratio)

        # Get elbow point if available
        elbow_point = None
        n_pcs_selected = None
        if 'pc_selection' in adata.uns:
            elbow_point = adata.uns['pc_selection'].get('elbow_point')
            n_pcs_selected = adata.uns['pc_selection'].get('n_pcs_selected')

        # Scree plot (individual variance)
        ax = axes[0]
        ax.bar(range(1, n_pcs + 1), variance_ratio * 100, alpha=0.7, color='steelblue')
        ax.set_xlabel('Principal Component')
        ax.set_ylabel('Variance Explained (%)')
        ax.set_title('Scree Plot')
        if elbow_point:
            ax.axvline(elbow_point, color='red', linestyle='--',
                       label=f'Elbow: PC{elbow_point}')
            ax.legend()

        # Cumulative variance plot
        ax = axes[1]
        ax.plot(range(1, n_pcs + 1), cum_var * 100, 'b-', linewidth=2)
        ax.scatter(range(1, n_pcs + 1), cum_var * 100, color='steelblue', s=20)
        ax.set_xlabel('Number of Principal Components')
        ax.set_ylabel('Cumulative Variance Explained (%)')
        ax.set_title('Cumulative Variance Explained')
        ax.axhline(90, color='gray', linestyle=':', alpha=0.5, label='90% threshold')
        if n_pcs_selected and n_pcs_selected <= len(cum_var):
            ax.axvline(n_pcs_selected, color='red', linestyle='--',
                       label=f'Selected: {n_pcs_selected} PCs ({cum_var[n_pcs_selected-1]*100:.1f}%)')
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        pca_path = os.path.join(output_dir, 'pca_variance.png')
        plt.savefig(pca_path, dpi=150, bbox_inches='tight')
        plt.close()
        plots_created.append(pca_path)

    logger.info(f"Created {len(plots_created)} plots")
    return plots_created


def generate_report(
    adata: sc.AnnData,
    qc_summary: pd.DataFrame,
    clusters: pd.DataFrame,
    markers: pd.DataFrame,
    output_dir: str,
    config: PipelineConfig,
    run_metadata: Optional[Dict[str, Any]] = None
) -> str:
    """Generate analysis report with statistical rationale."""
    logger.info("Generating report with statistical rationale...")

    cluster_key = 'leiden' if 'leiden' in adata.obs else 'louvain'
    n_clusters = adata.obs[cluster_key].nunique()

    # Get top markers per cluster
    top_markers_per_cluster = markers.groupby('cluster').head(5)

    # Get PC selection info
    n_pcs_used = run_metadata.get('n_pcs_used', config.n_pcs) if run_metadata else config.n_pcs
    n_hvgs = run_metadata.get('n_hvgs', config.n_top_genes) if run_metadata else config.n_top_genes

    report_lines = [
        "# scRNA-seq Analysis Report",
        "",
        "## Summary",
        "",
        f"- **Cells analyzed**: {adata.n_obs:,}",
        f"- **Genes (total)**: {adata.n_vars:,}",
        f"- **Highly variable genes**: {n_hvgs}",
        f"- **Principal components used**: {n_pcs_used}",
        f"- **Clusters found**: {n_clusters}",
        f"- **Clustering method**: {config.clustering_method}",
        f"- **Resolution**: {config.resolution}",
        "",
    ]

    # Add statistical rationale sections
    if run_metadata:
        # QC Rationale
        if 'qc_rationale' in run_metadata:
            report_lines.extend([
                "## Statistical Methods",
                "",
                run_metadata['qc_rationale'],
                "",
            ])

        # PC Selection Rationale
        if 'pc_rationale' in run_metadata:
            report_lines.extend([
                run_metadata['pc_rationale'],
                "",
            ])

    # QC Summary Table
    report_lines.extend([
        "## Quality Control Summary",
        "",
        "| Metric | Value | Description |",
        "|--------|-------|-------------|",
    ])

    for _, row in qc_summary.iterrows():
        val = row['value']
        desc = row.get('description', '')
        if isinstance(val, float):
            if abs(val) < 0.01 or abs(val) > 10000:
                val = f"{val:.2e}"
            else:
                val = f"{val:.2f}"
        report_lines.append(f"| {row['metric']} | {val} | {desc} |")

    # Cluster Distribution
    report_lines.extend([
        "",
        "## Cluster Distribution",
        "",
        "| Cluster | Cells | Percentage |",
        "|---------|-------|------------|",
    ])

    cluster_counts = adata.obs[cluster_key].value_counts().sort_index()
    for cluster, count in cluster_counts.items():
        pct = count / adata.n_obs * 100
        report_lines.append(f"| {cluster} | {count:,} | {pct:.1f}% |")

    # Enhanced marker table with more statistics
    report_lines.extend([
        "",
        "## Top Marker Genes",
        "",
        "| Cluster | Gene | Log2FC | Adj. P-value | % In | % Out |",
        "|---------|------|--------|--------------|------|-------|",
    ])

    for _, row in top_markers_per_cluster.iterrows():
        pval_str = f"{row['pval_adj']:.2e}" if row['pval_adj'] < 0.01 else f"{row['pval_adj']:.4f}"
        pct_in = row.get('pct_in', 'N/A')
        pct_out = row.get('pct_out', 'N/A')
        if isinstance(pct_in, (int, float)):
            pct_in = f"{pct_in:.1f}%"
        if isinstance(pct_out, (int, float)):
            pct_out = f"{pct_out:.1f}%"
        report_lines.append(
            f"| {row['cluster']} | {row['gene']} | {row['log2fc']:.2f} | {pval_str} | {pct_in} | {pct_out} |"
        )

    # DEG Analysis Statistics
    if run_metadata and 'n_significant_degs' in run_metadata:
        report_lines.extend([
            "",
            "### Differential Expression Statistics",
            "",
            f"- **Test method**: {run_metadata.get('deg_method', 'wilcoxon')} rank-sum test",
            f"- **Multiple testing correction**: {run_metadata.get('deg_correction', 'BH')}",
            f"- **Significance threshold**: FDR < 0.05, |log2FC| >= {run_metadata.get('min_fold_change', 0.25)}",
            f"- **Genes tested**: {run_metadata.get('n_genes_tested', 'N/A')}",
            f"- **Significant DEGs**: {run_metadata.get('n_significant_degs', 'N/A')}",
            "",
        ])

    # Add tissue annotation section if available
    if "celltype" in adata.obs.columns:
        try:
            from scrna_agent.tissue import generate_tissue_report_section
            tissue_section = generate_tissue_report_section(adata)
            report_lines.extend(["", tissue_section])
        except ImportError:
            pass

    # Output files
    report_lines.extend([
        "",
        "## Output Files",
        "",
        "| File | Description |",
        "|------|-------------|",
        "| `qc_summary.csv` | QC metrics and thresholds |",
        "| `clusters.csv` | Cell cluster assignments |",
        "| `markers.csv` | Top marker genes per cluster |",
        "| `markers_full.csv` | Complete DEG table (all genes, all clusters) |",
        "| `hvg_list.txt` | List of highly variable genes used |",
        "| `umap.png` | UMAP visualization |",
        "| `qc_plots.png` | QC metric distributions |",
        "| `pca_variance.png` | PCA explained variance (elbow plot) |",
        "| `run_metadata.json` | Complete analysis metadata |",
        "| `processed.h5ad` | Processed AnnData object |",
    ])

    # Add tissue-specific outputs if present
    if "celltype" in adata.obs.columns:
        report_lines.extend([
            "| `cell_metadata.csv` | Cell metadata with tissue and cell type |",
            "| `tissue_summary.csv` | Cell type counts per tissue |",
        ])

    # Reproducibility section
    report_lines.extend([
        "",
        "## Reproducibility",
        "",
        f"- **Random seed**: {config.seed}",
        f"- **Analysis timestamp**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "All modeling choices are recorded in `run_metadata.json` for full reproducibility.",
        "",
    ])

    # Parameters
    report_lines.extend([
        "## Parameters Used",
        "",
        "### QC Parameters",
        f"- Adaptive QC: {config.adaptive_qc}",
        f"- QC method: {config.qc_method}",
        f"- Number of MADs: {config.qc_nmads}",
        f"- Min genes floor: {config.min_genes}",
        f"- Min cells per gene: {config.min_cells}",
        f"- Max MT% ceiling: {config.max_mt_pct}",
        "",
        "### Normalization Parameters",
        f"- Target sum: {config.target_sum}",
        f"- HVG method: {config.hvg_flavor}",
        f"- Number of HVGs: {config.n_top_genes}",
        "",
        "### Dimensionality Reduction",
        f"- Max PCs computed: {config.n_pcs}",
        f"- PC selection method: {config.pc_selection_method}",
        f"- Variance threshold: {config.variance_threshold}",
        f"- k-neighbors: {config.n_neighbors}",
        "",
        "### Clustering",
        f"- Method: {config.clustering_method}",
        f"- Resolution: {config.resolution}",
        "",
        "### DEG Analysis",
        f"- Method: {config.deg_method}",
        f"- Correction: {config.deg_correction}",
        f"- Min fold change: {config.min_fold_change}",
    ])

    if config.tissue:
        report_lines.extend([
            "",
            "### Tissue/Annotation",
            f"- Tissue: {config.tissue}",
        ])

    report_lines.extend([
        "",
        "---",
        "*Generated by scrna-agent with statistical rigor*",
    ])

    report_content = "\n".join(report_lines)

    if config.report_format == "html":
        # Simple markdown to HTML conversion
        try:
            import markdown
            html_content = markdown.markdown(report_content, extensions=['tables'])
            html_full = f"""<!DOCTYPE html>
<html>
<head>
    <title>scRNA-seq Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; max-width: 900px; margin: 0 auto; padding: 20px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        tr:nth-child(even) {{ background-color: #f2f2f2; }}
    </style>
</head>
<body>
{html_content}
</body>
</html>"""
            report_path = os.path.join(output_dir, "report.html")
            with open(report_path, 'w') as f:
                f.write(html_full)
        except ImportError:
            # Fallback to markdown if markdown package not available
            report_path = os.path.join(output_dir, "report.md")
            with open(report_path, 'w') as f:
                f.write(report_content)
    else:
        report_path = os.path.join(output_dir, "report.md")
        with open(report_path, 'w') as f:
            f.write(report_content)

    logger.info(f"Report saved to {report_path}")
    return report_path


def run_pipeline(
    input_path: str,
    output_dir: str,
    config: Optional[PipelineConfig] = None,
    input_format: str = "auto"
) -> PipelineResults:
    """
    Run the complete scRNA-seq analysis pipeline with statistical rigor.

    This pipeline implements:
    1. Adaptive QC thresholds using robust statistics (MAD)
    2. Variance-aware normalization with size-factor correction
    3. Reproducible HVG selection (Seurat v3)
    4. Model-driven PC selection (elbow method)
    5. Graph-based clustering (Leiden/Louvain)
    6. Full differential expression with proper statistics
    7. Complete metadata recording for reproducibility

    Args:
        input_path: Path to input data (h5ad file or 10x directory)
        output_dir: Directory for output files
        config: Pipeline configuration (uses defaults if None)
        input_format: "h5ad", "10x", or "auto"

    Returns:
        PipelineResults with all outputs and metadata
    """
    if config is None:
        config = PipelineConfig()

    # Setup output directory
    os.makedirs(output_dir, exist_ok=True)

    results = PipelineResults(output_dir=output_dir)

    # Initialize run metadata for tracking all modeling choices
    run_metadata = {
        'pipeline_version': '2.0',
        'timestamp': datetime.now().isoformat(),
        'input_path': str(input_path),
        'input_format': input_format,
        'seed': config.seed,
    }

    logger.info("=" * 60)
    logger.info("scRNA-seq Analysis Pipeline (Statistical Mode)")
    logger.info("=" * 60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_dir}")
    logger.info("=" * 60)

    # 1. Load data
    adata = load_data(input_path, input_format)
    run_metadata['input_cells'] = adata.n_obs
    run_metadata['input_genes'] = adata.n_vars

    # 2. QC with adaptive thresholds
    qc_summary = run_qc(adata, config, run_metadata)
    qc_path = os.path.join(output_dir, "qc_summary.csv")
    qc_summary.to_csv(qc_path, index=False)
    results.qc_summary = qc_summary
    results.files_created.append(qc_path)

    # 3. Normalize and scale (with HVG selection)
    hvg_list = normalize_and_scale(adata, config, run_metadata)
    results.hvg_list = hvg_list

    # Save HVG list
    hvg_path = os.path.join(output_dir, "hvg_list.txt")
    with open(hvg_path, 'w') as f:
        f.write('\n'.join(hvg_list))
    results.files_created.append(hvg_path)

    # 4. Dimensionality reduction with model-driven PC selection
    n_pcs_used = run_dimensionality_reduction(adata, config, run_metadata)

    # 5. Clustering
    clusters = run_clustering(adata, config, run_metadata)
    clusters_path = os.path.join(output_dir, "clusters.csv")
    clusters.to_csv(clusters_path, index=False)
    results.clusters = clusters
    results.files_created.append(clusters_path)

    # 6. Find markers with full statistical output
    markers, markers_full = find_markers(adata, config, run_metadata)

    # Save top markers
    markers_path = os.path.join(output_dir, "markers.csv")
    markers.to_csv(markers_path, index=False)
    results.markers = markers
    results.files_created.append(markers_path)

    # Save full DEG table
    markers_full_path = os.path.join(output_dir, "markers_full.csv")
    markers_full.to_csv(markers_full_path, index=False)
    results.markers_full = markers_full
    results.files_created.append(markers_full_path)

    # 7. Tissue assignment and annotation
    if config.tissue is not None or config.tissue_map_path is not None or config.run_annotation:
        from scrna_agent.tissue import (
            load_tissue_map, assign_tissue, annotate_cells,
            generate_cell_metadata, generate_tissue_summary,
            generate_tissue_report_section
        )

        # Load tissue map if provided
        tissue_map = None
        if config.tissue_map_path:
            tissue_map = load_tissue_map(config.tissue_map_path, config.sample_col)

        # Assign tissue
        assign_tissue(
            adata,
            tissue=config.tissue,
            tissue_map=tissue_map,
            sample_col=config.sample_col,
            allow_missing=config.allow_missing_tissue
        )

        # Run annotation
        if config.run_annotation:
            if config.use_multi_method_annotation:
                # NEW: Multi-method annotation with biological validation
                # Key insight from JMV_ALI: organoid samples need context-aware annotation
                logger.info(f"Running multi-method annotation (sample_type={config.sample_type})...")
                try:
                    from scrna_agent.tissue import annotate_with_validation
                    adata, validation = annotate_with_validation(
                        adata,
                        sample_type=config.sample_type,
                        celltypist_model=config.celltypist_model,
                        output_dir=output_dir
                    )
                    # Store validation results
                    run_metadata['annotation_validation'] = {
                        'valid': validation['valid'],
                        'sample_type': validation['sample_type'],
                        'issues': len(validation['issues']),
                        'summary': validation['summary']
                    }
                    # Use consensus_celltype as primary celltype
                    if 'consensus_celltype' in adata.obs.columns:
                        adata.obs['celltype'] = adata.obs['consensus_celltype']
                except ImportError as e:
                    logger.warning(f"Multi-method annotation not available: {e}")
                    logger.warning("Falling back to standard marker-based annotation")
                    annotate_cells(adata, markers_path=config.markers_path)
            else:
                # Standard marker-based annotation
                logger.info("Running tissue-aware marker-based annotation...")
                annotate_cells(adata, markers_path=config.markers_path)

            # Generate and save cell metadata
            cell_metadata = generate_cell_metadata(adata, config.sample_col)
            cell_metadata_path = os.path.join(output_dir, "cell_metadata.csv")
            cell_metadata.to_csv(cell_metadata_path)
            results.cell_metadata = cell_metadata
            results.files_created.append(cell_metadata_path)

            # Generate and save tissue summary
            if "tissue" in adata.obs.columns and "celltype" in adata.obs.columns:
                tissue_summary = generate_tissue_summary(adata)
                tissue_summary_path = os.path.join(output_dir, "tissue_summary.csv")
                tissue_summary.to_csv(tissue_summary_path)
                results.tissue_summary = tissue_summary
                results.files_created.append(tissue_summary_path)

    # 8. Generate plots (including PCA variance plot)
    plot_paths = generate_plots(adata, output_dir, config)
    results.files_created.extend(plot_paths)

    # 9. Save processed data with all embeddings
    # Verify embeddings are stored in obsm
    logger.info("Stored embeddings in adata.obsm:")
    for key in adata.obsm.keys():
        logger.info(f"  - {key}: {adata.obsm[key].shape}")

    h5ad_path = os.path.join(output_dir, "processed.h5ad")
    adata.write_h5ad(h5ad_path)
    results.files_created.append(h5ad_path)

    # 10. Save run metadata (remove large arrays for JSON serialization)
    run_metadata_clean = {}
    for key, value in run_metadata.items():
        if key == 'hvg_list':
            # Store count instead of full list
            run_metadata_clean['hvg_count'] = len(value)
            run_metadata_clean['hvg_sample'] = value[:10] if len(value) > 10 else value
        elif key in ['qc_distributions', 'pc_selection']:
            # Keep these but ensure they're JSON serializable
            run_metadata_clean[key] = value
        elif key in ['qc_rationale', 'pc_rationale']:
            # Keep rationale strings
            run_metadata_clean[key] = value
        elif isinstance(value, (str, int, float, bool, type(None))):
            run_metadata_clean[key] = value
        elif isinstance(value, dict):
            # Recursively clean dict
            run_metadata_clean[key] = {
                k: v for k, v in value.items()
                if isinstance(v, (str, int, float, bool, type(None), list))
            }

    # Add summary statistics
    run_metadata_clean['summary'] = {
        'cells_analyzed': adata.n_obs,
        'genes_analyzed': adata.n_vars,
        'n_hvgs': len(hvg_list),
        'n_pcs_used': n_pcs_used,
        'n_clusters': run_metadata.get('n_clusters', 0),
        'clustering_method': config.clustering_method,
        'resolution': config.resolution,
    }

    metadata_path = os.path.join(output_dir, "run_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(run_metadata_clean, f, indent=2, default=str)
    results.files_created.append(metadata_path)
    results.run_metadata = run_metadata

    # 11. Generate report with statistical rationale
    report_path = generate_report(
        adata, qc_summary, clusters, markers, output_dir, config, run_metadata
    )
    results.files_created.append(report_path)

    results.adata = adata

    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info(f"Files created: {len(results.files_created)}")
    logger.info("=" * 60)

    return results


def create_toy_dataset(
    n_cells: int = 500,
    n_genes: int = 1000,
    n_clusters: int = 5,
    output_path: Optional[str] = None,
    seed: int = 42
) -> sc.AnnData:
    """
    Create a synthetic scRNA-seq dataset for testing.

    Args:
        n_cells: Number of cells
        n_genes: Number of genes
        n_clusters: Number of cell clusters
        output_path: If provided, save to this path
        seed: Random seed

    Returns:
        AnnData object with synthetic data
    """
    np.random.seed(seed)

    # Generate cluster assignments
    cluster_sizes = np.random.multinomial(n_cells, [1/n_clusters] * n_clusters)
    clusters = np.repeat(range(n_clusters), cluster_sizes)
    np.random.shuffle(clusters)

    # Generate base expression matrix
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)

    # Add cluster-specific patterns
    n_marker_genes = n_genes // n_clusters
    for c in range(n_clusters):
        cluster_mask = clusters == c
        start_gene = c * n_marker_genes
        end_gene = start_gene + n_marker_genes
        X[cluster_mask, start_gene:end_gene] += np.random.negative_binomial(
            10, 0.3, size=(cluster_mask.sum(), n_marker_genes)
        )

    # Create gene names
    gene_names = [f"Gene_{i}" for i in range(n_genes)]
    # Add some MT genes
    for i in range(min(20, n_genes)):
        gene_names[i] = f"MT-Gene_{i}"

    # Create cell barcodes
    cell_barcodes = [f"Cell_{i:04d}" for i in range(n_cells)]

    # Create AnnData
    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame({
            'true_cluster': clusters.astype(str),
            'batch': np.random.choice(['batch_A', 'batch_B'], n_cells)
        }, index=cell_barcodes),
        var=pd.DataFrame(index=gene_names)
    )

    if output_path:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        adata.write_h5ad(output_path)
        logger.info(f"Saved toy dataset to {output_path}")

    return adata
