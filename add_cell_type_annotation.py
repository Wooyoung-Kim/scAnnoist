#!/usr/bin/env python3
"""
Add Cell Type Annotation to scRNA-seq Analysis
This script helps annotate clusters with cell types and generates an annotation UMAP
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

def annotate_cell_types(adata_file, annotation_dict, cluster_key='louvain_0.8',
                        output_dir=None, annotation_key='cell_type'):
    """
    Annotate clusters with cell types and generate visualization

    Parameters:
    -----------
    adata_file : str or Path
        Path to the h5ad file
    annotation_dict : dict
        Dictionary mapping cluster IDs to cell type names
        Example: {'0': 'Basal cells', '1': 'Ciliated cells', ...}
    cluster_key : str
        Column name in adata.obs containing cluster IDs
    output_dir : str or Path
        Directory to save figures (if None, uses adata file directory)
    annotation_key : str
        Column name for the cell type annotation

    Returns:
    --------
    adata : AnnData
        Annotated AnnData object
    """

    print("="*80)
    print("Cell Type Annotation")
    print("="*80)

    # Load data
    print(f"\nLoading data from: {adata_file}")
    adata = sc.read_h5ad(adata_file)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Clusters: {adata.obs[cluster_key].nunique()}")

    # Set output directory
    if output_dir is None:
        output_dir = Path(adata_file).parent
    else:
        output_dir = Path(output_dir)

    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(exist_ok=True)

    # Annotate clusters
    print(f"\nAnnotating clusters with cell types...")
    adata.obs[annotation_key] = adata.obs[cluster_key].astype(str).map(annotation_dict)

    # Fill any unmapped clusters
    unmapped = adata.obs[annotation_key].isna()
    if unmapped.any():
        unmapped_clusters = adata.obs.loc[unmapped, cluster_key].unique()
        print(f"  Warning: {len(unmapped_clusters)} clusters not annotated: {unmapped_clusters}")
        adata.obs.loc[unmapped, annotation_key] = 'Unknown'

    # Print annotation summary
    print(f"\nAnnotation Summary:")
    annotation_counts = adata.obs[annotation_key].value_counts()
    for cell_type, count in annotation_counts.items():
        percentage = (count / adata.n_obs) * 100
        print(f"  {cell_type}: {count:,} cells ({percentage:.1f}%)")

    # Generate annotation UMAP
    print(f"\nGenerating annotation UMAP...")

    # Main annotation UMAP
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # UMAP colored by cell type
    sc.pl.umap(adata, color=annotation_key, ax=axes[0], show=False,
              title='Cell Type Annotation', legend_loc='right margin',
              frameon=False)

    # UMAP colored by cluster (for reference)
    sc.pl.umap(adata, color=cluster_key, ax=axes[1], show=False,
              title=f'Clusters ({cluster_key})', legend_loc='right margin',
              frameon=False)

    plt.tight_layout()
    output_file = figures_dir / 'umap_cell_type_annotation.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

    # Cell type by condition (if condition exists)
    if 'condition' in adata.obs.columns:
        print(f"\nGenerating cell type by condition plots...")

        # Stacked bar plot
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Cell type composition by condition
        ct = pd.crosstab(adata.obs['condition'], adata.obs[annotation_key])
        ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

        ax = axes[0]
        ct_pct.plot(kind='bar', stacked=True, ax=ax, colormap='tab20')
        ax.set_xlabel('Condition')
        ax.set_ylabel('Percentage of cells')
        ax.set_title('Cell Type Composition by Condition')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax.tick_params(axis='x', rotation=45)

        # Cell counts
        ax = axes[1]
        ct.plot(kind='bar', ax=ax, colormap='tab20')
        ax.set_xlabel('Condition')
        ax.set_ylabel('Number of cells')
        ax.set_title('Cell Counts by Condition')
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax.tick_params(axis='x', rotation=45)

        plt.tight_layout()
        output_file = figures_dir / 'cell_type_by_condition.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file}")

        # Faceted UMAP by condition
        conditions = adata.obs['condition'].unique()
        n_conditions = len(conditions)
        fig, axes = plt.subplots(1, n_conditions, figsize=(6*n_conditions, 5))
        if n_conditions == 1:
            axes = [axes]

        for idx, condition in enumerate(sorted(conditions)):
            adata_subset = adata[adata.obs['condition'] == condition]
            sc.pl.umap(adata_subset, color=annotation_key,
                      ax=axes[idx], show=False,
                      title=f'{condition} (n={adata_subset.n_obs:,})',
                      legend_loc='on data' if idx == 0 else None,
                      frameon=False)

        plt.tight_layout()
        output_file = figures_dir / 'umap_annotation_by_condition_faceted.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file}")

    # UMAP with cluster labels
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color=annotation_key, legend_loc='on data',
              ax=ax, show=False, title='Cell Type Annotation',
              frameon=False)
    output_file = figures_dir / 'umap_cell_type_with_labels.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

    # Save annotated data
    print(f"\nSaving annotated data...")
    output_h5ad = Path(adata_file).parent / f"{Path(adata_file).stem}_annotated.h5ad"
    adata.write(output_h5ad)
    print(f"  Saved: {output_h5ad}")

    # Save annotation table
    annotation_table = adata.obs[[cluster_key, annotation_key]].drop_duplicates().sort_values(cluster_key)
    output_csv = Path(adata_file).parent / 'cluster_annotation.csv'
    annotation_table.to_csv(output_csv, index=False)
    print(f"  Saved: {output_csv}")

    # Save cell type counts
    cell_type_summary = pd.DataFrame({
        'Cell Type': annotation_counts.index,
        'Number of Cells': annotation_counts.values,
        'Percentage': [f"{(v/adata.n_obs*100):.1f}%" for v in annotation_counts.values]
    })
    output_csv = Path(adata_file).parent / 'cell_type_summary.csv'
    cell_type_summary.to_csv(output_csv, index=False)
    print(f"  Saved: {output_csv}")

    print("\n" + "="*80)
    print("Annotation Complete!")
    print("="*80)
    print(f"\nNew files created:")
    print(f"  - {output_h5ad.name}")
    print(f"  - cluster_annotation.csv")
    print(f"  - cell_type_summary.csv")
    print(f"  - umap_cell_type_annotation.png")
    print(f"  - umap_cell_type_with_labels.png")
    if 'condition' in adata.obs.columns:
        print(f"  - cell_type_by_condition.png")
        print(f"  - umap_annotation_by_condition_faceted.png")
    print("="*80)

    return adata


# Example usage for JMV ALI analysis
if __name__ == "__main__":
    # Define cell type annotations based on marker genes
    # This is an example - adjust based on your marker gene analysis

    annotation_dict = {
        # Example annotations - MODIFY THESE based on your marker gene analysis
        '0': 'Basal cells',              # KRT5, TP63
        '1': 'Ciliated cells',           # FOXJ1, RSPH1
        '2': 'Secretory cells',          # SCGB1A1, SCGB3A2
        '3': 'Goblet cells',             # MUC5AC, MUC5B
        '4': 'Ionocytes',                # FOXI1, CFTR
        '5': 'Basal cells (cycling)',    # KRT5, MKI67
        '6': 'Club cells',               # SCGB1A1
        '7': 'Intermediate cells',       # KRT5, SCGB1A1
        '8': 'Secretory cells',          # SCGB3A2
        '9': 'Ciliated cells (mature)',  # FOXJ1, PIFO
        '10': 'Basal cells',             # KRT5
        '11': 'Ionocytes',               # FOXI1
        '12': 'Neuroendocrine',          # CHGA, CHGB
        '13': 'Secretory cells',         # SCGB1A1
        '14': 'Basal cells',             # KRT5
        '15': 'Ciliated cells',          # FOXJ1
        '16': 'Rare cell type 1',        # Check markers
        '17': 'Rare cell type 2',        # Check markers
        '18': 'Rare cell type 3',        # Check markers
        '19': 'Rare cell type 4',        # Check markers
        '20': 'Rare cell type 5',        # Check markers
    }

    # Path to your h5ad file
    adata_file = '/mnt2/kwy/scrna_agent/JMV_ALI_results/JMV_ALI_harmony_integrated.h5ad'

    # Run annotation
    adata = annotate_cell_types(
        adata_file=adata_file,
        annotation_dict=annotation_dict,
        cluster_key='louvain_0.8',
        annotation_key='cell_type'
    )

    print("\nTo regenerate PowerPoint with annotation:")
    print("  python generate_ppt_report.py")
