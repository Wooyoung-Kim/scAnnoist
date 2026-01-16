#!/usr/bin/env python
"""
Complete scRNA-seq Analysis Pipeline with Presentation Generation
Automates: preprocessing → annotation → visualization → presentation
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import json
import sys
import os
from pathlib import Path

# Set matplotlib backend
import matplotlib
matplotlib.use('Agg')
sc.settings.set_figure_params(dpi=150, facecolor='white', frameon=True)

sys.path.insert(0, '/mnt2/kwy/scrna_agent')

from scrna_agent.tools.model_management import (
    validate_and_add_tissue_metadata,
    get_models_for_adata,
    load_model_registry
)
from scrna_agent.tools.report_generator_v2 import generate_report_v2
from scrna_agent.tools.marker_annotation import run_marker_annotation, compare_annotations
import celltypist
from celltypist import models


def run_preprocessing(data_dir: str, output_dir: str, tissue_type: str = "Lung",
                      sample_type: str = "organoid") -> dict:
    """
    Step 1: Load and preprocess data

    Args:
        data_dir: Directory containing CellRanger outputs
        output_dir: Output directory for results
        tissue_type: Tissue type for annotation
        sample_type: Sample type description

    Returns:
        dict with preprocessing results
    """
    print("\n" + "=" * 80)
    print("STEP 1: PREPROCESSING")
    print("=" * 80)

    figures_dir = f"{output_dir}/figures"
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)

    # Find samples
    print("\nFinding samples...")
    sample_dirs = sorted([d for d in Path(data_dir).iterdir()
                         if d.is_dir() and (d / "outs" / "filtered_feature_bc_matrix.h5").exists()])
    print(f"Found {len(sample_dirs)} samples")

    # Load samples
    print("\nLoading samples...")
    adatas = []
    sample_names = []

    for sample_dir in sample_dirs:
        h5_file = sample_dir / "outs" / "filtered_feature_bc_matrix.h5"
        print(f"  Loading {sample_dir.name}...")
        adata = sc.read_10x_h5(h5_file)
        adata.var_names_make_unique()
        adata.obs['sample'] = sample_dir.name
        adatas.append(adata)
        sample_names.append(sample_dir.name)

    # Concatenate
    print(f"\nConcatenating {len(adatas)} samples...")
    adata = sc.concat(adatas, label="sample", keys=sample_names)
    print(f"  Total: {adata.n_obs:,} cells, {adata.n_vars:,} genes")

    # QC
    print("\nRunning QC...")
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # QC plots
    fig = plt.figure(figsize=(15, 4))
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                 jitter=0.4, multi_panel=True, show=False)
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/qc_metrics.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved QC plots")

    # Filter
    print("\nFiltering cells...")
    n_cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, min_counts=500)
    adata = adata[adata.obs['pct_counts_mt'] < 20, :]
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  {n_cells_before:,} → {adata.n_obs:,} cells")
    print(f"  {adata.n_vars:,} genes retained")

    # Preprocessing
    print("\nNormalizing and finding variable genes...")
    adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # Dimensionality reduction
    print("\nRunning PCA, UMAP, and clustering...")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    print(f"  Found {adata.obs['leiden'].nunique()} clusters")

    # Save UMAP plots
    print("\nSaving visualizations...")
    fig, ax = plt.subplots(figsize=(8, 6))
    sc.pl.umap(adata, color='leiden', legend_loc='on data',
               title='UMAP - Leiden Clusters', ax=ax, show=False, frameon=True)
    plt.savefig(f"{figures_dir}/umap_clusters.png", dpi=150, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(10, 6))
    sc.pl.umap(adata, color='sample', ax=ax, show=False, frameon=True,
               title='UMAP - By Sample')
    plt.savefig(f"{figures_dir}/umap_samples.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved UMAP plots")

    # Save preprocessed data
    print("\nSaving preprocessed data...")
    preprocessed_file = f"{output_dir}/preprocessed.h5ad"
    adata.write(preprocessed_file)

    # Save tissue mapping
    sample_info = [{'sample': s, 'tissue': tissue_type, 'sample_type': sample_type}
                   for s in sample_names]
    tissue_map_df = pd.DataFrame(sample_info)
    tissue_map_file = f"{output_dir}/tissue_map.csv"
    tissue_map_df.to_csv(tissue_map_file, index=False)

    print(f"✓ Preprocessing complete!")
    print(f"  Output: {preprocessed_file}")

    return {
        'adata': adata,
        'n_cells': adata.n_obs,
        'n_clusters': adata.obs['leiden'].nunique(),
        'preprocessed_file': preprocessed_file,
        'tissue_map_file': tissue_map_file
    }


def run_annotation(output_dir: str, tissue_type: str = "Lung") -> dict:
    """
    Step 2: Cell type annotation with CellTypist + Marker-based

    Args:
        output_dir: Output directory containing preprocessed data
        tissue_type: Tissue type for marker-based annotation

    Returns:
        dict with annotation results
    """
    print("\n" + "=" * 80)
    print("STEP 2: CELL TYPE ANNOTATION")
    print("=" * 80)

    celltypist_output = f"{output_dir}/celltypist_outputs"
    figures_dir = f"{output_dir}/figures"
    os.makedirs(celltypist_output, exist_ok=True)

    # Load data
    print("\nLoading preprocessed data...")
    adata = sc.read_h5ad(f"{output_dir}/preprocessed.h5ad")
    tissue_map = pd.read_csv(f"{output_dir}/tissue_map.csv")

    # Add tissue metadata
    print("\nAdding tissue metadata...")
    adata = validate_and_add_tissue_metadata(adata, tissue_map)

    # Get models
    print("\nDetermining models for each tissue...")
    tissue_info = get_models_for_adata(adata)

    # Run annotation
    print("\nRunning CellTypist annotation...")
    all_preds_list = []
    tissue_usage_list = []

    for tissue, info in tissue_info.items():
        model_name = info['model']
        print(f"\n  Annotating tissue: {tissue}")
        print(f"    Model: {model_name}")

        # Download model if needed
        models.download_models(model=[model_name])

        # Subset and normalize
        adata_tissue = adata[adata.obs['tissue'] == tissue].copy()

        if hasattr(adata_tissue, 'raw') and adata_tissue.raw is not None:
            adata_for_ct = adata_tissue.raw.to_adata()
            sc.pp.normalize_total(adata_for_ct, target_sum=1e4)
            sc.pp.log1p(adata_for_ct)
        else:
            adata_for_ct = adata_tissue.copy()

        # Annotate
        predictions = celltypist.annotate(
            adata_for_ct,
            model=model_name,
            majority_voting=True
        )

        # Extract results
        pred_df = predictions.predicted_labels
        pred_df['tissue'] = tissue
        pred_df['model_used'] = model_name
        pred_df['sample'] = adata_tissue.obs['sample'].values
        all_preds_list.append(pred_df)

        # Track usage
        tissue_usage_list.append({
            'tissue': tissue,
            'model_used': model_name,
            'selection_rule': info['rule'],
            'n_cells': adata_tissue.n_obs
        })

        print(f"    ✓ Annotated {adata_tissue.n_obs:,} cells")
        print(f"    Found {pred_df['majority_voting'].nunique()} cell types")

    # Merge annotations
    print("\nMerging annotations...")
    all_preds_df = pd.concat(all_preds_list)
    all_preds_df = all_preds_df.drop_duplicates(subset='cell_id', keep='first')

    pred_dict = {
        'label': dict(zip(all_preds_df.index, all_preds_df['majority_voting'])),
        'conf': dict(zip(all_preds_df.index, all_preds_df['conf_score'])),
        'model': dict(zip(all_preds_df.index, all_preds_df['model_used']))
    }

    adata.obs['celltypist_label'] = adata.obs.index.map(pred_dict['label'])
    adata.obs['celltypist_conf'] = adata.obs.index.map(pred_dict['conf'])
    adata.obs['celltypist_model'] = adata.obs.index.map(pred_dict['model'])

    # Save visualizations
    print("\nSaving visualizations...")

    # Cell type UMAP
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='celltypist_label', ax=ax, show=False, frameon=True,
               title='UMAP - Cell Types (CellTypist)', legend_fontsize=10)
    plt.savefig(f"{figures_dir}/umap_celltypes.png", dpi=150, bbox_inches='tight')
    plt.close()

    # Cell type distribution
    cell_counts = adata.obs['celltypist_label'].value_counts()
    fig, ax = plt.subplots(figsize=(10, 6))
    cell_counts.head(15).plot(kind='barh', ax=ax, color='#5EA8A7')
    ax.set_xlabel('Number of Cells')
    ax.set_ylabel('Cell Type')
    ax.set_title('Cell Type Distribution (Top 15)')
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/celltype_distribution.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved CellTypist visualizations")

    # Run marker-based annotation
    print(f"\nRunning marker-based annotation...")
    marker_results = run_marker_annotation(
        adata,
        tissue=tissue_type,
        cluster_key='leiden',
        output_col='marker_annotation'
    )

    # Marker-based UMAP
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='marker_annotation', ax=ax, show=False, frameon=True,
               title='UMAP - Cell Types (Marker-based)', legend_fontsize=10)
    plt.savefig(f"{figures_dir}/umap_marker.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved marker-based UMAP")

    # Compare annotations
    print(f"\nComparing annotation methods...")
    comparison_df = compare_annotations(
        adata,
        method1_col='celltypist_label',
        method2_col='marker_annotation',
        cluster_key='leiden'
    )
    comparison_df.to_csv(f"{celltypist_output}/annotation_comparison.csv", index=False)
    print(f"  ✓ Saved comparison results")

    # Create comparison visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color='celltypist_label', ax=ax1, show=False, frameon=True,
               title='CellTypist Annotation', legend_fontsize=8)
    sc.pl.umap(adata, color='marker_annotation', ax=ax2, show=False, frameon=True,
               title='Marker-based Annotation', legend_fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{figures_dir}/annotation_comparison.png", dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Saved comparison visualization")

    # Save outputs
    print("\nSaving results...")

    # CSV
    pred_output = all_preds_df.copy()
    pred_output['cell_id'] = pred_output.index
    pred_output = pred_output[['cell_id', 'sample', 'tissue', 'majority_voting', 'conf_score', 'model_used']]
    pred_output.columns = ['cell_id', 'sample', 'tissue', 'label', 'confidence', 'model_used']
    pred_output.to_csv(f"{celltypist_output}/celltypist_predictions.csv", index=False)

    # Tissue usage
    tissue_usage_df = pd.DataFrame(tissue_usage_list)
    tissue_usage_df.to_csv(f"{celltypist_output}/tissue_model_usage.csv", index=False)

    # Metadata
    metadata = {
        'total_cells': int(adata.n_obs),
        'total_tissues': len(tissue_info),
        'tissues': list(tissue_info.keys()),
        'tissue_model_mapping': {
            tissue: {'model': info['model'], 'rule': info['rule']}
            for tissue, info in tissue_info.items()
        },
        'celltypist': {
            'cell_types_identified': int(adata.obs['celltypist_label'].nunique()),
            'top_cell_types': adata.obs['celltypist_label'].value_counts().head(20).to_dict()
        },
        'marker_based': {
            'cell_types_identified': int(adata.obs['marker_annotation'].nunique()),
            'top_cell_types': adata.obs['marker_annotation'].value_counts().head(20).to_dict()
        },
        'comparison': {
            'agreement_rate': (adata.obs['celltypist_label'] == adata.obs['marker_annotation']).mean(),
            'disagreement_count': (adata.obs['celltypist_label'] != adata.obs['marker_annotation']).sum()
        },
        # Backward compatibility
        'cell_types_identified': int(adata.obs['celltypist_label'].nunique()),
        'top_cell_types': adata.obs['celltypist_label'].value_counts().head(20).to_dict()
    }
    with open(f"{celltypist_output}/run_metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"\n✓ Annotation complete!")
    print(f"  CellTypist: {metadata['celltypist']['cell_types_identified']} cell types")
    print(f"  Marker-based: {metadata['marker_based']['cell_types_identified']} cell types")
    print(f"  Agreement: {metadata['comparison']['agreement_rate']*100:.1f}%")

    return metadata


def run_report_generation(output_dir: str, title: str) -> str:
    """
    Step 3: Generate presentation

    Args:
        output_dir: Output directory containing analysis results
        title: Presentation title

    Returns:
        Path to generated presentation
    """
    print("\n" + "=" * 80)
    print("STEP 3: PRESENTATION GENERATION")
    print("=" * 80)

    celltypist_output = f"{output_dir}/celltypist_outputs"
    pptx_path = generate_report_v2(celltypist_output, title)

    print(f"\n✓ Presentation generated!")
    print(f"  Location: {pptx_path}")

    return str(pptx_path)


def run_full_pipeline(data_dir: str, output_dir: str, title: str,
                      tissue_type: str = "Lung", sample_type: str = "organoid") -> dict:
    """
    Run complete analysis pipeline

    Args:
        data_dir: Directory containing CellRanger outputs
        output_dir: Output directory for all results
        title: Presentation title
        tissue_type: Tissue type for annotation
        sample_type: Sample type description

    Returns:
        dict with all results
    """
    print("=" * 80)
    print("FULL scRNA-seq ANALYSIS PIPELINE")
    print("=" * 80)
    print(f"\nInput: {data_dir}")
    print(f"Output: {output_dir}")
    print(f"Tissue: {tissue_type}")
    print(f"Sample type: {sample_type}")

    # Step 1: Preprocessing
    preproc_results = run_preprocessing(data_dir, output_dir, tissue_type, sample_type)

    # Step 2: Annotation (with marker-based)
    annot_results = run_annotation(output_dir, tissue_type)

    # Step 3: Presentation
    pptx_path = run_report_generation(output_dir, title)

    print("\n" + "=" * 80)
    print("✓ PIPELINE COMPLETE!")
    print("=" * 80)
    print(f"\nResults:")
    print(f"  Cells analyzed: {annot_results['total_cells']:,}")
    print(f"  Cell types: {annot_results['cell_types_identified']}")
    print(f"  Presentation: {pptx_path}")

    return {
        'preprocessing': preproc_results,
        'annotation': annot_results,
        'presentation': pptx_path
    }


if __name__ == "__main__":
    # Example usage
    import argparse

    parser = argparse.ArgumentParser(description='Complete scRNA-seq analysis pipeline')
    parser.add_argument('data_dir', help='Directory containing CellRanger outputs')
    parser.add_argument('output_dir', help='Output directory')
    parser.add_argument('--title', default='scRNA-seq Analysis Results', help='Presentation title')
    parser.add_argument('--tissue', default='Lung', help='Tissue type')
    parser.add_argument('--sample-type', default='organoid', help='Sample type')

    args = parser.parse_args()

    results = run_full_pipeline(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        title=args.title,
        tissue_type=args.tissue,
        sample_type=args.sample_type
    )
