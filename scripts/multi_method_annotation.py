#!/usr/bin/env python
"""
Multi-method cell type annotation with consensus decision.
Combines celltypist (ML-based) with marker-based annotation for ALI organoids.
"""

import sys
from pathlib import Path
import logging

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import celltypist
from celltypist import models

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Paths
INPUT_H5AD = Path("/mnt2/kwy/scrna_agent/results/jmv_ali_epithelial/integrated.h5ad")
OUTPUT_DIR = Path("/mnt2/kwy/scrna_agent/results/jmv_ali_epithelial/multi_method_annotation")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Airway epithelial markers (for validation)
EPITHELIAL_MARKERS = {
    'Basal': ['KRT5', 'TP63', 'KRT15', 'KRT14', 'NGFR'],
    'Ciliated': ['FOXJ1', 'DNAH5', 'TPPP3', 'CAPS', 'SNTN'],
    'Secretory_Club': ['SCGB1A1', 'SCGB3A2', 'CYP2F1', 'LYPD2'],
    'Goblet': ['MUC5AC', 'MUC5B', 'TFF3', 'SPDEF'],
    'AT2_like': ['SFTPC', 'SFTPA1', 'SFTPB', 'ABCA3'],
    'Ionocyte': ['FOXI1', 'CFTR', 'ATP6V1C2']
}


def run_celltypist_annotation(adata, model_name='Cells_Lung_Airway.pkl'):
    """
    Run celltypist annotation using specified lung model.
    """
    logger.info("=" * 60)
    logger.info(f"METHOD 1: CellTypist with {model_name}")
    logger.info("=" * 60)

    # Download model if needed
    logger.info(f"Loading model: {model_name}")
    try:
        model = models.Model.load(model=model_name)
    except:
        logger.info(f"Downloading {model_name}...")
        models.download_models(model=[model_name.replace('.pkl', '')])
        model = models.Model.load(model=model_name)

    logger.info(f"Model loaded: {len(model.cell_types)} cell types")
    logger.info(f"Top cell types in model: {list(model.cell_types)[:20]}")

    # Run prediction
    logger.info("Running celltypist prediction...")
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=True  # Use majority voting for robustness
    )

    # Get results
    adata_ct = predictions.to_adata()

    # Summary
    celltypist_types = adata_ct.obs['majority_voting'].value_counts()
    logger.info(f"\nCellTypist results ({len(celltypist_types)} cell types):")
    for ct, count in celltypist_types.head(10).items():
        pct = count / len(adata_ct) * 100
        logger.info(f"  {ct}: {count} cells ({pct:.1f}%)")

    return adata_ct


def run_marker_based_annotation(adata):
    """
    Run marker-based annotation for epithelial subtypes.
    """
    logger.info("\n" + "=" * 60)
    logger.info("METHOD 2: Marker-Based Annotation")
    logger.info("=" * 60)

    # Score each cell type
    from scipy.sparse import issparse

    scores = {}
    for celltype, markers in EPITHELIAL_MARKERS.items():
        present_markers = [m for m in markers if m in adata.var_names]
        if not present_markers:
            scores[celltype] = np.zeros(adata.n_obs)
            continue

        logger.info(f"  {celltype}: {len(present_markers)}/{len(markers)} markers found")

        # Get expression
        expr = adata[:, present_markers].X
        if issparse(expr):
            expr = expr.toarray()

        # Mean expression across markers
        scores[celltype] = expr.mean(axis=1)

    # Assign cell types
    scores_df = pd.DataFrame(scores, index=adata.obs_names)
    marker_based_types = scores_df.idxmax(axis=1)
    marker_based_confidence = scores_df.max(axis=1)

    # Summary
    marker_counts = marker_based_types.value_counts()
    logger.info(f"\nMarker-based results ({len(marker_counts)} cell types):")
    for ct, count in marker_counts.items():
        pct = count / len(adata) * 100
        mean_conf = marker_based_confidence[marker_based_types == ct].mean()
        logger.info(f"  {ct}: {count} cells ({pct:.1f}%), avg score={mean_conf:.3f}")

    return marker_based_types, marker_based_confidence


def consensus_annotation(adata_ct, marker_types, marker_confidence):
    """
    Create consensus annotation combining celltypist and marker-based results.
    Implements "agent discussion" logic.
    """
    logger.info("\n" + "=" * 60)
    logger.info("CONSENSUS ANNOTATION (Agent Discussion)")
    logger.info("=" * 60)

    consensus_types = []
    consensus_confidence = []
    agreement = []

    celltypist_types = adata_ct.obs['majority_voting']
    celltypist_probs = adata_ct.obs['conf_score']  # Confidence from celltypist

    for idx in range(len(adata_ct)):
        ct_type = celltypist_types.iloc[idx]
        ct_conf = celltypist_probs.iloc[idx]
        mb_type = marker_types.iloc[idx]
        mb_conf = marker_confidence.iloc[idx]

        # Decision logic
        # 1. If both methods agree on epithelial subtype - use it
        if mb_type in ct_type or ct_type in mb_type:
            final_type = ct_type
            final_conf = "High"
            agree = "Agree"

        # 2. If marker-based has high confidence - prioritize it (organoid-specific)
        elif mb_conf > 0.5:
            final_type = mb_type
            final_conf = "Medium" if ct_conf < 0.5 else "High"
            agree = "Marker_Priority"

        # 3. If celltypist has high confidence and identifies epithelial - use it
        elif ct_conf > 0.7 and any(epi in ct_type.lower() for epi in ['basal', 'ciliated', 'secretory', 'club', 'goblet', 'epithelial']):
            final_type = ct_type
            final_conf = "High"
            agree = "CellTypist_Priority"

        # 4. Default to marker-based for organoids
        else:
            final_type = mb_type
            final_conf = "Low"
            agree = "Uncertain"

        consensus_types.append(final_type)
        consensus_confidence.append(final_conf)
        agreement.append(agree)

    # Summary
    logger.info("\nConsensus decision summary:")
    agreement_counts = pd.Series(agreement).value_counts()
    for decision, count in agreement_counts.items():
        pct = count / len(adata_ct) * 100
        logger.info(f"  {decision}: {count} cells ({pct:.1f}%)")

    logger.info("\nFinal consensus cell types:")
    consensus_counts = pd.Series(consensus_types).value_counts()
    for ct, count in consensus_counts.head(15).items():
        pct = count / len(adata_ct) * 100
        logger.info(f"  {ct}: {count} cells ({pct:.1f}%)")

    return consensus_types, consensus_confidence, agreement


def generate_comparison_plots(adata_ct, marker_types, consensus_types, output_dir):
    """
    Generate plots comparing all three annotation methods.
    """
    logger.info("\nGenerating comparison visualizations...")

    # Add all annotations to adata
    adata_ct.obs['marker_based'] = marker_types.values
    adata_ct.obs['consensus'] = consensus_types

    # 1. Three-method UMAP comparison
    fig, axes = plt.subplots(1, 3, figsize=(30, 8))

    sc.pl.umap(adata_ct, color='majority_voting', ax=axes[0], show=False,
               title='CellTypist (Cells_Lung_Airway)', legend_loc='right margin')
    sc.pl.umap(adata_ct, color='marker_based', ax=axes[1], show=False,
               title='Marker-Based', legend_loc='right margin')
    sc.pl.umap(adata_ct, color='consensus', ax=axes[2], show=False,
               title='Consensus', legend_loc='right margin')

    plt.tight_layout()
    plt.savefig(output_dir / 'three_method_comparison_umap.png', dpi=150, bbox_inches='tight')
    plt.close()

    # 2. Sankey/Alluvial-style comparison (bar plot version)
    fig, ax = plt.subplots(figsize=(16, 10))

    comparison_df = pd.DataFrame({
        'CellTypist': adata_ct.obs['majority_voting'].values,
        'Marker-Based': marker_types.values,
        'Consensus': consensus_types
    })

    # Count combinations
    combo_counts = comparison_df.groupby(['CellTypist', 'Marker-Based', 'Consensus']).size().reset_index(name='count')
    combo_counts = combo_counts.sort_values('count', ascending=False).head(20)

    ax.barh(range(len(combo_counts)), combo_counts['count'])
    labels = [f"CT:{row['CellTypist'][:15]}→MB:{row['Marker-Based'][:15]}→Con:{row['Consensus'][:15]}"
              for _, row in combo_counts.iterrows()]
    ax.set_yticks(range(len(combo_counts)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('Number of Cells', fontsize=12)
    ax.set_title('Top 20 Annotation Combinations\n(CellTypist → Marker-Based → Consensus)',
                 fontsize=14, fontweight='bold')
    ax.grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'annotation_method_flow.png', dpi=150, bbox_inches='tight')
    plt.close()

    # 3. Confusion matrix: CellTypist vs Consensus
    from sklearn.metrics import confusion_matrix
    import itertools

    ct_types_unique = adata_ct.obs['majority_voting'].unique()
    consensus_types_unique = pd.Series(consensus_types).unique()

    # Limit to top categories for readability
    top_ct = adata_ct.obs['majority_voting'].value_counts().head(15).index
    top_consensus = pd.Series(consensus_types).value_counts().head(15).index

    mask_ct = adata_ct.obs['majority_voting'].isin(top_ct)
    mask_consensus = pd.Series(consensus_types).isin(top_consensus)
    mask = mask_ct & mask_consensus

    if mask.sum() > 0:
        cm = confusion_matrix(
            adata_ct.obs.loc[mask, 'majority_voting'],
            pd.Series(consensus_types)[mask],
            labels=[t for t in top_ct if t in adata_ct.obs.loc[mask, 'majority_voting'].unique()]
        )

        fig, ax = plt.subplots(figsize=(16, 14))
        im = ax.imshow(cm, cmap='Blues', aspect='auto')

        labels_used = [t for t in top_ct if t in adata_ct.obs.loc[mask, 'majority_voting'].unique()]
        ax.set_xticks(range(len(labels_used)))
        ax.set_yticks(range(len(labels_used)))
        ax.set_xticklabels(labels_used, rotation=45, ha='right', fontsize=8)
        ax.set_yticklabels(labels_used, fontsize=8)
        ax.set_xlabel('Consensus Annotation', fontsize=12)
        ax.set_ylabel('CellTypist Annotation', fontsize=12)
        ax.set_title('Agreement: CellTypist vs Consensus (Top 15 Types)', fontsize=14, fontweight='bold')

        plt.colorbar(im, ax=ax, label='Number of Cells')
        plt.tight_layout()
        plt.savefig(output_dir / 'celltypist_consensus_confusion.png', dpi=150, bbox_inches='tight')
        plt.close()

    logger.info(f"Saved comparison plots to {output_dir}")


def main():
    """Main execution."""
    logger.info("=" * 70)
    logger.info("MULTI-METHOD CELL TYPE ANNOTATION WITH CONSENSUS")
    logger.info("=" * 70)

    # Load data
    logger.info(f"\nLoading data from {INPUT_H5AD}...")
    adata = sc.read_h5ad(INPUT_H5AD)
    logger.info(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Method 1: CellTypist with lung airway model
    adata_ct = run_celltypist_annotation(adata, model_name='Cells_Lung_Airway.pkl')

    # Method 2: Marker-based annotation
    marker_types, marker_confidence = run_marker_based_annotation(adata)

    # Consensus annotation (agent discussion)
    consensus_types, consensus_confidence, agreement = consensus_annotation(
        adata_ct, marker_types, marker_confidence
    )

    # Add all results to adata
    adata_ct.obs['marker_based'] = marker_types.values
    adata_ct.obs['marker_confidence'] = marker_confidence.values
    adata_ct.obs['consensus_celltype'] = consensus_types
    adata_ct.obs['consensus_confidence'] = consensus_confidence
    adata_ct.obs['agreement_type'] = agreement

    # Generate comparison plots
    generate_comparison_plots(adata_ct, marker_types, consensus_types, OUTPUT_DIR)

    # Save results
    output_h5ad = OUTPUT_DIR / "integrated_multimethod_annotated.h5ad"
    adata_ct.write_h5ad(output_h5ad)
    logger.info(f"\nSaved multi-method annotated data to {output_h5ad}")

    # Export comparison table
    comparison_df = pd.DataFrame({
        'cell_id': adata_ct.obs_names,
        'cluster': adata_ct.obs['leiden'],
        'celltypist': adata_ct.obs['majority_voting'],
        'celltypist_confidence': adata_ct.obs['conf_score'],
        'marker_based': marker_types.values,
        'marker_confidence': marker_confidence.values,
        'consensus': consensus_types,
        'consensus_confidence': consensus_confidence,
        'agreement_type': agreement,
        'sample': adata_ct.obs['sample'],
        'condition': adata_ct.obs['condition']
    })
    comparison_df.to_csv(OUTPUT_DIR / 'multi_method_comparison.csv', index=False)
    logger.info(f"Saved comparison table to {OUTPUT_DIR / 'multi_method_comparison.csv'}")

    # Summary statistics
    logger.info("\n" + "=" * 70)
    logger.info("MULTI-METHOD ANNOTATION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Output directory: {OUTPUT_DIR}")
    logger.info(f"\nKey files:")
    logger.info(f"  - integrated_multimethod_annotated.h5ad")
    logger.info(f"  - multi_method_comparison.csv")
    logger.info(f"  - three_method_comparison_umap.png")
    logger.info(f"  - annotation_method_flow.png")
    logger.info(f"  - celltypist_consensus_confusion.png")

    return adata_ct


if __name__ == "__main__":
    adata = main()
