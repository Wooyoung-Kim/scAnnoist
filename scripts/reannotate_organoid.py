#!/usr/bin/env python
"""
Re-annotate JMV_ALI dataset with biological context validation.
Organoid samples should NOT contain immune cells (neutrophils, monocytes, etc.)
"""

import sys
from pathlib import Path
import logging

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

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
OUTPUT_DIR = Path("/mnt2/kwy/scrna_agent/results/jmv_ali_epithelial/reannotation")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

# Airway epithelial markers (for organoids)
EPITHELIAL_MARKERS = {
    'Basal': ['KRT5', 'TP63', 'KRT15', 'KRT14', 'NGFR', 'KRT17'],
    'Ciliated': ['FOXJ1', 'DNAH5', 'TPPP3', 'CAPS', 'SNTN', 'PIFO'],
    'Secretory/Club': ['SCGB1A1', 'SCGB3A2', 'MUC5B', 'CYP2F1', 'LYPD2'],
    'Goblet': ['MUC5AC', 'MUC5B', 'TFF3', 'SPDEF', 'FCGBP'],
    'AT2-like': ['SFTPC', 'SFTPA1', 'SFTPB', 'LAMP3', 'ABCA3'],
    'Ionocyte': ['FOXI1', 'CFTR', 'ATP6V1C2', 'ATP6V0D2'],
    'General_Epithelial': ['EPCAM', 'KRT18', 'KRT8', 'KRT19', 'CDH1']
}

# Immune markers (should be ABSENT in organoids)
IMMUNE_MARKERS = {
    'Neutrophil': ['S100A8', 'S100A9', 'S100A12', 'CSF3R', 'CXCR2', 'FCGR3B'],
    'Monocyte/Macrophage': ['CD14', 'CD68', 'LYZ', 'FCGR1A', 'CD163'],
    'T_cell': ['CD3D', 'CD3E', 'CD3G', 'CD8A', 'CD4'],
    'B_cell': ['CD79A', 'CD79B', 'MS4A1', 'CD19'],
    'NK_cell': ['GNLY', 'NKG7', 'NCAM1'],
    'Pan_immune': ['PTPRC']  # CD45
}


def validate_organoid_annotation(adata):
    """
    Validate that organoid samples don't have inappropriate cell types.
    """
    logger.info("=" * 60)
    logger.info("BIOLOGICAL CONTEXT VALIDATION")
    logger.info("=" * 60)

    logger.info("Sample type: ALI Organoid")
    logger.info("Expected: Epithelial cells only")
    logger.info("Not expected: Immune cells (neutrophils, monocytes, etc.)")

    # Check current annotations
    celltype_counts = adata.obs['celltype'].value_counts()
    logger.info(f"\nCurrent annotation distribution:")
    for ct, count in celltype_counts.items():
        pct = count / len(adata) * 100
        logger.info(f"  {ct}: {count} cells ({pct:.1f}%)")

    # Flag suspicious annotations
    suspicious_types = ['Neutrophil', 'Monocyte', 'Macrophage', 'T_cell', 'B_cell',
                        'NK_cell', 'Plasma_cell', 'Dendritic_cell']

    issues = []
    for ct in suspicious_types:
        if ct in celltype_counts.index:
            count = celltype_counts[ct]
            pct = count / len(adata) * 100
            if pct > 5:  # More than 5% is suspicious
                issues.append(f"  ⚠️  {ct}: {count} cells ({pct:.1f}%) - UNEXPECTED in organoids!")
            elif pct > 0:
                issues.append(f"  ⚠️  {ct}: {count} cells ({pct:.1f}%) - unusual but possible contamination")

    if issues:
        logger.warning("\n⚠️ BIOLOGICAL VALIDATION ISSUES DETECTED:")
        for issue in issues:
            logger.warning(issue)
        logger.warning("\nOrganoid samples should NOT contain immune cells!")
        logger.warning("These cells are likely MISANNOTATED epithelial cells.")
        return False
    else:
        logger.info("\n✓ Annotation appears biologically consistent")
        return True


def check_marker_expression(adata, markers_dict, cell_mask=None):
    """
    Check expression of marker genes for cells.
    """
    if cell_mask is None:
        cell_mask = np.ones(adata.n_obs, dtype=bool)

    results = {}
    for celltype, markers in markers_dict.items():
        present_markers = [m for m in markers if m in adata.var_names]
        if not present_markers:
            results[celltype] = {'mean_expr': 0, 'pct_cells': 0, 'markers_found': 0}
            continue

        # Get expression
        expr = adata[cell_mask, present_markers].X
        if hasattr(expr, 'toarray'):
            expr = expr.toarray()

        mean_expr = expr.mean(axis=1).mean()  # Mean expression across cells
        pct_cells = (expr.max(axis=1) > 0).mean() * 100  # % cells expressing any marker

        results[celltype] = {
            'mean_expr': mean_expr,
            'pct_cells': pct_cells,
            'markers_found': len(present_markers)
        }

    return results


def reannotate_organoid_cells(adata):
    """
    Re-annotate cells using organoid-specific logic.
    """
    logger.info("\n" + "=" * 60)
    logger.info("RE-ANNOTATION WITH ORGANOID CONTEXT")
    logger.info("=" * 60)

    # For each cluster, check epithelial vs immune markers
    new_celltype = []
    new_confidence = []

    for cluster_id in adata.obs['leiden'].unique():
        cluster_mask = adata.obs['leiden'] == cluster_id
        n_cells = cluster_mask.sum()

        logger.info(f"\nCluster {cluster_id} ({n_cells} cells):")

        # Check epithelial markers
        epi_results = check_marker_expression(adata, EPITHELIAL_MARKERS, cluster_mask)

        # Check immune markers
        immune_results = check_marker_expression(adata, IMMUNE_MARKERS, cluster_mask)

        # Find best epithelial match
        best_epi_type = max(epi_results.items(), key=lambda x: x[1]['mean_expr'])

        # Find max immune score
        best_immune_type = max(immune_results.items(), key=lambda x: x[1]['mean_expr'])

        logger.info(f"  Epithelial markers: {best_epi_type[0]} (expr={best_epi_type[1]['mean_expr']:.3f}, pct={best_epi_type[1]['pct_cells']:.1f}%)")
        logger.info(f"  Immune markers: {best_immune_type[0]} (expr={best_immune_type[1]['mean_expr']:.3f}, pct={best_immune_type[1]['pct_cells']:.1f}%)")

        # Decision logic
        epi_score = best_epi_type[1]['mean_expr']
        immune_score = best_immune_type[1]['mean_expr']

        if immune_score > 0.5 and immune_score > epi_score * 2:
            # Strong immune signature - likely contamination, but flag it
            assigned_type = f"Contamination_{best_immune_type[0]}"
            confidence = "Low"
            logger.warning(f"  → {assigned_type} (contamination/doublet)")
        elif epi_score > 0.3:
            # Clear epithelial signature
            assigned_type = best_epi_type[0]
            confidence = "High" if epi_score > 0.7 else "Medium"
            logger.info(f"  → {assigned_type} (confidence={confidence})")
        else:
            # Unclear
            assigned_type = "Unknown_Epithelial"
            confidence = "Low"
            logger.info(f"  → {assigned_type}")

        # Assign to all cells in cluster
        for _ in range(n_cells):
            new_celltype.append(assigned_type)
            new_confidence.append(confidence)

    # Update annotations
    adata.obs['celltype_reannotated'] = new_celltype
    adata.obs['annotation_confidence'] = new_confidence

    # Summary
    logger.info("\n" + "=" * 60)
    logger.info("RE-ANNOTATION SUMMARY")
    logger.info("=" * 60)
    new_counts = adata.obs['celltype_reannotated'].value_counts()
    for ct, count in new_counts.items():
        pct = count / len(adata) * 100
        logger.info(f"  {ct}: {count} cells ({pct:.1f}%)")

    return adata


def generate_validation_plots(adata, output_dir):
    """
    Generate plots showing marker expression validation.
    """
    logger.info("\nGenerating validation plots...")

    # 1. UMAP comparison: old vs new annotation
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))

    sc.pl.umap(adata, color='celltype', ax=axes[0], show=False, title='Original Annotation')
    sc.pl.umap(adata, color='celltype_reannotated', ax=axes[1], show=False, title='Re-annotation (Organoid Context)')

    plt.tight_layout()
    plt.savefig(output_dir / 'annotation_comparison_umap.png', dpi=150, bbox_inches='tight')
    plt.close()

    # 2. Marker expression heatmap per cluster
    # Collect all markers
    all_epi_markers = []
    for markers in EPITHELIAL_MARKERS.values():
        all_epi_markers.extend(markers)

    all_immune_markers = []
    for markers in IMMUNE_MARKERS.values():
        all_immune_markers.extend(markers)

    # Filter to present markers
    present_epi = [m for m in all_epi_markers if m in adata.var_names][:20]  # Top 20
    present_immune = [m for m in all_immune_markers if m in adata.var_names][:10]  # Top 10

    markers_to_plot = present_epi + present_immune

    if len(markers_to_plot) > 5:
        try:
            fig, ax = plt.subplots(figsize=(12, 10))
            sc.pl.matrixplot(adata, markers_to_plot, groupby='leiden',
                            cmap='RdYlBu_r', ax=ax, show=False, swap_axes=True)
            plt.savefig(output_dir / 'marker_expression_by_cluster.png', dpi=150, bbox_inches='tight')
            plt.close()
        except Exception as e:
            logger.warning(f"Could not generate matrixplot: {e}")

    # 3. Composition comparison
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Old annotation
    old_counts = adata.obs['celltype'].value_counts()
    axes[0].barh(range(len(old_counts)), old_counts.values, color='lightcoral')
    axes[0].set_yticks(range(len(old_counts)))
    axes[0].set_yticklabels(old_counts.index)
    axes[0].set_xlabel('Number of Cells')
    axes[0].set_title('Original Annotation', fontweight='bold')
    axes[0].grid(axis='x', alpha=0.3)

    # New annotation
    new_counts = adata.obs['celltype_reannotated'].value_counts()
    axes[1].barh(range(len(new_counts)), new_counts.values, color='lightblue')
    axes[1].set_yticks(range(len(new_counts)))
    axes[1].set_yticklabels(new_counts.index)
    axes[1].set_xlabel('Number of Cells')
    axes[1].set_title('Re-annotation (Organoid Context)', fontweight='bold')
    axes[1].grid(axis='x', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_dir / 'composition_comparison.png', dpi=150, bbox_inches='tight')
    plt.close()

    logger.info(f"Saved validation plots to {output_dir}")


def main():
    """Main execution."""
    logger.info("=" * 70)
    logger.info("JMV_ALI ORGANOID RE-ANNOTATION WITH BIOLOGICAL VALIDATION")
    logger.info("=" * 70)

    # Load data
    logger.info(f"\nLoading data from {INPUT_H5AD}...")
    adata = sc.read_h5ad(INPUT_H5AD)
    logger.info(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Validate current annotation
    is_valid = validate_organoid_annotation(adata)

    if not is_valid:
        logger.info("\n" + "=" * 60)
        logger.info("PROCEEDING WITH RE-ANNOTATION...")
        logger.info("=" * 60)

        # Re-annotate
        adata = reannotate_organoid_cells(adata)

        # Generate validation plots
        generate_validation_plots(adata, OUTPUT_DIR)

        # Save results
        output_h5ad = OUTPUT_DIR / "integrated_reannotated.h5ad"
        adata.write_h5ad(output_h5ad)
        logger.info(f"\nSaved re-annotated data to {output_h5ad}")

        # Export annotation comparison
        comparison_df = pd.DataFrame({
            'cell_id': adata.obs_names,
            'cluster': adata.obs['leiden'],
            'original_annotation': adata.obs['celltype'],
            'reannotation': adata.obs['celltype_reannotated'],
            'confidence': adata.obs['annotation_confidence'],
            'sample': adata.obs['sample'],
            'condition': adata.obs['condition']
        })
        comparison_df.to_csv(OUTPUT_DIR / 'annotation_comparison.csv', index=False)
        logger.info(f"Saved annotation comparison to {OUTPUT_DIR / 'annotation_comparison.csv'}")

    logger.info("\n" + "=" * 70)
    logger.info("RE-ANNOTATION COMPLETE")
    logger.info("=" * 70)
    logger.info(f"Output directory: {OUTPUT_DIR}")

    return adata


if __name__ == "__main__":
    adata = main()
