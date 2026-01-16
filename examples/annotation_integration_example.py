#!/usr/bin/env python
"""
Example: Annotation Integration Workflow
==========================================

This example demonstrates how to:
1. Run tissue-aware CellTypist annotation
2. Run literature-based (PubMed + CellMarker) annotation
3. Compare the two annotation methods
4. Create optimal annotation by intelligently combining both
5. Validate the integrated results

Author: scRNA-agent
"""

import scanpy as sc
from scrna_agent.tools.annotation_tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware,
    annotate_with_literature_markers,
    compare_celltypist_with_literature,
    create_optimal_annotation,
    integrate_multiple_annotations,
)


def example_1_simple_integration():
    """
    Example 1: Simple integration workflow

    Run CellTypist + Literature annotation, then create optimal annotation
    """
    print("=" * 70)
    print("Example 1: Simple Annotation Integration")
    print("=" * 70)

    # Step 1: Load data
    print("\n[1/5] Loading data...")
    load_adata_for_annotation("your_data.h5ad")

    # Step 2: Add tissue metadata
    print("\n[2/5] Adding tissue metadata...")
    add_tissue_metadata(
        tissue_map="tissue_map.csv",
        sample_col="sample"
    )

    # Step 3: Run CellTypist annotation
    print("\n[3/5] Running CellTypist tissue-aware annotation...")
    result = celltypist_annotate_tissue_aware(
        sample_col="sample",
        majority_voting=True,
        cluster_key="leiden",
        output_dir="./integration_outputs"
    )
    print(result)

    # Step 4: Run literature-based annotation
    print("\n[4/5] Running literature-based annotation...")
    lit_result = annotate_with_literature_markers(
        tissue_type="Blood",  # or get from metadata
        cluster_key="leiden",
        top_n_genes=50,
        min_confidence=0.3,
        output_col="literature_annotation"
    )
    print(lit_result)

    # Step 5: Create optimal annotation
    print("\n[5/5] Creating optimal annotation...")
    optimal_result = create_optimal_annotation(
        celltypist_col="celltypist_label",
        celltypist_conf_col="celltypist_confidence",
        literature_col="literature_annotation",
        confidence_threshold=0.5,
        output_col="optimal_annotation"
    )
    print(optimal_result)

    print("\n✓ Integration complete!")
    print("  - CellTypist annotation: adata.obs['celltypist_label']")
    print("  - Literature annotation: adata.obs['literature_annotation']")
    print("  - Optimal annotation: adata.obs['optimal_annotation']")


def example_2_detailed_comparison():
    """
    Example 2: Detailed comparison and analysis

    Compare methods in detail before creating optimal annotation
    """
    print("\n" + "=" * 70)
    print("Example 2: Detailed Comparison and Integration")
    print("=" * 70)

    # Assume annotations already exist from Example 1

    # Step 1: Compare CellTypist vs Literature
    print("\n[1/3] Comparing CellTypist vs Literature annotations...")
    comparison = compare_celltypist_with_literature(
        celltypist_col="celltypist_label",
        literature_col="literature_annotation",
        cluster_key="leiden"
    )
    print(comparison)

    # Step 2: Create optimal annotation with custom threshold
    print("\n[2/3] Creating optimal annotation (high confidence threshold)...")
    optimal_result = create_optimal_annotation(
        celltypist_col="celltypist_label",
        celltypist_conf_col="celltypist_confidence",
        literature_col="literature_annotation",
        confidence_threshold=0.7,  # Higher threshold
        output_col="optimal_annotation_strict"
    )
    print(optimal_result)

    # Step 3: Create optimal annotation with lower threshold
    print("\n[3/3] Creating optimal annotation (lower confidence threshold)...")
    optimal_result_relaxed = create_optimal_annotation(
        celltypist_col="celltypist_label",
        celltypist_conf_col="celltypist_confidence",
        literature_col="literature_annotation",
        confidence_threshold=0.3,  # Lower threshold
        output_col="optimal_annotation_relaxed"
    )
    print(optimal_result_relaxed)

    print("\n✓ Detailed comparison complete!")
    print("  - Strict optimal: adata.obs['optimal_annotation_strict']")
    print("  - Relaxed optimal: adata.obs['optimal_annotation_relaxed']")


def example_3_three_way_integration():
    """
    Example 3: Three-way integration (CellTypist + Literature + ScType)

    If you also have ScType annotations, integrate all three methods
    """
    print("\n" + "=" * 70)
    print("Example 3: Three-Way Integration")
    print("=" * 70)

    # Assume we have run:
    # - CellTypist -> adata.obs['celltypist_label']
    # - Literature -> adata.obs['literature_annotation']
    # - ScType -> adata.obs['sctype_annotation']

    # Step 1: Compare all three methods
    print("\n[1/2] Integrating three annotation methods...")
    integrated_result = integrate_multiple_annotations(
        methods=["celltypist_label", "literature_annotation", "sctype_annotation"],
        strategy="voting",
        weights={
            "celltypist_label": 1.0,
            "literature_annotation": 0.8,
            "sctype_annotation": 0.7
        },
        confidence_cols={
            "celltypist_label": "celltypist_confidence"
        },
        output_col="three_way_integrated"
    )
    print(integrated_result)

    # Step 2: Use optimal annotation with all three
    print("\n[2/2] Creating optimal annotation from all three methods...")
    optimal_result = create_optimal_annotation(
        celltypist_col="celltypist_label",
        celltypist_conf_col="celltypist_confidence",
        literature_col="literature_annotation",
        sctype_col="sctype_annotation",
        confidence_threshold=0.5,
        output_col="optimal_three_way"
    )
    print(optimal_result)

    print("\n✓ Three-way integration complete!")
    print("  - Voting result: adata.obs['three_way_integrated']")
    print("  - Optimal result: adata.obs['optimal_three_way']")


def example_4_custom_integration():
    """
    Example 4: Custom integration strategy

    Use confidence-weighted integration with custom weights
    """
    print("\n" + "=" * 70)
    print("Example 4: Custom Integration Strategy")
    print("=" * 70)

    # Custom weighted integration
    print("\n[1/1] Custom weighted integration...")
    integrated_result = integrate_multiple_annotations(
        methods=["celltypist_label", "literature_annotation"],
        strategy="confidence",  # Confidence-weighted
        weights={
            "celltypist_label": 1.2,  # Give CellTypist more weight
            "literature_annotation": 0.8
        },
        confidence_cols={
            "celltypist_label": "celltypist_confidence"
        },
        output_col="custom_integrated"
    )
    print(integrated_result)

    print("\n✓ Custom integration complete!")
    print("  - Result: adata.obs['custom_integrated']")


def example_5_complete_workflow_with_analysis():
    """
    Example 5: Complete workflow with visualization and analysis

    Full pipeline from annotation to analysis
    """
    print("\n" + "=" * 70)
    print("Example 5: Complete Workflow with Analysis")
    print("=" * 70)

    from scrna_agent.tools.annotation_tools import (
        plot_annotation_umap,
        save_annotations
    )

    # Step 1-4: Run annotations (same as Example 1)
    print("\n[1/7] Loading and annotating...")
    load_adata_for_annotation("your_data.h5ad")
    add_tissue_metadata(tissue_map="tissue_map.csv", sample_col="sample")
    celltypist_annotate_tissue_aware(
        sample_col="sample",
        majority_voting=True,
        cluster_key="leiden",
        output_dir="./complete_outputs"
    )
    annotate_with_literature_markers(
        tissue_type="Blood",
        cluster_key="leiden",
        output_col="literature_annotation"
    )

    # Step 5: Compare
    print("\n[2/7] Comparing annotations...")
    comparison = compare_celltypist_with_literature(
        celltypist_col="celltypist_label",
        literature_col="literature_annotation",
        cluster_key="leiden"
    )
    print(comparison)

    # Step 6: Create optimal
    print("\n[3/7] Creating optimal annotation...")
    optimal_result = create_optimal_annotation(
        celltypist_col="celltypist_label",
        celltypist_conf_col="celltypist_confidence",
        literature_col="literature_annotation",
        confidence_threshold=0.5,
        output_col="optimal_annotation"
    )
    print(optimal_result)

    # Step 7: Visualize all annotations
    print("\n[4/7] Visualizing CellTypist annotation...")
    plot_annotation_umap(
        annotation_col="celltypist_label",
        output_file="./complete_outputs/umap_celltypist.png"
    )

    print("\n[5/7] Visualizing literature annotation...")
    plot_annotation_umap(
        annotation_col="literature_annotation",
        output_file="./complete_outputs/umap_literature.png"
    )

    print("\n[6/7] Visualizing optimal annotation...")
    plot_annotation_umap(
        annotation_col="optimal_annotation",
        output_file="./complete_outputs/umap_optimal.png"
    )

    # Step 8: Save all annotations
    print("\n[7/7] Saving annotations...")
    save_annotations(
        annotation_cols=[
            "celltypist_label",
            "literature_annotation",
            "optimal_annotation"
        ],
        output_file="./complete_outputs/all_annotations.csv"
    )

    print("\n✓ Complete workflow finished!")
    print("  Outputs saved to: ./complete_outputs/")
    print("  - all_annotations.csv: All annotation results")
    print("  - umap_*.png: UMAP visualizations")


if __name__ == "__main__":
    print("""
    Annotation Integration Examples
    ================================

    This script contains 5 examples demonstrating annotation integration:

    1. Simple integration: CellTypist + Literature → Optimal
    2. Detailed comparison: Compare methods before integration
    3. Three-way integration: CellTypist + Literature + ScType
    4. Custom integration: Custom weighted integration
    5. Complete workflow: Full pipeline with visualization

    Prerequisites:
    - Your data file (your_data.h5ad)
    - Tissue mapping file (tissue_map.csv)
    - Basic preprocessing done (normalization, clustering)

    To run an example, uncomment the corresponding function call below.
    """)

    # Uncomment the example you want to run:

    # example_1_simple_integration()
    # example_2_detailed_comparison()
    # example_3_three_way_integration()
    # example_4_custom_integration()
    # example_5_complete_workflow_with_analysis()

    print("\n" + "=" * 70)
    print("Ready to run! Uncomment an example function to start.")
    print("=" * 70)
