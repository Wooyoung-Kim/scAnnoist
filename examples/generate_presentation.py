#!/usr/bin/env python
"""
Example: Generate PPT Presentation from scRNA-seq Analysis Results

This example shows how to automatically generate a PowerPoint presentation
from your tissue-aware CellTypist annotation results.
"""

from scrna_agent.tools.report_generator import generate_report, PresentationGenerator


def example_1_simple():
    """Example 1: Simple presentation generation"""
    print("=" * 60)
    print("Example 1: Simple Presentation Generation")
    print("=" * 60)

    # Assuming you have run annotation and have results in ./celltypist_outputs
    output_dir = "./celltypist_outputs"

    # Generate presentation
    pdf_path = generate_report(
        output_dir=output_dir,
        title="My scRNA-seq Analysis Results"
    )

    print(f"\n✓ Presentation generated: {pdf_path}")


def example_2_custom():
    """Example 2: Custom presentation with more control"""
    print("\n" + "=" * 60)
    print("Example 2: Custom Presentation")
    print("=" * 60)

    output_dir = "./celltypist_outputs"

    # Create generator with custom title
    generator = PresentationGenerator(output_dir)

    # Load results
    results = generator.load_analysis_results()
    print(f"Loaded results with {len(results)} sections")

    # Create custom slide plan
    slides = generator.create_slide_plan(
        results,
        title="Tissue-Aware Cell Type Annotation Study"
    )
    print(f"Created plan with {len(slides)} slides")

    # Generate presentation
    pdf_path = generator.generate_presentation(
        title="Tissue-Aware Cell Type Annotation Study"
    )

    print(f"\n✓ Presentation generated: {pdf_path}")


def example_3_complete_workflow():
    """Example 3: Complete workflow from annotation to presentation"""
    print("\n" + "=" * 60)
    print("Example 3: Complete Workflow")
    print("=" * 60)

    from scrna_agent.tools.annotation_tools import (
        load_adata_for_annotation,
        add_tissue_metadata,
        celltypist_annotate_tissue_aware
    )

    # Step 1: Run annotation
    print("\nStep 1: Running annotation...")
    data_file = "your_data.h5ad"
    output_dir = "./complete_workflow_output"

    load_adata_for_annotation(data_file)

    add_tissue_metadata(
        tissue_map="tissue_map.csv",
        sample_col="sample"
    )

    result = celltypist_annotate_tissue_aware(
        sample_col="sample",
        majority_voting=True,
        cluster_key="leiden",
        output_dir=output_dir
    )

    print(result)

    # Step 2: Generate presentation
    print("\nStep 2: Generating presentation...")
    pdf_path = generate_report(
        output_dir=output_dir,
        title="Complete Analysis Results"
    )

    print(f"\n✓ Complete workflow finished!")
    print(f"  Results: {output_dir}")
    print(f"  Presentation: {pdf_path}")


if __name__ == "__main__":
    print("""
    PPT Presentation Generation Examples
    ====================================

    This script contains 3 examples:

    1. Simple presentation generation
    2. Custom presentation with more control
    3. Complete workflow from annotation to presentation

    Before running:
    - Make sure you have run annotation and have results
    - Update file paths in the examples
    - Ensure Nano Banana Pro is accessible

    To run an example, uncomment the corresponding function call below.
    """)

    # Uncomment the example you want to run:

    # example_1_simple()
    # example_2_custom()
    # example_3_complete_workflow()

    print("\nDone! Check the output directory for your presentation.")
