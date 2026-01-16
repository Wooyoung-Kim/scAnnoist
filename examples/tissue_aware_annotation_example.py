#!/usr/bin/env python
"""
Example: Tissue-Aware CellTypist Annotation

This example demonstrates how to use the tissue-aware CellTypist annotation
feature to automatically select and apply the appropriate model for each tissue
in your dataset.
"""

import scanpy as sc
from scrna_agent.tools.annotation_tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware
)


def example_single_tissue():
    """Example 1: Single tissue annotation"""
    print("=" * 60)
    print("Example 1: Single Tissue Annotation")
    print("=" * 60)

    # Load your data
    adata = sc.read_h5ad("your_data.h5ad")

    # Load into annotation system
    load_adata_for_annotation("your_data.h5ad")

    # Add tissue metadata (apply to all cells)
    result = add_tissue_metadata(
        tissue="Spleen",
        sample_col="sample"
    )
    print(result)

    # Run tissue-aware annotation
    result = celltypist_annotate_tissue_aware(
        sample_col="sample",
        majority_voting=True,
        cluster_key="leiden",
        output_dir="./output_single_tissue"
    )
    print(result)

    # Access results
    print("\nAnnotation columns added:")
    print("- celltypist_label: Predicted cell type")
    print("- celltypist_confidence: Confidence score")
    print("- celltypist_model_used: Model used for prediction")
    print("- celltypist_majority_vote: Majority vote label (if enabled)")


def example_multi_tissue():
    """Example 2: Multi-tissue annotation with mapping file"""
    print("\n" + "=" * 60)
    print("Example 2: Multi-Tissue Annotation")
    print("=" * 60)

    # Create a tissue mapping file
    import pandas as pd
    tissue_map = pd.DataFrame({
        "sample": ["sample1", "sample2", "sample3"],
        "tissue": ["Spleen", "Blood", "Liver"]
    })
    tissue_map.to_csv("tissue_map.csv", index=False)
    print("Created tissue_map.csv:")
    print(tissue_map)

    # Load your data
    load_adata_for_annotation("your_data.h5ad")

    # Add tissue metadata from mapping
    result = add_tissue_metadata(
        tissue_map="tissue_map.csv",
        sample_col="sample"
    )
    print("\n" + result)

    # Run tissue-aware annotation
    result = celltypist_annotate_tissue_aware(
        sample_col="sample",
        majority_voting=True,
        cluster_key="leiden",
        output_dir="./output_multi_tissue"
    )
    print(result)


def example_with_sample_type():
    """Example 3: Annotation with tissue and sample type"""
    print("\n" + "=" * 60)
    print("Example 3: Annotation with Tissue and Sample Type")
    print("=" * 60)

    # Create mapping files
    import pandas as pd

    tissue_map = pd.DataFrame({
        "sample": ["sample1", "sample2", "sample3"],
        "tissue": ["Blood", "Blood", "Lung"]
    })
    tissue_map.to_csv("tissue_map.csv", index=False)

    sample_type_map = pd.DataFrame({
        "sample": ["sample1", "sample2", "sample3"],
        "sample_type": ["PBMC", "tumor", "normal"]
    })
    sample_type_map.to_csv("sample_type_map.csv", index=False)

    print("Created mapping files:")
    print("\ntissue_map.csv:")
    print(tissue_map)
    print("\nsample_type_map.csv:")
    print(sample_type_map)

    # Load data
    load_adata_for_annotation("your_data.h5ad")

    # Add metadata
    result = add_tissue_metadata(
        tissue_map="tissue_map.csv",
        sample_type_map="sample_type_map.csv",
        sample_col="sample"
    )
    print("\n" + result)

    # Run annotation
    # Note: Blood|PBMC will use override model, Blood|tumor will use tissue model
    result = celltypist_annotate_tissue_aware(
        sample_col="sample",
        majority_voting=True,
        cluster_key="leiden",
        output_dir="./output_with_sample_type"
    )
    print(result)


def example_inspect_model_selection():
    """Example 4: Inspect model selection before annotation"""
    print("\n" + "=" * 60)
    print("Example 4: Inspect Model Selection")
    print("=" * 60)

    from scrna_agent.tools.model_management import (
        select_model_for_tissue,
        get_models_for_adata,
        load_model_registry
    )

    # Load registry
    registry = load_model_registry()

    # Test model selection for different scenarios
    scenarios = [
        ("Spleen", None),
        ("Blood", "PBMC"),
        ("Blood", "tumor"),
        ("Lung", None),
        ("Unknown", "PBMC"),
    ]

    print("\nModel Selection Examples:")
    print("-" * 60)
    for tissue, sample_type in scenarios:
        model, rule = select_model_for_tissue(tissue, sample_type, registry)
        print(f"Tissue: {tissue:15} Sample Type: {str(sample_type):10}")
        print(f"  → Model: {model}")
        print(f"  → Rule: {rule}")
        print()

    # For existing dataset, get model assignments
    load_adata_for_annotation("your_data.h5ad")
    add_tissue_metadata(tissue_map="tissue_map.csv", sample_col="sample")

    # Get models for all tissues in dataset
    adata = sc.read_h5ad("your_data.h5ad")
    tissue_info = get_models_for_adata(adata, sample_col="sample")

    print("\nModels for your dataset:")
    print("-" * 60)
    for tissue, info in tissue_info.items():
        print(f"Tissue: {tissue}")
        print(f"  Model: {info['model']}")
        print(f"  Selection rule: {info['rule']}")
        print(f"  Cells: {info['n_cells']:,}")
        print(f"  Samples: {', '.join(info['samples'])}")
        print()


def example_download_models():
    """Example 5: Download required models"""
    print("\n" + "=" * 60)
    print("Example 5: Download Models")
    print("=" * 60)

    from scrna_agent.tools.model_management import pull_celltypist_models

    # Download models for specific tissue
    print("Downloading models for Spleen...")
    pull_celltypist_models(tissue="Spleen", pull_all=False)

    # Or download all models in registry
    print("\nDownloading all models...")
    pull_celltypist_models(tissue=None, pull_all=True)


if __name__ == "__main__":
    print("""
    Tissue-Aware CellTypist Annotation Examples
    ============================================

    This script contains 5 examples:

    1. Single tissue annotation
    2. Multi-tissue annotation with mapping file
    3. Annotation with tissue and sample type
    4. Inspect model selection before annotation
    5. Download required models

    To run an example, uncomment the corresponding function call below.
    Make sure to replace 'your_data.h5ad' with your actual data file.
    """)

    # Uncomment the example you want to run:

    # example_single_tissue()
    # example_multi_tissue()
    # example_with_sample_type()
    # example_inspect_model_selection()
    # example_download_models()

    print("\nDone! Check the output directories for results.")
