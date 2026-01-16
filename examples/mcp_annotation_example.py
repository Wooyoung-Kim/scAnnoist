#!/usr/bin/env python3
"""
MCP-Based Cell Type Annotation Example

This example demonstrates how to use MCP (Model Context Protocol) tools
for literature-based cell type annotation and validation.

Requirements:
- scRNA-seq data with preliminary annotations
- Claude Code environment with PubMed MCP
"""

import scanpy as sc
from pathlib import Path

# Import MCP annotation tools
from scrna_agent.tools.mcp_annotation_tools import (
    mcp_search_marker_literature,
    mcp_validate_cell_type_annotation,
    mcp_get_marker_evidence,
    mcp_find_related_markers,
    hybrid_search_markers,
)

# Import MCP discovery tools
from scrna_agent.skills.mcp_tools import (
    list_mcp_tools,
    search_mcp_tools,
    get_mcp_tool_recommendations,
)


def example_1_discover_mcp_tools():
    """Example 1: Discover available MCP tools."""
    print("=" * 70)
    print("Example 1: Discover MCP Tools")
    print("=" * 70)

    # List all MCP tools
    print("\n1. Listing all MCP tools:")
    tools = list_mcp_tools()
    print(tools)

    # Search for specific tools
    print("\n2. Searching for literature-related tools:")
    search_result = search_mcp_tools("literature marker validation")
    print(search_result)

    # Get recommendations for a task
    print("\n3. Getting tool recommendations:")
    recommendations = get_mcp_tool_recommendations(
        "I need to validate cell type annotations with published research"
    )
    print(recommendations)


def example_2_search_marker_literature():
    """Example 2: Search PubMed for marker literature."""
    print("\n" + "=" * 70)
    print("Example 2: Search Marker Literature")
    print("=" * 70)

    # Search for T cell markers
    print("\n1. Searching for CD8 T cell markers:")
    result = mcp_search_marker_literature(
        cell_type="CD8 T cells",
        marker_genes=["CD8A", "CD8B", "GZMB", "PRF1"],
        species="human",
        tissue="blood",
        max_results=5
    )
    print(result)

    # Search for B cell markers
    print("\n2. Searching for B cell markers:")
    result = mcp_search_marker_literature(
        cell_type="B cells",
        marker_genes=["CD19", "MS4A1", "CD79A"],
        species="human",
        max_results=5
    )
    print(result)


def example_3_validate_annotation(adata_path: str = None):
    """Example 3: Validate cell type annotations."""
    print("\n" + "=" * 70)
    print("Example 3: Validate Annotations")
    print("=" * 70)

    # Define hypothetical cluster annotations
    cluster_annotations = {
        "0": {
            "cell_type": "CD8 T cells",
            "observed_markers": ["CD8A", "CD8B", "CD3D", "GZMB", "PRF1"]
        },
        "1": {
            "cell_type": "CD4 T cells",
            "observed_markers": ["CD4", "CD3D", "IL7R", "LDHB"]
        },
        "2": {
            "cell_type": "B cells",
            "observed_markers": ["CD19", "MS4A1", "CD79A", "CD79B"]
        },
        "3": {
            "cell_type": "NK cells",
            "observed_markers": ["NKG7", "GNLY", "NCAM1", "KLRD1", "GZMB"]
        },
    }

    # Validate each annotation
    for cluster_id, annotation in cluster_annotations.items():
        print(f"\nValidating Cluster {cluster_id}: {annotation['cell_type']}")
        print("-" * 70)

        validation = mcp_validate_cell_type_annotation(
            cell_type=annotation["cell_type"],
            observed_markers=annotation["observed_markers"],
            species="human",
            confidence_threshold=0.7
        )
        print(validation)


def example_4_find_related_markers():
    """Example 4: Discover related markers from literature."""
    print("\n" + "=" * 70)
    print("Example 4: Find Related Markers")
    print("=" * 70)

    # Find related markers for B cells
    print("\n1. Finding related markers for B cells:")
    related = mcp_find_related_markers(
        known_markers=["CD19", "MS4A1"],
        cell_type="B cells",
        species="human"
    )
    print(related)

    # Find related markers for regulatory T cells
    print("\n2. Finding related markers for Tregs:")
    related = mcp_find_related_markers(
        known_markers=["FOXP3", "IL2RA"],
        cell_type="regulatory T cells",
        species="human"
    )
    print(related)


def example_5_get_evidence():
    """Example 5: Get detailed evidence from papers."""
    print("\n" + "=" * 70)
    print("Example 5: Get Marker Evidence")
    print("=" * 70)

    # Example PMIDs (replace with actual PMIDs from search results)
    example_pmids = ["35486828", "34577062", "33264437"]

    print(f"\nRetrieving evidence from {len(example_pmids)} papers:")
    evidence = mcp_get_marker_evidence(
        pmids=example_pmids,
        extract_markers=True
    )
    print(evidence)


def example_6_hybrid_search():
    """Example 6: Use hybrid search with fallback."""
    print("\n" + "=" * 70)
    print("Example 6: Hybrid Literature Search")
    print("=" * 70)

    # Hybrid search tries MCP first, then Bio.Entrez
    print("\n1. Hybrid search for macrophage markers:")
    result = hybrid_search_markers(
        cell_type="macrophages",
        species="mouse",
        tissue="lung",
        use_mcp=True
    )
    print(result)


def example_7_complete_annotation_workflow(adata_path: str = None):
    """Example 7: Complete annotation workflow with MCP validation."""
    print("\n" + "=" * 70)
    print("Example 7: Complete Annotation Workflow")
    print("=" * 70)

    if adata_path is None:
        print("\nSkipping workflow example (no data provided)")
        print("To run this example, provide path to h5ad file:")
        print("  python mcp_annotation_example.py --data your_data.h5ad")
        return

    # Load data
    print(f"\n1. Loading data from {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"   Cells: {adata.n_obs:,}, Genes: {adata.n_vars:,}")

    # Get top marker genes for each cluster
    print("\n2. Computing marker genes...")
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # For each cluster, validate annotation
    print("\n3. Validating annotations with literature...")

    for cluster in adata.obs["leiden"].unique():
        print(f"\n   Cluster {cluster}:")

        # Get top markers
        markers_df = sc.get.rank_genes_groups_df(adata, group=cluster)
        top_markers = markers_df.head(5)["names"].tolist()

        print(f"   Top markers: {', '.join(top_markers)}")

        # Hypothetical annotation (replace with actual annotation logic)
        cell_type = "Unknown"  # Would come from automated annotation

        # Validate with MCP
        validation = mcp_validate_cell_type_annotation(
            cell_type=cell_type,
            observed_markers=top_markers,
            species="human",
            confidence_threshold=0.7
        )

        print(f"   Validation result:")
        print(f"   {validation[:200]}...")


def main():
    """Run all examples."""
    import argparse

    parser = argparse.ArgumentParser(
        description="MCP-based cell type annotation examples"
    )
    parser.add_argument(
        "--data",
        type=str,
        help="Path to h5ad file for complete workflow example"
    )
    parser.add_argument(
        "--example",
        type=int,
        choices=range(1, 8),
        help="Run specific example (1-7)"
    )

    args = parser.parse_args()

    print("=" * 70)
    print("MCP-Based Cell Type Annotation Examples")
    print("=" * 70)

    examples = {
        1: example_1_discover_mcp_tools,
        2: example_2_search_marker_literature,
        3: example_3_validate_annotation,
        4: example_4_find_related_markers,
        5: example_5_get_evidence,
        6: example_6_hybrid_search,
        7: lambda: example_7_complete_annotation_workflow(args.data),
    }

    if args.example:
        # Run specific example
        print(f"\nRunning Example {args.example}...\n")
        examples[args.example]()
    else:
        # Run all examples except workflow (needs data)
        for i in range(1, 7):
            try:
                examples[i]()
            except Exception as e:
                print(f"\n❌ Example {i} failed: {e}")
                print("   (This may be expected if MCP tools are not available)")

        # Optionally run workflow if data provided
        if args.data:
            try:
                examples[7]()
            except Exception as e:
                print(f"\n❌ Workflow example failed: {e}")

    print("\n" + "=" * 70)
    print("Examples Complete")
    print("=" * 70)
    print("\nFor more information, see:")
    print("  docs/MCP_INTEGRATION_GUIDE.md")


if __name__ == "__main__":
    main()
