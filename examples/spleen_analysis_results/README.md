# Spleen scRNA-seq Analysis Results

This directory contains example analysis results from a spleen single-cell RNA-seq dataset processed with scAnnoist.

## Contents

### Figures

**Quality Control & Processing:**
- `01_QC_metrics.png` - Quality control metrics
- `02_PCA_variance.png` - PCA variance explained
- `03_UMAP_before_integration.png` - UMAP visualization before batch correction
- `04_UMAP_after_integration.png` - UMAP visualization after Harmony integration
- `07_integration_before_after_comparison.png` - Before/after integration comparison

**Clustering & Resolution:**
- `05_clustering_resolutions.png` - Clustering at different resolutions
- `12_UMAP_louvain_clusters_labeled.png` - Louvain clustering results with labels
- `15_UMAP_cluster_numbers.png` - UMAP with cluster numbers

**Cell Type Annotation:**
- `07_annotation_comparison.png` - Comparison of annotation methods
- `08_major_cell_types.png` - Major cell type assignments
- `09_refined_cell_types.png` - Refined cell type annotations
- `10_UMAP_major_celltypes_labeled.png` - UMAP with major cell types
- `11_UMAP_refined_celltypes_labeled.png` - UMAP with refined cell types
- `13_UMAP_top10_celltypes_labeled.png` - Top 10 cell types visualization
- `14_UMAP_before_after_labeled.png` - Labeled before/after comparison

**Marker Genes:**
- `06_marker_genes_heatmap.png` - Heatmap of marker genes
- `09_marker_dotplot.png` - Dot plot of marker gene expression

### Data Files

- `cell_type_annotations.csv` - Cell type assignments for all cells
- `cluster_assignments.csv` - Cluster assignments
- `umap_coordinates.csv` - UMAP coordinates
- `umap_coordinates_annotated.csv` - UMAP coordinates with annotations
- `marker_genes_per_cluster.csv` - Top marker genes for each cluster
- `annotation_summary.txt` - Summary of annotation results

### Reports

- `ANALYSIS_REPORT.md` - Comprehensive analysis report
- `Spleen_Analysis_Presentation.pptx` - Automated presentation

### Workflow Diagrams

- `Hierarchical_Multi-Agent_System.pptx` - Multi-agent system architecture
- `Multi-Agent_Collaboration.pptx` - Agent collaboration workflow
- `scrna_agent_Workflow.pptx` - Overall analysis workflow

## Dataset Information

- **Tissue:** Spleen
- **Species:** (Check ANALYSIS_REPORT.md for details)
- **Processing:** Harmony integration for batch correction
- **Annotation:** Multi-method approach (CellTypist + marker-based)

## How to Use

These results demonstrate the complete scAnnoist pipeline output. You can use them as:

1. Reference for expected output format
2. Example for quality assessment
3. Template for your own analysis reports

For running your own analysis, see the parent `examples/` directory for code examples.
