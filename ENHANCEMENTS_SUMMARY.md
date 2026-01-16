# Enhanced JMV_ALI scRNA-seq Analysis Pipeline

## Overview
This document summarizes the enhancements made to the `scripts/run_jmv_ali.py` analysis pipeline for comprehensive analysis of the JMV_ALI dataset (24 ALI organoid samples across viral infection conditions).

## Dataset Information
- **Total Samples**: 24 (3 replicates × 8 conditions)
- **Tissue Type**: Airway epithelium (ALI organoid cultures)
- **Conditions**:
  - Control (3 samples)
  - HCoV_OC43 (3 samples)
  - MERS (3 samples)
  - SARS-CoV-1 (3 samples)
  - SARS-CoV-2_BA5 (3 samples)
  - SARS-CoV-2_BA275 (3 samples)
  - SARS-CoV-2_Delta (3 samples)
  - SARS-CoV-2_L (original strain, 3 samples)

## Key Enhancements

### 1. Enhanced Quality Control Visualization
**New Features:**
- Pre-QC and post-QC comprehensive visualizations
- Violin plots for QC metrics (genes, counts, MT%) per sample
- Scatter plots showing relationships between QC metrics
- Histograms with median/mean statistics
- Cell count distributions by sample and condition
- Per-sample QC statistics table exported to CSV

**Output Files:**
- `qc/qc_metrics_before.png` - Pre-filtering QC visualizations
- `qc/qc_metrics_after.png` - Post-filtering QC visualizations
- `qc/qc_statistics.csv` - Per-sample QC statistics table

### 2. Batch Integration Quality Assessment
**New Features:**
- Side-by-side comparison of before/after Harmony correction
- PCA plots colored by sample and condition (before and after)
- UMAP plots showing batch mixing (before and after)
- Visual assessment of integration quality

**Output Files:**
- `integration/integration_comparison.png` - 2×3 panel comparison

### 3. Enhanced Cell Type Annotation
**New Features:**
- Comprehensive marker gene visualizations
- Dotplot showing marker expression across clusters
- Heatmap of marker genes for airway cell types
- Cell type proportion analysis by condition

**Marker Genes Included:**
- **Ciliated**: FOXJ1, DNAH5, TPPP3, CAPS, SNTN
- **Basal**: KRT5, TP63, KRT15, KRT14, NGFR
- **Secretory/Club**: SCGB1A1, SCGB3A2, MUC5B, CYP2F1
- **Goblet**: MUC5AC, MUC5B, TFF3, SPDEF
- **AT2-like**: SFTPC, SFTPA1, SFTPB, LAMP3
- **Immune**: PTPRC, CD3D, CD79A, LYZ, CD68

**Output Files:**
- `annotation/marker_dotplot.png` - Dotplot of markers × clusters
- `annotation/celltype_markers_heatmap.png` - Heatmap of marker expression
- `annotation/celltype_proportions.csv` - Cell type percentages by condition
- `annotation/celltype_composition.png` - Stacked bar plot
- `annotation/cell_metadata.csv` - Cell-level annotations

### 4. Condition-Level Differential Expression Analysis
**New Features:**
- Pairwise comparisons: Each viral condition vs Control
- Full statistical results for each comparison
- Separate analysis for each of the 7 viral conditions:
  - HCoV_OC43 vs Control
  - MERS vs Control
  - SARS-CoV-1 vs Control
  - SARS-CoV-2_BA5 vs Control
  - SARS-CoV-2_BA275 vs Control
  - SARS-CoV-2_Delta vs Control
  - SARS-CoV-2_L vs Control

**Statistics Provided:**
- Log2 fold changes
- Raw and adjusted p-values (Benjamini-Hochberg FDR)
- Percentage of cells expressing gene in each group
- Wilcoxon rank-sum test scores

**Output Files:**
- `deg/condition_comparisons/{condition}_vs_Control_degs.csv` - Full DEG table per comparison
- `deg/condition_comparisons/volcano_{condition}_vs_Control.png` - Volcano plots

### 5. Enhanced Visualization Suite
**New Features:**

#### Extended UMAP Plots:
- UMAP colored by QC metrics (n_genes_by_counts)
- UMAP colored by key marker genes (FOXJ1, KRT5, SCGB1A1, MUC5AC)
- UMAP colored by viral response genes (ISG15, IFIT1, IFIT2, IFIT3)

#### Viral Response Genes Tracked:
ISG15, IFIT1, IFIT2, IFIT3, MX1, MX2, OAS1, IFI6, IFI44L, IFITM3

**Output Files:**
- `umap/umap_clusters.png` - UMAP by cluster
- `umap/umap_condition.png` - UMAP by condition
- `umap/umap_group.png` - UMAP by virus group
- `umap/umap_sample.png` - UMAP by sample
- `umap/umap_celltype.png` - UMAP by cell type annotation
- `umap/umap_n_genes.png` - UMAP colored by gene count
- `umap/umap_key_markers.png` - UMAP panel for key marker genes
- `umap/umap_viral_response_genes.png` - UMAP panel for viral response genes

### 6. Improved Output Organization
**New Directory Structure:**
```
results/jmv_ali_statistical/
├── integrated.h5ad                    # Main processed data
├── report.md                          # Analysis report
├── run_metadata.json                  # Analysis metadata
├── qc/                                # Quality control plots and stats
├── integration/                       # Batch integration assessment
├── umap/                              # All UMAP visualizations
├── annotation/                        # Cell type annotation results
├── clustering/                        # Clustering results and composition
├── deg/                               # Differential expression results
│   ├── markers.csv                    # Top cluster markers
│   ├── markers_full.csv               # All cluster markers
│   └── condition_comparisons/         # Virus vs Control DEGs
└── supplementary/                     # Supporting files (HVG list, etc.)
```

## Methodological Improvements

### Statistical Rigor
All existing statistical methods are preserved:
- **QC**: MAD-based adaptive thresholding (3 MADs from median)
- **Normalization**: Library size normalization (CPM) + log1p
- **HVG Selection**: Cell Ranger flavor (3000 genes)
- **Batch Correction**: Harmony integration (key='sample')
- **PC Selection**: Elbow method on scree plot
- **Clustering**: Leiden algorithm (resolution=0.8)
- **DEG Testing**: Wilcoxon rank-sum test with Benjamini-Hochberg FDR correction

### Reproducibility
- All analyses use fixed random seed = 42
- Complete metadata tracking
- Full parameter documentation

## Analysis Workflow

The enhanced pipeline follows these steps:

1. **Data Loading** - Load all 24 samples from 10X CellRanger output
2. **QC Metrics** - Calculate and visualize QC metrics (before filtering)
3. **Adaptive QC** - Apply MAD-based thresholds and visualize (after filtering)
4. **Normalization** - Library size normalization + log transformation
5. **HVG Selection** - Identify 3000 highly variable genes
6. **PCA** - Compute 50 principal components
7. **Pre-Integration UMAP** - Compute UMAP before Harmony correction
8. **Harmony Correction** - Batch integration across 24 samples
9. **Post-Integration UMAP** - Compute UMAP after Harmony correction
10. **Integration Visualization** - Compare before/after integration quality
11. **Clustering** - Leiden clustering (resolution=0.8)
12. **Cluster Markers** - Find marker genes per cluster
13. **Cell Type Annotation** - Marker-based annotation for airway cell types
14. **Marker Visualization** - Generate heatmaps and dotplots
15. **Condition DEG Analysis** - Compare each virus vs control
16. **DEG Visualization** - Generate volcano plots for each comparison
17. **Extended UMAP Suite** - Generate marker gene and viral response UMAPs
18. **Report Generation** - Create comprehensive markdown report

## Running the Enhanced Analysis

### Prerequisites
```bash
# Activate the scRNAseq conda environment
conda activate scRNAseq
```

### Execution
```bash
# Run the enhanced pipeline
cd /mnt2/kwy/scrna_agent
python scripts/run_jmv_ali.py
```

### Expected Runtime
- Full analysis: 15-30 minutes (depending on system)
- Dataset size: ~84,000 cells after QC, 44,700 genes

### Output Location
```bash
/mnt2/kwy/scrna_agent/results/jmv_ali_statistical/
```

## Key Findings to Expect

### Quality Control
- Total cells before QC: ~128,000
- Cells after QC: ~84,000 (65.2% retention rate)
- Median genes per cell: varies by sample
- Median UMI counts: varies by sample
- MT% threshold: 20%

### Cell Type Composition
- Expected cell types in ALI organoids:
  - Basal cells (progenitor population)
  - Ciliated cells (differentiated)
  - Secretory/Club cells
  - Goblet cells (mucus-producing)
  - Possible immune infiltrates

### Viral Response Signatures
- Interferon-stimulated genes (ISGs) expression
- Differential expression patterns across virus types
- Variant-specific transcriptional signatures
- Cell-type-specific viral responses

## Future Enhancement Opportunities

### Implemented in This Version:
- ✓ Enhanced QC visualization
- ✓ Integration quality assessment
- ✓ Marker-based annotation with visualization
- ✓ Condition-level DEG analysis
- ✓ Extended UMAP suite with key genes
- ✓ Organized output structure

### Potential Future Additions:
- Cell-type-specific condition DEG analysis
- SARS-CoV-2 variant comparison analysis
- Pathway enrichment analysis (GO/KEGG/Reactome)
- Trajectory analysis (differentiation pathways)
- Cell-cell interaction analysis
- HTML report generation with embedded plots

## Files Modified

**Primary Script:**
- `scripts/run_jmv_ali.py` - Enhanced from 451 to ~800+ lines

**New Functions Added:**
- `create_output_directories()` - Creates organized directory structure
- `generate_qc_visualizations()` - Comprehensive QC plotting
- `generate_integration_visualizations()` - Integration quality assessment
- `generate_marker_plots()` - Cell type marker visualizations
- `run_condition_deg_analysis()` - Condition-level DEG comparisons
- `generate_deg_visualizations()` - Volcano and MA plots
- `generate_extended_visualizations()` - Extended UMAP suite

## References

This enhanced pipeline builds upon:
- **Scanpy**: Wolf et al., Genome Biology 2018
- **Harmony**: Korsunsky et al., Nature Methods 2019
- **Statistical QC**: Lun et al., PLOS Computational Biology 2016
- **Leiden Clustering**: Traag et al., Scientific Reports 2019

## Contact & Support

For questions about the analysis pipeline:
- Pipeline location: `/mnt2/kwy/scrna_agent/scripts/run_jmv_ali.py`
- Output directory: `/mnt2/kwy/scrna_agent/results/jmv_ali_statistical/`
- Log files: Check console output or redirect to log file

---

*Enhanced pipeline version: 2.0*
*Last updated: 2026-01-15*
*Dataset: JMV_ALI (24 samples, ALI organoid viral infection study)*
