# JMV ALI Organoid scRNA-seq Analysis Summary

**Analysis Date:** 2026-01-15
**Sample Type:** ALI organoid (non-immune)
**Integration Method:** Harmony
**Clustering Method:** Louvain
**Number of PCs:** 50

---

## Dataset Overview

### Samples Analyzed
- **Total Samples:** 15 (3 replicates per condition)
- **Conditions:**
  - Control (3 replicates)
  - HCoV_OC43 (3 replicates)
  - MERS (3 replicates)
  - SARS-CoV-1 (3 replicates)
  - SARS-CoV-2 (3 replicates)

### Data Quality

**Before Filtering:**
- Total cells: 50,028
- Total genes: 78,733
- Median genes per cell: 642
- Median UMIs per cell: 1,193
- Median % mitochondrial: 11.29%

**After Quality Control:**
- Total cells: 35,245 (70.5% retained)
- Total genes: 49,285 (62.6% retained)
- Highly variable genes: 2,000

**QC Filtering Criteria:**
- Minimum genes per cell: 200
- Minimum cells per gene: 3
- Maximum genes per cell: 6,000
- Maximum % mitochondrial: 20%

---

## Cells per Condition

| Condition | Number of Cells | Percentage |
|-----------|----------------|------------|
| Control | 4,169 | 11.8% |
| HCoV_OC43 | 2,717 | 7.7% |
| MERS | 2,906 | 8.2% |
| SARS-CoV-1 | 2,260 | 6.4% |
| SARS-CoV-2 | 23,193 | 65.8% |

**Note:** SARS-CoV-2 condition has significantly more cells, primarily from sample JMV_SARS-CoV-2_L-3 (21,539 cells).

---

## Cells per Sample

| Sample | Cells | Condition |
|--------|-------|-----------|
| JMV_Control-1 | 592 | Control |
| JMV_Control-2 | 2,023 | Control |
| JMV_Control-3 | 1,554 | Control |
| JMV_HCoV_OC43-1 | 1,477 | HCoV_OC43 |
| JMV_HCoV_OC43-2 | 803 | HCoV_OC43 |
| JMV_HCoV_OC43-3 | 437 | HCoV_OC43 |
| JMV_MERS-1 | 455 | MERS |
| JMV_MERS-2 | 1,544 | MERS |
| JMV_MERS-3 | 907 | MERS |
| JMV_SARS-CoV-1-1 | 553 | SARS-CoV-1 |
| JMV_SARS-CoV-1-2 | 580 | SARS-CoV-1 |
| JMV_SARS-CoV-1-3 | 1,127 | SARS-CoV-1 |
| JMV_SARS-CoV-2_L-1 | 624 | SARS-CoV-2 |
| JMV_SARS-CoV-2_L-2 | 1,030 | SARS-CoV-2 |
| JMV_SARS-CoV-2_L-3 | 21,539 | SARS-CoV-2 |

---

## Harmony Integration

**Integration Strategy:**
- Batch correction performed across all 15 samples
- Used 50 principal components
- Harmony converged after 19 iterations

**Integration Results:**
- Successfully integrated samples while preserving biological variation
- Reduced batch effects between samples
- Maintained separation of cell types and conditions

---

## Clustering Analysis

**Louvain Clustering Results:**
- Resolution 0.5: 20 clusters
- Resolution 0.8: 21 clusters (default)
- Resolution 1.0: 23 clusters

### Cells per Cluster (Resolution 0.8)

| Cluster | Number of Cells | Percentage |
|---------|----------------|------------|
| 0 | 8,364 | 23.7% |
| 1 | 4,491 | 12.7% |
| 2 | 3,793 | 10.8% |
| 3 | 3,344 | 9.5% |
| 4 | 2,572 | 7.3% |
| 5 | 2,252 | 6.4% |
| 6 | 2,065 | 5.9% |
| 7 | 2,043 | 5.8% |
| 8 | 1,998 | 5.7% |
| 9 | 1,245 | 3.5% |
| 10 | 1,030 | 2.9% |
| 11 | 669 | 1.9% |
| 12 | 475 | 1.3% |
| 13 | 343 | 1.0% |
| 14 | 236 | 0.7% |
| 15 | 179 | 0.5% |
| 16 | 72 | 0.2% |
| 17 | 22 | 0.1% |
| 18 | 19 | 0.1% |
| 19 | 18 | 0.1% |
| 20 | 15 | 0.0% |

---

## Cell Type Markers Detected

The following airway epithelial cell markers were detected in the dataset:

### Basal Cells
- TP63
- KRT5

### Ciliated Cells
- FOXJ1
- RSPH1

### Club/Secretory Cells
- SCGB1A1
- SCGB3A2
- MUC5AC

### Goblet Cells
- MUC5AC
- MUC5B

### Ionocytes
- FOXI1
- CFTR

### Neuroendocrine
- CHGA
- CHGB

---

## Viral Receptor Expression

The following viral receptors were detected and visualized:

- **ACE2** - SARS-CoV-2 receptor
- **TMPRSS2** - Serine protease for viral entry
- **DPP4** - MERS-CoV receptor
- **ANPEP** - Coronavirus receptor (aminopeptidase N)

These markers can be used to identify cell populations susceptible to viral infection.

---

## Top Marker Genes (Selected Clusters)

### Cluster 0 (Largest cluster - 23.7% of cells)
Top markers include:
- SPRR3 (log2FC: 1.71)
- FTH1 (log2FC: 0.95)
- CSTB (log2FC: 1.10)
- ADIRF (log2FC: 1.54)
- KRT17 (log2FC: 1.48)

Suggests a basal/squamous cell type with stress response signatures.

*Note: Complete marker gene lists for all 21 clusters are available in separate CSV files.*

---

## Output Files

### Main Results Directory
`/mnt2/kwy/scrna_agent/JMV_ALI_results/`

### Key Files

1. **PowerPoint Report:** 📊
   - `JMV_ALI_Organoid_scRNA-seq_Analysis_report.pptx` - Comprehensive presentation (28 MB)
   - Contains all key figures, statistics, and marker gene tables
   - Ready for presentations and reports

2. **Processed Data:**
   - `JMV_ALI_harmony_integrated.h5ad` - Integrated AnnData object (7.3 GB)
   - `metadata.csv` - Cell metadata including cluster assignments

3. **Marker Genes:**
   - `markers_cluster_0.csv` to `markers_cluster_20.csv` - Differential expression results for each cluster

4. **Figures Directory:**
   `JMV_ALI_results/figures/`

   - `qc_metrics_before_filtering.png` - Quality control metrics
   - `highly_variable_genes.png` - Highly variable gene selection
   - `pca_variance_ratio.png` - PCA variance explained
   - `pca_before_harmony.png` - PCA plots before integration
   - `harmony_comparison.png` - Before/after Harmony integration
   - `umap_overview.png` - Comprehensive UMAP visualizations
   - `umap_by_condition_faceted.png` - UMAP split by condition
   - `cell_type_markers.png` - Expression of airway epithelial markers
   - `viral_receptors.png` - Expression of viral entry receptors
   - `cluster_composition_by_condition.png` - Cluster composition analysis
   - `cell_counts_per_cluster.png` - Cell counts by cluster and condition
   - `marker_genes_per_cluster.png` - Top marker genes for each cluster
   - `marker_genes_dotplot.png` - Dotplot of top markers

---

## Analysis Pipeline

1. **Data Loading:** Loaded 15 samples from 10X Genomics filtered matrices
2. **Quality Control:** Filtered low-quality cells and genes
3. **Normalization:** Log-normalization with target sum 10,000
4. **Feature Selection:** Identified 2,000 highly variable genes
5. **Dimensionality Reduction:** PCA with 50 components
6. **Batch Correction:** Harmony integration across samples
7. **Clustering:** Louvain clustering at multiple resolutions
8. **UMAP Visualization:** 2D embedding for visualization
9. **Marker Gene Detection:** Wilcoxon rank-sum test for differential expression
10. **Visualization:** Comprehensive plots for all aspects of the analysis

---

## Next Steps & Recommendations

1. **Cell Type Annotation:**
   - Use marker genes to annotate clusters with specific airway epithelial cell types
   - Consider automated annotation tools (e.g., SingleR, CellTypist)

2. **Differential Expression Analysis:**
   - Compare gene expression between conditions (Control vs. infected)
   - Identify condition-specific responses for each cell type

3. **Pathway Analysis:**
   - Perform Gene Ontology (GO) enrichment analysis on cluster markers
   - Identify biological pathways activated in response to viral infection

4. **Viral Tropism Analysis:**
   - Examine viral receptor expression across cell types
   - Correlate receptor expression with infection status

5. **Trajectory Analysis:**
   - Investigate differentiation trajectories in epithelial cells
   - Examine how infection affects cell state transitions

6. **Integration with Additional Data:**
   - If available, integrate with ATAC-seq, CITE-seq, or spatial data
   - Compare with published airway organoid datasets

---

## Technical Notes

- **Environment:** Python 3.9.18
- **Key Packages:**
  - scanpy (single-cell analysis)
  - harmonypy (batch correction)
  - python-igraph, louvain (clustering)
  - matplotlib, seaborn (visualization)

- **Computational Resources:**
  - Analysis completed in approximately 5-10 minutes
  - Peak memory usage for 35K cells dataset

---

## Contact & Support

For questions about this analysis or to request additional analyses:
- Analysis script: `JMV_ALI_analysis.py`
- Log file: `JMV_ALI_analysis.log`

---

**Analysis completed:** 2026-01-15 22:41 KST
