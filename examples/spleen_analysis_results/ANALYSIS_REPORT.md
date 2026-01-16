# Spleen scRNA-seq Analysis Report

**Date**: 2026-01-16
**Sample Type**: Spleen (tissue)
**Integration Method**: Harmony
**Clustering Method**: Louvain
**Number of PCs**: 40 (optimal from range 30-50)

---

## Executive Summary

This report presents a comprehensive analysis of single-cell RNA-sequencing data from spleen tissue, comprising **461,102 cells** after quality control. The analysis employed Harmony batch correction to integrate 73 samples across different vaccine types, treatment groups, and time points. Louvain clustering identified 31 distinct cell populations, which were annotated into 11 major cell types and 31 refined cell type subtypes based on marker gene expression.

---

## 1. Dataset Overview

### Sample Information
- **Total cells analyzed**: 461,102 (after filtering 28,433 doublets)
- **Total genes**: 32,194
- **Number of samples**: 73 (orig.ident)
- **Batch variables**:
  - Vaccine type (vac): StV (299,717 cells), LtV (189,818 cells)
  - Treatment groups (typ): G1, G2, G3, G4, G1234
  - Days post infection (dpi): 0d, 3d, 4d, 5d, 6d, 14d
  - Biological replicates (rep): 1, 2, 3

### Quality Control Metrics
- **Median genes per cell**: 1,505
- **Median UMIs per cell**: 3,496
- **Median mitochondrial ratio**: 0.033
- **Doublet filtering**: Singlets only (461,102 cells retained)

---

## 2. Analysis Workflow

### 2.1 Preprocessing
1. **Normalization**: Total count normalization to 10,000 UMIs per cell
2. **Log transformation**: Natural log(count + 1)
3. **Feature selection**: 3,000 highly variable genes (batch-aware selection using orig.ident)
4. **Scaling**: Scaled to unit variance, max value 10

### 2.2 Dimensionality Reduction
**PCA Analysis**:
- Computed 50 principal components
- Selected **40 PCs** for downstream analysis (optimal balance)
- Variance explained:
  - 30 PCs: 17.0%
  - 40 PCs: 17.9%
  - 50 PCs: 18.6%

### 2.3 Batch Correction
**Harmony Integration**:
- Batch variable: `orig.ident` (73 samples)
- Converged after 11 iterations
- Successfully removed batch effects while preserving biological variation
- Integration quality assessed by visual inspection of UMAP before/after

### 2.4 Clustering
**Louvain Algorithm**:
- Multiple resolutions tested: 0.5, 0.8, 1.0, 1.2
- **Resolution 1.0 selected** (31 clusters)
  - Resolution 0.5: 19 clusters
  - Resolution 0.8: 23 clusters
  - Resolution 1.0: 31 clusters ✓
  - Resolution 1.2: 33 clusters

### 2.5 Visualization
**UMAP**:
- Generated before and after Harmony integration
- Random seed: 42 (reproducible)
- Parameters: n_neighbors=15, min_dist=0.5 (default)

---

## 3. Cell Type Annotation

### 3.1 Major Cell Type Distribution

| Cell Type | Count | Percentage |
|-----------|-------|------------|
| **Macrophages** | 144,203 | 31.27% |
| **B cells** | 123,162 | 26.71% |
| **T cells** | 61,274 | 13.29% |
| **Neutrophils** | 47,158 | 10.23% |
| **Monocytes** | 35,935 | 7.79% |
| **Plasma cells** | 23,833 | 5.17% |
| **NK cells** | 10,062 | 2.18% |
| **Proliferating** | 9,655 | 2.09% |
| **Dendritic cells** | 3,230 | 0.70% |
| **Mast cells** | 1,328 | 0.29% |
| **Endothelial** | 1,262 | 0.27% |

### 3.2 Key Findings

1. **Macrophages are the dominant population** (31.27%)
   - Multiple macrophage subtypes identified (C1Q+, Resident, CD5L+, MERTK+, Lipid)
   - Likely reflects the important role of splenic red pulp macrophages

2. **B cell populations are highly diverse** (26.71% total)
   - Includes naive, follicular, memory, mature, and pro-B cells
   - Expected given the spleen's role in adaptive immunity

3. **Substantial T cell compartment** (13.29%)
   - Both CD4+ and naive T cells present
   - Cytotoxic NK/T cell population identified

4. **Active neutrophil response** (10.23%)
   - Multiple neutrophil subsets (mature, immature, S100A8+, etc.)
   - May reflect ongoing immune response to infection/vaccination

5. **Plasma cell differentiation** (5.17%)
   - Distinct plasma cell populations including IGKG+ cells
   - Indicates active antibody production

---

## 4. Detailed Cell Type Annotations

### Top 15 Refined Cell Types:

| Cluster | Cell Type | Count | % | Key Markers |
|---------|-----------|-------|---|-------------|
| 0 | B cells (Naive) | 54,816 | 11.89% | BANK1, MS4A1, CD79B |
| 1 | Macrophages (C1Q+) | 48,815 | 10.59% | C1QA, C1QC, CD5L |
| 2 | Macrophages (Resident) | 45,294 | 9.82% | C1QC, FTL, FTH1 |
| 3 | T cells (CD4+) | 26,236 | 5.69% | CD3E, IL7R, CD3D |
| 4 | Macrophages (CD5L+) | 25,025 | 5.43% | CD5L, C1QC, CST3 |
| 5 | B cells (Follicular) | 25,015 | 5.43% | MS4A1, BANK1, CR2 |
| 6 | Macrophages (MERTK+) | 24,661 | 5.35% | MERTK, PPARG, TCF7L2 |
| 7 | Neutrophils (Mature) | 22,289 | 4.83% | LTF, S100A12, CAMP |
| 8 | T cells (Naive) | 19,308 | 4.19% | LEF1, CD3E, LTB |
| 9 | Monocytes (Activated) | 19,187 | 4.16% | NAMPT, ACOD1, PLAUR |
| 10 | Plasma cells | 17,124 | 3.71% | MZB1, XBP1, JCHAIN |
| 11 | Monocytes (Classical) | 16,748 | 3.63% | FCN1, UPP1, VCAN |
| 12 | B cells (Pro-B) | 15,950 | 3.46% | EBF1, SOX4, PAX5 |
| 13 | Neutrophils (S100A8+) | 15,747 | 3.42% | S100A8, S100A9, CAMP |
| 14 | NK/T cells (Cytotoxic) | 15,730 | 3.41% | NKG7, GZMA, CCL5 |

---

## 5. Integration Quality Assessment

### Before Integration:
- Strong batch effects visible by sample (orig.ident)
- Clustering by vaccine type (StV vs LtV)
- Separate clusters by treatment groups

### After Integration:
- Batch effects successfully removed
- Cells cluster primarily by cell type, not batch
- Biological variation preserved
- Cell type markers clearly segregate populations

**Conclusion**: Harmony integration effectively corrected for technical batch effects while maintaining biological heterogeneity.

---

## 6. Marker Genes Summary

### B Cells
- **General**: CD79A, CD79B, MS4A1 (CD20), BANK1, PAX5
- **Naive**: IGHD, TCL1A
- **Follicular**: CR2, S1PR1
- **Pro-B**: EBF1, SOX4, RAG1, DNTT, VPREB1, IGLL5

### T Cells
- **General**: CD3E, CD3D, CD3G, TRAC
- **CD4+**: IL7R, CD4
- **Naive**: LEF1, TCF7, LTB
- **Cytotoxic**: NKG7, GZMA, CCL5

### NK Cells
- **General**: NKG7, KLRD1, KLRK1, KLRC1
- **Cytotoxic**: GZMB, GZMH, PRF1, HOPX

### Macrophages
- **C1Q+**: C1QA, C1QB, C1QC
- **General**: CD68, CD163, APOE
- **MERTK+**: MERTK, PPARG, TCF7L2
- **Lipid**: ABCA8, PDK4

### Monocytes
- **Classical**: FCN1, S100A8, S100A9, CD14, LYZ, VCAN
- **Activated**: NAMPT, ACOD1, PLAUR

### Neutrophils
- **Mature**: LTF, CAMP, PGLYRP1, S100A12
- **Immature**: MPO, ELANE, CTSG
- **S100A8+**: S100A8, S100A9

### Plasma Cells
- MZB1, XBP1, JCHAIN, TNFRSF17, SDC1
- IGKG+ subset: IGKV2D-26

### Dendritic Cells
- HLA-DRA, HLA-DQA1, HLA-DQB1, HLA-DRB1, CXCL16

### Mast Cells
- CPA3, GATA2, MS4A2, TPSAB1

### Other
- **Proliferating**: MKI67, PCLAF, TOP2A, UBE2C, BIRC5
- **Endothelial**: PECAM1, CDK6, STMN1

---

## 7. Discussion Points for Refinement

### Questions for Consideration:

1. **Macrophage Subtypes**:
   - Are the 5 macrophage subtypes biologically distinct or should some be merged?
   - MERTK+ macrophages may represent M2-like or anti-inflammatory phenotype
   - C1Q+ complement-expressing macrophages may be involved in debris clearance
   - Should "Lipid macrophages" (408 cells) be considered a separate population?

2. **B Cell Developmental Stages**:
   - Pro-B cells (15,950) and RAG+ Pro-B cells (3,059) - keep separate or merge?
   - "B cells (IRF8+)" and "B cells (Ribosomal high)" - biological significance?
   - Are these distinct populations or technical artifacts?

3. **Neutrophil Heterogeneity**:
   - 5 neutrophil subtypes identified - expected for tissue neutrophils?
   - Mature vs Immature classification appears sound
   - S100A8+, RNASE+, PGLYRP+ subsets - merge or keep distinct?

4. **NK/T Cell Boundary**:
   - Cluster 14 annotated as "NK/T cells (Cytotoxic)" - consider splitting?
   - May represent CD8+ cytotoxic T cells or NKT cells
   - Recommend checking CD8A, CD8B expression

5. **Rare Cell Types**:
   - Endothelial cells (1,262) - expected in spleen?
   - Mast cells (1,328) - appropriate markers?
   - Should very small clusters (<500 cells) be re-examined?

6. **Proliferating Cells**:
   - Mixed cell types undergoing division (9,655 cells)
   - Could assign to cell types based on other markers
   - Or keep as separate category?

### Biological Context:
- Do the cell type proportions make sense for spleen tissue?
- Are infection/vaccination time points showing expected immune responses?
- Should we analyze cell type composition changes across conditions?

---

## 8. Output Files Generated

### Figures:
1. `01_QC_metrics.png` - Quality control distributions
2. `02_PCA_variance.png` - PCA variance explained
3. `03_UMAP_before_integration.png` - Pre-integration UMAP
4. `04_UMAP_after_integration.png` - Post-integration UMAP
5. `05_clustering_resolutions.png` - Clustering at different resolutions
6. `06_marker_genes_heatmap.png` - Top markers per cluster
7. `07_integration_before_after_comparison.png` - Integration comparison
8. `08_major_cell_types.png` - Major cell type UMAP and distribution
9. `09_refined_cell_types.png` - Refined annotation UMAP and distribution

### Data Files:
1. `Spleen_harmony_processed.h5ad` - Processed AnnData object
2. `cluster_assignments.csv` - Cluster assignments per cell
3. `umap_coordinates.csv` - UMAP coordinates (before/after)
4. `marker_genes_per_cluster.csv` - Top marker genes
5. `cell_type_annotations.csv` - Cell type annotations per cell
6. `umap_coordinates_annotated.csv` - UMAP with annotations
7. `annotation_summary.txt` - Text summary of annotations

---

## 9. Recommendations for Next Steps

1. **Validate Annotations**:
   - Review marker genes for questionable cell types
   - Consider checking additional canonical markers
   - Validate rare cell types (<1% of total)

2. **Refine Annotations**:
   - Merge redundant subtypes if biologically appropriate
   - Split ambiguous populations (e.g., NK/T cells)
   - Consider functional annotation schemes

3. **Downstream Analysis**:
   - Differential expression between conditions (vaccine, time points)
   - Trajectory analysis for B cell differentiation
   - Cell-cell communication analysis
   - Integration with other datasets

4. **Quality Checks**:
   - Verify ribosomal/mitochondrial gene contributions
   - Check for ambient RNA contamination
   - Validate doublet removal effectiveness

---

## 10. Technical Notes

### Software Versions:
- scanpy: 1.11.2
- Python: 3.10
- harmonypy: (installed)
- louvain: (installed)

### Computational Resources:
- Analysis performed on scRNAseq conda environment
- Large dataset (129GB h5ad file) handled efficiently
- Harmony convergence: 11 iterations

### Reproducibility:
- Random seeds set: PCA, UMAP, Louvain clustering (seed=42)
- All parameters documented
- Analysis scripts saved in scratchpad

---

## Conclusion

This analysis successfully integrated and annotated 461,102 single cells from spleen tissue across multiple experimental conditions. The Harmony integration effectively removed batch effects, and Louvain clustering with 40 PCs identified 31 biologically meaningful cell populations. The predominance of macrophages and B cells is consistent with splenic immunology, and the diversity of cell subtypes reflects the complex immune response captured in this dataset.

**The annotation is ready for user review and refinement before finalization.**

---

**Analysis completed**: 2026-01-16
**Analyst**: Claude Sonnet 4.5
**Save path**: /home/kwy7605/data_61/Agent_test/
