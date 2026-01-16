# Cell Type Annotation Performance Evaluation

## Dataset Information

- Total cells: 461,102
- Number of clusters: 31

## Annotation Methods Compared

**Original Method:**
- Coarse level: 8 cell types (normalized from 8 abbreviated forms)
- Fine level: 41 cell types

**Agent-based Method (scAnnoist):**
- Major cell types: 11 types
- Refined cell types: 31 types

## Performance Metrics

| Method | ARI | NMI | Homogeneity | Completeness | V-measure | N_cells | N_true_types | N_pred_types |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Original Coarse vs Agent Major | 0.8528614834732602 | 0.802198054259288 | 0.8337820480246741 | 0.7729195557733527 | 0.8021980542592881 | 461102 | 8 | 11 |
| Original Fine vs Agent Refined | 0.5899258187775831 | 0.7716200399309158 | 0.7613415143265497 | 0.7821798950380958 | 0.7716200399309159 | 461102 | 41 | 31 |
| Original Coarse vs Agent Refined | 0.3381998037477464 | 0.6303295113219296 | 0.8837678766647438 | 0.4898539537003394 | 0.6303295113219296 | 461102 | 8 | 31 |
| Original Fine vs Agent Major | 0.4262240093582391 | 0.6549217858297882 | 0.5180416092962151 | 0.8901129724362484 | 0.6549217858297881 | 461102 | 41 | 11 |

### Metric Interpretation

- **ARI (Adjusted Rand Index)**: Measures similarity between two clusterings (1.0 = perfect match, 0.0 = random)
- **NMI (Normalized Mutual Information)**: Measures mutual dependence (1.0 = perfect match)
- **Homogeneity**: Whether all clusters contain only members of a single class
- **Completeness**: Whether all members of a class are in the same cluster
- **V-measure**: Harmonic mean of homogeneity and completeness

## Cluster-Level Agreement

**Agreement Rate** (how many clusters have matching annotations):

- **Coarse Level**: 23/31 clusters (74.2%)
- **Fine Level**: 0/31 clusters (0.0%)

### All Clusters Comparison

| Cluster | N_cells | Original_Coarse | Agent_Major | Match | Original_Fine | Agent_Refined | Match |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0 | 54816 | B cells | B cells | ✓ | NB1 | B cells (Naive) | ✗ |
| 1 | 48815 | Macrophages | Macrophages | ✓ | MQ1 | Macrophages (C1Q+) | ✗ |
| 2 | 45294 | Macrophages | Macrophages | ✓ | MQ1 | Macrophages (Resident) | ✗ |
| 3 | 26236 | T cells | T cells | ✓ | CD4_TEM | T cells (CD4+) | ✗ |
| 4 | 25025 | Macrophages | Macrophages | ✓ | MQ21 | Macrophages (CD5L+) | ✗ |
| 5 | 25015 | B cells | B cells | ✓ | MB2 | B cells (Follicular) | ✗ |
| 6 | 24661 | Macrophages | Macrophages | ✓ | MQ22 | Macrophages (MERTK+) | ✗ |
| 7 | 22289 | Neutrophils | Neutrophils | ✓ | Neu2 | Neutrophils (Mature) | ✗ |
| 8 | 19308 | T cells | T cells | ✓ | CD4_N | T cells (Naive) | ✗ |
| 9 | 19187 | Neutrophils | Monocytes | ✗ | Neu1 | Monocytes (Activated) | ✗ |
| 10 | 17124 | Plasma cells | Plasma cells | ✓ | PC1 | Plasma cells | ✗ |
| 11 | 16748 | Monocytes | Monocytes | ✓ | IM | Monocytes (Classical) | ✗ |
| 12 | 15950 | B cells | B cells | ✓ | NB3 | B cells (Pro-B) | ✗ |
| 13 | 15747 | Neutrophils | Neutrophils | ✓ | Neu3 | Neutrophils (S100A8+) | ✗ |
| 14 | 15730 | T cells | T cells | ✓ | CD8_CTL | NK/T cells (Cytotoxic) | ✗ |
| 15 | 13753 | B cells | B cells | ✓ | NB2 | B cells (Mature) | ✗ |
| 16 | 9655 | Plasmablast | Proliferating | ✗ | PB1 | Proliferating cells | ✗ |
| 17 | 6709 | Plasma cells | Plasma cells | ✓ | PC1 | Plasma cells (IGKG+) | ✗ |
| 18 | 6622 | T cells | NK cells | ✗ | CD8_TEM | NK cells | ✗ |
| 19 | 6237 | Neutrophils | Neutrophils | ✓ | GMP | Neutrophils (Immature) | ✗ |
| 20 | 4677 | B cells | B cells | ✓ | MB1 | B cells (Memory) | ✗ |
| 21 | 3440 | T cells | NK cells | ✗ | NK | NK cells (GZMB+) | ✗ |
| 22 | 3230 | Monocytes | Dendritic cells | ✗ | cDC2 | Dendritic cells | ✗ |
| 23 | 3059 | B cells | B cells | ✓ | ImmB | Pro-B cells (RAG+) | ✗ |
| 24 | 2971 | Dendritic cells | B cells | ✗ | pDC | B cells (IRF8+) | ✗ |
| 25 | 2921 | B cells | B cells | ✓ | PC2 | B cells (Ribosomal high) | ✗ |
| 26 | 1818 | Neutrophils | Neutrophils | ✓ | Platelet | Neutrophils (RNASE+) | ✗ |
| 27 | 1328 | Monocytes | Mast cells | ✗ | Mast | Mast cells | ✗ |
| 28 | 1262 | Monocytes | Endothelial | ✗ | HSC | Endothelial cells | ✗ |
| 29 | 1067 | Neutrophils | Neutrophils | ✓ | Neu2 | Neutrophils (PGLYRP+) | ✗ |
| 30 | 408 | Macrophages | Macrophages | ✓ | MQ22 | Macrophages (Lipid) | ✗ |

## Cluster-level Purity

| Statistic | Original_Coarse_Purity | Original_Fine_Purity | Agent_Major_Purity | Agent_Refined_Purity |
| --- | --- | --- | --- | --- |
| count | 31.0000 | 31.0000 | 31.0000 | 31.0000 |
| mean | 0.9069 | 0.8006 | 1.0000 | 1.0000 |
| std | 0.1345 | 0.1931 | 0.0000 | 0.0000 |
| min | 0.4307 | 0.2889 | 1.0000 | 1.0000 |
| 25% | 0.8824 | 0.6264 | 1.0000 | 1.0000 |
| 50% | 0.9593 | 0.8949 | 1.0000 | 1.0000 |
| 75% | 0.9947 | 0.9369 | 1.0000 | 1.0000 |
| max | 1.0000 | 0.9963 | 1.0000 | 1.0000 |

## Top 10 Largest Clusters

| Cluster | N_cells | Original_Coarse | Agent_Major | Original_Coarse_Purity | Agent_Major_Purity |
| --- | --- | --- | --- | --- | --- |
| 0 | 54816 | B cells | B cells | 0.9261711908931699 | 1.0 |
| 1 | 48815 | Macrophages | Macrophages | 0.9991191232203216 | 1.0 |
| 2 | 45294 | Macrophages | Macrophages | 0.9825363182761514 | 1.0 |
| 3 | 26236 | T cells | T cells | 0.9993520353712456 | 1.0 |
| 4 | 25025 | Macrophages | Macrophages | 0.998001998001998 | 1.0 |
| 5 | 25015 | B cells | B cells | 0.998360983409954 | 1.0 |
| 6 | 24661 | Macrophages | Macrophages | 0.9545841612262277 | 1.0 |
| 7 | 22289 | Neutrophils | Neutrophils | 0.9380860514154964 | 1.0 |
| 8 | 19308 | T cells | T cells | 0.9804226227470478 | 1.0 |
| 9 | 19187 | Neutrophils | Monocytes | 0.9163496117162663 | 1.0 |

## Files Generated

- `annotation_metrics.csv`: Detailed performance metrics (ARI, NMI, etc.)
- `cluster_match_analysis.csv`: Cluster-by-cluster agreement analysis
- `cluster_level_agreement.csv`: Per-cluster purity metrics
- `confusion_matrix_coarse.png`: Coarse-level confusion matrix
- `confusion_matrix_fine.png`: Fine-level confusion matrix
- `purity_comparison.png`: Cluster purity comparison plots
- `EVALUATION_REPORT.md`: This report
