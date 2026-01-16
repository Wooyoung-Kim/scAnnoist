# Cell Type Annotation Performance Evaluation

## Dataset Information

- Total cells: 461,102
- Number of clusters: 31

## Annotation Methods Compared

**Original Method:**
- Coarse level: 8 cell types
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
| 0 | 54816 | B | B cells | 0.9261711908931699 | 1.0 |
| 1 | 48815 | MQ | Macrophages | 0.9991191232203216 | 1.0 |
| 2 | 45294 | MQ | Macrophages | 0.9825363182761514 | 1.0 |
| 3 | 26236 | T | T cells | 0.9993520353712456 | 1.0 |
| 4 | 25025 | MQ | Macrophages | 0.998001998001998 | 1.0 |
| 5 | 25015 | B | B cells | 0.998360983409954 | 1.0 |
| 6 | 24661 | MQ | Macrophages | 0.9545841612262277 | 1.0 |
| 7 | 22289 | Neutrophil | Neutrophils | 0.9380860514154964 | 1.0 |
| 8 | 19308 | T | T cells | 0.9804226227470478 | 1.0 |
| 9 | 19187 | Neutrophil | Monocytes | 0.9163496117162663 | 1.0 |

## Files Generated

- `annotation_metrics.csv`: Detailed performance metrics
- `cluster_level_agreement.csv`: Per-cluster annotation comparison
- `confusion_matrix_coarse.png`: Coarse-level confusion matrix
- `confusion_matrix_fine.png`: Fine-level confusion matrix
- `purity_comparison.png`: Cluster purity comparison plots
- `EVALUATION_REPORT.md`: This report
