# Quick Start Guide: Enhanced JMV_ALI Analysis

## Running the Enhanced Analysis

### 1. Activate the Environment
```bash
conda activate scRNAseq
```

### 2. Navigate to Project Directory
```bash
cd /mnt2/kwy/scrna_agent
```

### 3. Run the Analysis
```bash
python scripts/run_jmv_ali.py
```

**Expected Runtime:** 15-30 minutes

### 4. Monitor Progress
The script will print progress messages to the console. Watch for these key steps:
- ✓ Loading 24 samples
- ✓ QC metrics calculation and visualization
- ✓ Adaptive QC filtering
- ✓ Normalization and HVG selection
- ✓ PCA and Harmony batch correction
- ✓ UMAP computation and clustering
- ✓ Marker gene identification
- ✓ Cell type annotation
- ✓ Condition-level DEG analysis
- ✓ Visualization generation

## Output Location

All results will be saved to:
```
/mnt2/kwy/scrna_agent/results/jmv_ali_statistical/
```

## Key Output Files

### Main Results
- `integrated.h5ad` - Complete processed dataset
- `report.md` - Analysis report with summary statistics

### Quality Control
- `qc/qc_metrics_before.png` - Pre-filtering QC plots
- `qc/qc_metrics_after.png` - Post-filtering QC plots
- `qc/qc_statistics.csv` - Per-sample QC statistics

### Batch Integration
- `integration/integration_comparison.png` - Before/after Harmony correction

### Cell Type Annotation
- `annotation/marker_dotplot.png` - Marker gene dotplot
- `annotation/celltype_markers_heatmap.png` - Marker expression heatmap
- `annotation/celltype_proportions.csv` - Cell type percentages
- `annotation/cell_metadata.csv` - Cell-level annotations

### Clustering
- `clustering/clusters.csv` - Cell cluster assignments
- `clustering/cluster_composition.png` - Composition by virus group
- `clustering/pca_variance.png` - PC variance explained

### Differential Expression
- `deg/markers.csv` - Top 25 marker genes per cluster
- `deg/markers_full.csv` - All cluster markers with statistics
- `deg/condition_comparisons/*.csv` - DEGs for each virus vs control
- `deg/condition_comparisons/volcano_*.png` - Volcano plots

### Visualizations
- `umap/umap_clusters.png` - UMAP by cluster
- `umap/umap_condition.png` - UMAP by viral condition
- `umap/umap_celltype.png` - UMAP by cell type
- `umap/umap_sample.png` - UMAP by sample
- `umap/umap_n_genes.png` - UMAP colored by gene count
- `umap/umap_key_markers.png` - UMAPs for key marker genes
- `umap/umap_viral_response_genes.png` - UMAPs for ISG genes

## What to Look For

### 1. Quality Control Results
Check `qc/qc_statistics.csv` to see:
- Number of cells per sample (should be ~3,000-5,000 per sample)
- Mean genes per cell (typically 1,000-2,000)
- Mean UMI counts (typically 2,000-5,000)
- MT% (should be <20%)

### 2. Batch Effects
Check `integration/integration_comparison.png`:
- **Top row**: Before Harmony - samples should cluster separately
- **Bottom row**: After Harmony - samples should be well-mixed

### 3. Cell Types
Check `annotation/celltype_markers_heatmap.png`:
- Look for clear expression patterns of:
  - **FOXJ1** (Ciliated cells)
  - **KRT5** (Basal cells)
  - **SCGB1A1** (Secretory cells)
  - **MUC5AC** (Goblet cells)

### 4. Viral Response
Check `deg/condition_comparisons/`:
- Each CSV file shows genes differentially expressed in viral vs control
- Look for interferon-stimulated genes (ISGs):
  - ISG15, IFIT1, IFIT2, IFIT3
  - MX1, MX2, OAS1
  - IFI6, IFI44L, IFITM3

### 5. Condition-Specific Signatures
Compare volcano plots to identify:
- **SARS-CoV-2 variants**: Do they show similar or distinct responses?
- **MERS vs SARS-CoV-1**: Are there unique signatures?
- **HCoV_OC43**: Milder seasonal coronavirus response?

## Common Questions

### Q: How do I load the results in Python?
```python
import scanpy as sc

# Load processed data
adata = sc.read_h5ad('results/jmv_ali_statistical/integrated.h5ad')

# View UMAP
sc.pl.umap(adata, color='celltype')

# Access cluster assignments
clusters = adata.obs['leiden']

# Access DEG results
import pandas as pd
degs = pd.read_csv('results/jmv_ali_statistical/deg/markers.csv')
```

### Q: How do I compare specific conditions?
```python
# Load condition-specific DEGs
mers_degs = pd.read_csv('results/jmv_ali_statistical/deg/condition_comparisons/MERS_vs_Control_degs.csv')

# Filter for significant genes
sig_degs = mers_degs[(mers_degs['pvals_adj'] < 0.05) & (mers_degs['logfoldchanges'].abs() > 0.5)]

# Top upregulated genes
top_up = sig_degs.nlargest(20, 'logfoldchanges')
print(top_up[['names', 'logfoldchanges', 'pvals_adj']])
```

### Q: How do I visualize specific genes?
```python
# Visualize a gene of interest
sc.pl.umap(adata, color=['ISG15', 'FOXJ1', 'KRT5'], cmap='Reds')

# Violin plot by condition
sc.pl.violin(adata, keys=['ISG15'], groupby='condition', rotation=90)

# Dotplot for multiple genes
genes = ['ISG15', 'IFIT1', 'MX1', 'OAS1']
sc.pl.dotplot(adata, genes, groupby='condition')
```

### Q: Where are cell type proportions?
Check two files:
1. `annotation/celltype_proportions.csv` - Percentages by condition
2. `annotation/celltype_composition.png` - Visual stacked bar chart

## Troubleshooting

### Script fails with "ModuleNotFoundError"
Make sure you activated the conda environment:
```bash
conda activate scRNAseq
```

### "Out of memory" error
The dataset has ~84,000 cells. If running on a machine with <32GB RAM, you may need to:
1. Close other applications
2. Process fewer samples
3. Use a high-memory node

### Plots look crowded
This is expected with 24 samples. To make plots clearer:
- Zoom in on saved PNG files
- Regenerate specific plots with larger figure sizes
- Focus on summary plots (by condition/group rather than individual samples)

## Next Steps

### 1. Explore Top DEGs
```bash
head -50 results/jmv_ali_statistical/deg/condition_comparisons/MERS_vs_Control_degs.csv
```

### 2. Check Viral Response Genes
Look for ISGs in the DEG tables:
```bash
grep -i "ISG15\|IFIT\|MX1\|MX2\|OAS" results/jmv_ali_statistical/deg/condition_comparisons/*.csv
```

### 3. Compare Cell Type Proportions
```bash
cat results/jmv_ali_statistical/annotation/celltype_proportions.csv
```

### 4. Generate Custom Plots
Use the integrated.h5ad file to create additional visualizations in Python/R.

## Support

For questions or issues:
- Script location: `/mnt2/kwy/scrna_agent/scripts/run_jmv_ali.py`
- Documentation: `/mnt2/kwy/scrna_agent/ENHANCEMENTS_SUMMARY.md`
- Log issues: Check console output for error messages

---

*Enhanced analysis pipeline version 2.0*
*Dataset: 24 ALI organoid samples, 8 viral conditions*
