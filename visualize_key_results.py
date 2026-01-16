#!/usr/bin/env python3
"""
Quick visualization script to display key results from JMV ALI analysis
"""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load the processed data
print("Loading processed data...")
adata = sc.read_h5ad('/mnt2/kwy/scrna_agent/JMV_ALI_results/JMV_ALI_harmony_integrated.h5ad')

output_dir = Path('/mnt2/kwy/scrna_agent/JMV_ALI_results/figures')

print(f"\nDataset: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"Clusters: {adata.obs['louvain_0.8'].nunique()}")
print(f"Conditions: {adata.obs['condition'].nunique()}")

# Create a summary figure
print("\nGenerating summary figure...")
fig = plt.figure(figsize=(20, 16))
gs = fig.add_gridspec(4, 3, hspace=0.3, wspace=0.3)

# Row 1: UMAP plots
ax1 = fig.add_subplot(gs[0, 0])
sc.pl.umap(adata, color='louvain_0.8', legend_loc='on data', ax=ax1, show=False, title='Clusters (Louvain)')

ax2 = fig.add_subplot(gs[0, 1])
sc.pl.umap(adata, color='condition', ax=ax2, show=False, title='Condition')

ax3 = fig.add_subplot(gs[0, 2])
sc.pl.umap(adata, color='sample', ax=ax3, show=False, title='Sample', legend_fontsize=6)

# Row 2: QC metrics on UMAP
ax4 = fig.add_subplot(gs[1, 0])
sc.pl.umap(adata, color='n_genes_by_counts', ax=ax4, show=False, title='Number of genes', cmap='viridis')

ax5 = fig.add_subplot(gs[1, 1])
sc.pl.umap(adata, color='total_counts', ax=ax5, show=False, title='Total UMI counts', cmap='viridis')

ax6 = fig.add_subplot(gs[1, 2])
sc.pl.umap(adata, color='pct_counts_mt', ax=ax6, show=False, title='% Mitochondrial', cmap='viridis')

# Row 3: Cell type markers
markers_to_plot = ['TP63', 'FOXJ1', 'SCGB1A1']
for idx, marker in enumerate(markers_to_plot):
    if marker in adata.var_names:
        ax = fig.add_subplot(gs[2, idx])
        sc.pl.umap(adata, color=marker, ax=ax, show=False, title=f'{marker} expression',
                   vmax='p99', cmap='Reds')

# Row 4: Viral receptors
viral_receptors = ['ACE2', 'DPP4', 'TMPRSS2']
for idx, receptor in enumerate(viral_receptors):
    if receptor in adata.var_names:
        ax = fig.add_subplot(gs[3, idx])
        sc.pl.umap(adata, color=receptor, ax=ax, show=False, title=f'{receptor} expression',
                   vmax='p99', cmap='plasma')

plt.savefig(output_dir / 'summary_overview.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: {output_dir / 'summary_overview.png'}")

# Create condition-specific analysis
print("\nGenerating condition comparison figure...")
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Cell count bar plot
ax = axes[0, 0]
condition_counts = adata.obs['condition'].value_counts().sort_index()
condition_counts.plot(kind='bar', ax=ax, color='steelblue')
ax.set_xlabel('Condition')
ax.set_ylabel('Number of cells')
ax.set_title('Cells per condition')
ax.tick_params(axis='x', rotation=45)

# Cluster distribution per condition - stacked bar
ax = axes[0, 1]
ct = pd.crosstab(adata.obs['condition'], adata.obs['louvain_0.8'])
ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100
ct_pct.T.plot(kind='bar', stacked=True, ax=ax, legend=False)
ax.set_xlabel('Cluster')
ax.set_ylabel('Percentage')
ax.set_title('Cluster composition by condition')
ax.tick_params(axis='x', rotation=0)

# Add legend outside
ax.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

# QC metrics by condition - violin plots
ax = axes[0, 2]
df = adata.obs[['condition', 'n_genes_by_counts']].copy()
sns.violinplot(data=df, x='condition', y='n_genes_by_counts', ax=ax)
ax.set_xlabel('Condition')
ax.set_ylabel('Genes per cell')
ax.set_title('Gene counts by condition')
ax.tick_params(axis='x', rotation=45)

ax = axes[1, 0]
df = adata.obs[['condition', 'total_counts']].copy()
sns.violinplot(data=df, x='condition', y='total_counts', ax=ax)
ax.set_xlabel('Condition')
ax.set_ylabel('Total counts')
ax.set_title('UMI counts by condition')
ax.tick_params(axis='x', rotation=45)

ax = axes[1, 1]
df = adata.obs[['condition', 'pct_counts_mt']].copy()
sns.violinplot(data=df, x='condition', y='pct_counts_mt', ax=ax)
ax.set_xlabel('Condition')
ax.set_ylabel('% Mitochondrial')
ax.set_title('Mitochondrial % by condition')
ax.tick_params(axis='x', rotation=45)

# Cluster size distribution
ax = axes[1, 2]
cluster_counts = adata.obs['louvain_0.8'].value_counts().sort_index()
ax.bar(range(len(cluster_counts)), cluster_counts.values, color='coral')
ax.set_xlabel('Cluster')
ax.set_ylabel('Number of cells')
ax.set_title('Cells per cluster')
ax.set_xticks(range(len(cluster_counts)))
ax.set_xticklabels(cluster_counts.index, rotation=0)

plt.tight_layout()
plt.savefig(output_dir / 'condition_comparison.png', dpi=300, bbox_inches='tight')
plt.close()

print(f"Saved: {output_dir / 'condition_comparison.png'}")

# Print summary statistics
print("\n" + "="*80)
print("SUMMARY STATISTICS")
print("="*80)

print("\nCells per condition:")
print(condition_counts)

print("\nTop 5 clusters by size:")
print(cluster_counts.head())

print("\nQC metrics summary:")
for condition in sorted(adata.obs['condition'].unique()):
    subset = adata[adata.obs['condition'] == condition]
    print(f"\n{condition}:")
    print(f"  Cells: {subset.n_obs}")
    print(f"  Median genes: {subset.obs['n_genes_by_counts'].median():.0f}")
    print(f"  Median UMIs: {subset.obs['total_counts'].median():.0f}")
    print(f"  Median % MT: {subset.obs['pct_counts_mt'].median():.2f}%")

print("\n" + "="*80)
print("Visualization complete!")
print("="*80)
