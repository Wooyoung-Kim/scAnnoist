#!/usr/bin/env python3
"""
Evaluate cell type annotation performance between original and agent-based methods.

Compares:
1. Original annotations (original_cl_ct, original_sbcl_ct)
2. Agent-based annotations (major_cell_type, cell_type_refined)

Metrics:
- Adjusted Rand Index (ARI)
- Normalized Mutual Information (NMI)
- Homogeneity, Completeness, V-measure
- Cluster-level agreement
"""

import pandas as pd
import numpy as np
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    homogeneity_completeness_v_measure,
    confusion_matrix
)
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


def normalize_cell_type_names(label):
    """
    Normalize cell type names to a consistent format.
    Maps abbreviated forms to full names for fair comparison.
    """
    # Mapping dictionary for original abbreviated names
    name_mapping = {
        'B': 'B cells',
        'T': 'T cells',
        'MQ': 'Macrophages',
        'Neutrophil': 'Neutrophils',
        'PC': 'Plasma cells',
        'NK': 'NK cells',
        'DC': 'Dendritic cells',
        'Mono': 'Monocytes',
        'Monocyte': 'Monocytes',  # Alternative form
        'Mast': 'Mast cells',
        'Endo': 'Endothelial',
        'PB': 'Proliferating',  # Plasmablasts (proliferating B cells)
    }

    # Apply mapping if found
    if label in name_mapping:
        return name_mapping[label]

    # Return as-is if no mapping needed
    return label


def load_annotations(csv_path):
    """Load cell type annotations."""
    df = pd.read_csv(csv_path, index_col=0)
    print(f"Loaded {len(df):,} cells")
    print(f"Columns: {df.columns.tolist()}")

    # Normalize original annotation names
    print("\nNormalizing cell type names...")
    df['original_cl_ct_normalized'] = df['original_cl_ct'].apply(normalize_cell_type_names)

    # Check mapping results
    original_types = df['original_cl_ct'].unique()
    normalized_types = df['original_cl_ct_normalized'].unique()
    print(f"Original types ({len(original_types)}): {sorted(original_types)}")
    print(f"Normalized types ({len(normalized_types)}): {sorted(normalized_types)}")

    return df


def calculate_metrics(true_labels, pred_labels, method_name):
    """Calculate clustering comparison metrics."""
    # Remove any missing values
    mask = ~(pd.isna(true_labels) | pd.isna(pred_labels))
    true_labels = true_labels[mask]
    pred_labels = pred_labels[mask]

    ari = adjusted_rand_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)
    homo, comp, v_measure = homogeneity_completeness_v_measure(true_labels, pred_labels)

    metrics = {
        'Method': method_name,
        'ARI': ari,
        'NMI': nmi,
        'Homogeneity': homo,
        'Completeness': comp,
        'V-measure': v_measure,
        'N_cells': len(true_labels),
        'N_true_types': len(np.unique(true_labels)),
        'N_pred_types': len(np.unique(pred_labels))
    }

    return metrics


def cluster_level_agreement(df, cluster_col='louvain'):
    """Calculate cluster-level agreement between annotation methods."""
    clusters = df[cluster_col].unique()

    results = []
    for cluster in sorted(clusters):
        cluster_cells = df[df[cluster_col] == cluster]
        n_cells = len(cluster_cells)

        # Get most common cell type for each annotation method (using normalized names)
        original_coarse = cluster_cells['original_cl_ct_normalized'].mode()[0] if len(cluster_cells['original_cl_ct_normalized'].mode()) > 0 else 'NA'
        original_fine = cluster_cells['original_sbcl_ct'].mode()[0] if len(cluster_cells['original_sbcl_ct'].mode()) > 0 else 'NA'
        agent_major = cluster_cells['major_cell_type'].mode()[0] if len(cluster_cells['major_cell_type'].mode()) > 0 else 'NA'
        agent_refined = cluster_cells['cell_type_refined'].mode()[0] if len(cluster_cells['cell_type_refined'].mode()) > 0 else 'NA'

        # Calculate purity for each annotation (using normalized names)
        original_coarse_purity = (cluster_cells['original_cl_ct_normalized'] == original_coarse).sum() / n_cells
        original_fine_purity = (cluster_cells['original_sbcl_ct'] == original_fine).sum() / n_cells
        agent_major_purity = (cluster_cells['major_cell_type'] == agent_major).sum() / n_cells
        agent_refined_purity = (cluster_cells['cell_type_refined'] == agent_refined).sum() / n_cells

        results.append({
            'Cluster': cluster,
            'N_cells': n_cells,
            'Original_Coarse': original_coarse,
            'Original_Coarse_Purity': original_coarse_purity,
            'Original_Fine': original_fine,
            'Original_Fine_Purity': original_fine_purity,
            'Agent_Major': agent_major,
            'Agent_Major_Purity': agent_major_purity,
            'Agent_Refined': agent_refined,
            'Agent_Refined_Purity': agent_refined_purity
        })

    return pd.DataFrame(results)


def plot_confusion_matrix(true_labels, pred_labels, title, output_path):
    """Plot confusion matrix between two annotation methods."""
    # Get unique labels
    all_labels = sorted(set(list(true_labels.unique()) + list(pred_labels.unique())))

    # Calculate confusion matrix
    cm = confusion_matrix(true_labels, pred_labels, labels=all_labels)

    # Normalize by true labels (rows)
    cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]

    # Create figure
    plt.figure(figsize=(max(12, len(all_labels) * 0.6), max(10, len(all_labels) * 0.5)))
    sns.heatmap(cm_normalized,
                xticklabels=all_labels,
                yticklabels=all_labels,
                cmap='Blues',
                fmt='.2f',
                cbar_kws={'label': 'Proportion'},
                annot=False)

    plt.title(title, fontsize=14, pad=20)
    plt.xlabel('Agent-based Annotation', fontsize=12)
    plt.ylabel('Original Annotation', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved confusion matrix: {output_path}")
    plt.close()


def plot_purity_comparison(cluster_agreement_df, output_path):
    """Plot cluster purity comparison."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    metrics = [
        ('Original_Coarse_Purity', 'Agent_Major_Purity', 'Coarse Level', axes[0, 0]),
        ('Original_Fine_Purity', 'Agent_Refined_Purity', 'Fine Level', axes[0, 1]),
        ('Original_Coarse_Purity', 'Original_Fine_Purity', 'Original: Coarse vs Fine', axes[1, 0]),
        ('Agent_Major_Purity', 'Agent_Refined_Purity', 'Agent: Major vs Refined', axes[1, 1])
    ]

    for x_col, y_col, title, ax in metrics:
        ax.scatter(cluster_agreement_df[x_col],
                  cluster_agreement_df[y_col],
                  s=cluster_agreement_df['N_cells'] / 100,
                  alpha=0.6)

        # Add diagonal line
        ax.plot([0, 1], [0, 1], 'r--', alpha=0.5)

        ax.set_xlabel(x_col.replace('_', ' '), fontsize=11)
        ax.set_ylabel(y_col.replace('_', ' '), fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.set_xlim(0, 1.05)
        ax.set_ylim(0, 1.05)
        ax.grid(alpha=0.3)

        # Add correlation
        corr = cluster_agreement_df[[x_col, y_col]].corr().iloc[0, 1]
        ax.text(0.05, 0.95, f'Correlation: {corr:.3f}',
                transform=ax.transAxes,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved purity comparison: {output_path}")
    plt.close()


def main():
    # Configuration
    input_csv = "/home/kwy7605/data_61/Agent_test/cell_type_annotations.csv"
    output_dir = Path("/mnt2/kwy/scrna_agent/examples/spleen_analysis_results/performance_evaluation")
    output_dir.mkdir(exist_ok=True, parents=True)

    print("="*80)
    print("Cell Type Annotation Performance Evaluation")
    print("="*80)

    # Load data
    print("\n1. Loading annotations...")
    df = load_annotations(input_csv)

    # Calculate metrics
    print("\n2. Calculating performance metrics...")
    print("\nComparing annotation levels (with normalized names):")

    comparisons = [
        ('original_cl_ct_normalized', 'major_cell_type', 'Original Coarse vs Agent Major'),
        ('original_sbcl_ct', 'cell_type_refined', 'Original Fine vs Agent Refined'),
        ('original_cl_ct_normalized', 'cell_type_refined', 'Original Coarse vs Agent Refined'),
        ('original_sbcl_ct', 'major_cell_type', 'Original Fine vs Agent Major'),
    ]

    all_metrics = []
    for true_col, pred_col, name in comparisons:
        metrics = calculate_metrics(df[true_col], df[pred_col], name)
        all_metrics.append(metrics)
        print(f"\n{name}:")
        print(f"  ARI: {metrics['ARI']:.4f}")
        print(f"  NMI: {metrics['NMI']:.4f}")
        print(f"  V-measure: {metrics['V-measure']:.4f}")
        print(f"  True types: {metrics['N_true_types']}, Pred types: {metrics['N_pred_types']}")

    # Save metrics
    metrics_df = pd.DataFrame(all_metrics)
    metrics_path = output_dir / "annotation_metrics.csv"
    metrics_df.to_csv(metrics_path, index=False)
    print(f"\nSaved metrics: {metrics_path}")

    # Cluster-level agreement
    print("\n3. Calculating cluster-level agreement...")
    cluster_agreement = cluster_level_agreement(df)
    cluster_path = output_dir / "cluster_level_agreement.csv"
    cluster_agreement.to_csv(cluster_path, index=False)
    print(f"Saved cluster agreement: {cluster_path}")

    # Summary statistics
    print("\n" + "="*80)
    print("Cluster-level Purity Summary:")
    print("="*80)
    for col in ['Original_Coarse_Purity', 'Original_Fine_Purity',
                'Agent_Major_Purity', 'Agent_Refined_Purity']:
        mean_purity = cluster_agreement[col].mean()
        median_purity = cluster_agreement[col].median()
        print(f"{col:30s}: Mean={mean_purity:.3f}, Median={median_purity:.3f}")

    # Plot confusion matrices
    print("\n4. Generating confusion matrices...")
    plot_confusion_matrix(
        df['original_cl_ct_normalized'],
        df['major_cell_type'],
        'Original Coarse vs Agent Major Cell Types (Normalized)',
        output_dir / "confusion_matrix_coarse.png"
    )

    plot_confusion_matrix(
        df['original_sbcl_ct'],
        df['cell_type_refined'],
        'Original Fine vs Agent Refined Cell Types',
        output_dir / "confusion_matrix_fine.png"
    )

    # Plot purity comparison
    print("\n5. Generating purity comparison plots...")
    plot_purity_comparison(cluster_agreement, output_dir / "purity_comparison.png")

    # Create summary report
    print("\n6. Creating summary report...")
    report_path = output_dir / "EVALUATION_REPORT.md"

    with open(report_path, 'w') as f:
        f.write("# Cell Type Annotation Performance Evaluation\n\n")
        f.write("## Dataset Information\n\n")
        f.write(f"- Total cells: {len(df):,}\n")
        f.write(f"- Number of clusters: {df['louvain'].nunique()}\n\n")

        f.write("## Annotation Methods Compared\n\n")
        f.write("**Original Method:**\n")
        f.write(f"- Coarse level: {df['original_cl_ct_normalized'].nunique()} cell types (normalized from {df['original_cl_ct'].nunique()} abbreviated forms)\n")
        f.write(f"- Fine level: {df['original_sbcl_ct'].nunique()} cell types\n\n")

        f.write("**Agent-based Method (scAnnoist):**\n")
        f.write(f"- Major cell types: {df['major_cell_type'].nunique()} types\n")
        f.write(f"- Refined cell types: {df['cell_type_refined'].nunique()} types\n\n")

        f.write("## Performance Metrics\n\n")
        f.write("| " + " | ".join(metrics_df.columns) + " |\n")
        f.write("| " + " | ".join(["---"] * len(metrics_df.columns)) + " |\n")
        for _, row in metrics_df.iterrows():
            f.write("| " + " | ".join([str(v) for v in row.values]) + " |\n")
        f.write("\n")

        f.write("### Metric Interpretation\n\n")
        f.write("- **ARI (Adjusted Rand Index)**: Measures similarity between two clusterings (1.0 = perfect match, 0.0 = random)\n")
        f.write("- **NMI (Normalized Mutual Information)**: Measures mutual dependence (1.0 = perfect match)\n")
        f.write("- **Homogeneity**: Whether all clusters contain only members of a single class\n")
        f.write("- **Completeness**: Whether all members of a class are in the same cluster\n")
        f.write("- **V-measure**: Harmonic mean of homogeneity and completeness\n\n")

        f.write("## Cluster-level Purity\n\n")
        purity_summary = cluster_agreement[['Original_Coarse_Purity', 'Original_Fine_Purity',
                                           'Agent_Major_Purity', 'Agent_Refined_Purity']].describe()
        f.write("| Statistic | " + " | ".join(purity_summary.columns) + " |\n")
        f.write("| --- | " + " | ".join(["---"] * len(purity_summary.columns)) + " |\n")
        for idx, row in purity_summary.iterrows():
            f.write(f"| {idx} | " + " | ".join([f"{v:.4f}" for v in row.values]) + " |\n")
        f.write("\n")

        f.write("## Top 10 Largest Clusters\n\n")
        top_clusters = cluster_agreement.nlargest(10, 'N_cells')[
            ['Cluster', 'N_cells', 'Original_Coarse', 'Agent_Major',
             'Original_Coarse_Purity', 'Agent_Major_Purity']
        ]
        f.write("| " + " | ".join(top_clusters.columns) + " |\n")
        f.write("| " + " | ".join(["---"] * len(top_clusters.columns)) + " |\n")
        for _, row in top_clusters.iterrows():
            f.write("| " + " | ".join([str(v) for v in row.values]) + " |\n")
        f.write("\n")

        f.write("## Files Generated\n\n")
        f.write("- `annotation_metrics.csv`: Detailed performance metrics\n")
        f.write("- `cluster_level_agreement.csv`: Per-cluster annotation comparison\n")
        f.write("- `confusion_matrix_coarse.png`: Coarse-level confusion matrix\n")
        f.write("- `confusion_matrix_fine.png`: Fine-level confusion matrix\n")
        f.write("- `purity_comparison.png`: Cluster purity comparison plots\n")
        f.write("- `EVALUATION_REPORT.md`: This report\n")

    print(f"Saved evaluation report: {report_path}")

    print("\n" + "="*80)
    print("Evaluation Complete!")
    print("="*80)
    print(f"\nResults saved to: {output_dir}")


if __name__ == "__main__":
    main()
