"""
Marker-based Cell Type Annotation
Simple implementation inspired by ScType
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
import scanpy as sc


# Marker gene database for major tissues
MARKER_DATABASE = {
    "Lung": {
        "AT1 cells": ["AGER", "PDPN", "HOPX", "CAV1", "CLIC5"],
        "AT2 cells": ["SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "LAMP3"],
        "Club cells": ["SCGB1A1", "SCGB3A2", "CYP2F1", "LYPD2"],
        "Ciliated cells": ["FOXJ1", "PIFO", "SNTN", "TUBA1A", "DNAI2"],
        "Basal cells": ["TP63", "KRT5", "KRT14", "KRT15", "NGFR"],
        "Goblet cells": ["MUC5AC", "MUC5B", "TFF3", "SPDEF"],
        "Endothelial cells": ["PECAM1", "CDH5", "VWF", "ENG", "CD34"],
        "Fibroblasts": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"],
        "Smooth muscle": ["ACTA2", "TAGLN", "MYH11", "CNN1"],
        "Macrophages": ["CD68", "CD14", "FCGR3A", "CSF1R", "CD163"],
        "T cells": ["CD3D", "CD3E", "CD3G", "TRAC", "IL7R"],
        "B cells": ["CD19", "MS4A1", "CD79A", "CD79B"],
        "NK cells": ["NKG7", "GNLY", "NCAM1", "KLRD1"],
        "Dendritic cells": ["CD1C", "CLEC9A", "FCER1A", "ITGAX"],
        "Mast cells": ["TPSAB1", "TPSB2", "CPA3", "KIT"],
        "Neutrophils": ["S100A8", "S100A9", "CSF3R", "FCGR3B"],
    },
    "Blood": {
        "CD14 Monocytes": ["CD14", "LYZ", "S100A8", "S100A9", "FCGR3A"],
        "CD16 Monocytes": ["FCGR3A", "MS4A7", "CD14", "LYZ"],
        "B cells": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5"],
        "Naive B cells": ["CD19", "MS4A1", "IGHD", "IGHM", "TCL1A"],
        "Memory B cells": ["CD19", "MS4A1", "CD27", "IGHG1"],
        "Plasma cells": ["JCHAIN", "IGKC", "IGLC2", "MZB1", "SDC1"],
        "CD4 T cells": ["CD3D", "CD3E", "CD4", "IL7R", "CCR7"],
        "CD8 T cells": ["CD3D", "CD3E", "CD8A", "CD8B"],
        "NK cells": ["NKG7", "GNLY", "NCAM1", "KLRD1", "GZMB"],
        "Dendritic cells": ["FCER1A", "CD1C", "CLEC10A"],
        "pDCs": ["LILRA4", "IL3RA", "CLEC4C", "IRF7"],
        "Platelets": ["PPBP", "PF4", "GP9", "ITGA2B"],
    },
    "PBMC": {
        "CD14_mono": ["CD14", "LYZ", "S100A8", "S100A9"],
        "CD4.Th17": ["CD4", "IL17A", "IL17F", "RORC", "CCR6"],
        "CD4.Naive": ["CD4", "CCR7", "SELL", "IL7R", "TCF7"],
        "B_immature": ["CD19", "MS4A1", "CD24", "CD38"],
        "B_switched_memory": ["CD19", "MS4A1", "CD27", "IGHG1"],
        "B_malignant": ["CD19", "MS4A1", "MKI67", "TOP2A"],
        "HSC_CD38neg": ["CD34", "KIT", "THY1", "CD38"],
        "NKT": ["CD3D", "CD3E", "NCAM1", "NKG7"],
    }
}


def run_marker_annotation(
    adata,
    tissue: str = "Lung",
    cluster_key: str = "leiden",
    min_markers: int = 2,
    output_col: str = "marker_annotation"
) -> Dict:
    """
    Run marker-based cell type annotation.

    Args:
        adata: AnnData object (should be log-normalized)
        tissue: Tissue type for marker database
        cluster_key: Column containing cluster labels
        min_markers: Minimum markers required for assignment
        output_col: Column name for output annotations

    Returns:
        dict with annotation results
    """
    print(f"\n{'='*60}")
    print(f"Marker-Based Annotation ({tissue})")
    print(f"{'='*60}")

    # Get markers for tissue
    if tissue not in MARKER_DATABASE:
        available = list(MARKER_DATABASE.keys())
        raise ValueError(f"Tissue '{tissue}' not in database. Available: {available}")

    markers = MARKER_DATABASE[tissue]

    # Check cluster key
    if cluster_key not in adata.obs:
        raise ValueError(f"Cluster column '{cluster_key}' not found")

    print(f"\nDatabase: {len(markers)} cell types")
    print(f"Clusters: {adata.obs[cluster_key].nunique()}")

    # Score each cluster for each cell type
    cluster_annotations = {}
    cluster_scores_all = {}

    for cluster in sorted(adata.obs[cluster_key].unique()):
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_data = adata[cluster_mask].copy()

        # Make obs names unique to avoid indexing issues
        cluster_data.obs_names_make_unique()

        scores = {}
        for cell_type, marker_genes in markers.items():
            # Find markers present in data (case-insensitive)
            var_names_upper = [g.upper() for g in adata.var_names]
            present_markers = []
            for marker in marker_genes:
                marker_upper = marker.upper()
                if marker_upper in var_names_upper:
                    idx = var_names_upper.index(marker_upper)
                    present_markers.append(adata.var_names[idx])

            if len(present_markers) < min_markers:
                scores[cell_type] = 0.0
                continue

            # Calculate mean expression
            if hasattr(cluster_data.X, 'toarray'):
                expr = cluster_data[:, present_markers].X.toarray()
            else:
                expr = cluster_data[:, present_markers].X

            # Mean expression score
            mean_expr = np.mean(expr)
            # Fraction of cells expressing (log-normalized > 0)
            frac_expr = np.mean(expr > 0)

            # Combined score: mean expression * fraction expressing
            scores[cell_type] = mean_expr * frac_expr

        # Get best cell type
        if scores:
            best_type = max(scores, key=scores.get)
            best_score = scores[best_type]
            cluster_annotations[cluster] = best_type
            cluster_scores_all[cluster] = {
                "annotation": best_type,
                "score": best_score,
                "all_scores": scores
            }
        else:
            cluster_annotations[cluster] = "Unknown"
            cluster_scores_all[cluster] = {
                "annotation": "Unknown",
                "score": 0.0,
                "all_scores": {}
            }

    # Add to adata
    adata.obs[output_col] = adata.obs[cluster_key].map(cluster_annotations).astype(str)

    # Summary
    print(f"\n{'='*60}")
    print("Annotation Results:")
    print(f"{'='*60}")

    for cluster in sorted(cluster_annotations.keys()):
        ann = cluster_scores_all[cluster]
        n_cells = (adata.obs[cluster_key] == cluster).sum()
        print(f"Cluster {cluster}: {ann['annotation']} (score={ann['score']:.3f}, n={n_cells})")

    # Cell type distribution
    print(f"\n{'='*60}")
    print("Cell Type Distribution:")
    print(f"{'='*60}")
    type_counts = adata.obs[output_col].value_counts()
    for cell_type, count in type_counts.items():
        pct = (count / len(adata)) * 100
        print(f"  {cell_type}: {count:,} ({pct:.1f}%)")

    return {
        "cluster_annotations": cluster_annotations,
        "cluster_scores": cluster_scores_all,
        "cell_type_counts": type_counts.to_dict(),
        "total_cell_types": len(type_counts)
    }


def compare_annotations(
    adata,
    method1_col: str,
    method2_col: str,
    cluster_key: str = "leiden"
) -> pd.DataFrame:
    """
    Compare two annotation methods.

    Args:
        adata: AnnData object
        method1_col: First annotation column
        method2_col: Second annotation column
        cluster_key: Cluster column for grouping

    Returns:
        DataFrame with comparison results
    """
    print(f"\n{'='*60}")
    print(f"Annotation Comparison")
    print(f"{'='*60}")
    print(f"Method 1: {method1_col}")
    print(f"Method 2: {method2_col}")

    # Overall agreement
    agreement = (adata.obs[method1_col] == adata.obs[method2_col]).mean()
    print(f"\nOverall agreement: {agreement*100:.1f}%")

    # Cluster-level comparison
    comparison_data = []
    for cluster in sorted(adata.obs[cluster_key].unique()):
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_obs = adata.obs[cluster_mask]

        method1_ann = cluster_obs[method1_col].mode()[0] if len(cluster_obs[method1_col].mode()) > 0 else "Unknown"
        method2_ann = cluster_obs[method2_col].mode()[0] if len(cluster_obs[method2_col].mode()) > 0 else "Unknown"

        method1_pct = (cluster_obs[method1_col] == method1_ann).mean() * 100
        method2_pct = (cluster_obs[method2_col] == method2_ann).mean() * 100

        agree = method1_ann == method2_ann

        comparison_data.append({
            "cluster": cluster,
            "n_cells": cluster_mask.sum(),
            f"{method1_col}_annotation": method1_ann,
            f"{method1_col}_consistency": f"{method1_pct:.1f}%",
            f"{method2_col}_annotation": method2_ann,
            f"{method2_col}_consistency": f"{method2_pct:.1f}%",
            "agreement": "✓" if agree else "✗"
        })

    comparison_df = pd.DataFrame(comparison_data)

    print(f"\n{'='*60}")
    print("Cluster-Level Comparison:")
    print(f"{'='*60}")
    print(comparison_df.to_string(index=False))

    return comparison_df
