"""
Annotation Integration Module
Compare and integrate multiple annotation methods (CellTypist, Literature-based, ScType)
to produce optimal cell type annotations
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Any
from collections import Counter


def compare_annotation_methods(
    adata,
    methods: List[str],
    cluster_key: str = "leiden"
) -> Dict[str, Any]:
    """
    Compare multiple annotation methods.

    Args:
        adata: AnnData object with multiple annotation columns
        methods: List of annotation column names to compare
        cluster_key: Cluster column for cluster-level comparison

    Returns:
        Dictionary with comparison statistics
    """
    comparison = {
        "methods": methods,
        "available": [],
        "missing": [],
        "agreement": {},
        "disagreement": {},
        "cluster_consistency": {}
    }

    # Check which methods are available
    for method in methods:
        if method in adata.obs.columns:
            comparison["available"].append(method)
        else:
            comparison["missing"].append(method)

    if len(comparison["available"]) < 2:
        return comparison

    # Pairwise agreement
    available = comparison["available"]
    for i, method1 in enumerate(available):
        for method2 in available[i+1:]:
            pair_key = f"{method1} vs {method2}"
            agreement_rate = (adata.obs[method1] == adata.obs[method2]).mean()
            comparison["agreement"][pair_key] = agreement_rate

            # Find disagreement cells
            disagreement_mask = adata.obs[method1] != adata.obs[method2]
            disagreement_cells = adata.obs.index[disagreement_mask].tolist()
            comparison["disagreement"][pair_key] = {
                "count": disagreement_mask.sum(),
                "percentage": (disagreement_mask.sum() / len(adata)) * 100,
                "cells": disagreement_cells[:100]  # First 100 for summary
            }

    # Cluster-level consistency
    if cluster_key in adata.obs.columns:
        for cluster in adata.obs[cluster_key].unique():
            cluster_mask = adata.obs[cluster_key] == cluster
            cluster_data = adata.obs[cluster_mask]

            cluster_comparison = {}
            for method in available:
                # Most common annotation in this cluster
                mode_annotation = cluster_data[method].mode()
                if len(mode_annotation) > 0:
                    mode_annotation = mode_annotation.iloc[0]
                    mode_count = (cluster_data[method] == mode_annotation).sum()
                    mode_percentage = (mode_count / len(cluster_data)) * 100

                    cluster_comparison[method] = {
                        "annotation": mode_annotation,
                        "consistency": mode_percentage,
                        "count": mode_count
                    }

            comparison["cluster_consistency"][str(cluster)] = cluster_comparison

    return comparison


def calculate_annotation_confidence(
    adata,
    annotation_col: str,
    confidence_col: Optional[str] = None,
    cluster_key: str = "leiden"
) -> pd.Series:
    """
    Calculate confidence scores for annotations.

    Args:
        adata: AnnData object
        annotation_col: Annotation column name
        confidence_col: Optional confidence column (e.g., from CellTypist)
        cluster_key: Cluster column for cluster-based confidence

    Returns:
        Series of confidence scores per cell
    """
    if confidence_col and confidence_col in adata.obs.columns:
        # Use provided confidence scores
        return adata.obs[confidence_col]

    # Calculate confidence based on cluster consistency
    confidence_scores = np.zeros(len(adata))

    if cluster_key in adata.obs.columns:
        for cluster in adata.obs[cluster_key].unique():
            cluster_mask = adata.obs[cluster_key] == cluster
            cluster_annotations = adata.obs.loc[cluster_mask, annotation_col]

            # Calculate consistency as confidence
            for annotation in cluster_annotations.unique():
                annotation_mask = cluster_annotations == annotation
                consistency = annotation_mask.sum() / len(cluster_annotations)

                # Assign confidence to cells with this annotation
                cell_indices = cluster_annotations[annotation_mask].index
                idx_positions = [adata.obs.index.get_loc(idx) for idx in cell_indices]
                confidence_scores[idx_positions] = consistency

    return pd.Series(confidence_scores, index=adata.obs.index)


def integrate_annotations_voting(
    adata,
    methods: List[str],
    weights: Optional[Dict[str, float]] = None,
    confidence_cols: Optional[Dict[str, str]] = None,
    output_col: str = "integrated_annotation"
) -> Tuple[pd.Series, pd.Series]:
    """
    Integrate multiple annotations using weighted voting.

    Args:
        adata: AnnData object
        methods: List of annotation column names
        weights: Optional weights for each method (default: equal weights)
        confidence_cols: Optional confidence columns for each method
        output_col: Name for output column

    Returns:
        Tuple of (integrated annotations, confidence scores)
    """
    # Default equal weights
    if weights is None:
        weights = {method: 1.0 for method in methods}

    # Normalize weights
    total_weight = sum(weights.values())
    weights = {k: v/total_weight for k, v in weights.items()}

    integrated_annotations = []
    integrated_confidences = []

    for cell_idx in range(len(adata)):
        votes = {}
        total_confidence = {}

        for method in methods:
            if method not in adata.obs.columns:
                continue

            annotation = adata.obs[method].iloc[cell_idx]

            # Get confidence for this annotation
            if confidence_cols and method in confidence_cols:
                conf_col = confidence_cols[method]
                if conf_col in adata.obs.columns:
                    confidence = adata.obs[conf_col].iloc[cell_idx]
                else:
                    confidence = 1.0
            else:
                confidence = 1.0

            # Weighted vote
            vote_weight = weights[method] * confidence

            if annotation not in votes:
                votes[annotation] = 0
                total_confidence[annotation] = 0

            votes[annotation] += vote_weight
            total_confidence[annotation] += confidence

        # Select annotation with highest weighted vote
        if votes:
            best_annotation = max(votes, key=votes.get)
            best_confidence = total_confidence[best_annotation] / len(methods)
        else:
            best_annotation = "Unknown"
            best_confidence = 0.0

        integrated_annotations.append(best_annotation)
        integrated_confidences.append(best_confidence)

    return (
        pd.Series(integrated_annotations, index=adata.obs.index),
        pd.Series(integrated_confidences, index=adata.obs.index)
    )


def resolve_conflicts(
    adata,
    primary_method: str,
    secondary_method: str,
    primary_confidence_col: Optional[str] = None,
    secondary_confidence_col: Optional[str] = None,
    confidence_threshold: float = 0.5,
    output_col: str = "resolved_annotation"
) -> pd.Series:
    """
    Resolve conflicts between two annotation methods using confidence scores.

    Args:
        adata: AnnData object
        primary_method: Primary annotation method (e.g., CellTypist)
        secondary_method: Secondary annotation method (e.g., literature)
        primary_confidence_col: Confidence column for primary method
        secondary_confidence_col: Confidence column for secondary method
        confidence_threshold: Threshold for using primary method
        output_col: Name for output column

    Returns:
        Series of resolved annotations
    """
    resolved = []

    for cell_idx in range(len(adata)):
        primary_annotation = adata.obs[primary_method].iloc[cell_idx]
        secondary_annotation = adata.obs[secondary_method].iloc[cell_idx]

        # Get confidence scores
        if primary_confidence_col and primary_confidence_col in adata.obs.columns:
            primary_conf = adata.obs[primary_confidence_col].iloc[cell_idx]
        else:
            primary_conf = 1.0

        if secondary_confidence_col and secondary_confidence_col in adata.obs.columns:
            secondary_conf = adata.obs[secondary_confidence_col].iloc[cell_idx]
        else:
            secondary_conf = 1.0

        # Resolution logic
        if primary_annotation == secondary_annotation:
            # Agreement: use either
            resolved.append(primary_annotation)
        elif primary_conf >= confidence_threshold:
            # High confidence in primary: use primary
            resolved.append(primary_annotation)
        elif secondary_conf > primary_conf:
            # Secondary has higher confidence: use secondary
            resolved.append(secondary_annotation)
        else:
            # Default to primary
            resolved.append(primary_annotation)

    return pd.Series(resolved, index=adata.obs.index)


def validate_integrated_annotation(
    adata,
    annotation_col: str,
    marker_dict: Dict[str, List[str]],
    min_expression_threshold: float = 0.5
) -> Dict[str, Any]:
    """
    Validate integrated annotations against known markers.

    Args:
        adata: AnnData object
        annotation_col: Annotation column to validate
        marker_dict: Dictionary of cell_type -> marker genes
        min_expression_threshold: Minimum expression level for validation

    Returns:
        Dictionary with validation results
    """
    validation_results = {}

    for cell_type, markers in marker_dict.items():
        # Find cells with this annotation
        type_mask = adata.obs[annotation_col] == cell_type
        if type_mask.sum() == 0:
            continue

        # Check marker expression
        present_markers = [m for m in markers if m in adata.var_names]
        if not present_markers:
            validation_results[cell_type] = {
                "status": "no_markers",
                "n_cells": type_mask.sum()
            }
            continue

        # Get expression for these cells
        type_cells = adata[type_mask, present_markers]

        if hasattr(type_cells.X, 'toarray'):
            expr = type_cells.X.toarray()
        else:
            expr = type_cells.X

        # Calculate mean expression per marker
        mean_expr = np.mean(expr, axis=0)

        # Validation: how many markers are expressed above threshold
        expressed_markers = sum(mean_expr > min_expression_threshold)
        validation_score = expressed_markers / len(present_markers)

        validation_results[cell_type] = {
            "status": "validated" if validation_score > 0.5 else "questionable",
            "n_cells": type_mask.sum(),
            "n_markers": len(present_markers),
            "expressed_markers": expressed_markers,
            "validation_score": validation_score,
            "marker_expression": {
                present_markers[i]: float(mean_expr[i])
                for i in range(len(present_markers))
            }
        }

    return validation_results


def create_integrated_annotation_report(
    adata,
    methods: List[str],
    integrated_col: str,
    output_file: Optional[str] = None
) -> str:
    """
    Create a detailed report of annotation integration.

    Args:
        adata: AnnData object
        methods: List of methods that were integrated
        integrated_col: Column with integrated annotations
        output_file: Optional file path to save report

    Returns:
        Report string
    """
    report_lines = []

    report_lines.append("=" * 70)
    report_lines.append("ANNOTATION INTEGRATION REPORT")
    report_lines.append("=" * 70)

    # Methods used
    report_lines.append(f"\nMethods integrated: {', '.join(methods)}")

    # Available methods
    available = [m for m in methods if m in adata.obs.columns]
    report_lines.append(f"Available methods: {', '.join(available)}")

    # Integrated annotation summary
    report_lines.append(f"\nIntegrated annotation: {integrated_col}")
    integrated_counts = adata.obs[integrated_col].value_counts()
    report_lines.append(f"Total cell types identified: {len(integrated_counts)}")

    report_lines.append("\nTop 15 cell types:")
    for cell_type, count in integrated_counts.head(15).items():
        percentage = (count / len(adata)) * 100
        report_lines.append(f"  {cell_type}: {count:,} ({percentage:.1f}%)")

    # Agreement between methods
    if len(available) >= 2:
        report_lines.append("\nPairwise Agreement:")
        for i, method1 in enumerate(available):
            for method2 in available[i+1:]:
                agreement = (adata.obs[method1] == adata.obs[method2]).mean()
                report_lines.append(f"  {method1} vs {method2}: {agreement:.1%}")

    # Confidence distribution (if available)
    if f"{integrated_col}_confidence" in adata.obs.columns:
        conf_col = f"{integrated_col}_confidence"
        report_lines.append(f"\nConfidence Statistics:")
        report_lines.append(f"  Mean: {adata.obs[conf_col].mean():.3f}")
        report_lines.append(f"  Median: {adata.obs[conf_col].median():.3f}")
        report_lines.append(f"  Low confidence (<0.5): {(adata.obs[conf_col] < 0.5).sum():,} cells")

    report_lines.append("\n" + "=" * 70)

    report = "\n".join(report_lines)

    # Save to file if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report)

    return report
