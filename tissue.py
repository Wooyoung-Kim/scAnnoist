"""
Tissue handling and marker-based cell type annotation module.

Provides:
- Tissue metadata loading and validation
- Tissue-aware marker-based cell type annotation
- Summary generation
"""

import os
import logging
from pathlib import Path
from typing import Optional, Dict, List, Any, Union
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import yaml
import scanpy as sc

logger = logging.getLogger(__name__)

# Default markers file path
DEFAULT_MARKERS_PATH = Path(__file__).parent / "data" / "markers" / "default_markers.yaml"


@dataclass
class TissueConfig:
    """Configuration for tissue handling."""
    tissue: Optional[str] = None  # Single tissue for all cells
    tissue_map_path: Optional[str] = None  # Path to tissue-map CSV
    sample_col: str = "sample"  # Column name for sample IDs
    allow_missing_tissue: bool = False  # Allow samples without tissue assignment


@dataclass
class AnnotationConfig:
    """Configuration for marker-based annotation."""
    markers_path: Optional[str] = None  # Path to markers YAML
    min_markers_expressed: int = 2  # Minimum markers that must be expressed
    min_expression_threshold: float = 0.5  # Minimum mean expression
    score_method: str = "mean"  # "mean" or "sum"


class TissueMapError(Exception):
    """Error in tissue map processing."""
    pass


def load_tissue_map(
    tissue_map_path: str,
    sample_col: str = "sample"
) -> Dict[str, str]:
    """
    Load tissue map from CSV file.

    Args:
        tissue_map_path: Path to CSV with sample,tissue columns
        sample_col: Name of sample column (default: "sample")

    Returns:
        Dictionary mapping sample -> tissue

    Raises:
        TissueMapError: If file format is invalid
    """
    logger.info(f"Loading tissue map from {tissue_map_path}")

    df = pd.read_csv(tissue_map_path)

    # Check required columns
    if sample_col not in df.columns:
        raise TissueMapError(
            f"Tissue map must have '{sample_col}' column. "
            f"Found columns: {list(df.columns)}"
        )
    if "tissue" not in df.columns:
        raise TissueMapError(
            f"Tissue map must have 'tissue' column. "
            f"Found columns: {list(df.columns)}"
        )

    # Build mapping
    tissue_map = dict(zip(df[sample_col].astype(str), df["tissue"].astype(str)))

    logger.info(f"Loaded tissue map with {len(tissue_map)} entries")
    return tissue_map


def assign_tissue(
    adata: sc.AnnData,
    tissue: Optional[str] = None,
    tissue_map: Optional[Dict[str, str]] = None,
    sample_col: str = "sample",
    allow_missing: bool = False
) -> sc.AnnData:
    """
    Assign tissue metadata to cells in AnnData.

    Args:
        adata: AnnData object
        tissue: Single tissue value for all cells (used if tissue_map is None)
        tissue_map: Dictionary mapping sample -> tissue
        sample_col: Column containing sample IDs
        allow_missing: If True, fill missing tissues with "unknown"

    Returns:
        AnnData with tissue column added to obs

    Raises:
        TissueMapError: If validation fails
    """
    if tissue_map is not None:
        # Use tissue map
        if sample_col not in adata.obs.columns:
            raise TissueMapError(
                f"Sample column '{sample_col}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )

        samples = adata.obs[sample_col].astype(str).unique()
        missing_samples = [s for s in samples if s not in tissue_map]
        unknown_samples = [s for s in tissue_map if s not in samples]

        if unknown_samples:
            logger.warning(
                f"Tissue map contains {len(unknown_samples)} samples not in data: "
                f"{unknown_samples[:5]}{'...' if len(unknown_samples) > 5 else ''}"
            )

        if missing_samples:
            if allow_missing:
                logger.warning(
                    f"Filling {len(missing_samples)} samples with 'unknown' tissue: "
                    f"{missing_samples[:5]}{'...' if len(missing_samples) > 5 else ''}"
                )
                for s in missing_samples:
                    tissue_map[s] = "unknown"
            else:
                raise TissueMapError(
                    f"Samples missing from tissue map: {missing_samples}. "
                    "Use --allow-missing-tissue to fill with 'unknown'."
                )

        # Assign tissue to each cell
        adata.obs["tissue"] = adata.obs[sample_col].astype(str).map(tissue_map)

    elif tissue is not None:
        # Single tissue for all cells
        adata.obs["tissue"] = tissue
        logger.info(f"Assigned tissue '{tissue}' to all {adata.n_obs} cells")

    else:
        # No tissue specified
        adata.obs["tissue"] = "unknown"
        logger.warning("No tissue specified, using 'unknown' for all cells")

    return adata


def load_markers(markers_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load markers from YAML file.

    Args:
        markers_path: Path to markers YAML file (uses default if None)

    Returns:
        Dictionary with markers and settings
    """
    if markers_path is None:
        markers_path = DEFAULT_MARKERS_PATH

    logger.info(f"Loading markers from {markers_path}")

    with open(markers_path, 'r') as f:
        data = yaml.safe_load(f)

    return data


def get_tissue_markers(
    markers_data: Dict[str, Any],
    tissue: str
) -> Dict[str, List[str]]:
    """
    Get markers for a specific tissue, falling back to default.

    Args:
        markers_data: Full markers dictionary
        tissue: Tissue name

    Returns:
        Dictionary of cell_type -> marker_genes
    """
    markers = markers_data.get("markers", {})

    # Normalize tissue name
    tissue_lower = tissue.lower().replace(" ", "_")

    # Check for tissue-specific markers
    tissue_markers = markers.get(tissue_lower, {})
    default_markers = markers.get("default", {})

    # Merge: tissue-specific takes precedence
    combined = {**default_markers, **tissue_markers}

    if tissue_markers:
        logger.info(
            f"Using tissue-specific markers for '{tissue}' "
            f"({len(tissue_markers)} cell types) + default ({len(default_markers)} cell types)"
        )
    else:
        logger.info(
            f"No tissue-specific markers for '{tissue}', using default markers "
            f"({len(default_markers)} cell types)"
        )

    return combined


def score_cell_types(
    adata: sc.AnnData,
    markers: Dict[str, List[str]],
    min_markers: int = 2,
    min_expression: float = 0.5,
    score_method: str = "mean"
) -> pd.DataFrame:
    """
    Score cells for each cell type based on marker expression.

    Args:
        adata: AnnData with normalized expression
        markers: Dictionary of cell_type -> marker_genes
        min_markers: Minimum markers that must be expressed
        min_expression: Minimum mean expression threshold
        score_method: "mean" or "sum"

    Returns:
        DataFrame with cell type scores (cells x cell_types)
    """
    # Get expression matrix (use raw if available for normalized values)
    if adata.raw is not None:
        expr_data = adata.raw.to_adata()
    else:
        expr_data = adata

    scores = {}

    for cell_type, marker_genes in markers.items():
        # Filter to genes present in data
        available_genes = [g for g in marker_genes if g in expr_data.var_names]

        if len(available_genes) < min_markers:
            logger.debug(
                f"Skipping {cell_type}: only {len(available_genes)}/{len(marker_genes)} "
                f"markers available (min: {min_markers})"
            )
            continue

        # Get expression for marker genes
        marker_expr = expr_data[:, available_genes].X

        # Handle sparse matrices
        if hasattr(marker_expr, 'toarray'):
            marker_expr = marker_expr.toarray()

        # Calculate score
        if score_method == "sum":
            score = np.sum(marker_expr, axis=1)
        else:  # mean
            score = np.mean(marker_expr, axis=1)

        # Count expressed markers per cell
        n_expressed = np.sum(marker_expr > min_expression, axis=1)

        # Only score cells with minimum markers expressed
        score[n_expressed < min_markers] = 0

        scores[cell_type] = score.flatten()

    return pd.DataFrame(scores, index=adata.obs_names)


def annotate_cells(
    adata: sc.AnnData,
    markers_path: Optional[str] = None,
    min_markers: int = 2,
    min_expression: float = 0.5,
    score_method: str = "mean"
) -> sc.AnnData:
    """
    Annotate cells using tissue-aware marker-based approach.

    Args:
        adata: AnnData with tissue column in obs
        markers_path: Path to markers YAML
        min_markers: Minimum markers that must be expressed
        min_expression: Minimum mean expression threshold
        score_method: "mean" or "sum"

    Returns:
        AnnData with celltype and celltype_score columns
    """
    if "tissue" not in adata.obs.columns:
        logger.warning("No tissue column found, using 'unknown'")
        adata.obs["tissue"] = "unknown"

    markers_data = load_markers(markers_path)
    settings = markers_data.get("settings", {})

    # Use settings from YAML if not overridden
    min_markers = settings.get("min_markers_expressed", min_markers)
    min_expression = settings.get("min_expression_threshold", min_expression)
    score_method = settings.get("score_method", score_method)

    # Get unique tissues
    tissues = adata.obs["tissue"].unique()
    logger.info(f"Annotating cells for {len(tissues)} tissue(s): {list(tissues)}")

    # Initialize annotation columns
    adata.obs["celltype"] = "Unknown"
    adata.obs["celltype_score"] = 0.0

    # Annotate each tissue separately
    for tissue in tissues:
        mask = adata.obs["tissue"] == tissue
        adata_tissue = adata[mask].copy()

        logger.info(f"Annotating {mask.sum()} cells for tissue '{tissue}'")

        # Get markers for this tissue
        tissue_markers = get_tissue_markers(markers_data, tissue)

        if not tissue_markers:
            logger.warning(f"No markers available for tissue '{tissue}'")
            continue

        # Score cells
        scores_df = score_cell_types(
            adata_tissue,
            tissue_markers,
            min_markers=min_markers,
            min_expression=min_expression,
            score_method=score_method
        )

        if scores_df.empty:
            logger.warning(f"No valid scores for tissue '{tissue}'")
            continue

        # Assign best cell type
        best_celltype = scores_df.idxmax(axis=1)
        best_score = scores_df.max(axis=1)

        # Mark low-scoring cells as Unknown
        best_celltype[best_score < min_expression] = "Unknown"

        # Update main adata
        adata.obs.loc[mask, "celltype"] = best_celltype.values
        adata.obs.loc[mask, "celltype_score"] = best_score.values

    # Log summary
    celltype_counts = adata.obs["celltype"].value_counts()
    logger.info(f"Annotation complete. Cell types found: {len(celltype_counts)}")
    logger.info(f"Top cell types:\n{celltype_counts.head(10)}")

    return adata


def generate_tissue_summary(adata: sc.AnnData) -> pd.DataFrame:
    """
    Generate summary of cell types per tissue.

    Args:
        adata: Annotated AnnData

    Returns:
        DataFrame with tissue x celltype counts
    """
    if "tissue" not in adata.obs.columns or "celltype" not in adata.obs.columns:
        raise ValueError("AnnData must have 'tissue' and 'celltype' columns")

    # Cross-tabulation
    summary = pd.crosstab(
        adata.obs["tissue"],
        adata.obs["celltype"],
        margins=True,
        margins_name="Total"
    )

    return summary


def generate_cell_metadata(adata: sc.AnnData, sample_col: str = "sample") -> pd.DataFrame:
    """
    Generate cell metadata DataFrame for export.

    Args:
        adata: Annotated AnnData
        sample_col: Sample column name

    Returns:
        DataFrame with cell metadata
    """
    cols = ["tissue", "celltype", "celltype_score"]

    # Add sample column if present
    if sample_col in adata.obs.columns:
        cols = [sample_col] + cols

    # Add cluster if present
    for cluster_col in ["leiden", "louvain", "cluster"]:
        if cluster_col in adata.obs.columns:
            cols.append(cluster_col)
            break

    # Select available columns
    available_cols = [c for c in cols if c in adata.obs.columns]

    metadata = adata.obs[available_cols].copy()
    metadata.index.name = "cell_id"

    return metadata


def generate_tissue_report_section(adata: sc.AnnData) -> str:
    """
    Generate report section for tissue-aware annotation.

    Args:
        adata: Annotated AnnData

    Returns:
        Markdown string for report
    """
    lines = [
        "## Tissue-Aware Annotation",
        "",
    ]

    # Tissue distribution
    if "tissue" in adata.obs.columns:
        tissue_counts = adata.obs["tissue"].value_counts()
        lines.extend([
            "### Tissue Distribution",
            "",
            "| Tissue | Cells | Percentage |",
            "|--------|-------|------------|",
        ])
        for tissue, count in tissue_counts.items():
            pct = count / adata.n_obs * 100
            lines.append(f"| {tissue} | {count:,} | {pct:.1f}% |")
        lines.append("")

    # Cell type distribution
    if "celltype" in adata.obs.columns:
        celltype_counts = adata.obs["celltype"].value_counts()
        lines.extend([
            "### Cell Type Distribution",
            "",
            "| Cell Type | Cells | Percentage |",
            "|-----------|-------|------------|",
        ])
        for celltype, count in celltype_counts.head(15).items():
            pct = count / adata.n_obs * 100
            lines.append(f"| {celltype} | {count:,} | {pct:.1f}% |")
        if len(celltype_counts) > 15:
            lines.append(f"| ... | ... | ... |")
        lines.append("")

    # Cell type by tissue
    if "tissue" in adata.obs.columns and "celltype" in adata.obs.columns:
        lines.extend([
            "### Cell Types by Tissue",
            "",
        ])
        for tissue in adata.obs["tissue"].unique():
            mask = adata.obs["tissue"] == tissue
            tissue_celltypes = adata.obs.loc[mask, "celltype"].value_counts().head(5)
            lines.append(f"**{tissue}** ({mask.sum()} cells):")
            for ct, count in tissue_celltypes.items():
                lines.append(f"  - {ct}: {count}")
            lines.append("")

    # Annotation quality
    if "celltype_score" in adata.obs.columns:
        scores = adata.obs["celltype_score"]
        lines.extend([
            "### Annotation Quality",
            "",
            f"- Mean score: {scores.mean():.3f}",
            f"- Median score: {scores.median():.3f}",
            f"- Cells with score > 1.0: {(scores > 1.0).sum()} ({(scores > 1.0).mean()*100:.1f}%)",
            f"- Unknown cells: {(adata.obs['celltype'] == 'Unknown').sum()}",
            "",
        ])

    return "\n".join(lines)


# =============================================================================
# Multi-Method Annotation with Biological Validation
# =============================================================================

def annotate_with_validation(
    adata: sc.AnnData,
    sample_type: str = 'primary_tissue',
    markers_dict: Optional[Dict[str, List[str]]] = None,
    celltypist_model: Optional[str] = None,
    output_dir: Optional[Path] = None
) -> tuple:
    """
    Annotate cells with multi-method consensus and biological context validation.

    This function combines CellTypist (ML-based) and marker-based annotation
    with a consensus decision logic that respects biological constraints based
    on sample type.

    Key insight: Organoid samples should NOT have immune cells. If the annotation
    produces immune cells in organoids, it indicates misannotation (e.g., S100A8/S100A9
    stress markers being incorrectly interpreted as neutrophil markers).

    Args:
        adata: AnnData object (should be log-normalized)
        sample_type: Type of sample:
            - 'organoid': Only epithelial/stem cells allowed (no immune)
            - 'primary_tissue': All cell types allowed
            - 'cell_line': Homogeneous population expected
            - 'pbmc': Only immune cells expected
        markers_dict: Optional custom markers dictionary (cell_type -> gene_list)
        celltypist_model: CellTypist model name (auto-selected if None)
        output_dir: Optional output directory for saving comparison CSV

    Returns:
        Tuple of (adata, validation_result)
            - adata: AnnData with annotations added to obs:
                - 'celltypist': CellTypist prediction
                - 'marker_based': Marker-based prediction
                - 'consensus_celltype': Final consensus annotation
                - 'annotation_confidence': 'High', 'Medium', or 'Low'
                - 'annotation_agreement': Agreement type between methods
            - validation_result: Dict with validation details

    Example:
        >>> adata, validation = annotate_with_validation(
        ...     adata, sample_type='organoid'
        ... )
        >>> if not validation['valid']:
        ...     print("Warning:", validation['summary'])
    """
    try:
        from tools.multi_method_annotation import MultiMethodAnnotator
        from tools.biological_validator import BiologicalValidator
    except ImportError:
        # Try relative import
        try:
            from .tools.multi_method_annotation import MultiMethodAnnotator
            from .tools.biological_validator import BiologicalValidator
        except ImportError:
            logger.error(
                "Multi-method annotation tools not found. "
                "Please ensure tools/multi_method_annotation.py and "
                "tools/biological_validator.py exist."
            )
            raise

    logger.info("=" * 60)
    logger.info(f"MULTI-METHOD ANNOTATION WITH VALIDATION")
    logger.info(f"Sample type: {sample_type}")
    logger.info("=" * 60)

    # Initialize annotator and validator
    annotator = MultiMethodAnnotator(sample_type=sample_type)
    validator = BiologicalValidator()

    # Run full annotation pipeline
    results = annotator.annotate(
        adata,
        markers_dict=markers_dict,
        celltypist_model=celltypist_model,
        add_to_adata=True
    )

    # Validate results
    validation = validator.validate_annotation(adata, sample_type, celltype_col='consensus_celltype')

    if not validation['valid']:
        logger.warning("=" * 60)
        logger.warning("BIOLOGICAL VALIDATION ISSUES DETECTED")
        logger.warning("=" * 60)
        for issue in validation['issues']:
            logger.warning(f"  - {issue['message']}")
        logger.warning("")
        logger.warning(f"Summary: {validation['summary']}")

        # Suggest reannotation if needed
        suggestion = validator.suggest_reannotation(adata, sample_type, celltype_col='consensus_celltype')
        if suggestion['needs_reannotation']:
            logger.warning(f"Severity: {suggestion['severity']}")
            logger.warning(f"Action: {suggestion['action']}")
    else:
        logger.info("Biological validation PASSED")

    # Save comparison CSV if output_dir provided
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

        comparison_df = pd.DataFrame({
            'cell_id': adata.obs_names,
            'celltypist': adata.obs.get('celltypist', 'N/A'),
            'marker_based': adata.obs.get('marker_based', 'N/A'),
            'consensus': adata.obs.get('consensus_celltype', 'N/A'),
            'confidence': adata.obs.get('annotation_confidence', 'N/A'),
            'agreement': adata.obs.get('annotation_agreement', 'N/A')
        })
        csv_path = output_dir / 'multi_method_comparison.csv'
        comparison_df.to_csv(csv_path, index=False)
        logger.info(f"Saved annotation comparison to {csv_path}")

    return adata, validation
