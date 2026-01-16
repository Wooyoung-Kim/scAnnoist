# Annotation Tools
# Actual implementations for automatic cell type annotation

import os
import pickle
from typing import Annotated, Optional, List, Dict, Any
import numpy as np
import pandas as pd

## 3rd party
from langchain_core.tools import tool

# Global state to store loaded data and models
_annotation_state: Dict[str, Any] = {
    "adata": None,
    "celltypist_model": None,
    "celltypist_predictions": None,
    "scanvi_model": None,
    "sctype_results": None,
    "scibet_model": None,
    "all_annotations": {},
}


def _get_adata():
    """Get the current AnnData object from state."""
    if _annotation_state["adata"] is None:
        raise ValueError("No data loaded. Use load_data tool first.")
    return _annotation_state["adata"]


# ============================================================
# Data Loading
# ============================================================
@tool
def load_adata_for_annotation(
    file_path: Annotated[str, "Path to the h5ad file"],
) -> Annotated[str, "Summary of loaded data"]:
    """Load an h5ad file for annotation. The data should already be QC'd and normalized."""
    import scanpy as sc
    
    adata = sc.read_h5ad(file_path)
    _annotation_state["adata"] = adata
    
    # Check for required preprocessing
    has_umap = "X_umap" in adata.obsm
    has_pca = "X_pca" in adata.obsm
    has_clusters = "leiden" in adata.obs or "louvain" in adata.obs
    
    summary = f"""
Loaded: {file_path}
- Cells: {adata.n_obs:,}
- Genes: {adata.n_vars:,}
- Has UMAP: {has_umap}
- Has PCA: {has_pca}
- Has clusters: {has_clusters}
- Obs columns: {list(adata.obs.columns)[:10]}...
"""
    return summary


# ============================================================
# CellTypist Tools
# ============================================================
@tool
def celltypist_list_models() -> Annotated[str, "Available CellTypist models"]:
    """List all available pre-trained CellTypist models."""
    try:
        import celltypist
        from celltypist import models
        
        # Get model descriptions
        model_info = models.models_description()
        
        # Format output
        output = "Available CellTypist Models:\n"
        output += "=" * 50 + "\n"
        for idx, row in model_info.iterrows():
            output += f"\n{row['model']}\n"
            output += f"  - {row.get('description', 'No description')[:100]}...\n"
        
        return output
    except ImportError:
        return "CellTypist not installed. Run: pip install celltypist"


@tool
def celltypist_download_model(
    model_name: Annotated[str, "Name of the model to download (e.g., 'Immune_All_High.pkl')"],
) -> Annotated[str, "Download result"]:
    """Download a CellTypist pre-trained model."""
    try:
        import celltypist
        from celltypist import models
        
        models.download_models(model=model_name)
        return f"Successfully downloaded CellTypist model: {model_name}"
    except Exception as e:
        return f"Error downloading model: {str(e)}"


@tool
def celltypist_annotate(
    model: Annotated[str, "CellTypist model name (e.g., 'Immune_All_High.pkl')"],
    majority_voting: Annotated[bool, "Use cluster-level majority voting for consensus"] = True,
    cluster_key: Annotated[Optional[str], "Column for majority voting (e.g., 'leiden')"] = None,
) -> Annotated[str, "CellTypist annotation results"]:
    """
    Run CellTypist automatic annotation.
    Uses pre-trained deep learning models for cell type prediction.
    Requires log-normalized data.
    """
    try:
        import celltypist
        from celltypist import models
        
        adata = _get_adata()
        
        # Load model
        ct_model = models.Model.load(model=model)
        _annotation_state["celltypist_model"] = ct_model
        
        # Run prediction
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=majority_voting
        )
        
        # Store predictions
        _annotation_state["celltypist_predictions"] = predictions
        
        # Add to adata
        adata.obs["celltypist_prediction"] = predictions.predicted_labels.predicted_labels
        if majority_voting:
            adata.obs["celltypist_majority"] = predictions.predicted_labels.majority_voting
        
        # Get confidence scores
        adata.obs["celltypist_confidence"] = predictions.probability_matrix.max(axis=1).values
        
        # Summary statistics
        pred_counts = adata.obs["celltypist_prediction"].value_counts()
        
        result = f"""
CellTypist Annotation Complete
==============================
Model: {model}
Majority voting: {majority_voting}

Cell Type Distribution:
{pred_counts.head(15).to_string()}

Confidence Statistics:
- Mean: {adata.obs['celltypist_confidence'].mean():.3f}
- Median: {adata.obs['celltypist_confidence'].median():.3f}
- Min: {adata.obs['celltypist_confidence'].min():.3f}
- Low confidence (<0.5): {(adata.obs['celltypist_confidence'] < 0.5).sum():,} cells
"""
        return result
        
    except ImportError:
        return "CellTypist not installed. Run: pip install celltypist"
    except Exception as e:
        return f"Error running CellTypist: {str(e)}"


@tool
def celltypist_get_probabilities(
    top_n: Annotated[int, "Number of top cell types to show per cell"] = 3,
) -> Annotated[str, "Probability information"]:
    """Get CellTypist probability matrix and top predictions per cell."""
    predictions = _annotation_state.get("celltypist_predictions")
    if predictions is None:
        return "No CellTypist predictions found. Run celltypist_annotate first."

    prob_matrix = predictions.probability_matrix

    # Get top N predictions for each cell
    top_types = []
    for idx in range(min(5, len(prob_matrix))):  # Show first 5 cells as example
        cell_probs = prob_matrix.iloc[idx].sort_values(ascending=False)
        top = cell_probs.head(top_n)
        top_types.append({
            "cell": prob_matrix.index[idx],
            "predictions": [(t, f"{p:.3f}") for t, p in zip(top.index, top.values)]
        })

    result = f"""
CellTypist Probability Matrix
=============================
Shape: {prob_matrix.shape}
Cell types: {list(prob_matrix.columns)[:10]}...

Example predictions (top {top_n} per cell):
"""
    for cell_info in top_types:
        result += f"\n{cell_info['cell']}: {cell_info['predictions']}"

    return result


@tool
def celltypist_annotate_tissue_aware(
    sample_col: Annotated[str, "Column containing sample IDs"] = "sample",
    majority_voting: Annotated[bool, "Use cluster-level majority voting for consensus"] = True,
    cluster_key: Annotated[Optional[str], "Column for majority voting (e.g., 'leiden')"] = None,
    output_dir: Annotated[Optional[str], "Directory to save outputs"] = None,
) -> Annotated[str, "Tissue-aware CellTypist annotation results"]:
    """
    Run tissue-aware CellTypist annotation.

    Automatically selects the appropriate model for each tissue in the dataset
    based on the model registry. If multiple tissues are present, runs CellTypist
    separately per tissue subset and merges results.

    Requires tissue metadata to be present in adata.obs["tissue"].
    Use validate_and_add_tissue_metadata first to add tissue annotations.

    Args:
        sample_col: Column containing sample IDs
        majority_voting: Whether to use majority voting
        cluster_key: Column for majority voting (if enabled)
        output_dir: Directory to save prediction CSVs

    Returns:
        Summary of annotation results
    """
    try:
        import celltypist
        from celltypist import models
        from scrna_agent.tools.model_management import (
            get_models_for_adata,
            save_run_metadata
        )
    except ImportError as e:
        return f"Import error: {str(e)}. Install celltypist and check model_management.py"

    adata = _get_adata()

    # Validate tissue metadata exists
    if "tissue" not in adata.obs.columns:
        return """
Error: Tissue metadata not found in adata.obs["tissue"].
Please run validate_and_add_tissue_metadata first to add tissue annotations.
"""

    # Get model assignments for each tissue
    try:
        tissue_info = get_models_for_adata(adata, sample_col=sample_col)
    except Exception as e:
        return f"Error determining models for tissues: {str(e)}"

    # Store tissue info for later use
    _annotation_state["tissue_model_info"] = tissue_info

    # Initialize result columns
    adata.obs["celltypist_label"] = ""
    adata.obs["celltypist_confidence"] = 0.0
    if majority_voting and cluster_key:
        adata.obs["celltypist_majority_vote"] = ""
    adata.obs["celltypist_model_used"] = ""

    # Run annotation per tissue
    all_predictions = []
    summary_lines = []

    summary_lines.append("Tissue-Aware CellTypist Annotation")
    summary_lines.append("=" * 50)
    summary_lines.append(f"Total tissues: {len(tissue_info)}")
    summary_lines.append(f"Sample column: {sample_col}")
    summary_lines.append(f"Majority voting: {majority_voting}")
    summary_lines.append("")

    for tissue, info in tissue_info.items():
        model_name = info["model"]
        n_cells = info["n_cells"]
        rule = info["rule"]

        summary_lines.append(f"Tissue: {tissue}")
        summary_lines.append(f"  Model: {model_name} (selected by {rule})")
        summary_lines.append(f"  Cells: {n_cells:,}")
        summary_lines.append(f"  Samples: {', '.join(info['samples'])}")

        # Subset data for this tissue
        tissue_mask = adata.obs["tissue"] == tissue
        adata_tissue = adata[tissue_mask].copy()

        try:
            # Load model
            ct_model = models.Model.load(model=model_name)

            # Run prediction
            if majority_voting and cluster_key and cluster_key in adata_tissue.obs:
                predictions = celltypist.annotate(
                    adata_tissue,
                    model=model_name,
                    majority_voting=True
                )

                # Extract predictions
                pred_labels = predictions.predicted_labels.predicted_labels
                majority_labels = predictions.predicted_labels.majority_voting
                confidences = predictions.probability_matrix.max(axis=1).values

                # Add to main adata
                cell_indices = adata.obs.index[tissue_mask]
                adata.obs.loc[cell_indices, "celltypist_label"] = pred_labels.values
                adata.obs.loc[cell_indices, "celltypist_majority_vote"] = majority_labels.values
                adata.obs.loc[cell_indices, "celltypist_confidence"] = confidences
                adata.obs.loc[cell_indices, "celltypist_model_used"] = model_name

            else:
                predictions = celltypist.annotate(
                    adata_tissue,
                    model=model_name,
                    majority_voting=False
                )

                # Extract predictions
                pred_labels = predictions.predicted_labels.predicted_labels
                confidences = predictions.probability_matrix.max(axis=1).values

                # Add to main adata
                cell_indices = adata.obs.index[tissue_mask]
                adata.obs.loc[cell_indices, "celltypist_label"] = pred_labels.values
                adata.obs.loc[cell_indices, "celltypist_confidence"] = confidences
                adata.obs.loc[cell_indices, "celltypist_model_used"] = model_name

            # Store predictions for this tissue
            all_predictions.append({
                "tissue": tissue,
                "model": model_name,
                "predictions": predictions
            })

            # Summary stats
            label_counts = pd.Series(pred_labels).value_counts()
            summary_lines.append(f"  Top cell types: {', '.join(label_counts.head(3).index.tolist())}")
            summary_lines.append(f"  Mean confidence: {confidences.mean():.3f}")
            summary_lines.append("")

        except Exception as e:
            summary_lines.append(f"  ERROR: {str(e)}")
            summary_lines.append("")
            continue

    # Store all predictions in state
    _annotation_state["celltypist_tissue_predictions"] = all_predictions

    # Generate outputs if output_dir specified
    if output_dir:
        try:
            save_celltypist_outputs(adata, tissue_info, output_dir, sample_col)
            summary_lines.append(f"Outputs saved to: {output_dir}")
        except Exception as e:
            summary_lines.append(f"Warning: Could not save outputs: {str(e)}")

    # Overall statistics
    summary_lines.append("")
    summary_lines.append("Overall Statistics:")
    summary_lines.append("=" * 50)
    overall_counts = adata.obs["celltypist_label"].value_counts()
    summary_lines.append(f"Total unique cell types: {len(overall_counts)}")
    summary_lines.append(f"\nTop 15 cell types:")
    summary_lines.append(overall_counts.head(15).to_string())
    summary_lines.append(f"\nConfidence statistics:")
    summary_lines.append(f"  Mean: {adata.obs['celltypist_confidence'].mean():.3f}")
    summary_lines.append(f"  Median: {adata.obs['celltypist_confidence'].median():.3f}")
    summary_lines.append(f"  Low confidence (<0.5): {(adata.obs['celltypist_confidence'] < 0.5).sum():,} cells")

    return "\n".join(summary_lines)


def save_celltypist_outputs(
    adata,
    tissue_info: Dict[str, Dict[str, Any]],
    output_dir: str,
    sample_col: str
) -> None:
    """
    Save CellTypist prediction outputs.

    Generates:
    - celltypist_predictions.csv: Per-cell predictions with metadata
    - tissue_model_usage.csv: Summary of models used per tissue

    Args:
        adata: Annotated AnnData object
        tissue_info: Tissue information from get_models_for_adata
        output_dir: Output directory
        sample_col: Sample column name
    """
    import os
    from scrna_agent.tools.model_management import save_run_metadata

    os.makedirs(output_dir, exist_ok=True)

    # 1. Save per-cell predictions
    prediction_cols = ["celltypist_label", "celltypist_confidence", "celltypist_model_used"]
    if "celltypist_majority_vote" in adata.obs.columns:
        prediction_cols.append("celltypist_majority_vote")

    metadata_cols = []
    if sample_col in adata.obs.columns:
        metadata_cols.append(sample_col)
    if "tissue" in adata.obs.columns:
        metadata_cols.append("tissue")
    if "sample_type" in adata.obs.columns:
        metadata_cols.append("sample_type")

    # Build output dataframe
    output_df = pd.DataFrame(index=adata.obs.index)
    output_df["cell_id"] = adata.obs.index

    for col in metadata_cols:
        output_df[col] = adata.obs[col].values

    output_df["label"] = adata.obs["celltypist_label"].values
    output_df["confidence"] = adata.obs["celltypist_confidence"].values
    output_df["model_used"] = adata.obs["celltypist_model_used"].values

    if "celltypist_majority_vote" in adata.obs.columns:
        output_df["majority_vote"] = adata.obs["celltypist_majority_vote"].values

    predictions_path = os.path.join(output_dir, "celltypist_predictions.csv")
    output_df.to_csv(predictions_path, index=False)
    print(f"Saved predictions to {predictions_path}")

    # 2. Save tissue model usage summary
    tissue_summary = []
    for tissue, info in tissue_info.items():
        tissue_summary.append({
            "tissue": tissue,
            "model_used": info["model"],
            "selection_rule": info["rule"],
            "n_cells": info["n_cells"],
            "samples": ";".join(info["samples"]),
            "sample_type": info.get("sample_type", "")
        })

    tissue_df = pd.DataFrame(tissue_summary)
    tissue_usage_path = os.path.join(output_dir, "tissue_model_usage.csv")
    tissue_df.to_csv(tissue_usage_path, index=False)
    print(f"Saved tissue model usage to {tissue_usage_path}")

    # 3. Save run metadata JSON
    save_run_metadata(output_dir, tissue_info)


@tool
def add_tissue_metadata(
    tissue: Annotated[Optional[str], "Tissue type to apply to all cells"] = None,
    tissue_map: Annotated[Optional[str], "Path to CSV with sample->tissue mapping"] = None,
    sample_type: Annotated[Optional[str], "Sample type to apply to all cells"] = None,
    sample_type_map: Annotated[Optional[str], "Path to CSV with sample->sample_type mapping"] = None,
    sample_col: Annotated[str, "Column name for sample IDs"] = "sample",
    allow_missing_tissue: Annotated[bool, "Fill missing tissues with 'unknown'"] = False,
) -> Annotated[str, "Tissue metadata validation result"]:
    """
    Validate and add tissue/sample_type metadata to AnnData.

    This must be called before celltypist_annotate_tissue_aware.

    Args:
        tissue: Tissue type to apply to all cells
        tissue_map: Path to CSV with sample->tissue mapping
        sample_type: Sample type to apply to all cells
        sample_type_map: Path to CSV with sample->sample_type mapping
        sample_col: Column name for sample IDs
        allow_missing_tissue: If True, fill missing with "unknown"

    Returns:
        Summary of added metadata
    """
    try:
        from scrna_agent.tools.model_management import validate_and_add_tissue_metadata
    except ImportError as e:
        return f"Import error: {str(e)}. Check model_management.py"

    adata = _get_adata()

    try:
        validate_and_add_tissue_metadata(
            adata,
            tissue=tissue,
            tissue_map=tissue_map,
            sample_type=sample_type,
            sample_type_map=sample_type_map,
            sample_col=sample_col,
            allow_missing_tissue=allow_missing_tissue
        )
    except Exception as e:
        return f"Error adding tissue metadata: {str(e)}"

    # Summary
    result = ["Tissue Metadata Added Successfully", "=" * 50]

    if "tissue" in adata.obs.columns:
        tissue_counts = adata.obs["tissue"].value_counts()
        result.append(f"\nTissue distribution:")
        result.append(tissue_counts.to_string())

    if "sample_type" in adata.obs.columns:
        sample_type_counts = adata.obs["sample_type"].value_counts()
        result.append(f"\nSample type distribution:")
        result.append(sample_type_counts.to_string())

    result.append(f"\nTotal cells: {adata.n_obs:,}")
    result.append(f"Sample column: {sample_col}")

    return "\n".join(result)


# ============================================================
# scANVI Annotation Tools
# ============================================================
@tool
def scanvi_setup_reference(
    reference_path: Annotated[str, "Path to reference h5ad file with labels"],
    labels_key: Annotated[str, "Column containing cell type labels"],
    batch_key: Annotated[Optional[str], "Column for batch correction"] = None,
    layer: Annotated[str, "Layer with raw counts"] = "counts",
) -> Annotated[str, "Reference setup result"]:
    """
    Prepare reference data for scANVI annotation.
    Reference must have cell type labels.
    """
    import scvi
    import scanpy as sc
    
    # Load reference
    adata_ref = sc.read_h5ad(reference_path)
    
    # Verify labels exist
    if labels_key not in adata_ref.obs:
        return f"Error: Labels column '{labels_key}' not found in reference"
    
    # Setup for scVI
    setup_kwargs = {"layer": layer, "labels_key": labels_key}
    if batch_key:
        setup_kwargs["batch_key"] = batch_key
    
    scvi.model.SCVI.setup_anndata(adata_ref, **setup_kwargs)
    
    _annotation_state["reference_adata"] = adata_ref
    _annotation_state["reference_labels_key"] = labels_key
    
    label_counts = adata_ref.obs[labels_key].value_counts()
    
    return f"""
Reference Setup Complete
========================
Reference: {reference_path}
- Cells: {adata_ref.n_obs:,}
- Genes: {adata_ref.n_vars:,}
- Labels column: {labels_key}
- Batch column: {batch_key}

Cell Type Distribution:
{label_counts.to_string()}
"""


@tool
def scanvi_train_and_transfer(
    max_epochs_scvi: Annotated[int, "Max epochs for scVI pre-training"] = 400,
    max_epochs_scanvi: Annotated[int, "Max epochs for scANVI training"] = 100,
    n_latent: Annotated[int, "Latent dimension"] = 30,
) -> Annotated[str, "scANVI training and transfer results"]:
    """
    Train scVI + scANVI on reference and transfer labels to query.
    Assumes reference is already setup and query is in _annotation_state['adata'].
    """
    import scvi
    import scanpy as sc
    
    adata_ref = _annotation_state.get("reference_adata")
    adata_query = _get_adata()
    labels_key = _annotation_state.get("reference_labels_key")
    
    if adata_ref is None:
        return "Error: Reference not setup. Run scanvi_setup_reference first."
    
    # Find common genes
    common_genes = adata_ref.var_names.intersection(adata_query.var_names)
    adata_ref = adata_ref[:, common_genes].copy()
    adata_query = adata_query[:, common_genes].copy()
    
    # Mark query as unlabeled
    adata_query.obs[labels_key] = "Unknown"
    
    # Concatenate
    adata_combined = sc.concat([adata_ref, adata_query], label="dataset", keys=["reference", "query"])
    
    # Setup
    scvi.model.SCVI.setup_anndata(
        adata_combined,
        layer="counts" if "counts" in adata_combined.layers else None,
        batch_key="dataset",
        labels_key=labels_key,
    )
    
    # Train scVI
    scvi_model = scvi.model.SCVI(adata_combined, n_latent=n_latent, n_layers=2)
    scvi_model.train(max_epochs=max_epochs_scvi, early_stopping=True)
    
    # Convert to scANVI
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_key,
        unlabeled_category="Unknown"
    )
    scanvi_model.train(max_epochs=max_epochs_scanvi, early_stopping=True)
    
    _annotation_state["scanvi_model"] = scanvi_model
    _annotation_state["combined_adata"] = adata_combined
    
    # Get predictions
    predictions = scanvi_model.predict()
    adata_combined.obs["scanvi_prediction"] = predictions
    
    # Get soft predictions
    soft_preds = scanvi_model.predict(soft=True)
    adata_combined.obs["scanvi_confidence"] = soft_preds.max(axis=1)
    
    # Extract query predictions
    query_mask = adata_combined.obs["dataset"] == "query"
    query_preds = adata_combined.obs.loc[query_mask, "scanvi_prediction"]
    query_conf = adata_combined.obs.loc[query_mask, "scanvi_confidence"]
    
    # Update original query adata
    adata_query = _get_adata()
    adata_query.obs["scanvi_prediction"] = query_preds.values
    adata_query.obs["scanvi_confidence"] = query_conf.values
    
    pred_counts = query_preds.value_counts()
    
    return f"""
scANVI Training and Transfer Complete
=====================================
scVI epochs: {max_epochs_scvi}
scANVI epochs: {max_epochs_scanvi}
Latent dimension: {n_latent}
Common genes: {len(common_genes):,}

Query Cell Type Predictions:
{pred_counts.head(15).to_string()}

Confidence Statistics:
- Mean: {query_conf.mean():.3f}
- Median: {query_conf.median():.3f}
- Low confidence (<0.7): {(query_conf < 0.7).sum():,} cells
"""


# ============================================================
# ScType Tools (R-based functions called via rpy2)
# ============================================================
@tool
def sctype_annotate(
    tissue: Annotated[str, "Tissue type for marker database"],
    cluster_key: Annotated[str, "Column containing cluster labels"] = "leiden",
) -> Annotated[str, "ScType annotation results"]:
    """
    Run ScType marker-based annotation.
    Available tissues: Immune system, Brain, Lung, Liver, Spleen, etc.
    Note: This is a simplified Python implementation inspired by ScType.
    """
    adata = _get_adata()
    
    # ScType marker database (simplified for key tissues)
    sctype_markers = {
        "Immune system": {
            "B cells": ["Cd19", "Cd79a", "Ms4a1", "Pax5"],
            "T cells": ["Cd3d", "Cd3e", "Cd3g", "Trac"],
            "CD4 T cells": ["Cd4", "Cd3d", "Il7r"],
            "CD8 T cells": ["Cd8a", "Cd8b1", "Gzmb"],
            "NK cells": ["Nkg7", "Klrb1c", "Ncr1", "Klrd1"],
            "Monocytes": ["Cd14", "Lyz2", "Csf1r"],
            "Macrophages": ["Adgre1", "Cd68", "Mrc1"],
            "Dendritic cells": ["Itgax", "Cd74", "H2-Ab1"],
            "Neutrophils": ["S100a8", "S100a9", "Ly6g"],
            "Plasma cells": ["Sdc1", "Xbp1", "Prdm1"],
        },
        "Spleen": {
            "ImmB": ["Cd19", "Cd93", "Ighm"],
            "NB1": ["Cd19", "Ighd", "Cd21"],
            "NB2": ["Cd19", "Ighd", "Cd23"],
            "MZB": ["Cd19", "Cd21", "Cd1d"],
            "PB": ["Sdc1", "Xbp1", "Prdm1"],
            "T cells": ["Cd3d", "Cd3e"],
            "Macrophages": ["Adgre1", "Cd68"],
            "Neutrophils": ["S100a8", "S100a9"],
        },
    }
    
    if tissue not in sctype_markers:
        available = list(sctype_markers.keys())
        return f"Tissue '{tissue}' not found. Available: {available}"
    
    markers = sctype_markers[tissue]
    
    # Check cluster key
    if cluster_key not in adata.obs:
        return f"Cluster column '{cluster_key}' not found in adata.obs"
    
    # Score each cluster for each cell type
    cluster_scores = {}
    for cluster in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_data = adata[cluster_mask]
        
        scores = {}
        for cell_type, marker_genes in markers.items():
            # Find markers present in data
            present_markers = [g for g in marker_genes if g in adata.var_names]
            if not present_markers:
                scores[cell_type] = 0
                continue
            
            # Calculate mean expression
            if hasattr(cluster_data.X, 'toarray'):
                expr = cluster_data[:, present_markers].X.toarray()
            else:
                expr = cluster_data[:, present_markers].X
            
            scores[cell_type] = np.mean(expr)
        
        # Get best cell type
        best_type = max(scores, key=scores.get)
        best_score = scores[best_type]
        cluster_scores[cluster] = {"type": best_type, "score": best_score}
    
    # Add to adata
    adata.obs["sctype_annotation"] = adata.obs[cluster_key].map(
        lambda x: cluster_scores[x]["type"]
    )
    
    _annotation_state["sctype_results"] = cluster_scores
    
    # Format results
    result = f"""
ScType Annotation Complete
==========================
Tissue: {tissue}
Cluster column: {cluster_key}

Cluster Annotations:
"""
    for cluster, info in sorted(cluster_scores.items(), key=lambda x: str(x[0])):
        n_cells = (adata.obs[cluster_key] == cluster).sum()
        result += f"  Cluster {cluster}: {info['type']} (score: {info['score']:.3f}, n={n_cells:,})\n"
    
    return result


# ============================================================
# Comparison and Consensus Tools
# ============================================================
@tool
def compare_annotations(
    methods: Annotated[List[str], "List of annotation columns to compare"],
) -> Annotated[str, "Comparison table"]:
    """Compare annotations from multiple methods."""
    adata = _get_adata()
    
    # Check which methods are available
    available = [m for m in methods if m in adata.obs.columns]
    if not available:
        return f"None of the specified columns found. Available: {list(adata.obs.columns)}"
    
    missing = [m for m in methods if m not in adata.obs.columns]
    
    result = f"""
Annotation Comparison
=====================
Comparing: {available}
Missing: {missing}

"""
    
    # Per-method cell type counts
    for method in available:
        counts = adata.obs[method].value_counts()
        result += f"\n{method}:\n{counts.head(10).to_string()}\n"
    
    # Pairwise agreement if multiple methods
    if len(available) >= 2:
        result += "\nPairwise Agreement:\n"
        for i, m1 in enumerate(available):
            for m2 in available[i+1:]:
                agreement = (adata.obs[m1] == adata.obs[m2]).mean()
                result += f"  {m1} vs {m2}: {agreement:.1%}\n"
    
    return result


@tool
def create_consensus_annotation(
    methods: Annotated[List[str], "Annotation columns to combine"],
    strategy: Annotated[str, "Strategy: 'majority', 'first_available', or 'highest_confidence'"] = "majority",
    confidence_cols: Annotated[Optional[List[str]], "Confidence columns for 'highest_confidence' strategy"] = None,
    output_col: Annotated[str, "Output column name"] = "consensus_annotation",
) -> Annotated[str, "Consensus annotation result"]:
    """Create consensus annotation by combining multiple methods."""
    adata = _get_adata()
    
    available = [m for m in methods if m in adata.obs.columns]
    if not available:
        return f"None of the specified columns found."
    
    if strategy == "majority":
        # Majority voting
        def get_majority(row):
            votes = [row[m] for m in available if pd.notna(row[m])]
            if not votes:
                return "Unknown"
            return pd.Series(votes).mode().iloc[0]
        
        adata.obs[output_col] = adata.obs.apply(get_majority, axis=1)
        
    elif strategy == "first_available":
        # Use first non-null value
        adata.obs[output_col] = adata.obs[available[0]]
        for method in available[1:]:
            mask = adata.obs[output_col].isna() | (adata.obs[output_col] == "Unknown")
            adata.obs.loc[mask, output_col] = adata.obs.loc[mask, method]
            
    elif strategy == "highest_confidence":
        if not confidence_cols or len(confidence_cols) != len(available):
            return "For 'highest_confidence', provide matching confidence columns."
        
        # Select annotation with highest confidence
        conf_df = adata.obs[confidence_cols]
        best_idx = conf_df.values.argmax(axis=1)
        adata.obs[output_col] = [
            adata.obs[available[idx]].iloc[i] 
            for i, idx in enumerate(best_idx)
        ]
    
    # Summary
    counts = adata.obs[output_col].value_counts()
    
    return f"""
Consensus Annotation Created
============================
Strategy: {strategy}
Methods used: {available}
Output column: {output_col}

Cell Type Distribution:
{counts.head(15).to_string()}
"""


@tool
def validate_with_markers(
    annotation_col: Annotated[str, "Annotation column to validate"],
    markers: Annotated[Dict[str, List[str]], "Dict of cell_type: [marker_genes]"],
) -> Annotated[str, "Marker validation result"]:
    """Validate cell type annotation by checking expected marker expression."""
    import scanpy as sc
    
    adata = _get_adata()
    
    if annotation_col not in adata.obs:
        return f"Annotation column '{annotation_col}' not found."
    
    result = f"""
Marker Validation
=================
Annotation column: {annotation_col}

"""
    
    for cell_type, marker_genes in markers.items():
        # Check which markers are in data
        present = [g for g in marker_genes if g in adata.var_names]
        missing = [g for g in marker_genes if g not in adata.var_names]
        
        if not present:
            result += f"\n{cell_type}: No markers found in data. Missing: {missing}\n"
            continue
        
        # Get cells of this type
        type_mask = adata.obs[annotation_col] == cell_type
        if type_mask.sum() == 0:
            result += f"\n{cell_type}: No cells found with this annotation.\n"
            continue
        
        # Calculate mean expression
        type_cells = adata[type_mask, present]
        other_cells = adata[~type_mask, present]
        
        if hasattr(type_cells.X, 'toarray'):
            type_expr = type_cells.X.toarray().mean(axis=0)
            other_expr = other_cells.X.toarray().mean(axis=0)
        else:
            type_expr = type_cells.X.mean(axis=0)
            other_expr = other_cells.X.mean(axis=0)
        
        result += f"\n{cell_type} (n={type_mask.sum():,}):\n"
        for i, gene in enumerate(present):
            fold = type_expr[i] / (other_expr[i] + 0.01)
            status = "✓" if fold > 1.5 else "⚠" if fold > 1.0 else "✗"
            result += f"  {status} {gene}: {type_expr[i]:.2f} vs {other_expr[i]:.2f} (FC={fold:.1f})\n"
        
        if missing:
            result += f"  Missing: {missing}\n"
    
    return result


@tool
def save_annotations(
    output_path: Annotated[str, "Path to save annotated h5ad file"],
) -> Annotated[str, "Save result"]:
    """Save the annotated AnnData to file."""
    adata = _get_adata()
    adata.write_h5ad(output_path)
    
    # Also save annotations as CSV
    csv_path = output_path.replace(".h5ad", "_annotations.csv")
    annotation_cols = [c for c in adata.obs.columns if "annotation" in c.lower() or 
                       "prediction" in c.lower() or "celltypist" in c.lower() or
                       "scanvi" in c.lower() or "sctype" in c.lower() or
                       "consensus" in c.lower()]
    
    if annotation_cols:
        adata.obs[annotation_cols].to_csv(csv_path)
    
    return f"""
Saved Annotations
=================
H5AD: {output_path}
CSV: {csv_path}
Annotation columns: {annotation_cols}
"""


@tool
def plot_annotation_umap(
    annotation_col: Annotated[str, "Annotation column to plot"],
    output_path: Annotated[str, "Path to save the plot"],
    figsize: Annotated[tuple, "Figure size (width, height)"] = (10, 8),
) -> Annotated[str, "Path to saved plot"]:
    """Generate UMAP colored by annotation."""
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    adata = _get_adata()
    
    if "X_umap" not in adata.obsm:
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
    
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(adata, color=annotation_col, ax=ax, show=False, legend_loc="on data")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    return f"Saved annotation UMAP to {output_path}"


# ============================================================
# PubMed RAG-Based Marker Annotation Tools
# ============================================================

# Cache for literature-derived markers
_literature_markers_cache: Dict[str, Dict[str, List[str]]] = {}


@tool
def search_pubmed_markers(
    cell_type: Annotated[str, "Cell type to search markers for (e.g., 'B cells', 'T helper cells')"],
    species: Annotated[str, "Species: 'mouse' or 'human'"] = "mouse",
    tissue: Annotated[Optional[str], "Specific tissue context (e.g., 'spleen', 'bone marrow')"] = None,
    max_results: Annotated[int, "Maximum number of papers to search"] = 10,
) -> Annotated[str, "Marker genes found in literature"]:
    """
    Search PubMed for marker genes of a specific cell type.
    Uses NCBI E-utilities to find relevant papers and extract marker information.
    
    This tool queries PubMed for papers about the specified cell type and extracts
    commonly mentioned marker genes from abstracts.
    """
    from Bio import Entrez
    import os
    import re
    
    # Setup Entrez (requires NCBI_API_KEY and EMAIL in environment)
    email = os.getenv("EMAIL1") or os.getenv("ENTREZ_EMAIL") or "user@example.com"
    api_key = os.getenv("NCBI_API_KEY1") or os.getenv("NCBI_API_KEY")
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key
    
    # Build search query
    query_parts = [f'"{cell_type}"[Title/Abstract]', f"{species}[Title/Abstract]"]
    if tissue:
        query_parts.append(f'"{tissue}"[Title/Abstract]')
    query_parts.append("(marker OR markers OR expression)[Title/Abstract]")
    query_parts.append("(single cell OR scRNA-seq OR RNA-seq)[Title/Abstract]")
    
    query = " AND ".join(query_parts)
    
    try:
        # Search PubMed
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        
        pmids = record.get("IdList", [])
        
        if not pmids:
            return f"No papers found for '{cell_type}' in {species}. Try broader search terms."
        
        # Fetch abstracts
        handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
        abstracts = handle.read()
        handle.close()
        
        # Extract potential marker genes (gene symbols are typically uppercase or Title case)
        # Common patterns: Cd19, CD19, Pax5, PAX5, etc.
        gene_pattern = r'\b([A-Z][a-z0-9]{1,6}|[A-Z]{2,6}[0-9]*)\b'
        potential_genes = re.findall(gene_pattern, abstracts)
        
        # Filter to likely gene names (common patterns)
        # Exclude common words that match the pattern
        exclude_words = {'The', 'This', 'That', 'With', 'From', 'For', 'And', 'Not', 
                        'Are', 'Was', 'Were', 'Has', 'Have', 'Had', 'May', 'Can',
                        'Cell', 'Cells', 'RNA', 'DNA', 'PCR', 'FACS', 'Fig', 'Table'}
        
        gene_counts = {}
        for gene in potential_genes:
            if gene not in exclude_words and len(gene) >= 2:
                gene_counts[gene] = gene_counts.get(gene, 0) + 1
        
        # Get top mentioned genes
        top_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)[:20]
        
        # Store in cache
        cache_key = f"{cell_type}_{species}_{tissue or 'any'}"
        _literature_markers_cache[cache_key] = {
            "markers": [g[0] for g in top_genes],
            "counts": dict(top_genes),
            "pmids": pmids,
        }
        
        result = f"""
PubMed Marker Search Results
============================
Query: {cell_type} ({species}{f', {tissue}' if tissue else ''})
Papers found: {len(pmids)}
PubMed IDs: {', '.join(pmids[:5])}{'...' if len(pmids) > 5 else ''}

Top Mentioned Genes (potential markers):
"""
        for gene, count in top_genes[:15]:
            result += f"  - {gene}: mentioned {count} times\n"
        
        result += f"""
Note: These are candidate markers extracted from abstracts. 
Validate with known marker databases before using for annotation.
"""
        return result
        
    except Exception as e:
        return f"Error searching PubMed: {str(e)}. Check NCBI_API_KEY and EMAIL environment variables."


@tool
def get_cellmarker_database(
    species: Annotated[str, "Species: 'mouse' or 'human'"] = "mouse",
    tissue: Annotated[Optional[str], "Tissue type (e.g., 'Spleen', 'Bone marrow')"] = None,
) -> Annotated[str, "Cell markers from CellMarker database"]:
    """
    Get curated cell type markers from the CellMarker database.
    CellMarker is a manually curated resource of cell markers from literature.
    
    This provides a baseline of known markers to complement PubMed search.
    """
    # Curated markers from CellMarker database (subset)
    # Source: http://biocc.hrbmu.edu.cn/CellMarker/
    cellmarker_db = {
        "mouse": {
            "B cells": ["Cd19", "Cd79a", "Cd79b", "Ms4a1", "Pax5", "Ebf1"],
            "Naive B cells": ["Cd19", "Ighd", "Fcer2a", "Cd38"],
            "Memory B cells": ["Cd19", "Cd27", "Nt5e"],
            "Plasma cells": ["Sdc1", "Xbp1", "Prdm1", "Ighg1"],
            "Plasmablasts": ["Sdc1", "Mki67", "Xbp1"],
            "Germinal center B cells": ["Bcl6", "Aicda", "Fas", "Cd38"],
            "Marginal zone B cells": ["Cd21", "Cd1d", "Cd9"],
            "Immature B cells": ["Cd19", "Cd93", "Ighm", "Kit"],
            "T cells": ["Cd3d", "Cd3e", "Cd3g", "Trac", "Trbc1"],
            "CD4 T cells": ["Cd4", "Cd3d", "Il7r"],
            "CD8 T cells": ["Cd8a", "Cd8b1", "Cd3d", "Gzmb"],
            "Tfh cells": ["Cxcr5", "Bcl6", "Pdcd1", "Icos", "Cd4"],
            "Treg cells": ["Foxp3", "Il2ra", "Cd4", "Ctla4"],
            "NK cells": ["Ncr1", "Klrb1c", "Nkg7", "Klrd1", "Prf1"],
            "Monocytes": ["Cd14", "Lyz2", "Csf1r", "Ccr2"],
            "Macrophages": ["Adgre1", "Cd68", "Mrc1", "Csf1r"],
            "Dendritic cells": ["Itgax", "Cd74", "H2-Ab1", "Flt3"],
            "Neutrophils": ["S100a8", "S100a9", "Ly6g", "Cxcr2"],
        },
        "human": {
            "B cells": ["CD19", "CD79A", "CD79B", "MS4A1", "PAX5"],
            "Naive B cells": ["CD19", "IGHD", "FCER2", "TCL1A"],
            "Memory B cells": ["CD19", "CD27", "NT5E"],
            "Plasma cells": ["SDC1", "XBP1", "PRDM1", "MZB1"],
            "Germinal center B cells": ["BCL6", "AICDA", "FAS"],
            "T cells": ["CD3D", "CD3E", "CD3G", "TRAC"],
            "CD4 T cells": ["CD4", "CD3D", "IL7R"],
            "CD8 T cells": ["CD8A", "CD8B", "CD3D", "GZMB"],
            "Tfh cells": ["CXCR5", "BCL6", "PDCD1", "ICOS"],
            "Treg cells": ["FOXP3", "IL2RA", "CD4", "CTLA4"],
            "NK cells": ["NCR1", "KLRB1", "NKG7", "GNLY"],
            "Monocytes": ["CD14", "LYZ", "CSF1R", "FCGR3A"],
            "Macrophages": ["CD68", "MRC1", "MARCO"],
            "Dendritic cells": ["ITGAX", "FLT3", "HLA-DRA"],
        }
    }
    
    species_markers = cellmarker_db.get(species.lower(), {})
    
    if not species_markers:
        return f"Species '{species}' not in database. Available: mouse, human"
    
    result = f"""
CellMarker Database ({species})
{'=' * 30}
"""
    
    for cell_type, markers in species_markers.items():
        if tissue is None or tissue.lower() in cell_type.lower():
            result += f"\n{cell_type}:\n  {', '.join(markers)}\n"
    
    return result


@tool
def annotate_with_literature_markers(
    cell_type_markers: Annotated[Dict[str, List[str]], "Dict mapping cell types to marker genes from literature"],
    cluster_key: Annotated[str, "Column containing cluster labels"] = "leiden",
    min_score: Annotated[float, "Minimum score threshold for assignment"] = 0.1,
) -> Annotated[str, "Literature-based annotation results"]:
    """
    Annotate clusters using marker genes derived from literature search.
    Uses the markers obtained from search_pubmed_markers or get_cellmarker_database.
    
    Args:
        cell_type_markers: Dictionary like {"B cells": ["Cd19", "Cd79a"], "T cells": ["Cd3d", "Cd3e"]}
        cluster_key: Column with cluster assignments
        min_score: Minimum expression score for cell type assignment
    """
    adata = _get_adata()
    
    if cluster_key not in adata.obs:
        return f"Cluster column '{cluster_key}' not found in adata.obs"
    
    # Score each cluster for each cell type
    cluster_annotations = {}
    cluster_scores = {}
    
    for cluster in adata.obs[cluster_key].unique():
        cluster_mask = adata.obs[cluster_key] == cluster
        cluster_data = adata[cluster_mask]
        
        scores = {}
        for cell_type, markers in cell_type_markers.items():
            # Find markers present in data
            present_markers = [g for g in markers if g in adata.var_names]
            if not present_markers:
                scores[cell_type] = 0
                continue
            
            # Calculate mean expression
            if hasattr(cluster_data.X, 'toarray'):
                expr = cluster_data[:, present_markers].X.toarray()
            else:
                expr = cluster_data[:, present_markers].X
            
            # Score = mean expression * fraction of markers present
            coverage = len(present_markers) / len(markers)
            scores[cell_type] = np.mean(expr) * coverage
        
        # Get best cell type
        best_type = max(scores, key=scores.get)
        best_score = scores[best_type]
        
        if best_score >= min_score:
            cluster_annotations[cluster] = best_type
        else:
            cluster_annotations[cluster] = "Unknown"
        
        cluster_scores[cluster] = scores
    
    # Add to adata
    adata.obs["literature_annotation"] = adata.obs[cluster_key].map(cluster_annotations)
    _annotation_state["all_annotations"]["literature"] = cluster_annotations
    
    # Format results
    result = f"""
Literature-Based Annotation Complete
====================================
Cluster column: {cluster_key}
Min score threshold: {min_score}
Cell types used: {list(cell_type_markers.keys())}

Cluster Annotations:
"""
    for cluster in sorted(cluster_annotations.keys(), key=str):
        n_cells = (adata.obs[cluster_key] == cluster).sum()
        annotation = cluster_annotations[cluster]
        score = cluster_scores[cluster].get(annotation, 0)
        result += f"  Cluster {cluster}: {annotation} (score: {score:.3f}, n={n_cells:,})\n"
    
    # Summary
    annotation_counts = adata.obs["literature_annotation"].value_counts()
    result += f"""
Cell Type Distribution:
{annotation_counts.to_string()}
"""
    return result


@tool
def combine_pubmed_and_database_markers(
    cell_types: Annotated[List[str], "List of cell types to get markers for"],
    species: Annotated[str, "Species: 'mouse' or 'human'"] = "mouse",
    tissue: Annotated[Optional[str], "Tissue context"] = None,
) -> Annotated[str, "Combined marker gene dictionary"]:
    """
    Combine markers from PubMed search and CellMarker database.
    Creates a comprehensive marker set for annotation.

    This should be called after search_pubmed_markers to merge results.
    """
    # Get database markers
    db_markers = {
        "mouse": {
            "B cells": ["Cd19", "Cd79a", "Cd79b", "Ms4a1", "Pax5"],
            "Naive B cells": ["Cd19", "Ighd", "Fcer2a"],
            "Memory B cells": ["Cd19", "Cd27"],
            "Plasma cells": ["Sdc1", "Xbp1", "Prdm1"],
            "Immature B cells": ["Cd19", "Cd93", "Ighm"],
            "T cells": ["Cd3d", "Cd3e", "Trac"],
            "CD4 T cells": ["Cd4", "Cd3d"],
            "CD8 T cells": ["Cd8a", "Cd3d", "Gzmb"],
        },
        "human": {
            "B cells": ["CD19", "CD79A", "MS4A1", "PAX5"],
            "Naive B cells": ["CD19", "IGHD", "TCL1A"],
            "Memory B cells": ["CD19", "CD27"],
            "Plasma cells": ["SDC1", "XBP1", "PRDM1"],
            "T cells": ["CD3D", "CD3E", "TRAC"],
        }
    }

    species_db = db_markers.get(species.lower(), {})

    combined = {}
    for cell_type in cell_types:
        markers = set()

        # Add from database
        if cell_type in species_db:
            markers.update(species_db[cell_type])

        # Add from literature cache
        cache_key = f"{cell_type}_{species}_{tissue or 'any'}"
        if cache_key in _literature_markers_cache:
            # Add top 5 from literature
            lit_markers = _literature_markers_cache[cache_key].get("markers", [])[:5]
            markers.update(lit_markers)

        combined[cell_type] = list(markers)

    result = f"""
Combined Marker Dictionary
==========================
Species: {species}
Tissue: {tissue or 'any'}

Markers by Cell Type:
"""
    for cell_type, markers in combined.items():
        result += f"\n{cell_type}:\n  {', '.join(markers)}\n"

    # Store for later use
    _annotation_state["combined_markers"] = combined

    result += """
Use these markers with annotate_with_literature_markers() tool.
"""
    return result


# ============================================================
# Integrated Annotation Tools
# ============================================================

@tool
def integrate_multiple_annotations(
    methods: Annotated[List[str], "List of annotation column names to integrate"],
    strategy: Annotated[str, "Integration strategy: 'voting', 'confidence', or 'resolve'"] = "voting",
    weights: Annotated[Optional[Dict[str, float]], "Weights for each method (voting strategy)"] = None,
    confidence_cols: Annotated[Optional[Dict[str, str]], "Confidence columns for each method"] = None,
    output_col: Annotated[str, "Output column name"] = "integrated_annotation",
) -> Annotated[str, "Integration result"]:
    """
    Integrate multiple annotation methods to produce optimal cell type annotations.

    Combines CellTypist, literature-based, ScType, and other methods intelligently.

    Strategies:
    - 'voting': Weighted voting based on confidence
    - 'confidence': Use method with highest confidence per cell
    - 'resolve': Resolve conflicts between primary and secondary methods

    Args:
        methods: List of annotation columns (e.g., ['celltypist_label', 'literature_annotation', 'sctype_annotation'])
        strategy: Integration approach
        weights: Optional weights for voting strategy
        confidence_cols: Optional confidence columns for each method
        output_col: Name for integrated annotation column

    Returns:
        Summary of integration results
    """
    from scrna_agent.tools.annotation_integration import (
        integrate_annotations_voting,
        compare_annotation_methods,
        create_integrated_annotation_report
    )

    adata = _get_adata()

    # Validate methods exist
    available_methods = [m for m in methods if m in adata.obs.columns]
    missing_methods = [m for m in methods if m not in adata.obs.columns]

    if len(available_methods) < 2:
        return f"Error: Need at least 2 annotation methods. Available: {available_methods}, Missing: {missing_methods}"

    # Integrate annotations
    if strategy == "voting":
        integrated, confidence = integrate_annotations_voting(
            adata,
            available_methods,
            weights=weights,
            confidence_cols=confidence_cols,
            output_col=output_col
        )

        adata.obs[output_col] = integrated
        adata.obs[f"{output_col}_confidence"] = confidence

    elif strategy == "confidence":
        # Use method with highest confidence per cell
        if not confidence_cols:
            return "Error: 'confidence' strategy requires confidence_cols parameter"

        integrated = []
        max_confidences = []

        for idx in range(len(adata)):
            best_method = None
            best_conf = -1

            for method in available_methods:
                if method in confidence_cols and confidence_cols[method] in adata.obs.columns:
                    conf = adata.obs[confidence_cols[method]].iloc[idx]
                    if conf > best_conf:
                        best_conf = conf
                        best_method = method

            if best_method:
                integrated.append(adata.obs[best_method].iloc[idx])
                max_confidences.append(best_conf)
            else:
                integrated.append("Unknown")
                max_confidences.append(0.0)

        adata.obs[output_col] = integrated
        adata.obs[f"{output_col}_confidence"] = max_confidences

    else:
        return f"Error: Unknown strategy '{strategy}'. Use 'voting' or 'confidence'."

    # Generate report
    report = create_integrated_annotation_report(
        adata,
        available_methods,
        output_col
    )

    # Store in state
    _annotation_state["integrated_annotation"] = output_col

    return report


@tool
def compare_celltypist_with_literature(
    celltypist_col: Annotated[str, "CellTypist annotation column"] = "celltypist_label",
    literature_col: Annotated[str, "Literature annotation column"] = "literature_annotation",
    cluster_key: Annotated[str, "Cluster column"] = "leiden",
) -> Annotated[str, "Comparison report"]:
    """
    Compare CellTypist and literature-based annotation results.

    Provides detailed comparison including agreement rates, disagreements,
    and cluster-level consistency.

    Args:
        celltypist_col: Column with CellTypist predictions
        literature_col: Column with literature-based predictions
        cluster_key: Cluster column for cluster-level comparison

    Returns:
        Detailed comparison report
    """
    from scrna_agent.tools.annotation_integration import compare_annotation_methods

    adata = _get_adata()

    # Check columns exist
    if celltypist_col not in adata.obs.columns:
        return f"Error: CellTypist column '{celltypist_col}' not found"
    if literature_col not in adata.obs.columns:
        return f"Error: Literature column '{literature_col}' not found"

    # Run comparison
    comparison = compare_annotation_methods(
        adata,
        [celltypist_col, literature_col],
        cluster_key=cluster_key
    )

    # Format report
    result_lines = []
    result_lines.append("=" * 70)
    result_lines.append("CELLTYPIST vs LITERATURE ANNOTATION COMPARISON")
    result_lines.append("=" * 70)

    # Overall agreement
    pair_key = f"{celltypist_col} vs {literature_col}"
    if pair_key in comparison["agreement"]:
        agreement = comparison["agreement"][pair_key]
        result_lines.append(f"\nOverall Agreement: {agreement:.1%}")

    # Disagreements
    if pair_key in comparison["disagreement"]:
        dis_info = comparison["disagreement"][pair_key]
        result_lines.append(f"Disagreement Count: {dis_info['count']:,} cells ({dis_info['percentage']:.1f}%)")

    # Cluster-level consistency
    if comparison["cluster_consistency"]:
        result_lines.append(f"\nCluster-Level Consistency:")
        result_lines.append("-" * 70)

        for cluster, cluster_info in list(comparison["cluster_consistency"].items())[:10]:
            result_lines.append(f"\nCluster {cluster}:")

            if celltypist_col in cluster_info:
                ct_info = cluster_info[celltypist_col]
                result_lines.append(
                    f"  CellTypist: {ct_info['annotation']} "
                    f"({ct_info['consistency']:.1f}% consistent, n={ct_info['count']:,})"
                )

            if literature_col in cluster_info:
                lit_info = cluster_info[literature_col]
                result_lines.append(
                    f"  Literature: {lit_info['annotation']} "
                    f"({lit_info['consistency']:.1f}% consistent, n={lit_info['count']:,})"
                )

            # Agreement in this cluster
            if celltypist_col in cluster_info and literature_col in cluster_info:
                if cluster_info[celltypist_col]['annotation'] == cluster_info[literature_col]['annotation']:
                    result_lines.append("  ✓ AGREEMENT")
                else:
                    result_lines.append("  ✗ DISAGREEMENT")

    result_lines.append("\n" + "=" * 70)

    return "\n".join(result_lines)


@tool
def create_optimal_annotation(
    celltypist_col: Annotated[str, "CellTypist annotation column"] = "celltypist_label",
    celltypist_conf_col: Annotated[str, "CellTypist confidence column"] = "celltypist_confidence",
    literature_col: Annotated[Optional[str], "Literature annotation column"] = None,
    sctype_col: Annotated[Optional[str], "ScType annotation column"] = None,
    confidence_threshold: Annotated[float, "Confidence threshold for CellTypist"] = 0.5,
    output_col: Annotated[str, "Output column name"] = "optimal_annotation",
) -> Annotated[str, "Optimal annotation result"]:
    """
    Create optimal cell type annotation by intelligently combining multiple methods.

    Strategy:
    1. Use CellTypist if confidence > threshold
    2. If CellTypist has low confidence, check agreement with other methods
    3. If methods agree, use consensus
    4. If methods disagree, use voting with confidence weighting

    Args:
        celltypist_col: CellTypist prediction column
        celltypist_conf_col: CellTypist confidence column
        literature_col: Optional literature-based annotation
        sctype_col: Optional ScType annotation
        confidence_threshold: Threshold for trusting CellTypist alone
        output_col: Output column name

    Returns:
        Summary of optimal annotation
    """
    from scrna_agent.tools.annotation_integration import (
        integrate_annotations_voting,
        create_integrated_annotation_report
    )

    adata = _get_adata()

    # Check required columns
    if celltypist_col not in adata.obs.columns:
        return f"Error: CellTypist column '{celltypist_col}' not found"
    if celltypist_conf_col not in adata.obs.columns:
        return f"Error: Confidence column '{celltypist_conf_col}' not found"

    # Collect available methods
    methods = [celltypist_col]
    confidence_cols = {celltypist_col: celltypist_conf_col}

    if literature_col and literature_col in adata.obs.columns:
        methods.append(literature_col)
    if sctype_col and sctype_col in adata.obs.columns:
        methods.append(sctype_col)

    if len(methods) == 1:
        # Only CellTypist available
        adata.obs[output_col] = adata.obs[celltypist_col]
        adata.obs[f"{output_col}_confidence"] = adata.obs[celltypist_conf_col]
        adata.obs[f"{output_col}_method"] = "celltypist_only"

        summary = f"""
Optimal Annotation Created (CellTypist Only)
============================================
Method: CellTypist only (no other methods available)
Output column: {output_col}

Cell Type Distribution:
{adata.obs[output_col].value_counts().head(15).to_string()}
"""
        return summary

    # Multiple methods available - integrate
    optimal_annotations = []
    optimal_confidences = []
    methods_used = []

    for idx in range(len(adata)):
        ct_annotation = adata.obs[celltypist_col].iloc[idx]
        ct_confidence = adata.obs[celltypist_conf_col].iloc[idx]

        # High confidence in CellTypist
        if ct_confidence >= confidence_threshold:
            optimal_annotations.append(ct_annotation)
            optimal_confidences.append(ct_confidence)
            methods_used.append("celltypist")
            continue

        # Low confidence - check other methods
        annotations_this_cell = [ct_annotation]
        confidences_this_cell = [ct_confidence * 0.5]  # Reduce weight for low confidence

        if literature_col and literature_col in adata.obs.columns:
            lit_annotation = adata.obs[literature_col].iloc[idx]
            annotations_this_cell.append(lit_annotation)
            confidences_this_cell.append(0.8)  # Literature gets fixed confidence

        if sctype_col and sctype_col in adata.obs.columns:
            st_annotation = adata.obs[sctype_col].iloc[idx]
            annotations_this_cell.append(st_annotation)
            confidences_this_cell.append(0.7)  # ScType gets fixed confidence

        # Check for consensus
        if len(set(annotations_this_cell)) == 1:
            # All methods agree
            optimal_annotations.append(annotations_this_cell[0])
            optimal_confidences.append(0.9)  # High confidence for consensus
            methods_used.append("consensus")
        else:
            # Weighted voting
            vote_counts = {}
            for ann, conf in zip(annotations_this_cell, confidences_this_cell):
                if ann not in vote_counts:
                    vote_counts[ann] = 0
                vote_counts[ann] += conf

            best_annotation = max(vote_counts, key=vote_counts.get)
            best_confidence = vote_counts[best_annotation] / sum(confidences_this_cell)

            optimal_annotations.append(best_annotation)
            optimal_confidences.append(best_confidence)
            methods_used.append("voting")

    # Add to adata
    adata.obs[output_col] = optimal_annotations
    adata.obs[f"{output_col}_confidence"] = optimal_confidences
    adata.obs[f"{output_col}_method"] = methods_used

    # Generate report
    report_lines = []
    report_lines.append("=" * 70)
    report_lines.append("OPTIMAL ANNOTATION CREATED")
    report_lines.append("=" * 70)
    report_lines.append(f"\nMethods integrated: {', '.join(methods)}")
    report_lines.append(f"Confidence threshold: {confidence_threshold}")
    report_lines.append(f"Output column: {output_col}")

    # Method usage statistics
    method_counts = pd.Series(methods_used).value_counts()
    report_lines.append(f"\nMethod Usage:")
    for method, count in method_counts.items():
        percentage = (count / len(adata)) * 100
        report_lines.append(f"  {method}: {count:,} cells ({percentage:.1f}%)")

    # Cell type distribution
    optimal_counts = pd.Series(optimal_annotations).value_counts()
    report_lines.append(f"\nTop 15 Cell Types:")
    for cell_type, count in optimal_counts.head(15).items():
        percentage = (count / len(adata)) * 100
        report_lines.append(f"  {cell_type}: {count:,} ({percentage:.1f}%)")

    # Confidence statistics
    report_lines.append(f"\nConfidence Statistics:")
    report_lines.append(f"  Mean: {pd.Series(optimal_confidences).mean():.3f}")
    report_lines.append(f"  Median: {pd.Series(optimal_confidences).median():.3f}")
    report_lines.append(f"  Low confidence (<0.5): {sum([c < 0.5 for c in optimal_confidences]):,} cells")

    report_lines.append("\n" + "=" * 70)

    return "\n".join(report_lines)


# ============================================================
# Biological Validation Tools (NEW - from JMV_ALI session)
# ============================================================
# Key insight: Organoid samples were incorrectly annotated as 80.6% neutrophils
# because S100A8/S100A9 stress markers were misinterpreted as immune markers.
# These tools add biological context awareness to prevent such errors.

@tool
def validate_annotation_biology(
    sample_type: Annotated[str, "Sample type: 'organoid', 'primary_tissue', 'cell_line', 'pbmc'"],
    celltype_col: Annotated[str, "Column containing cell type annotations"] = "celltype",
) -> Annotated[str, "Biological validation result"]:
    """
    Validate cell type annotations against biological context constraints.

    CRITICAL: Organoid samples should NOT contain immune cells (neutrophils, monocytes, etc.)
    This tool detects biologically impossible annotations based on sample type.

    Sample type constraints:
    - organoid: Only epithelial/stem cells (immune cells = misannotation)
    - primary_tissue: All cell types allowed
    - cell_line: Homogeneous population expected
    - pbmc: Only immune cells expected

    Example usage:
        validate_annotation_biology(sample_type="organoid")
        # Will flag if >5% of cells are annotated as immune cells
    """
    try:
        from scrna_agent.tools.biological_validator import BiologicalValidator
    except ImportError:
        return "Error: BiologicalValidator not found. Check tools/biological_validator.py"

    adata = _get_adata()

    if celltype_col not in adata.obs.columns:
        return f"Error: Cell type column '{celltype_col}' not found in adata.obs"

    validator = BiologicalValidator()
    validation = validator.validate_annotation(adata, sample_type, celltype_col=celltype_col)

    # Format result
    result_lines = []
    result_lines.append("=" * 70)
    result_lines.append("BIOLOGICAL VALIDATION RESULT")
    result_lines.append("=" * 70)
    result_lines.append(f"\nSample type: {sample_type}")
    result_lines.append(f"Cell type column: {celltype_col}")
    result_lines.append(f"Validation result: {'PASSED' if validation['valid'] else 'FAILED'}")
    result_lines.append(f"\nSummary: {validation['summary']}")

    if validation['issues']:
        result_lines.append(f"\nCRITICAL ISSUES ({len(validation['issues'])}):")
        for issue in validation['issues']:
            result_lines.append(f"  - {issue['message']}")

    if validation['warnings']:
        result_lines.append(f"\nWarnings ({len(validation['warnings'])}):")
        for warning in validation['warnings']:
            result_lines.append(f"  - {warning}")

    # Show current cell type distribution
    celltype_counts = adata.obs[celltype_col].value_counts()
    result_lines.append(f"\nCurrent cell type distribution:")
    for ct, count in celltype_counts.head(10).items():
        pct = count / len(adata) * 100
        result_lines.append(f"  {ct}: {count:,} ({pct:.1f}%)")

    result_lines.append("\n" + "=" * 70)

    return "\n".join(result_lines)


@tool
def check_ambiguous_marker(
    marker: Annotated[str, "Gene symbol to check (e.g., 'S100A8')"],
    sample_type: Annotated[str, "Sample type for context"],
) -> Annotated[str, "Marker interpretation based on context"]:
    """
    Check if a marker gene has context-dependent interpretation.

    Key example: S100A8/S100A9 are commonly interpreted as neutrophil markers,
    but in organoid samples they indicate epithelial stress response, NOT immune cells.

    This tool helps avoid misannotation by providing context-aware interpretation.

    Example:
        check_ambiguous_marker(marker="S100A8", sample_type="organoid")
        # Returns: "In organoids, S100A8 indicates stress response, not neutrophils"
    """
    try:
        from scrna_agent.tools.biological_validator import BiologicalValidator
    except ImportError:
        return "Error: BiologicalValidator not found"

    validator = BiologicalValidator()
    result = validator.check_ambiguous_markers(marker, sample_type)

    if result['is_ambiguous']:
        return f"""
Ambiguous Marker Alert: {marker}
================================
Common interpretation: {result['common_interpretation']}
Context ({sample_type}): {result['context_interpretation']}
Recommendation: {result['recommendation']}

WARNING: {result['warning']}
"""
    else:
        return f"Marker '{marker}' is not in the ambiguous markers database."


@tool
def get_recommended_celltypist_model(
    tissue_type: Annotated[str, "Tissue type: 'lung_airway', 'pbmc', 'gut', etc."],
) -> Annotated[str, "Recommended CellTypist model"]:
    """
    Get the recommended CellTypist model for a specific tissue type.

    Proper model selection is critical for accurate annotation.
    For example, lung airway organoids should use 'Cells_Lung_Airway.pkl' (78 cell types).

    Available tissue types:
    - lung_airway: Cells_Lung_Airway.pkl (78 types)
    - lung_immune: Immune_All_Low.pkl (98 types)
    - pbmc: Immune_All_High.pkl (32 types)
    - gut: Cells_Intestinal_Tract.pkl (50 types)
    """
    try:
        from scrna_agent.tools.biological_validator import BiologicalValidator
    except ImportError:
        return "Error: BiologicalValidator not found"

    validator = BiologicalValidator()
    model = validator.get_recommended_celltypist_model(tissue_type)

    return f"""
Recommended CellTypist Model
============================
Tissue type: {tissue_type}
Model: {model}

Use this model with celltypist_annotate(model="{model}")
"""


@tool
def run_multi_method_annotation(
    sample_type: Annotated[str, "Sample type: 'organoid', 'primary_tissue', 'cell_line', 'pbmc'"],
    celltypist_model: Annotated[Optional[str], "CellTypist model (auto-selected if None)"] = None,
    output_col: Annotated[str, "Output column for consensus annotation"] = "consensus_celltype",
) -> Annotated[str, "Multi-method annotation result"]:
    """
    Run multi-method cell type annotation with biological validation.

    Combines:
    1. CellTypist (ML-based): Pre-trained deep learning model
    2. Marker-based: Canonical marker gene expression scoring
    3. Consensus: Context-aware decision logic

    For organoid samples, immune cell annotations from CellTypist are automatically
    overridden with marker-based epithelial annotations.

    This approach solved the JMV_ALI issue where 80.6% cells were incorrectly
    annotated as neutrophils. After multi-method annotation:
    - Basal: 26.2%
    - General_Epithelial: 62.1%
    - Ciliated: 5.7%
    - Secretory_Club: 5.7%

    Example:
        run_multi_method_annotation(sample_type="organoid")
    """
    try:
        from scrna_agent.tools.multi_method_annotation import MultiMethodAnnotator
        from scrna_agent.tools.biological_validator import BiologicalValidator
    except ImportError as e:
        return f"Error importing modules: {e}"

    adata = _get_adata()

    result_lines = []
    result_lines.append("=" * 70)
    result_lines.append("MULTI-METHOD ANNOTATION WITH BIOLOGICAL VALIDATION")
    result_lines.append("=" * 70)
    result_lines.append(f"\nSample type: {sample_type}")
    result_lines.append(f"Cells: {adata.n_obs:,}")

    # Initialize annotator
    annotator = MultiMethodAnnotator(sample_type=sample_type)
    validator = BiologicalValidator()

    # Auto-select model if not provided
    if celltypist_model is None:
        # Infer tissue type from sample type
        if sample_type == 'organoid':
            tissue_type = 'lung_airway'  # Default for organoids
        elif sample_type == 'pbmc':
            tissue_type = 'pbmc'
        else:
            tissue_type = 'lung_airway'
        celltypist_model = validator.get_recommended_celltypist_model(tissue_type)

    result_lines.append(f"CellTypist model: {celltypist_model}")

    try:
        # Run CellTypist
        result_lines.append("\n--- Running CellTypist ---")
        adata_ct = annotator.run_celltypist(adata, model_name=celltypist_model)

        ct_counts = adata_ct.obs['majority_voting'].value_counts()
        result_lines.append(f"CellTypist results ({len(ct_counts)} types):")
        for ct, count in ct_counts.head(5).items():
            pct = count / len(adata_ct) * 100
            result_lines.append(f"  {ct}: {count:,} ({pct:.1f}%)")

        # Run marker-based
        result_lines.append("\n--- Running Marker-Based ---")
        marker_types, marker_conf = annotator.run_marker_based(adata)

        mb_counts = marker_types.value_counts()
        result_lines.append(f"Marker-based results ({len(mb_counts)} types):")
        for ct, count in mb_counts.items():
            pct = count / len(adata) * 100
            result_lines.append(f"  {ct}: {count:,} ({pct:.1f}%)")

        # Consensus
        result_lines.append("\n--- Creating Consensus ---")
        consensus, conf, agreement = annotator.consensus_annotation(
            adata, adata_ct, marker_types, marker_conf
        )

        # Add to adata
        adata.obs['celltypist'] = adata_ct.obs['majority_voting'].values
        adata.obs['celltypist_conf'] = adata_ct.obs['conf_score'].values
        adata.obs['marker_based'] = marker_types.values
        adata.obs['marker_conf'] = marker_conf.values
        adata.obs[output_col] = consensus
        adata.obs[f'{output_col}_confidence'] = conf
        adata.obs[f'{output_col}_agreement'] = agreement

        # Consensus summary
        consensus_counts = pd.Series(consensus).value_counts()
        result_lines.append(f"\nConsensus results ({len(consensus_counts)} types):")
        for ct, count in consensus_counts.items():
            pct = count / len(adata) * 100
            result_lines.append(f"  {ct}: {count:,} ({pct:.1f}%)")

        # Agreement summary
        agree_counts = pd.Series(agreement).value_counts()
        result_lines.append(f"\nDecision method distribution:")
        for method, count in agree_counts.items():
            pct = count / len(adata) * 100
            result_lines.append(f"  {method}: {count:,} ({pct:.1f}%)")

        # Validate
        result_lines.append("\n--- Biological Validation ---")
        adata.obs['celltype'] = consensus  # Set for validation
        validation = validator.validate_annotation(adata, sample_type)
        result_lines.append(f"Validation: {'PASSED' if validation['valid'] else 'FAILED'}")
        result_lines.append(f"Summary: {validation['summary']}")

        if validation['issues']:
            for issue in validation['issues']:
                result_lines.append(f"  WARNING: {issue['message']}")

    except Exception as e:
        result_lines.append(f"\nError during annotation: {str(e)}")
        import traceback
        result_lines.append(traceback.format_exc())

    result_lines.append("\n" + "=" * 70)
    result_lines.append(f"Results saved to columns: {output_col}, celltypist, marker_based")
    result_lines.append("=" * 70)

    return "\n".join(result_lines)


@tool
def list_sample_types() -> Annotated[str, "Available sample types and constraints"]:
    """
    List all available sample types and their biological constraints.

    Use this to understand what cell types are expected/forbidden for each sample type.
    """
    try:
        from scrna_agent.tools.biological_validator import BiologicalValidator
    except ImportError:
        return "Error: BiologicalValidator not found"

    validator = BiologicalValidator()
    sample_types = validator.list_sample_types()

    result_lines = []
    result_lines.append("Available Sample Types")
    result_lines.append("=" * 50)

    for st in sample_types:
        info = validator.get_sample_type_info(st)
        result_lines.append(f"\n{st}:")
        result_lines.append(f"  Description: {info.get('description', 'N/A')}")
        result_lines.append(f"  Allowed: {info.get('allowed_cell_types', [])}")
        result_lines.append(f"  Forbidden: {info.get('forbidden_cell_types', [])}")
        if info.get('validation_warning'):
            result_lines.append(f"  Warning: {info['validation_warning']}")

    return "\n".join(result_lines)

