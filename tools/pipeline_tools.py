# Pipeline Tools
# Shared state management and scANVI fine-tuning for end-to-end pipeline

import os
from typing import Annotated, Optional, List, Dict, Any
import numpy as np
import pandas as pd

## 3rd party
from langchain_core.tools import tool

# ============================================================
# Global Pipeline State (Shared across annotation & integration)
# ============================================================
_pipeline_state: Dict[str, Any] = {
    "adata": None,
    "original_path": None,
    "batch_key": None,
    "labels_key": None,
    "species": None,
    # Annotation results
    "annotation_method": None,
    "annotation_col": None,
    # Integration results
    "integration_method": None,
    "scvi_model": None,
    "scanvi_model": None,
    # Fine-tuning
    "scanvi_finetuned_model": None,
    # Workflow tracking
    "completed_steps": [],
}


def get_pipeline_state() -> Dict[str, Any]:
    """Get current pipeline state."""
    return _pipeline_state


def set_pipeline_state(key: str, value: Any):
    """Set a value in pipeline state."""
    _pipeline_state[key] = value


def _get_adata():
    """Get the current AnnData object from pipeline state."""
    if _pipeline_state["adata"] is None:
        raise ValueError("No data loaded. Use load_data_for_pipeline first.")
    return _pipeline_state["adata"]


# ============================================================
# Data Loading for Pipeline
# ============================================================
@tool
def load_data_for_pipeline(
    file_path: Annotated[str, "Path to the h5ad file"],
    batch_key: Annotated[str, "Column containing batch/sample information"],
    species: Annotated[str, "Species: 'mouse' or 'human'"] = "mouse",
) -> Annotated[str, "Summary of loaded data"]:
    """
    Load h5ad file and initialize the pipeline state.
    This is the entry point for the integrated workflow.
    Data should already be QC'd and normalized.
    """
    import scanpy as sc
    
    adata = sc.read_h5ad(file_path)
    
    # Store in state
    _pipeline_state["adata"] = adata
    _pipeline_state["original_path"] = file_path
    _pipeline_state["batch_key"] = batch_key
    _pipeline_state["species"] = species
    _pipeline_state["completed_steps"] = ["load_data"]
    
    # Check data properties
    has_raw = adata.raw is not None
    has_counts = "counts" in adata.layers
    has_umap = "X_umap" in adata.obsm
    has_pca = "X_pca" in adata.obsm
    has_clusters = "leiden" in adata.obs or "louvain" in adata.obs
    
    # Check batch key
    if batch_key in adata.obs:
        batch_counts = adata.obs[batch_key].value_counts()
        n_batches = len(batch_counts)
    else:
        n_batches = 0
        batch_counts = "Batch key not found!"
    
    summary = f"""
Pipeline Data Loaded
====================
File: {file_path}
Species: {species}
Batch key: {batch_key}

Data Overview:
- Cells: {adata.n_obs:,}
- Genes: {adata.n_vars:,}
- Has .raw: {has_raw}
- Has counts layer: {has_counts}
- Has PCA: {has_pca}
- Has UMAP: {has_umap}
- Has clusters: {has_clusters}
- Number of batches: {n_batches}

Batch Distribution:
{batch_counts if isinstance(batch_counts, str) else batch_counts.head(10).to_string()}

Obs columns: {list(adata.obs.columns)[:15]}...

Next Step: Run annotation (CellTypist recommended for {species} immune cells)
"""
    return summary


# ============================================================
# scANVI Fine-Tuning Tools
# ============================================================
@tool
def setup_scanvi_finetuning(
    annotation_col: Annotated[str, "Column containing cell type annotations to use as labels"],
    unknown_labels: Annotated[Optional[List[str]], "Labels to treat as unknown (e.g., ['Unknown', 'Unassigned'])"] = None,
) -> Annotated[str, "Setup result"]:
    """
    Setup scANVI fine-tuning using annotations from CellTypist/ScType as labels.
    This refines integration by using cell type labels from annotation step.
    
    Workflow: Initial Annotation → scANVI fine-tuning → Refined annotation + integration
    """
    import scvi
    
    adata = _get_adata()
    
    if annotation_col not in adata.obs:
        available = [c for c in adata.obs.columns if 
                    'annotation' in c.lower() or 'prediction' in c.lower() or 
                    'celltypist' in c.lower() or 'sctype' in c.lower()]
        return f"Annotation column '{annotation_col}' not found. Available annotation columns: {available}"
    
    # Mark low-confidence or unknown cells
    if unknown_labels is None:
        unknown_labels = ['Unknown', 'Unassigned', 'Low_confidence']
    
    # Create labels column for scANVI
    labels = adata.obs[annotation_col].copy()
    
    # Check for confidence column
    conf_col = None
    for possible_conf in ['celltypist_confidence', 'scanvi_confidence', 'conf_score']:
        if possible_conf in adata.obs:
            conf_col = possible_conf
            break
    
    # Mark low confidence cells as Unknown (if confidence available)
    if conf_col:
        low_conf_mask = adata.obs[conf_col] < 0.5
        labels[low_conf_mask] = "Unknown"
        n_low_conf = low_conf_mask.sum()
    else:
        n_low_conf = 0
    
    # Mark specified unknown labels
    for ul in unknown_labels:
        labels[labels == ul] = "Unknown"
    
    adata.obs["scanvi_labels"] = labels
    _pipeline_state["labels_key"] = "scanvi_labels"
    
    label_counts = labels.value_counts()
    n_unknown = (labels == "Unknown").sum()
    n_known = len(labels) - n_unknown
    
    return f"""
scANVI Fine-Tuning Setup Complete
=================================
Source annotation: {annotation_col}
Labels column created: scanvi_labels
Confidence column used: {conf_col}

Label Statistics:
- Known labels: {n_known:,} ({n_known/len(labels)*100:.1f}%)
- Unknown (to refine): {n_unknown:,} ({n_unknown/len(labels)*100:.1f}%)
- Low confidence flagged: {n_low_conf:,}

Label Distribution:
{label_counts.head(15).to_string()}

Next Step: Run scANVI fine-tuning with run_scanvi_finetuning()
"""


@tool
def run_scanvi_finetuning(
    n_latent: Annotated[int, "Latent dimension"] = 30,
    n_layers: Annotated[int, "Number of neural network layers"] = 2,
    max_epochs_scvi: Annotated[int, "Max epochs for scVI pre-training"] = 400,
    max_epochs_scanvi: Annotated[int, "Max epochs for scANVI training"] = 100,
    layer: Annotated[str, "Layer with raw counts"] = "counts",
) -> Annotated[str, "scANVI fine-tuning results"]:
    """
    Run scANVI fine-tuning: scVI → scANVI with annotation labels.
    This simultaneously integrates data AND refines cell type predictions.
    
    Uses labels from setup_scanvi_finetuning().
    """
    import scvi
    import scanpy as sc
    
    adata = _get_adata()
    batch_key = _pipeline_state.get("batch_key")
    labels_key = _pipeline_state.get("labels_key", "scanvi_labels")
    
    if batch_key is None:
        return "Error: batch_key not set. Run load_data_for_pipeline first."
    
    if labels_key not in adata.obs:
        return f"Error: Labels column '{labels_key}' not found. Run setup_scanvi_finetuning first."
    
    # Check for counts
    if layer not in adata.layers:
        if adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()
            layer = "counts"
        else:
            return "Error: No counts layer found. Need raw counts for scVI/scANVI."
    
    # Setup anndata
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
        labels_key=labels_key,
    )
    
    # ======== Phase 1: scVI Pre-training ========
    scvi_model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood="nb",
    )
    
    scvi_model.train(
        max_epochs=max_epochs_scvi,
        early_stopping=True,
        early_stopping_patience=20,
    )
    
    _pipeline_state["scvi_model"] = scvi_model
    
    # Get scVI latent
    adata.obsm["X_scVI"] = scvi_model.get_latent_representation()
    
    # ======== Phase 2: scANVI Fine-tuning ========
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_key,
        unlabeled_category="Unknown",
    )
    
    scanvi_model.train(
        max_epochs=max_epochs_scanvi,
        early_stopping=True,
    )
    
    _pipeline_state["scanvi_model"] = scanvi_model
    _pipeline_state["scanvi_finetuned_model"] = scanvi_model
    
    # Get scANVI latent
    latent = scanvi_model.get_latent_representation()
    adata.obsm["X_scANVI"] = latent
    
    # Get refined predictions
    predictions = scanvi_model.predict()
    adata.obs["scanvi_refined"] = predictions
    
    # Get prediction probabilities
    soft_preds = scanvi_model.predict(soft=True)
    adata.obs["scanvi_refined_confidence"] = soft_preds.max(axis=1)
    
    # Track step
    _pipeline_state["completed_steps"].append("scanvi_finetuning")
    _pipeline_state["integration_method"] = "scANVI"
    _pipeline_state["annotation_col"] = "scanvi_refined"
    
    # Summary statistics
    pred_counts = adata.obs["scanvi_refined"].value_counts()
    mean_conf = adata.obs["scanvi_refined_confidence"].mean()
    
    # Compare with original annotation
    original_labels = adata.obs[labels_key]
    known_mask = original_labels != "Unknown"
    if known_mask.sum() > 0:
        agreement = (adata.obs.loc[known_mask, "scanvi_refined"] == original_labels[known_mask]).mean()
    else:
        agreement = 0.0
    
    return f"""
scANVI Fine-Tuning Complete
===========================
Parameters:
- Batch key: {batch_key}
- Labels key: {labels_key}
- Latent dimension: {n_latent}
- Layers: {n_layers}
- scVI epochs: {max_epochs_scvi}
- scANVI epochs: {max_epochs_scanvi}

Integration Results:
- scVI latent: adata.obsm['X_scVI'] {adata.obsm['X_scVI'].shape}
- scANVI latent: adata.obsm['X_scANVI'] {adata.obsm['X_scANVI'].shape}

Refined Annotation Results:
- Predictions stored in: adata.obs['scanvi_refined']
- Confidence stored in: adata.obs['scanvi_refined_confidence']
- Mean confidence: {mean_conf:.3f}
- Agreement with original labels: {agreement:.1%}

Cell Type Distribution (Refined):
{pred_counts.head(15).to_string()}

Next Step: Compute UMAP using X_scANVI embedding and validate results
"""


@tool
def compute_integrated_umap(
    use_rep: Annotated[str, "Embedding to use: 'X_scVI', 'X_scANVI', or 'X_harmony'"] = "X_scANVI",
    n_neighbors: Annotated[int, "Number of neighbors"] = 15,
) -> Annotated[str, "UMAP computation result"]:
    """Compute neighbors and UMAP from integrated embedding."""
    import scanpy as sc
    
    adata = _get_adata()
    
    if use_rep not in adata.obsm:
        available = [k for k in adata.obsm.keys() if k.startswith("X_")]
        return f"Embedding '{use_rep}' not found. Available: {available}"
    
    # Compute neighbors
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
    
    # Compute UMAP
    sc.tl.umap(adata)
    
    _pipeline_state["completed_steps"].append("umap")
    
    return f"""
UMAP Computed
=============
Using: {use_rep}
Neighbors: {n_neighbors}
UMAP shape: {adata.obsm['X_umap'].shape}

Stored in: adata.obsm['X_umap']

Next Step: Plot results and save final data
"""


@tool
def plot_pipeline_results(
    output_dir: Annotated[str, "Directory to save plots"],
    annotation_col: Annotated[str, "Annotation column to plot"] = "scanvi_refined",
) -> Annotated[str, "Paths to saved plots"]:
    """Generate summary plots for the pipeline results."""
    import scanpy as sc
    import matplotlib.pyplot as plt
    import os
    
    adata = _get_adata()
    batch_key = _pipeline_state.get("batch_key", "batch")
    
    os.makedirs(output_dir, exist_ok=True)
    
    plots_saved = []
    
    # 1. UMAP colored by batch
    fig, ax = plt.subplots(figsize=(10, 8))
    if batch_key in adata.obs:
        sc.pl.umap(adata, color=batch_key, ax=ax, show=False, title="Batches (Integration Check)")
    plt.tight_layout()
    batch_plot = os.path.join(output_dir, "umap_batch.png")
    plt.savefig(batch_plot, dpi=150, bbox_inches="tight")
    plt.close()
    plots_saved.append(batch_plot)
    
    # 2. UMAP colored by cell type
    fig, ax = plt.subplots(figsize=(12, 10))
    if annotation_col in adata.obs:
        sc.pl.umap(adata, color=annotation_col, ax=ax, show=False, 
                   title="Cell Type Annotation", legend_loc="on data", legend_fontsize=6)
    plt.tight_layout()
    celltype_plot = os.path.join(output_dir, "umap_celltype.png")
    plt.savefig(celltype_plot, dpi=150, bbox_inches="tight")
    plt.close()
    plots_saved.append(celltype_plot)
    
    # 3. Confidence distribution
    conf_col = "scanvi_refined_confidence" if "scanvi_refined_confidence" in adata.obs else None
    if conf_col:
        fig, ax = plt.subplots(figsize=(8, 5))
        adata.obs[conf_col].hist(bins=50, ax=ax, edgecolor='black')
        ax.set_xlabel("Confidence Score")
        ax.set_ylabel("Number of Cells")
        ax.set_title("Annotation Confidence Distribution")
        ax.axvline(0.5, color='red', linestyle='--', label='Threshold (0.5)')
        ax.legend()
        plt.tight_layout()
        conf_plot = os.path.join(output_dir, "confidence_distribution.png")
        plt.savefig(conf_plot, dpi=150, bbox_inches="tight")
        plt.close()
        plots_saved.append(conf_plot)
    
    return f"""
Pipeline Plots Saved
====================
Directory: {output_dir}
Plots:
{chr(10).join(['- ' + p for p in plots_saved])}

Completed pipeline steps: {_pipeline_state.get('completed_steps', [])}
"""


@tool
def save_pipeline_results(
    output_path: Annotated[str, "Path to save the final h5ad file"],
    model_dir: Annotated[Optional[str], "Directory to save scVI/scANVI models"] = None,
) -> Annotated[str, "Save result"]:
    """Save final results: annotated and integrated AnnData + models."""
    import os
    
    adata = _get_adata()
    adata.write_h5ad(output_path)
    
    saved_items = [f"H5AD: {output_path}"]
    
    # Save models
    if model_dir:
        os.makedirs(model_dir, exist_ok=True)
        
        scvi_model = _pipeline_state.get("scvi_model")
        if scvi_model:
            scvi_path = os.path.join(model_dir, "scvi_model")
            scvi_model.save(scvi_path, overwrite=True)
            saved_items.append(f"scVI model: {scvi_path}")
        
        scanvi_model = _pipeline_state.get("scanvi_model")
        if scanvi_model:
            scanvi_path = os.path.join(model_dir, "scanvi_model")
            scanvi_model.save(scanvi_path, overwrite=True)
            saved_items.append(f"scANVI model: {scanvi_path}")
    
    # Save annotations CSV
    csv_path = output_path.replace(".h5ad", "_annotations.csv")
    annotation_cols = [c for c in adata.obs.columns if 
                      any(k in c.lower() for k in ['annotation', 'prediction', 'celltypist', 
                                                     'sctype', 'scanvi', 'refined', 'consensus'])]
    if annotation_cols:
        adata.obs[annotation_cols].to_csv(csv_path)
        saved_items.append(f"Annotations CSV: {csv_path}")
    
    return f"""
Pipeline Results Saved
======================
{chr(10).join(saved_items)}

Summary:
- Cells: {adata.n_obs:,}
- Genes: {adata.n_vars:,}
- Embeddings: {[k for k in adata.obsm.keys() if k.startswith('X_')]}
- Annotation columns: {annotation_cols}
- Completed steps: {_pipeline_state.get('completed_steps', [])}
"""


@tool
def get_pipeline_status() -> Annotated[str, "Current pipeline status"]:
    """Get the current status of the pipeline."""
    adata = _pipeline_state.get("adata")
    
    if adata is None:
        return "Pipeline not started. Use load_data_for_pipeline to begin."
    
    status = f"""
Pipeline Status
===============
Data: {_pipeline_state.get('original_path', 'Unknown')}
Species: {_pipeline_state.get('species', 'Unknown')}
Batch key: {_pipeline_state.get('batch_key', 'Not set')}
Labels key: {_pipeline_state.get('labels_key', 'Not set')}

Data Shape: {adata.n_obs:,} cells × {adata.n_vars:,} genes

Completed Steps:
{chr(10).join(['  ✓ ' + s for s in _pipeline_state.get('completed_steps', [])])}

Available Embeddings:
{[k for k in adata.obsm.keys() if k.startswith('X_')]}

Models Trained:
- scVI: {'✓' if _pipeline_state.get('scvi_model') else '✗'}
- scANVI: {'✓' if _pipeline_state.get('scanvi_model') else '✗'}

Integration method: {_pipeline_state.get('integration_method', 'Not set')}
Annotation column: {_pipeline_state.get('annotation_col', 'Not set')}
"""
    return status


# ============================================================
# Comprehensive Visualization and Reporting Orchestration
# ============================================================

@tool
def generate_qc_visualizations(
    adata_before_path: Annotated[str, "Path to h5ad file before QC"],
    adata_after_path: Annotated[str, "Path to h5ad file after QC"],
    output_dir: Annotated[str, "Output directory for QC plots"] = "./qc_plots",
) -> Annotated[str, "Summary of QC visualizations generated"]:
    """
    Generate comprehensive QC visualizations including before/after comparison.

    This orchestrates:
    - QC metrics before filtering
    - QC metrics after filtering
    - Before/after comparison plots

    Returns paths to all generated plots.
    """
    from scrna_agent.tools.qc_visualization import (
        plot_qc_metrics_before,
        plot_qc_metrics_after,
        plot_qc_comparison,
        get_qc_plot_paths,
    )

    results = []

    # Generate before QC plots
    try:
        before_result = plot_qc_metrics_before.invoke({"adata_path": adata_before_path, "output_dir": output_dir})
        results.append(f"✓ Before QC plots: {before_result}")
    except Exception as e:
        results.append(f"✗ Before QC plots failed: {e}")

    # Generate after QC plots
    try:
        after_result = plot_qc_metrics_after.invoke({"adata_path": adata_after_path, "output_dir": output_dir})
        results.append(f"✓ After QC plots: {after_result}")
    except Exception as e:
        results.append(f"✗ After QC plots failed: {e}")

    # Generate comparison plots
    try:
        comp_result = plot_qc_comparison.invoke({
            "adata_before_path": adata_before_path,
            "adata_after_path": adata_after_path,
            "output_dir": output_dir
        })
        results.append(f"✓ QC comparison plots: {comp_result}")
    except Exception as e:
        results.append(f"✗ QC comparison plots failed: {e}")

    # Get all plot paths
    plot_paths = get_qc_plot_paths.invoke({})

    summary = f"""
QC Visualizations Generated
============================
{chr(10).join(results)}

All plots saved to: {output_dir}

Generated plots:
{chr(10).join(f'  - {k}: {v}' for k, v in plot_paths.items())}
"""

    return summary


@tool
def generate_integration_visualizations(
    adata_before_path: Annotated[str, "Path to h5ad file before integration"],
    adata_after_path: Annotated[str, "Path to h5ad file after integration"],
    batch_key: Annotated[str, "Column name for batch information"] = "batch",
    celltype_key: Annotated[str, "Column name for cell type annotation"] = "cell_type",
    embedding_key: Annotated[str, "Embedding to use (X_scVI, X_harmony, etc.)"] = "X_scVI",
    output_dir: Annotated[str, "Output directory for integration plots"] = "./integration_plots",
) -> Annotated[str, "Summary of integration visualizations generated"]:
    """
    Generate comprehensive integration visualizations including before/after comparison.

    This orchestrates:
    - Integration metrics before correction
    - Integration metrics after correction
    - Before/after comparison plots
    - Batch mixing analysis

    Returns paths to all generated plots.
    """
    from scrna_agent.tools.integration_visualization import (
        plot_integration_before,
        plot_integration_after,
        plot_integration_comparison,
        get_integration_plot_paths,
    )

    results = []

    # Generate before integration plots
    try:
        before_result = plot_integration_before.invoke({
            "adata_path": adata_before_path,
            "batch_key": batch_key,
            "output_dir": output_dir
        })
        results.append(f"✓ Before integration plots: {before_result}")
    except Exception as e:
        results.append(f"✗ Before integration plots failed: {e}")

    # Generate after integration plots
    try:
        after_result = plot_integration_after.invoke({
            "adata_path": adata_after_path,
            "batch_key": batch_key,
            "celltype_key": celltype_key,
            "embedding_key": embedding_key,
            "output_dir": output_dir
        })
        results.append(f"✓ After integration plots: {after_result}")
    except Exception as e:
        results.append(f"✗ After integration plots failed: {e}")

    # Generate comparison plots
    try:
        comp_result = plot_integration_comparison.invoke({
            "adata_before_path": adata_before_path,
            "adata_after_path": adata_after_path,
            "batch_key": batch_key,
            "output_dir": output_dir
        })
        results.append(f"✓ Integration comparison plots: {comp_result}")
    except Exception as e:
        results.append(f"✗ Integration comparison plots failed: {e}")

    # Get all plot paths
    plot_paths = get_integration_plot_paths.invoke({})

    summary = f"""
Integration Visualizations Generated
=====================================
{chr(10).join(results)}

All plots saved to: {output_dir}

Generated plots:
{chr(10).join(f'  - {k}: {v}' for k, v in plot_paths.items())}
"""

    return summary


@tool
def generate_condition_analysis(
    adata_path: Annotated[str, "Path to h5ad file with annotations"],
    condition_keys: Annotated[List[str], "List of condition variables to analyze (e.g., ['age', 'treatment'])"],
    celltype_key: Annotated[str, "Column name for cell type annotation"] = "cell_type",
    output_dir: Annotated[str, "Output directory for condition plots"] = "./condition_analysis",
) -> Annotated[str, "Summary of condition-wise analysis generated"]:
    """
    Generate comprehensive condition-wise cell type distribution analysis.

    Analyzes how cell type composition changes across experimental conditions.
    Useful for:
    - Young vs Aged comparisons
    - Treatment vs Control
    - Disease vs Healthy
    - Time-course experiments

    Returns paths to all generated plots.
    """
    from scrna_agent.tools.integration_visualization import (
        plot_condition_celltype_distribution,
    )

    results = []

    for condition_key in condition_keys:
        try:
            result = plot_condition_celltype_distribution.invoke({
                "adata_path": adata_path,
                "condition_key": condition_key,
                "celltype_key": celltype_key,
                "output_dir": output_dir
            })
            results.append(f"✓ Condition '{condition_key}' analysis: Success")
        except Exception as e:
            results.append(f"✗ Condition '{condition_key}' analysis failed: {e}")

    summary = f"""
Condition-wise Analysis Generated
==================================
Analyzed {len(condition_keys)} condition(s): {', '.join(condition_keys)}

{chr(10).join(results)}

All plots saved to: {output_dir}
"""

    return summary


@tool
def generate_comprehensive_report(
    project_name: Annotated[str, "Name of the project"],
    qc_metrics: Annotated[Optional[Dict[str, Any]], "QC metrics dictionary"] = None,
    annotation_results: Annotated[Optional[Dict[str, Any]], "Annotation results dictionary"] = None,
    integration_results: Annotated[Optional[Dict[str, Any]], "Integration results dictionary"] = None,
    figure_paths: Annotated[Optional[Dict[str, List[str]]], "Dictionary of figure paths by category"] = None,
    output_dir: Annotated[str, "Output directory for report"] = "./reports",
) -> Annotated[str, "Path to generated report"]:
    """
    Generate a comprehensive analysis report with all results and figures.

    This orchestrates the complete report generation including:
    - QC summary
    - Annotation summary
    - Integration summary
    - All figures with captions
    - Methods section

    Returns path to generated markdown report.
    """
    from scrna_agent.tools.comprehensive_report import (
        initialize_report,
        add_qc_metrics,
        add_annotation_results,
        add_integration_results,
        add_figure,
        generate_markdown_report,
    )

    # Initialize report
    initialize_report.invoke({
        "project_name": project_name,
        "output_dir": output_dir
    })

    # Add QC metrics if provided
    if qc_metrics:
        add_qc_metrics.invoke(qc_metrics)

    # Add annotation results if provided
    if annotation_results:
        add_annotation_results.invoke(annotation_results)

    # Add integration results if provided
    if integration_results:
        add_integration_results.invoke(integration_results)

    # Add figures if provided
    if figure_paths:
        for category, paths in figure_paths.items():
            for path in paths:
                add_figure.invoke({
                    "figure_name": Path(path).stem,
                    "figure_path": path,
                    "figure_caption": f"{category.title()} analysis figure",
                    "figure_category": category
                })

    # Generate final report
    report_path = os.path.join(output_dir, f"{project_name}_comprehensive_report.md")
    result = generate_markdown_report.invoke({"output_file": report_path})

    return result


@tool
def run_complete_visualization_pipeline(
    project_name: Annotated[str, "Name of the project"],
    adata_raw_path: Annotated[str, "Path to raw h5ad (before QC)"],
    adata_qc_path: Annotated[str, "Path to QC'd h5ad"],
    adata_before_integration_path: Annotated[str, "Path to h5ad before integration"],
    adata_final_path: Annotated[str, "Path to final integrated h5ad"],
    batch_key: Annotated[str, "Batch key column name"] = "batch",
    celltype_key: Annotated[str, "Cell type column name"] = "cell_type",
    condition_keys: Annotated[Optional[List[str]], "Condition variables to analyze"] = None,
    embedding_key: Annotated[str, "Embedding to use"] = "X_scVI",
    base_output_dir: Annotated[str, "Base output directory"] = "./analysis_output",
) -> Annotated[str, "Summary of complete visualization pipeline"]:
    """
    Run the complete end-to-end visualization and reporting pipeline.

    This is the master orchestration tool that runs:
    1. QC visualizations (before/after comparison)
    2. Integration visualizations (before/after comparison)
    3. Condition-wise analysis (if condition_keys provided)
    4. Comprehensive report generation

    Returns summary of all generated outputs.
    """
    import os
    from pathlib import Path

    base_dir = Path(base_output_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    results = []
    all_figures = {}

    # 1. QC Visualizations
    qc_dir = base_dir / "qc_plots"
    try:
        qc_result = generate_qc_visualizations.invoke({
            "adata_before_path": adata_raw_path,
            "adata_after_path": adata_qc_path,
            "output_dir": str(qc_dir)
        })
        results.append("✓ QC visualizations completed")

        # Collect QC figure paths
        from scrna_agent.tools.qc_visualization import get_qc_plot_paths
        all_figures["qc"] = list(get_qc_plot_paths.invoke({}).values())
    except Exception as e:
        results.append(f"✗ QC visualizations failed: {e}")

    # 2. Integration Visualizations
    integ_dir = base_dir / "integration_plots"
    try:
        integ_result = generate_integration_visualizations.invoke({
            "adata_before_path": adata_before_integration_path,
            "adata_after_path": adata_final_path,
            "batch_key": batch_key,
            "celltype_key": celltype_key,
            "embedding_key": embedding_key,
            "output_dir": str(integ_dir)
        })
        results.append("✓ Integration visualizations completed")

        # Collect integration figure paths
        from scrna_agent.tools.integration_visualization import get_integration_plot_paths
        all_figures["integration"] = list(get_integration_plot_paths.invoke({}).values())
    except Exception as e:
        results.append(f"✗ Integration visualizations failed: {e}")

    # 3. Condition Analysis (if requested)
    if condition_keys:
        cond_dir = base_dir / "condition_analysis"
        try:
            cond_result = generate_condition_analysis.invoke({
                "adata_path": adata_final_path,
                "condition_keys": condition_keys,
                "celltype_key": celltype_key,
                "output_dir": str(cond_dir)
            })
            results.append(f"✓ Condition analysis completed for {len(condition_keys)} condition(s)")

            # Collect condition figure paths
            all_figures["condition"] = [
                str(cond_dir / f"condition_{cond}_celltype_distribution.png")
                for cond in condition_keys
            ]
        except Exception as e:
            results.append(f"✗ Condition analysis failed: {e}")

    # 4. Generate Comprehensive Report
    report_dir = base_dir / "reports"
    try:
        report_result = generate_comprehensive_report.invoke({
            "project_name": project_name,
            "figure_paths": all_figures,
            "output_dir": str(report_dir)
        })
        results.append("✓ Comprehensive report generated")
    except Exception as e:
        results.append(f"✗ Report generation failed: {e}")

    summary = f"""
Complete Visualization Pipeline Completed
==========================================
Project: {project_name}
Output directory: {base_output_dir}

Pipeline Steps:
{chr(10).join(results)}

Generated Outputs:
- QC plots: {qc_dir}
- Integration plots: {integ_dir}
{f'- Condition analysis: {cond_dir}' if condition_keys else ''}
- Reports: {report_dir}

Total figures generated: {sum(len(v) for v in all_figures.values())}

Next steps:
- Review QC plots to ensure filtering was appropriate
- Review integration plots to verify batch mixing
- Review condition analysis to understand biological differences
- Read comprehensive report for full summary
"""

    return summary
