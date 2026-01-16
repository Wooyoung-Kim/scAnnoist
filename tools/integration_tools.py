# Integration Tools
# Actual implementations for data integration and batch correction

import os
import pickle
from typing import Annotated, Optional, List, Dict, Any
import numpy as np
import pandas as pd

## 3rd party
from langchain_core.tools import tool

# Global state to store loaded data and models
_integration_state: Dict[str, Any] = {
    "adata": None,
    "scvi_model": None,
    "scanvi_model": None,
    "benchmark_results": None,
    "embeddings": {},
}


def _get_adata():
    """Get the current AnnData object from state."""
    if _integration_state["adata"] is None:
        raise ValueError("No data loaded. Use load_data tool first.")
    return _integration_state["adata"]


# ============================================================
# Data Loading and Preparation
# ============================================================
@tool
def load_adata_for_integration(
    file_path: Annotated[str, "Path to the h5ad file"],
) -> Annotated[str, "Summary of loaded data"]:
    """Load an h5ad file for integration. Data should be QC'd."""
    import scanpy as sc
    
    adata = sc.read_h5ad(file_path)
    _integration_state["adata"] = adata
    
    # Check for raw counts
    has_counts = "counts" in adata.layers
    has_raw = adata.raw is not None
    
    # Check batch info
    batch_cols = [c for c in adata.obs.columns if 
                  any(k in c.lower() for k in ["batch", "sample", "donor", "library"])]
    
    summary = f"""
Loaded: {file_path}
- Cells: {adata.n_obs:,}
- Genes: {adata.n_vars:,}
- Has counts layer: {has_counts}
- Has .raw: {has_raw}
- Potential batch columns: {batch_cols}
- Obs columns: {list(adata.obs.columns)[:15]}
"""
    return summary


@tool
def prepare_for_integration(
    batch_key: Annotated[str, "Column containing batch/sample information"],
    n_top_genes: Annotated[int, "Number of highly variable genes to select"] = 3000,
    layer: Annotated[Optional[str], "Layer containing raw counts (if None, use .X)"] = "counts",
) -> Annotated[str, "Preparation summary"]:
    """
    Prepare data for integration: select HVGs and ensure raw counts are available.
    This is required before running scVI/scANVI.
    """
    import scanpy as sc
    
    adata = _get_adata()
    
    # Verify batch key
    if batch_key not in adata.obs:
        available = [c for c in adata.obs.columns if adata.obs[c].dtype == "category" or 
                     adata.obs[c].dtype == "object"]
        return f"Batch key '{batch_key}' not found. Categorical columns: {available}"
    
    # Check/create counts layer
    if layer and layer not in adata.layers:
        if adata.raw is not None:
            adata.layers["counts"] = adata.raw.X.copy()
            layer = "counts"
        else:
            return f"Layer '{layer}' not found and no .raw available. Provide raw counts."
    
    # Select HVGs
    if layer:
        # Use counts for HVG selection
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            flavor="seurat_v3",
            layer=layer,
            batch_key=batch_key,
        )
    else:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_genes,
            batch_key=batch_key,
        )
    
    # Subset to HVGs
    n_hvg = adata.var["highly_variable"].sum()
    
    # Get batch distribution
    batch_counts = adata.obs[batch_key].value_counts()
    
    return f"""
Data Prepared for Integration
=============================
Batch key: {batch_key}
Counts layer: {layer}
Highly variable genes: {n_hvg:,}

Batch Distribution:
{batch_counts.to_string()}
"""


# ============================================================
# scVI Integration
# ============================================================
@tool
def run_scvi(
    batch_key: Annotated[str, "Column containing batch information"],
    n_latent: Annotated[int, "Latent dimension"] = 30,
    n_layers: Annotated[int, "Number of neural network layers"] = 2,
    max_epochs: Annotated[int, "Maximum training epochs"] = 400,
    gene_likelihood: Annotated[str, "Gene likelihood: 'nb' (negative binomial) or 'zinb'"] = "nb",
    layer: Annotated[Optional[str], "Layer with raw counts"] = "counts",
    early_stopping: Annotated[bool, "Use early stopping"] = True,
) -> Annotated[str, "scVI training results"]:
    """
    Train scVI model for unsupervised batch correction and integration.
    scVI uses a VAE to learn batch-corrected latent representations.
    Requires raw counts.
    """
    import scvi
    
    adata = _get_adata()
    
    # Setup anndata
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=layer,
        batch_key=batch_key,
    )
    
    # Create model
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood=gene_likelihood,
    )
    
    # Train
    model.train(
        max_epochs=max_epochs,
        early_stopping=early_stopping,
        early_stopping_patience=20,
    )
    
    # Store model
    _integration_state["scvi_model"] = model
    
    # Get latent representation
    latent = model.get_latent_representation()
    adata.obsm["X_scVI"] = latent
    _integration_state["embeddings"]["X_scVI"] = latent
    
    # Get training history
    history = model.history
    final_loss = history["elbo_train"].iloc[-1] if "elbo_train" in history else "N/A"
    
    return f"""
scVI Training Complete
======================
Parameters:
- Batch key: {batch_key}
- Latent dimension: {n_latent}
- Layers: {n_layers}
- Gene likelihood: {gene_likelihood}
- Max epochs: {max_epochs}
- Early stopping: {early_stopping}

Results:
- Epochs trained: {len(history)}
- Final ELBO: {final_loss}
- Latent shape: {latent.shape}
- Stored in: adata.obsm['X_scVI']
"""


@tool
def run_scanvi_from_scvi(
    labels_key: Annotated[str, "Column containing cell type labels"],
    unlabeled_category: Annotated[str, "Label for unlabeled cells"] = "Unknown",
    max_epochs: Annotated[int, "Maximum training epochs"] = 100,
) -> Annotated[str, "scANVI training results"]:
    """
    Train scANVI model from pre-trained scVI for semi-supervised integration.
    Use when partial labels are available.
    Must run run_scvi first.
    """
    import scvi
    
    adata = _get_adata()
    scvi_model = _integration_state.get("scvi_model")
    
    if scvi_model is None:
        return "Error: Must run run_scvi first before scANVI."
    
    if labels_key not in adata.obs:
        return f"Labels key '{labels_key}' not found."
    
    # Convert to scANVI
    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )
    
    # Train
    scanvi_model.train(
        max_epochs=max_epochs,
        early_stopping=True,
    )
    
    # Store model
    _integration_state["scanvi_model"] = scanvi_model
    
    # Get latent representation
    latent = scanvi_model.get_latent_representation()
    adata.obsm["X_scANVI"] = latent
    _integration_state["embeddings"]["X_scANVI"] = latent
    
    # Get predictions
    predictions = scanvi_model.predict()
    adata.obs["scanvi_prediction"] = predictions
    
    # Label stats
    known_cells = (adata.obs[labels_key] != unlabeled_category).sum()
    unknown_cells = (adata.obs[labels_key] == unlabeled_category).sum()
    
    return f"""
scANVI Training Complete
========================
Parameters:
- Labels key: {labels_key}
- Unlabeled category: {unlabeled_category}
- Max epochs: {max_epochs}

Label Statistics:
- Known labels: {known_cells:,}
- Unknown (to predict): {unknown_cells:,}

Results:
- Latent shape: {latent.shape}
- Stored in: adata.obsm['X_scANVI']
- Predictions stored in: adata.obs['scanvi_prediction']
"""


# ============================================================
# Harmony Integration
# ============================================================
@tool
def run_harmony(
    batch_key: Annotated[str, "Column containing batch information"],
    n_pcs: Annotated[int, "Number of PCs to use"] = 50,
    max_iter: Annotated[int, "Maximum iterations"] = 30,
) -> Annotated[str, "Harmony integration results"]:
    """
    Run Harmony integration on PCA space.
    Fast and effective for simpler batch effects.
    Requires PCA to be computed first.
    """
    import scanpy as sc
    
    try:
        import harmonypy as hm
    except ImportError:
        return "Harmony not installed. Run: pip install harmonypy"
    
    adata = _get_adata()
    
    # Compute PCA if not present
    if "X_pca" not in adata.obsm:
        sc.pp.pca(adata, n_comps=n_pcs)
    
    # Run Harmony
    pca_data = adata.obsm["X_pca"][:, :n_pcs]
    batch_labels = adata.obs[batch_key].values
    
    ho = hm.run_harmony(
        pca_data,
        adata.obs,
        batch_key,
        max_iter_harmony=max_iter,
    )
    
    # Store corrected embedding
    adata.obsm["X_harmony"] = ho.Z_corr.T
    _integration_state["embeddings"]["X_harmony"] = ho.Z_corr.T
    
    return f"""
Harmony Integration Complete
============================
Parameters:
- Batch key: {batch_key}
- PCs used: {n_pcs}
- Max iterations: {max_iter}

Results:
- Harmony embedding shape: {ho.Z_corr.T.shape}
- Stored in: adata.obsm['X_harmony']
"""


# ============================================================
# Scanorama Integration
# ============================================================
@tool
def run_scanorama(
    batch_key: Annotated[str, "Column containing batch information"],
) -> Annotated[str, "Scanorama integration results"]:
    """
    Run Scanorama integration using mutual nearest neighbors.
    Good for panoramic stitching of batches.
    """
    import scanpy as sc
    
    try:
        import scanorama
    except ImportError:
        return "Scanorama not installed. Run: pip install scanorama"
    
    adata = _get_adata()
    
    # Split by batch
    batches = adata.obs[batch_key].unique()
    adatas = [adata[adata.obs[batch_key] == b].copy() for b in batches]
    
    # Get gene expression matrices
    datasets = [ad.X if not hasattr(ad.X, 'toarray') else ad.X.toarray() for ad in adatas]
    genes_list = [ad.var_names.tolist() for ad in adatas]
    
    # Run Scanorama
    corrected, _ = scanorama.correct(datasets, genes_list, return_dimred=True)
    
    # Combine corrected embeddings
    corrected_embedding = np.vstack([c.X for c in corrected])
    
    # Reorder to match original adata
    adata.obsm["X_scanorama"] = corrected_embedding
    _integration_state["embeddings"]["X_scanorama"] = corrected_embedding
    
    return f"""
Scanorama Integration Complete
==============================
Parameters:
- Batch key: {batch_key}
- Number of batches: {len(batches)}

Results:
- Scanorama embedding shape: {corrected_embedding.shape}
- Stored in: adata.obsm['X_scanorama']
"""


# ============================================================
# Benchmarking
# ============================================================
@tool
def benchmark_integration(
    embeddings: Annotated[List[str], "List of embedding keys in obsm to compare"],
    batch_key: Annotated[str, "Column containing batch information"],
    label_key: Annotated[str, "Column containing cell type labels for bio conservation"],
) -> Annotated[str, "Benchmarking results"]:
    """
    Benchmark integration methods using scib-metrics.
    Computes batch correction and bio-conservation metrics.
    
    Metrics computed:
    - Batch correction: silhouette_batch, iLISI
    - Bio conservation: silhouette_label, cLISI, NMI, ARI
    """
    import scanpy as sc
    
    try:
        from scib_metrics.benchmark import Benchmarker
    except ImportError:
        return "scib-metrics not installed. Run: pip install scib-metrics"
    
    adata = _get_adata()
    
    # Verify embeddings exist
    available = [e for e in embeddings if e in adata.obsm]
    missing = [e for e in embeddings if e not in adata.obsm]
    
    if not available:
        return f"No embeddings found. Available: {list(adata.obsm.keys())}"
    
    if missing:
        print(f"Warning: Missing embeddings: {missing}")
    
    # Add PCA if not in list
    if "X_pca" not in available and "X_pca" in adata.obsm:
        available = ["X_pca"] + available
    
    # Compute neighbors for each embedding (required for some metrics)
    for emb in available:
        key_added = f"neighbors_{emb.replace('X_', '')}"
        sc.pp.neighbors(adata, use_rep=emb, key_added=key_added)
    
    # Run benchmarker
    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        embedding_obsm_keys=available,
        n_jobs=-1,
    )
    bm.benchmark()
    
    # Get results
    results_df = bm.get_results(min_max_scale=False)
    _integration_state["benchmark_results"] = results_df
    
    # Format output
    result = f"""
Integration Benchmark Results
=============================
Embeddings compared: {available}
Missing: {missing}
Batch key: {batch_key}
Label key: {label_key}

Metrics (higher is better):
{results_df.to_string()}

Best Methods:
- Overall: {results_df['Total'].idxmax()}
- Batch correction: {results_df['Batch correction'].idxmax()}
- Bio conservation: {results_df['Bio conservation'].idxmax()}
"""
    return result


@tool
def select_best_integration(
    prioritize: Annotated[str, "Priority: 'batch', 'bio', or 'balanced'"] = "balanced",
) -> Annotated[str, "Selected best method"]:
    """
    Select the best integration method based on benchmark results.
    """
    results = _integration_state.get("benchmark_results")
    if results is None:
        return "No benchmark results. Run benchmark_integration first."
    
    if prioritize == "batch":
        best = results["Batch correction"].idxmax()
        score = results.loc[best, "Batch correction"]
    elif prioritize == "bio":
        best = results["Bio conservation"].idxmax()
        score = results.loc[best, "Bio conservation"]
    else:  # balanced
        best = results["Total"].idxmax()
        score = results.loc[best, "Total"]
    
    # Store selection
    _integration_state["selected_method"] = best
    
    return f"""
Best Integration Method
=======================
Priority: {prioritize}
Selected: {best}
Score: {score:.4f}

All Scores:
{results.loc[best].to_string()}

Recommendation: Use adata.obsm['{best}'] for downstream analysis.
"""


# ============================================================
# Post-Integration Processing
# ============================================================
@tool
def compute_umap_from_embedding(
    use_rep: Annotated[str, "Key in obsm for the integrated representation"],
    n_neighbors: Annotated[int, "Number of neighbors"] = 15,
    min_dist: Annotated[float, "UMAP min_dist parameter"] = 0.3,
) -> Annotated[str, "UMAP computation result"]:
    """Compute neighbors and UMAP from integrated representation."""
    import scanpy as sc
    
    adata = _get_adata()
    
    if use_rep not in adata.obsm:
        return f"Embedding '{use_rep}' not found. Available: {list(adata.obsm.keys())}"
    
    # Compute neighbors
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=n_neighbors)
    
    # Compute UMAP
    sc.tl.umap(adata, min_dist=min_dist)
    
    return f"""
UMAP Computed
=============
Using representation: {use_rep}
- n_neighbors: {n_neighbors}
- min_dist: {min_dist}

UMAP embedding stored in: adata.obsm['X_umap']
"""


@tool
def plot_integration_comparison(
    embeddings: Annotated[List[str], "List of embedding keys to compare"],
    color_by: Annotated[str, "Column to color points by (e.g., batch or cell type)"],
    output_path: Annotated[str, "Path to save the comparison plot"],
) -> Annotated[str, "Path to saved plot"]:
    """Generate side-by-side UMAP comparison of integration methods."""
    import scanpy as sc
    import matplotlib.pyplot as plt
    
    adata = _get_adata()
    
    available = [e for e in embeddings if e in adata.obsm]
    n_plots = len(available)
    
    if n_plots == 0:
        return "No valid embeddings to plot."
    
    # Compute UMAP for each embedding
    fig, axes = plt.subplots(1, n_plots, figsize=(5*n_plots, 4))
    if n_plots == 1:
        axes = [axes]
    
    for i, emb in enumerate(available):
        # Compute neighbors and UMAP for this embedding
        temp_neighbors_key = f"neighbors_temp_{i}"
        sc.pp.neighbors(adata, use_rep=emb, key_added=temp_neighbors_key)
        sc.tl.umap(adata, neighbors_key=temp_neighbors_key)
        
        # Store temporarily
        umap_key = f"X_umap_{emb.replace('X_', '')}"
        adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        
        # Plot
        sc.pl.embedding(adata, basis=umap_key, color=color_by, ax=axes[i], 
                       show=False, title=emb.replace("X_", ""))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    
    return f"Saved integration comparison to {output_path}"


@tool
def save_integration_model(
    method: Annotated[str, "Method name: 'scvi' or 'scanvi'"],
    output_dir: Annotated[str, "Directory to save the model"],
) -> Annotated[str, "Save result"]:
    """Save trained integration model for reproducibility."""
    os.makedirs(output_dir, exist_ok=True)
    
    if method.lower() == "scvi":
        model = _integration_state.get("scvi_model")
        if model is None:
            return "No scVI model found."
        model.save(output_dir, overwrite=True)
    elif method.lower() == "scanvi":
        model = _integration_state.get("scanvi_model")
        if model is None:
            return "No scANVI model found."
        model.save(output_dir, overwrite=True)
    else:
        return f"Unknown method: {method}. Use 'scvi' or 'scanvi'."
    
    return f"Saved {method} model to {output_dir}"


@tool
def save_integrated_adata(
    output_path: Annotated[str, "Path to save the integrated h5ad file"],
) -> Annotated[str, "Save result"]:
    """Save the integrated AnnData with all embeddings."""
    adata = _get_adata()
    adata.write_h5ad(output_path)
    
    embeddings_saved = [k for k in adata.obsm.keys() if k.startswith("X_")]
    
    return f"""
Saved Integrated Data
=====================
Path: {output_path}
Embeddings included: {embeddings_saved}
Cells: {adata.n_obs:,}
Genes: {adata.n_vars:,}
"""
