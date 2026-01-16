# Detailed Workflow Documentation

## Overview

The scRNA-seq Analysis Agent follows a structured workflow that combines multiple annotation methods and integration techniques to produce high-quality cell type annotations and batch-corrected data.

## Workflow Diagram

```
┌──────────────────────────────────────────────────────────────────────┐
│                        INPUT: QC'd h5ad file                          │
│                    (normalized, with raw counts)                      │
└────────────────────────────────┬─────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────┐
│                   PHASE 1: Data Initialization                        │
│  • Load data with load_data_for_pipeline()                           │
│  • Verify batch_key, species, raw counts                             │
│  • Check for existing embeddings (PCA, UMAP)                         │
└────────────────────────────────┬─────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────┐
│                  PHASE 2: Collaborative Annotation                    │
│                                                                       │
│  ┌─────────────────────────────────────────────────────────────────┐ │
│  │              Annotation Coordinator                              │ │
│  │                                                                  │ │
│  │  ┌──────────────┐ ┌──────────────┐ ┌──────────────────────────┐│ │
│  │  │  CellTypist  │ │   ScType     │ │   Literature RAG         ││ │
│  │  │  Specialist  │ │  Specialist  │ │    Specialist            ││ │
│  │  │              │ │              │ │                          ││ │
│  │  │ • Pre-trained│ │ • Marker DB  │ │ • PubMed search          ││ │
│  │  │   DL models  │ │ • Cluster    │ │ • CellMarker DB          ││ │
│  │  │ • Per-cell   │ │   scoring    │ │ • Marker validation      ││ │
│  │  │   confidence │ │              │ │                          ││ │
│  │  └──────┬───────┘ └──────┬───────┘ └────────────┬─────────────┘│ │
│  │         │                │                      │              │ │
│  │         └────────────────┼──────────────────────┘              │ │
│  │                          ▼                                     │ │
│  │              ┌───────────────────────┐                         │ │
│  │              │  Compare Annotations  │                         │ │
│  │              │  • Agreement rate     │                         │ │
│  │              │  • Identify conflicts │                         │ │
│  │              └───────────┬───────────┘                         │ │
│  │                          ▼                                     │ │
│  │              ┌───────────────────────┐                         │ │
│  │              │  Resolve Conflicts    │                         │ │
│  │              │  • Literature check   │                         │ │
│  │              │  • Marker validation  │                         │ │
│  │              └───────────┬───────────┘                         │ │
│  │                          ▼                                     │ │
│  │              ┌───────────────────────┐                         │ │
│  │              │  Consensus Annotation │                         │ │
│  │              │  (stored in adata.obs)│                         │ │
│  │              └───────────────────────┘                         │ │
│  └─────────────────────────────────────────────────────────────────┘ │
└────────────────────────────────┬─────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────┐
│                      PHASE 3: Integration                             │
│                                                                       │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐       │
│  │      scVI       │  │    Harmony      │  │   (Scanorama)   │       │
│  │  • VAE-based    │  │  • PCA-based    │  │   • MNN-based   │       │
│  │  • X_scVI       │  │  • X_harmony    │  │   • Optional    │       │
│  └────────┬────────┘  └────────┬────────┘  └────────┬────────┘       │
│           │                    │                    │                 │
│           └────────────────────┼────────────────────┘                 │
│                                ▼                                      │
│                    ┌───────────────────────┐                          │
│                    │  Benchmark (scib)     │                          │
│                    │  • Batch correction   │                          │
│                    │  • Bio conservation   │                          │
│                    └───────────┬───────────┘                          │
│                                ▼                                      │
│                    ┌───────────────────────┐                          │
│                    │  Select Best Method   │                          │
│                    └───────────────────────┘                          │
└────────────────────────────────┬─────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────┐
│                   PHASE 4: scANVI Fine-Tuning                         │
│                                                                       │
│  ┌─────────────────────────────────────────────────────────────────┐ │
│  │  Setup: Use consensus annotations as labels                     │ │
│  │  • Mark low-confidence cells as "Unknown"                       │ │
│  │  • Prepare for semi-supervised learning                         │ │
│  └────────────────────────────────────────────────────────────────┘ │
│                                │                                      │
│                                ▼                                      │
│  ┌─────────────────────────────────────────────────────────────────┐ │
│  │  scVI Pre-training → scANVI Fine-tuning                         │ │
│  │  • Uses known labels for supervision                            │ │
│  │  • Predicts labels for "Unknown" cells                          │ │
│  │  • Produces refined integration + annotations                   │ │
│  └─────────────────────────────────────────────────────────────────┘ │
│                                │                                      │
│                                ▼                                      │
│  ┌─────────────────┐  ┌─────────────────┐                            │
│  │   X_scANVI      │  │ scanvi_refined  │                            │
│  │  (embedding)    │  │  (predictions)  │                            │
│  └─────────────────┘  └─────────────────┘                            │
└────────────────────────────────┬─────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────┐
│                     PHASE 5: Finalization                             │
│                                                                       │
│  • Compute UMAP from best embedding                                  │
│  • Generate summary plots (batch, cell type, confidence)             │
│  • Save annotated h5ad file                                          │
│  • Save trained models (scVI, scANVI)                                │
│  • Export annotations as CSV                                         │
└────────────────────────────────┬─────────────────────────────────────┘
                                 │
                                 ▼
┌──────────────────────────────────────────────────────────────────────┐
│                            OUTPUTS                                    │
│                                                                       │
│  📁 output/                                                          │
│  ├── annotated_integrated.h5ad    # Final AnnData                    │
│  ├── annotations.csv              # Cell type annotations            │
│  ├── models/                                                         │
│  │   ├── scvi_model/              # Trained scVI                     │
│  │   └── scanvi_model/            # Trained scANVI                   │
│  └── plots/                                                          │
│      ├── umap_batch.png           # Integration check                │
│      ├── umap_celltype.png        # Cell type distribution           │
│      └── confidence_distribution.png                                 │
└──────────────────────────────────────────────────────────────────────┘
```

## Data Requirements

### Input h5ad Requirements

| Requirement | Description |
|-------------|-------------|
| `.X` | Normalized expression matrix (log1p) |
| `.layers['counts']` | Raw counts (required for scVI/scANVI) |
| `.obs[batch_key]` | Batch/sample information |
| `.obsm['X_pca']` | PCA embedding (optional, will compute if missing) |
| `.obs['leiden']` | Cluster assignments (optional, for ScType) |

### Output Structure

```python
adata.obs['celltypist_prediction']     # CellTypist per-cell predictions
adata.obs['celltypist_confidence']     # CellTypist confidence scores
adata.obs['sctype_annotation']         # ScType cluster annotations
adata.obs['literature_annotation']     # Literature-based annotations
adata.obs['consensus_annotation']      # Final consensus
adata.obs['scanvi_refined']            # scANVI-refined predictions
adata.obs['scanvi_refined_confidence'] # scANVI confidence

adata.obsm['X_scVI']                   # scVI latent space
adata.obsm['X_scANVI']                 # scANVI latent space
adata.obsm['X_harmony']                # Harmony-corrected PCA
adata.obsm['X_umap']                   # Final UMAP
```

## Agent Collaboration Protocol

### How Specialists Collaborate

1. **Independent Annotation**: Each specialist runs their method independently
2. **Results Collection**: Coordinator collects all predictions
3. **Comparison**: `compare_annotations()` finds agreements and conflicts
4. **Discussion**: For conflicts:
   - Check Literature RAG for published markers
   - Consider confidence scores
   - Evaluate marker expression
5. **Resolution**: Coordinator decides based on:
   - Method reliability for specific cell types
   - Literature support
   - Expression validation
6. **Consensus**: `create_consensus_annotation()` produces final labels

### Conflict Resolution Hierarchy

1. **High agreement (≥2/3 methods)**: Use majority vote
2. **Literature support**: Prefer method with literature backing
3. **Confidence scores**: Use highest confidence prediction
4. **Expression validation**: Check actual marker expression
5. **Manual flag**: Mark for user review if unresolved

## Customization

### Using Custom Markers

```python
from scrna_agent.tools import annotate_with_literature_markers

custom_markers = {
    "My_CellType_1": ["Gene1", "Gene2", "Gene3"],
    "My_CellType_2": ["Gene4", "Gene5"],
}

annotate_with_literature_markers(
    cell_type_markers=custom_markers,
    cluster_key="leiden"
)
```

### Using a Specific Model

```python
agent = create_scrna_pipeline_agent(model_name="gpt-4o")
# or
agent = create_scrna_pipeline_agent(model_name="claude-3-sonnet-20240229")
```

### Skipping Phases

Use individual components if you only need certain parts:

```python
# Only annotation
from scrna_agent import create_annotation_coordinator
coordinator = create_annotation_coordinator()

# Only integration
from scrna_agent.tools import run_scvi, run_harmony
```
