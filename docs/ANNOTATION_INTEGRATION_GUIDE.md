# Annotation Integration Guide

## Overview

scRNA-agent supports integrating multiple cell type annotation methods to create more reliable and accurate annotations. This guide explains how to combine CellTypist (deep learning), literature-based (PubMed + CellMarker), and ScType annotations for optimal results.

## Why Integrate Multiple Methods?

Each annotation method has strengths and weaknesses:

| Method | Strengths | Weaknesses |
|--------|-----------|------------|
| **CellTypist** | Fast, consistent, trained on large datasets | May miss novel cell types, tissue-specific limitations |
| **Literature** | Custom markers, recent discoveries, tissue-specific | Requires API keys, slower, marker availability varies |
| **ScType** | Database-driven, no training needed | Limited to known markers, may be outdated |

**Integration Benefits:**
- ✅ Higher confidence through consensus
- ✅ Better handling of edge cases
- ✅ Validation across methods
- ✅ Discovery of annotation conflicts

---

## Quick Start

### Basic Integration Workflow

```python
from scrna_agent.tools.annotation_tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware,
    annotate_with_literature_markers,
    create_optimal_annotation,
)

# 1. Load data and add tissue metadata
load_adata_for_annotation("data.h5ad")
add_tissue_metadata(tissue_map="tissue_map.csv", sample_col="sample")

# 2. Run CellTypist annotation
celltypist_annotate_tissue_aware(
    sample_col="sample",
    majority_voting=True,
    output_dir="./outputs"
)

# 3. Run literature-based annotation
annotate_with_literature_markers(
    tissue_type="Blood",
    cluster_key="leiden",
    output_col="literature_annotation"
)

# 4. Create optimal annotation
result = create_optimal_annotation(
    celltypist_col="celltypist_label",
    celltypist_conf_col="celltypist_confidence",
    literature_col="literature_annotation",
    confidence_threshold=0.5,
    output_col="optimal_annotation"
)

print(result)
# adata.obs now has 'optimal_annotation' column
```

---

## Integration Strategies

### 1. Optimal Annotation (Recommended)

**Function:** `create_optimal_annotation()`

**Strategy:** Intelligent combination based on confidence and agreement

**Logic:**
```
For each cell:
  if celltypist_confidence >= threshold:
      use CellTypist (high confidence)
  else:
      collect predictions from all methods
      if all methods agree:
          use consensus (confidence = 0.9)
      else:
          weighted voting with method-specific weights
```

**Example:**
```python
result = create_optimal_annotation(
    celltypist_col="celltypist_label",
    celltypist_conf_col="celltypist_confidence",
    literature_col="literature_annotation",
    sctype_col="sctype_annotation",  # Optional
    confidence_threshold=0.5,
    output_col="optimal_annotation"
)
```

**Parameters:**
- `celltypist_col`: CellTypist annotation column (required)
- `celltypist_conf_col`: CellTypist confidence column
- `literature_col`: Literature annotation column (optional)
- `sctype_col`: ScType annotation column (optional)
- `confidence_threshold`: Threshold for using CellTypist directly (default: 0.5)
- `output_col`: Output column name (default: "optimal_annotation")

**Output:**
- Creates `optimal_annotation` column
- Creates `optimal_annotation_confidence` column
- Creates `optimal_annotation_source` column (shows which method was used)

---

### 2. Weighted Voting

**Function:** `integrate_multiple_annotations(strategy="voting")`

**Strategy:** Combine methods using weighted votes

**Example:**
```python
result = integrate_multiple_annotations(
    methods=["celltypist_label", "literature_annotation", "sctype_annotation"],
    strategy="voting",
    weights={
        "celltypist_label": 1.0,
        "literature_annotation": 0.8,
        "sctype_annotation": 0.7
    },
    output_col="voting_result"
)
```

**When to use:**
- You have 3+ methods
- You want to give different weights to methods
- You don't have confidence scores

---

### 3. Confidence-Weighted Integration

**Function:** `integrate_multiple_annotations(strategy="confidence")`

**Strategy:** Weight votes by method confidence scores

**Example:**
```python
result = integrate_multiple_annotations(
    methods=["celltypist_label", "literature_annotation"],
    strategy="confidence",
    confidence_cols={
        "celltypist_label": "celltypist_confidence",
        "literature_annotation": "literature_confidence"
    },
    output_col="confidence_integrated"
)
```

**When to use:**
- You have confidence scores for all methods
- You want data-driven weighting
- Methods have varying reliability per cell

---

## Comparison and Analysis

### Compare Two Methods

```python
from scrna_agent.tools.annotation_tools import compare_celltypist_with_literature

comparison = compare_celltypist_with_literature(
    celltypist_col="celltypist_label",
    literature_col="literature_annotation",
    cluster_key="leiden"
)

print(comparison)
```

**Output includes:**
- Overall agreement rate
- Per-cluster agreement rates
- Cells with disagreements
- Cluster-level consistency
- Statistics per cell type

**Example output:**
```
================================================================================
ANNOTATION COMPARISON REPORT
================================================================================

Methods Compared: celltypist_label vs literature_annotation
Total cells: 10,000

OVERALL AGREEMENT
-----------------
Agreement rate: 78.5%
Disagreement: 2,150 cells (21.5%)

CLUSTER-LEVEL AGREEMENT
-----------------------
Cluster 0: 92.3% agreement (450/487 cells)
  CellTypist: B cells (95%)
  Literature: B cells (92%)

Cluster 1: 45.2% agreement (123/272 cells)
  CellTypist: T cells (60%)
  Literature: NK cells (40%)
  ⚠ High disagreement - manual review recommended
...
```

---

## Validation

### Validate Against Known Markers

```python
from scrna_agent.tools.annotation_integration import validate_integrated_annotation

# Define marker genes for each cell type
markers = {
    "B cells": ["CD19", "MS4A1", "CD79A"],
    "T cells": ["CD3D", "CD3E", "CD8A"],
    "NK cells": ["GNLY", "NKG7", "NCAM1"]
}

validation = validate_integrated_annotation(
    adata,
    annotation_col="optimal_annotation",
    marker_dict=markers,
    min_expression_threshold=0.5
)

print(validation)
```

**Output:**
```python
{
    "B cells": {
        "status": "validated",
        "n_cells": 2500,
        "validation_score": 0.85,
        "expressed_markers": 2,
        "marker_expression": {
            "CD19": 1.2,
            "MS4A1": 0.8,
            "CD79A": 0.9
        }
    },
    ...
}
```

---

## Complete Workflow Examples

### Example 1: Two-Method Integration

```python
# Load data
load_adata_for_annotation("data.h5ad")
add_tissue_metadata(tissue="Blood")

# Method 1: CellTypist
celltypist_annotate_tissue_aware(output_dir="./outputs")

# Method 2: Literature
annotate_with_literature_markers(
    tissue_type="Blood",
    cluster_key="leiden",
    output_col="literature_annotation"
)

# Compare
comparison = compare_celltypist_with_literature()
print(comparison)

# Integrate
result = create_optimal_annotation(
    celltypist_col="celltypist_label",
    celltypist_conf_col="celltypist_confidence",
    literature_col="literature_annotation",
    output_col="final_annotation"
)
```

### Example 2: Three-Method Integration

```python
# Assuming all three methods have been run:
# - CellTypist -> celltypist_label
# - Literature -> literature_annotation
# - ScType -> sctype_annotation

# Weighted voting
result = integrate_multiple_annotations(
    methods=["celltypist_label", "literature_annotation", "sctype_annotation"],
    strategy="voting",
    weights={
        "celltypist_label": 1.0,
        "literature_annotation": 0.8,
        "sctype_annotation": 0.7
    },
    output_col="three_way_voting"
)

# Or use optimal strategy
result = create_optimal_annotation(
    celltypist_col="celltypist_label",
    celltypist_conf_col="celltypist_confidence",
    literature_col="literature_annotation",
    sctype_col="sctype_annotation",
    output_col="three_way_optimal"
)
```

---

## Best Practices

### 1. Start with CellTypist
CellTypist is fast and provides good baseline annotations with confidence scores.

### 2. Add Literature for Specific Tissues
Use literature-based annotation for tissue-specific or rare cell types.

### 3. Check Agreement First
Always compare methods before integration to identify conflicts.

### 4. Use Optimal Strategy
The `create_optimal_annotation()` function provides the best balance of speed and accuracy.

### 5. Validate Results
Use marker gene validation to verify integrated annotations.

### 6. Adjust Confidence Threshold
- **High threshold (0.7+)**: More conservative, relies more on CellTypist
- **Low threshold (0.3)**: More liberal, integrates more methods
- **Default (0.5)**: Good balance

### 7. Manual Review for Conflicts
When methods disagree significantly (>30% disagreement), manually review those clusters.

---

## Understanding Confidence Scores

### CellTypist Confidence
- Based on prediction probability from the neural network
- Range: 0-1
- >0.7: High confidence
- 0.4-0.7: Medium confidence
- <0.4: Low confidence

### Optimal Annotation Confidence
- **High confidence (0.8-1.0)**: Methods agree or CellTypist very confident
- **Medium confidence (0.5-0.8)**: Weighted voting used, partial agreement
- **Low confidence (<0.5)**: Methods disagree, uncertain

---

## Troubleshooting

### Issue: Low Agreement Between Methods

**Symptoms:** Agreement <60%

**Possible causes:**
- Different tissues/sample types
- Different granularity (broad vs specific cell types)
- One method using wrong tissue model

**Solutions:**
```python
# 1. Check tissue models
from scrna_agent.tools.model_management import get_models_for_adata
models = get_models_for_adata(adata)
print(models)

# 2. Use more permissive integration
result = create_optimal_annotation(
    confidence_threshold=0.3,  # Lower threshold
    output_col="permissive_annotation"
)
```

### Issue: Many "Unknown" Annotations

**Symptoms:** >20% cells labeled "Unknown"

**Possible causes:**
- Novel cell types not in training data
- Poor quality data
- Wrong tissue model

**Solutions:**
```python
# 1. Run literature search with more markers
annotate_with_literature_markers(
    top_n_genes=100,  # Use more markers
    min_confidence=0.2  # Lower threshold
)

# 2. Use ScType as additional method
sctype_annotate(tissue="Blood", output_col="sctype_annotation")
```

### Issue: Contradictory Annotations

**Symptoms:** Same cluster assigned different cell types

**Possible causes:**
- Mixed populations in cluster
- Cluster resolution too low
- Transitional states

**Solutions:**
```python
# 1. Increase clustering resolution
import scanpy as sc
sc.tl.leiden(adata, resolution=1.0)  # Increase from 0.5

# 2. Check cluster purity
comparison = compare_celltypist_with_literature(cluster_key="leiden")
# Review clusters with low agreement

# 3. Use optimal annotation which handles this
result = create_optimal_annotation()
```

---

## API Reference

### Main Functions

#### `create_optimal_annotation()`
Create optimal annotation by intelligently combining methods.

**Returns:** Dictionary with statistics

**Creates columns:**
- `{output_col}`: Optimal annotation
- `{output_col}_confidence`: Confidence scores
- `{output_col}_source`: Source method used

---

#### `integrate_multiple_annotations()`
Integrate multiple methods using voting or confidence.

**Parameters:**
- `methods`: List of column names
- `strategy`: "voting" or "confidence"
- `weights`: Optional weight dictionary
- `confidence_cols`: Optional confidence column mapping
- `output_col`: Output column name

---

#### `compare_celltypist_with_literature()`
Compare CellTypist and literature annotations.

**Returns:** Detailed comparison report string

---

## CLI Usage (Future)

Future versions may support CLI integration:

```bash
# Run CellTypist
scrna-agent annotate \
  --file data.h5ad \
  --tissue Blood \
  --output annotated.h5ad

# Integrate with literature
scrna-agent integrate \
  --file annotated.h5ad \
  --methods celltypist,literature \
  --strategy optimal \
  --output integrated.h5ad
```

---

## Examples

See `examples/annotation_integration_example.py` for:
- Simple integration workflow
- Detailed comparison analysis
- Three-way integration
- Custom integration strategies
- Complete workflow with visualization

---

## Summary

**Installation:** Already included in scRNA-agent

**Basic usage:**
```python
from scrna_agent.tools import (
    celltypist_annotate_tissue_aware,
    annotate_with_literature_markers,
    create_optimal_annotation
)

# Run both methods
celltypist_annotate_tissue_aware(output_dir="./outputs")
annotate_with_literature_markers(tissue_type="Blood")

# Integrate
create_optimal_annotation(
    celltypist_col="celltypist_label",
    literature_col="literature_annotation"
)
```

**Key benefits:**
- Higher accuracy through consensus
- Better handling of edge cases
- Confidence scoring
- Conflict detection and resolution

---

**Ready to create optimal annotations!** 🎯
