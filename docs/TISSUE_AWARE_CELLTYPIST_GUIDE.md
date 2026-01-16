# Tissue-Aware CellTypist Annotation Guide

## Overview

The tissue-aware CellTypist annotation automatically selects the best CellTypist model for each tissue type in your dataset. If your dataset contains multiple tissues, the tool will run CellTypist separately for each tissue with the appropriate model and merge the results.

## Features

- **Automatic model selection** based on tissue and sample type
- **Per-tissue annotation** for multi-tissue datasets
- **Priority-based model selection**: overrides > tissue > sample_type > default
- **Comprehensive outputs**: predictions CSV, tissue usage summary, run metadata
- **Model management**: Easy downloading and caching of models

---

## Installation

```bash
# Install the package
pip install -e .

# Or with full dependencies
pip install -e ".[full]"

# Install CellTypist
pip install celltypist

# Install PyYAML for model registry
pip install pyyaml
```

---

## Quick Start

### 1. Basic Usage (Single Tissue)

```bash
# Apply tissue to all cells
scrna-agent annotate \
  --file data.h5ad \
  --tissue "Spleen" \
  --species mouse \
  --output annotated.h5ad
```

### 2. Multi-Tissue Dataset (with mapping file)

Create a tissue mapping CSV (`tissue_map.csv`):
```csv
sample,tissue
sample1,Spleen
sample2,Blood
sample3,Liver
```

Run annotation:
```bash
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --sample-col sample \
  --species mouse \
  --output annotated.h5ad
```

### 3. With Sample Type Information

Create a sample type mapping CSV (`sample_type_map.csv`):
```csv
sample,sample_type
sample1,PBMC
sample2,tumor
sample3,organoid
```

Run annotation:
```bash
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --sample-type-map sample_type_map.csv \
  --sample-col sample \
  --output annotated.h5ad
```

---

## CLI Commands

### Models Management

List available models in registry:
```bash
scrna-agent models list
```

Download models for a specific tissue:
```bash
scrna-agent models pull --tissue Spleen
```

Download all models in registry:
```bash
scrna-agent models pull --all
```

### Annotation Options

```bash
scrna-agent annotate --help
```

**Key Options:**
- `--tissue <name>`: Apply tissue type to all cells
- `--tissue-map <csv>`: CSV file with sample→tissue mapping
- `--sample-type <name>`: Apply sample type to all cells
- `--sample-type-map <csv>`: CSV file with sample→sample_type mapping
- `--sample-col <column>`: Column name for sample IDs (default: "sample")
- `--allow-missing-tissue`: Fill missing tissues with "unknown"
- `--species <species>`: Species (mouse/human)
- `--output <path>`: Output file path

---

## Python API Usage

### Example 1: Single Tissue

```python
import scanpy as sc
from scrna_agent.tools.annotation_tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware
)

# Load data
adata = sc.read_h5ad("data.h5ad")
load_adata_for_annotation("data.h5ad")

# Add tissue metadata
add_tissue_metadata(tissue="Spleen", sample_col="sample")

# Run tissue-aware annotation
result = celltypist_annotate_tissue_aware(
    sample_col="sample",
    majority_voting=True,
    cluster_key="leiden",
    output_dir="./outputs"
)

print(result)
```

### Example 2: Multi-Tissue with Mapping

```python
import scanpy as sc
from scrna_agent.tools.annotation_tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware
)

# Load data
load_adata_for_annotation("data.h5ad")

# Add tissue metadata from mapping file
add_tissue_metadata(
    tissue_map="tissue_map.csv",
    sample_type_map="sample_type_map.csv",
    sample_col="sample",
    allow_missing_tissue=False
)

# Run tissue-aware annotation
result = celltypist_annotate_tissue_aware(
    sample_col="sample",
    majority_voting=True,
    cluster_key="leiden",
    output_dir="./outputs"
)

print(result)
```

### Example 3: Direct Model Management API

```python
from scrna_agent.tools.model_management import (
    select_model_for_tissue,
    load_model_registry,
    pull_celltypist_models
)

# Load model registry
registry = load_model_registry()

# Select model for specific tissue
model, rule = select_model_for_tissue(
    tissue="Spleen",
    sample_type="PBMC",
    registry=registry
)
print(f"Selected model: {model} (rule: {rule})")

# Download models
pull_celltypist_models(tissue="Spleen")
```

---

## Model Selection Priority

The tool uses the following priority order to select models:

1. **Override** (`tissue|sample_type`): Exact combination match
   - Example: `Blood|PBMC` → `Healthy_COVID19_PBMC.pkl`

2. **Tissue**: Tissue-specific model
   - Example: `Spleen` → `Immune_All_High.pkl`

3. **Sample Type**: Sample type-specific model
   - Example: `PBMC` → `Healthy_COVID19_PBMC.pkl`

4. **Default**: Fallback model
   - Example: `Immune_All_High.pkl`

---

## Model Registry Configuration

The model registry is stored at: `scrna_agent/resources/celltypist_models.yaml`

### Example Registry Structure

```yaml
default:
  model: "Immune_All_High.pkl"
  description: "Pan-immune model for general immune cell annotation"

tissue:
  "Spleen":
    model: "Immune_All_High.pkl"
    description: "Spleen immune populations"

  "Blood":
    model: "Immune_All_High.pkl"
    description: "Peripheral blood mononuclear cells"

sample_type:
  "PBMC":
    model: "Healthy_COVID19_PBMC.pkl"
    description: "Peripheral blood mononuclear cells"

overrides:
  "Blood|PBMC":
    model: "Healthy_COVID19_PBMC.pkl"
    description: "Blood PBMC samples"
```

### Adding Custom Models

You can add custom models to the registry:

```yaml
tissue:
  "MyTissue":
    model: "/path/to/custom_model.pkl"
    description: "My custom tissue model"
```

---

## Output Files

The tool generates several output files:

### 1. `celltypist_predictions.csv`
Per-cell predictions with metadata:
- `cell_id`: Cell barcode
- `sample`: Sample ID
- `tissue`: Tissue type
- `sample_type`: Sample type (if provided)
- `label`: Predicted cell type
- `confidence`: Prediction confidence score
- `model_used`: CellTypist model used
- `majority_vote`: Majority vote label (if enabled)

### 2. `tissue_model_usage.csv`
Summary of models used per tissue:
- `tissue`: Tissue name
- `model_used`: Model name
- `selection_rule`: How model was selected
- `n_cells`: Number of cells
- `samples`: Sample IDs (semicolon-separated)
- `sample_type`: Sample type

### 3. `run_metadata.json`
Detailed run metadata including:
- Tissue-model mapping
- Total tissues
- Total cells
- Model selection rules

---

## Common Workflows

### Workflow 1: Multi-Sample Study

```bash
# 1. Create tissue mapping
cat > tissue_map.csv << EOF
sample,tissue
ctrl_1,Spleen
ctrl_2,Spleen
treat_1,Spleen
treat_2,Spleen
EOF

# 2. Download models
scrna-agent models pull --tissue Spleen

# 3. Run annotation
scrna-agent annotate \
  --file samples.h5ad \
  --tissue-map tissue_map.csv \
  --sample-col sample \
  --output annotated_samples.h5ad
```

### Workflow 2: Multi-Tissue Atlas

```bash
# 1. Create tissue mapping
cat > atlas_tissues.csv << EOF
sample,tissue,sample_type
spleen_1,Spleen,normal
blood_1,Blood,PBMC
liver_1,Liver,normal
lung_1,Lung,normal
EOF

# 2. Download all models
scrna-agent models pull --all

# 3. Run annotation
scrna-agent annotate \
  --file atlas.h5ad \
  --tissue-map atlas_tissues.csv \
  --sample-col sample \
  --output annotated_atlas.h5ad
```

### Workflow 3: Tumor Samples

```bash
# Create mapping
cat > tumor_samples.csv << EOF
sample,tissue,sample_type
tumor_1,Lung,tumor
tumor_2,Lung,tumor
tumor_3,Breast,tumor
EOF

# Run annotation
scrna-agent annotate \
  --file tumor_data.h5ad \
  --tissue-map tumor_samples.csv \
  --sample-col sample \
  --output annotated_tumor.h5ad
```

---

## Troubleshooting

### Error: "Tissue metadata not found"

**Solution**: Make sure to run `add_tissue_metadata()` before `celltypist_annotate_tissue_aware()`:

```python
# Add this first
add_tissue_metadata(tissue="Spleen", sample_col="sample")

# Then run annotation
celltypist_annotate_tissue_aware(sample_col="sample")
```

### Error: "Unknown sample IDs in tissue mapping"

**Solution**: Either:
1. Fix the tissue mapping CSV to include all samples
2. Use `--allow-missing-tissue` flag:

```bash
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --allow-missing-tissue
```

### Error: "Model not found"

**Solution**: Download the required model:

```bash
scrna-agent models pull --tissue <tissue_name>
```

### Error: "Sample column not found"

**Solution**: Specify the correct sample column name:

```bash
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --sample-col "sample_id"  # Use your actual column name
```

---

## Best Practices

1. **Always download models first**:
   ```bash
   scrna-agent models pull --all
   ```

2. **Validate tissue mapping** before running annotation:
   - Check that all sample IDs in your h5ad are in the mapping CSV
   - Use consistent tissue names (check registry for available tissues)

3. **Use majority voting** for cleaner annotations:
   ```python
   celltypist_annotate_tissue_aware(
       majority_voting=True,
       cluster_key="leiden"
   )
   ```

4. **Specify output directory** to save results:
   ```python
   celltypist_annotate_tissue_aware(
       output_dir="./celltypist_outputs"
   )
   ```

5. **Check confidence scores** in the output CSV:
   - Low confidence (<0.5) may indicate uncertain predictions
   - Review and validate low-confidence cells manually

---

## Advanced Usage

### Custom Model Registry

Create a custom registry file:

```python
import yaml

custom_registry = {
    "default": {"model": "Immune_All_High.pkl"},
    "tissue": {
        "MyCustomTissue": {
            "model": "/path/to/my_model.pkl",
            "description": "Custom trained model"
        }
    }
}

with open("custom_registry.yaml", "w") as f:
    yaml.dump(custom_registry, f)
```

### Programmatic Model Selection

```python
from scrna_agent.tools.model_management import get_models_for_adata

# Get model assignments for all tissues in dataset
tissue_info = get_models_for_adata(adata, sample_col="sample")

# Inspect model selections
for tissue, info in tissue_info.items():
    print(f"Tissue: {tissue}")
    print(f"  Model: {info['model']}")
    print(f"  Selection rule: {info['rule']}")
    print(f"  Cells: {info['n_cells']}")
```

---

## Citation

If you use this tissue-aware CellTypist annotation feature, please cite:

1. **CellTypist**:
   Domínguez Conde C, et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science, 376(6594), eabl5197.

2. **This Tool**:
   [Your citation information]

---

## Support

For issues, questions, or feature requests:
- GitHub Issues: [repository URL]
- Documentation: [docs URL]
- Email: [contact email]

---

## License

[Your license information]
