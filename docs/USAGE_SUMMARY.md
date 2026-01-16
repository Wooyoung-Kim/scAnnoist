# Tissue-Aware CellTypist Annotation - Usage Summary

## What Has Been Implemented

The scRNA-seq agent now includes **tissue-aware CellTypist annotation** that automatically selects the best model for each tissue type in your dataset.

### Key Features

✅ **Automatic model selection** based on tissue/sample_type
✅ **Per-tissue annotation** for multi-tissue datasets
✅ **Priority-based selection** (override > tissue > sample_type > default)
✅ **Comprehensive outputs** (predictions CSV, tissue usage, metadata)
✅ **Model management** (download, cache, list)
✅ **Full test coverage** (unit, integration, end-to-end)

---

## Quick Start for Users

### 1. Installation

```bash
cd /mnt2/kwy/scrna_agent
pip install -e .
pip install celltypist pyyaml
```

### 2. Download Models

```bash
# List available models
scrna-agent models list

# Download all models
scrna-agent models pull --all

# Or download for specific tissue
scrna-agent models pull --tissue Spleen
```

### 3. Run Annotation

**Option A: Single tissue (all cells)**
```bash
scrna-agent annotate \
  --file your_data.h5ad \
  --tissue "Spleen" \
  --species mouse \
  --output annotated.h5ad
```

**Option B: Multi-tissue (with mapping file)**
```bash
# Create mapping file
cat > tissue_map.csv << EOF
sample,tissue
sample1,Spleen
sample2,Blood
sample3,Liver
EOF

# Run annotation
scrna-agent annotate \
  --file your_data.h5ad \
  --tissue-map tissue_map.csv \
  --sample-col sample \
  --output annotated.h5ad
```

### 4. Check Outputs

The tool generates:
- `celltypist_predictions.csv` - Per-cell predictions with metadata
- `tissue_model_usage.csv` - Summary of models used
- `run_metadata.json` - Detailed run information
- Updated h5ad file with annotation columns

---

## Python API Example

```python
import scanpy as sc
from scrna_agent.tools.annotation_tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware
)

# Load data
load_adata_for_annotation("data.h5ad")

# Add tissue metadata
add_tissue_metadata(
    tissue_map="tissue_map.csv",
    sample_col="sample"
)

# Run annotation
result = celltypist_annotate_tissue_aware(
    sample_col="sample",
    majority_voting=True,
    cluster_key="leiden",
    output_dir="./outputs"
)

print(result)
```

---

## CLI Options Reference

```bash
scrna-agent annotate [OPTIONS]

Required:
  --file, -f PATH              Input h5ad file

Tissue/Sample Metadata:
  --tissue TEXT                Tissue type (apply to all cells)
  --tissue-map PATH            CSV with sample→tissue mapping
  --sample-type TEXT           Sample type (apply to all cells)
  --sample-type-map PATH       CSV with sample→sample_type mapping
  --sample-col TEXT            Sample ID column [default: sample]
  --allow-missing-tissue       Fill missing with 'unknown'

Other:
  --species, -s TEXT           Species [default: mouse]
  --output, -o PATH            Output file
  --methods LIST               Annotation methods to use
```

---

## Model Selection Logic

The tool selects models using this priority:

```
1. Override (tissue|sample_type) → Exact match
   Example: Blood|PBMC → Healthy_COVID19_PBMC.pkl

2. Tissue → Tissue-specific
   Example: Spleen → Immune_All_High.pkl

3. Sample Type → Sample type-specific
   Example: PBMC → Healthy_COVID19_PBMC.pkl

4. Default → Fallback
   Example: Any → Immune_All_High.pkl
```

### Example Scenarios

| Tissue | Sample Type | Model Used | Rule |
|--------|-------------|------------|------|
| Blood | PBMC | Healthy_COVID19_PBMC.pkl | override |
| Spleen | - | Immune_All_High.pkl | tissue |
| Unknown | PBMC | Healthy_COVID19_PBMC.pkl | sample_type |
| Unknown | - | Immune_All_High.pkl | default |

---

## File Structure

```
scrna_agent/
├── tools/
│   ├── annotation_tools.py          # Annotation functions
│   ├── model_management.py          # Model selection & metadata
│   └── __init__.py                  # Exports
├── resources/
│   └── celltypist_models.yaml       # Model registry
├── tests/
│   └── test_tissue_aware_celltypist.py  # Tests
├── examples/
│   ├── tissue_aware_annotation_example.py  # Python examples
│   └── run_tissue_aware_annotation.sh      # Shell examples
└── docs/
    ├── TISSUE_AWARE_CELLTYPIST_GUIDE.md   # Full guide
    └── USAGE_SUMMARY.md                    # This file
```

---

## Example Workflows

### Workflow 1: Basic Multi-Sample Study

```bash
# 1. Prepare mapping
echo "sample,tissue" > map.csv
echo "ctrl_1,Spleen" >> map.csv
echo "ctrl_2,Spleen" >> map.csv
echo "treat_1,Spleen" >> map.csv

# 2. Download models
scrna-agent models pull --tissue Spleen

# 3. Annotate
scrna-agent annotate --file data.h5ad --tissue-map map.csv --output out.h5ad
```

### Workflow 2: Multi-Tissue Atlas

```bash
# 1. Create complete mapping
cat > atlas_map.csv << EOF
sample,tissue,sample_type
s1,Spleen,normal
s2,Blood,PBMC
s3,Liver,normal
s4,Lung,normal
EOF

# 2. Download all models
scrna-agent models pull --all

# 3. Annotate
scrna-agent annotate \
  --file atlas.h5ad \
  --tissue-map atlas_map.csv \
  --output atlas_annotated.h5ad
```

### Workflow 3: Tumor Samples

```bash
# 1. Create mapping with sample types
cat > tumor_map.csv << EOF
sample,tissue,sample_type
t1,Lung,tumor
t2,Lung,tumor
t3,Breast,tumor
n1,Lung,normal
EOF

# 2. Annotate
scrna-agent annotate \
  --file tumor_data.h5ad \
  --tissue-map tumor_map.csv \
  --output tumor_annotated.h5ad
```

---

## Output Files Explained

### 1. celltypist_predictions.csv

```csv
cell_id,sample,tissue,sample_type,label,confidence,model_used
cell_001,s1,Spleen,,B cells,0.95,Immune_All_High.pkl
cell_002,s1,Spleen,,T cells,0.88,Immune_All_High.pkl
cell_003,s2,Blood,PBMC,NK cells,0.92,Healthy_COVID19_PBMC.pkl
```

### 2. tissue_model_usage.csv

```csv
tissue,model_used,selection_rule,n_cells,samples
Spleen,Immune_All_High.pkl,tissue(Spleen),5000,s1;s2
Blood,Healthy_COVID19_PBMC.pkl,override(Blood|PBMC),3000,s3
```

### 3. run_metadata.json

```json
{
  "tissue_model_mapping": {
    "Spleen": {
      "model": "Immune_All_High.pkl",
      "rule": "tissue(Spleen)",
      "n_cells": 5000,
      "samples": ["s1", "s2"]
    }
  },
  "total_tissues": 2,
  "total_cells": 8000
}
```

---

## Troubleshooting

### Issue: "Tissue metadata not found"
**Fix**: Add tissue metadata before annotation
```python
add_tissue_metadata(tissue="Spleen", sample_col="sample")
celltypist_annotate_tissue_aware(sample_col="sample")
```

### Issue: "Unknown sample IDs"
**Fix**: Use `--allow-missing-tissue` or fix mapping
```bash
scrna-agent annotate --file data.h5ad --tissue-map map.csv --allow-missing-tissue
```

### Issue: "Model not found"
**Fix**: Download the model first
```bash
scrna-agent models pull --tissue Spleen
```

### Issue: "Sample column not found"
**Fix**: Specify correct column name
```bash
scrna-agent annotate --file data.h5ad --tissue-map map.csv --sample-col "sample_id"
```

---

## Advanced: Custom Model Registry

Edit `scrna_agent/resources/celltypist_models.yaml`:

```yaml
# Add your custom tissue
tissue:
  "MyTissue":
    model: "/path/to/custom_model.pkl"
    description: "Custom tissue model"

# Add custom override
overrides:
  "MyTissue|tumor":
    model: "/path/to/tumor_model.pkl"
    description: "Tumor-specific model"
```

---

## Testing

Run the test suite:

```bash
cd /mnt2/kwy/scrna_agent
python -m pytest tests/test_tissue_aware_celltypist.py -v
```

Tests cover:
- ✅ Model selection priority (override > tissue > sample_type > default)
- ✅ Tissue metadata validation
- ✅ CSV parsing and error handling
- ✅ Per-tissue annotation and merging
- ✅ Output file generation
- ✅ End-to-end workflows

---

## Documentation

- **Full Guide**: `docs/TISSUE_AWARE_CELLTYPIST_GUIDE.md`
- **Python Examples**: `examples/tissue_aware_annotation_example.py`
- **Shell Examples**: `examples/run_tissue_aware_annotation.sh`
- **Usage Summary**: `docs/USAGE_SUMMARY.md` (this file)

---

## Support

For help:
1. Check the full guide: `docs/TISSUE_AWARE_CELLTYPIST_GUIDE.md`
2. Run examples: `python examples/tissue_aware_annotation_example.py`
3. Check CLI help: `scrna-agent annotate --help`

---

## Next Steps

1. **Try the examples**:
   ```bash
   python examples/tissue_aware_annotation_example.py
   # or
   bash examples/run_tissue_aware_annotation.sh
   ```

2. **Read the full guide**:
   ```bash
   cat docs/TISSUE_AWARE_CELLTYPIST_GUIDE.md
   ```

3. **Run on your data**:
   - Create tissue mapping CSV
   - Download models
   - Run annotation
   - Check outputs

---

## Summary of Changes

### New CLI Commands
- `scrna-agent models list` - List models
- `scrna-agent models pull` - Download models
- `scrna-agent annotate --tissue-map` - Multi-tissue annotation
- `scrna-agent annotate --sample-type-map` - Sample type metadata

### New Python Functions
- `add_tissue_metadata()` - Add/validate tissue metadata
- `celltypist_annotate_tissue_aware()` - Run tissue-aware annotation
- `select_model_for_tissue()` - Select model programmatically
- `get_models_for_adata()` - Get model assignments
- `pull_celltypist_models()` - Download models

### New Files
- `tools/model_management.py` - Model management
- `resources/celltypist_models.yaml` - Model registry
- `tests/test_tissue_aware_celltypist.py` - Test suite
- `docs/TISSUE_AWARE_CELLTYPIST_GUIDE.md` - Full documentation
- `examples/tissue_aware_annotation_example.py` - Python examples
- `examples/run_tissue_aware_annotation.sh` - Shell examples

---

**Ready to use!** Start with the quick start section above or run the example scripts.
