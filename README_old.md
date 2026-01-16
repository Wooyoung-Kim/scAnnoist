# scrna-agent

A command-line tool for single-cell RNA sequencing (scRNA-seq) analysis. Runs a complete analysis pipeline from raw data to annotated clusters with marker genes.

## Features

### Core Pipeline
- **QC filtering**: Filter cells by gene count, UMI count, and mitochondrial percentage
- **Normalization**: Library size normalization, log transformation, HVG selection
- **Dimensionality reduction**: PCA and UMAP
- **Clustering**: Leiden or Louvain clustering
- **Marker detection**: Differential expression analysis to find cluster markers
- **Report generation**: Markdown/HTML report with summary statistics

### Advanced Annotation
- **Tissue-aware CellTypist**: Automatic tissue-specific model selection and annotation
- **Literature-based annotation**: PubMed + CellMarker database integration
- **Annotation integration**: Combine multiple methods for optimal cell type identification
- **Presentation generation**: Auto-create professional slide decks from results

## Installation

### Prerequisites

Create a conda environment with required dependencies:

```bash
conda env create -f environment.yml
conda activate scrna-agent
pip install -e .
```

### Quick Install (if dependencies already available)

```bash
pip install -e .
```

## Quick Start (30-second demo)

```bash
# Generate toy dataset
scrna-agent generate-toy --output ./data/toy/toy_data.h5ad

# Run the analysis pipeline
scrna-agent run --input ./data/toy/toy_data.h5ad --out ./demo_out

# Check the outputs
ls ./demo_out/
```

## Usage

### Basic Analysis

```bash
# Run on h5ad file
scrna-agent run --input data.h5ad --out ./results

# Run with custom parameters
scrna-agent run --input data.h5ad --out ./results \
    --min-genes 300 \
    --max-mt-pct 15 \
    --resolution 0.8
```

### CellRanger Input

```bash
# Run on 10x CellRanger output directory
scrna-agent prerun --input /path/to/cellranger_output --out ./results

# Example with real data path:
scrna-agent prerun --input /home/kwy7605/data_61/SARS/Count/JMV_ALI --out ./sars_results
```

### Commands

| Command | Description |
|---------|-------------|
| `run` | Run analysis pipeline on h5ad or 10x data |
| `prerun` | Run pipeline on CellRanger output directory |
| `generate-toy` | Generate synthetic dataset for testing |
| `config` | Generate default YAML config file |

### Options for `run` command

| Option | Description | Default |
|--------|-------------|---------|
| `--input`, `-i` | Input file (.h5ad) or directory (10x) | Required |
| `--out`, `-o` | Output directory | Required |
| `--config`, `-c` | YAML config file | - |
| `--min-genes` | Min genes per cell | 200 |
| `--min-cells` | Min cells per gene | 3 |
| `--max-mt-pct` | Max mitochondrial % | 20 |
| `--resolution` | Clustering resolution | 0.5 |
| `--clustering` | Method: leiden/louvain | leiden |
| `--seed` | Random seed | 42 |
| `--verbose`, `-v` | Verbose output | False |

## Input Formats

### Primary: AnnData (.h5ad)

Standard format for single-cell data. Load with:
```python
import scanpy as sc
adata = sc.read_h5ad("data.h5ad")
```

### 10x CellRanger Export

Directory containing:
- `matrix.mtx` or `matrix.mtx.gz`
- `barcodes.tsv` or `barcodes.tsv.gz`
- `features.tsv` or `genes.tsv` (with `.gz` variants)

### Seurat (TODO)

Seurat RDS import is not yet implemented.

## Output Files

| File | Description |
|------|-------------|
| `qc_summary.csv` | QC metrics and filtering statistics |
| `umap.png` | UMAP visualization colored by clusters |
| `clusters.csv` | Cell-to-cluster assignments |
| `markers.csv` | Marker genes per cluster (DEG results) |
| `report.md` | Analysis report with summary |
| `processed.h5ad` | Processed AnnData with all results |
| `qc_plots.png` | QC metric distributions |

## Configuration

Generate a default config file:

```bash
scrna-agent config --output my_config.yaml
```

Example `config.yaml`:

```yaml
min_genes: 200
min_cells: 3
max_mt_pct: 20.0
target_sum: 10000.0
n_top_genes: 2000
n_pcs: 50
n_neighbors: 15
resolution: 0.5
clustering_method: leiden
deg_method: wilcoxon
n_genes_per_cluster: 25
seed: 42
report_format: md
```

Run with config:

```bash
scrna-agent run --input data.h5ad --out ./results --config my_config.yaml
```

## Pipeline Steps

1. **Load data**: Read h5ad or 10x matrix files
2. **QC filtering**: Remove low-quality cells and genes
3. **Normalization**: Library size normalization + log1p
4. **HVG selection**: Identify highly variable genes
5. **PCA**: Dimensionality reduction
6. **Neighbors**: Build k-NN graph
7. **UMAP**: 2D embedding for visualization
8. **Clustering**: Leiden/Louvain community detection
9. **Markers**: Differential expression for cluster markers
10. **Report**: Generate summary report

## Testing

Run the test suite:

```bash
pytest tests/ -v
```

## Examples

### Full workflow example

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate scrna-agent

# 2. Install package
pip install -e .

# 3. Generate test data
scrna-agent generate-toy --output ./data/toy/toy_data.h5ad --n-cells 1000

# 4. Run analysis
scrna-agent run --input ./data/toy/toy_data.h5ad --out ./results -v

# 5. Check outputs
cat ./results/report.md
```

### Using with Python

```python
from scrna_agent.pipeline import run_pipeline, PipelineConfig

config = PipelineConfig(
    min_genes=300,
    resolution=0.8,
    seed=123
)

results = run_pipeline(
    input_path="data.h5ad",
    output_dir="./output",
    config=config
)

# Access results
print(f"Found {results.clusters['cluster'].nunique()} clusters")
print(f"Top markers: {results.markers.head()}")
```

## FAQ

**Q: How do I adjust clustering resolution?**

A: Use `--resolution` flag. Higher values = more clusters. Try 0.3-1.5.

**Q: Pipeline fails with "No counts layer found"**

A: Ensure your h5ad has raw counts. The tool expects `adata.X` to contain counts.

**Q: How do I run on multiple samples?**

A: Currently processes one sample at a time. For batch integration, use the advanced agent mode (requires `pip install -e ".[full]"`).

**Q: Can I use GPU acceleration?**

A: The MVP pipeline uses CPU. GPU support is available in full mode with scVI.

## Advanced Features

### Tissue-Aware Cell Type Annotation

Automatic tissue-specific CellTypist annotation with intelligent model selection:

```bash
# CLI usage
scrna-agent annotate \
  --file data.h5ad \
  --tissue-map tissue_map.csv \
  --output annotated.h5ad

# Python API
from scrna_agent.tools import (
    load_adata_for_annotation,
    add_tissue_metadata,
    celltypist_annotate_tissue_aware
)

load_adata_for_annotation("data.h5ad")
add_tissue_metadata(tissue_map="tissue_map.csv")
celltypist_annotate_tissue_aware(output_dir="./outputs")
```

**Features:**
- Automatic tissue-specific model selection
- Per-tissue annotation with merging
- Comprehensive output reports
- Model caching and management

See `docs/TISSUE_AWARE_CELLTYPIST_GUIDE.md` for details.

### Annotation Integration

Combine multiple annotation methods (CellTypist, Literature, ScType) for optimal results:

```python
from scrna_agent.tools import (
    celltypist_annotate_tissue_aware,
    annotate_with_literature_markers,
    create_optimal_annotation
)

# Run CellTypist
celltypist_annotate_tissue_aware(output_dir="./outputs")

# Run literature-based annotation
annotate_with_literature_markers(tissue_type="Blood")

# Create optimal annotation
create_optimal_annotation(
    celltypist_col="celltypist_label",
    literature_col="literature_annotation",
    output_col="optimal_annotation"
)
```

**Features:**
- Intelligent integration using confidence scores
- Consensus-based annotation
- Weighted voting strategies
- Validation and comparison tools

See `docs/ANNOTATION_INTEGRATION_GUIDE.md` for details.

### Presentation Generation

Auto-generate professional presentations from analysis results:

```bash
# Run annotation
scrna-agent annotate --file data.h5ad --tissue Blood

# Generate presentation
scrna-agent report -o ./outputs -t "My Study"
```

Creates a complete slide deck with visualizations and statistics.

See `docs/PPT_GENERATION_GUIDE.md` for details.

### Multi-Agent Mode

For advanced LangChain-based features (collaborative annotation, scVI integration, literature RAG), install with:

```bash
pip install -e ".[full]"
```

This enables the full multi-agent system with additional capabilities.

## License

MIT License

## Contributing

Contributions welcome! Please open an issue or PR.
