# New Features Guide - scRNA-agent Enhanced Capabilities

## Overview

This guide covers the new features added to scRNA-agent:

1. **Domain Specialist Agents** - Expert agents for specific biological contexts
2. **Comprehensive QC Visualization** - Before/after QC comparison plots
3. **Integration Visualization** - Batch mixing and integration quality plots
4. **Condition-wise Analysis** - Cell type distribution across experimental conditions
5. **Comprehensive Reporting** - Automated markdown reports with all figures
6. **Pipeline Orchestration** - End-to-end visualization and reporting

---

## 1. Domain Specialist Agents

### Overview

Six new domain expert agents have been added to provide biological interpretation during annotation:

- **Hematopoiesis Specialist**: Blood/BM/spleen, HSC, progenitors, immune lineages
- **Infection/Inflammation Specialist**: ISG-high states, cytokine responses, stress signatures
- **Aging Specialist**: Age-related changes, inflammaging, exhaustion, senescence
- **Cell State/Plasticity Specialist**: Activation, effector, memory, cycling, trajectories
- **Epithelial Specialist**: Epithelial lineages (basal, luminal, secretory, ciliated)
- **Organoid/Stemness Specialist**: Organoid cultures, stem cells, culture artifacts

### How They Work

These specialists are automatically consulted by the **Annotation Coordinator** when:
- Methods disagree on cell type labels
- Annotations have low confidence
- Clusters show ambiguous marker expression
- Domain-specific expertise is needed

### Usage

The domain specialists are integrated into the annotation workflow automatically. They work behind the scenes during the annotation phase:

```python
from scrna_agent.agents import create_annotation_coordinator

# Create coordinator (includes all specialists)
coordinator = create_annotation_coordinator(model_name="gpt-4o")

# Run annotation - specialists are consulted automatically
result = coordinator.invoke({
    "messages": [{
        "role": "user",
        "content": "Annotate this PBMC dataset with consensus from multiple methods"
    }]
})
```

### When Each Specialist is Consulted

**Hematopoiesis Specialist:**
- Tissue: BM, PBMC, Spleen
- Questions: Progenitor vs mature, lineage ambiguity, emergency myelopoiesis

**Infection/Inflammation Specialist:**
- Context: Infection experiments, viral/bacterial challenge
- Questions: ISG-high clusters, inflammatory vs resting states, stress vs inflammation

**Aging Specialist:**
- Context: Young vs aged comparisons
- Questions: Age-associated programs, exhaustion vs activation, inflammaging

**Cell State/Plasticity Specialist:**
- Context: Time-course, activation experiments
- Questions: Cycling cells, effector vs memory, activation states

**Epithelial Specialist:**
- Context: Epithelial tissues (lung, gut, skin)
- Questions: Epithelial subtypes, stress vs immune contamination

**Organoid/Stemness Specialist:**
- Context: Organoid cultures, stem cell datasets
- Questions: Culture artifacts, stemness vs differentiation, WNT/YAP/Notch states

---

## 2. Comprehensive QC Visualization

### Overview

Generate publication-quality QC plots showing data quality before and after filtering.

### Features

- **Before QC plots**: Distribution of n_genes, n_counts, pct_counts_mt
- **After QC plots**: Same metrics after filtering
- **Comparison plots**: Side-by-side before/after visualization
- **Statistical summaries**: Median, mean, filtering statistics

### Usage

#### Individual Plot Generation

```python
from scrna_agent.tools import (
    plot_qc_metrics_before,
    plot_qc_metrics_after,
    plot_qc_comparison
)

# Generate before QC plots
plot_qc_metrics_before.invoke({
    "adata_path": "data_raw.h5ad",
    "output_dir": "./qc_plots"
})

# Generate after QC plots
plot_qc_metrics_after.invoke({
    "adata_path": "data_qc.h5ad",
    "output_dir": "./qc_plots"
})

# Generate comparison plots
plot_qc_comparison.invoke({
    "adata_before_path": "data_raw.h5ad",
    "adata_after_path": "data_qc.h5ad",
    "output_dir": "./qc_plots"
})
```

#### Orchestrated QC Visualization

```python
from scrna_agent.tools import generate_qc_visualizations

# Generate all QC plots in one call
result = generate_qc_visualizations.invoke({
    "adata_before_path": "data_raw.h5ad",
    "adata_after_path": "data_qc.h5ad",
    "output_dir": "./qc_plots"
})

print(result)
```

### Output Files

- `qc_before_filtering.png` - QC metrics before filtering
- `qc_after_filtering.png` - QC metrics after filtering
- `qc_before_after_comparison.png` - Side-by-side comparison

---

## 3. Integration Visualization

### Overview

Visualize batch integration quality and batch mixing.

### Features

- **Before integration plots**: Batch clustering (showing batch effect)
- **After integration plots**: Batch mixing (showing correction)
- **Comparison plots**: Before/after side-by-side
- **Batch mixing score**: Quantitative mixing metric

### Usage

#### Individual Plot Generation

```python
from scrna_agent.tools import (
    plot_integration_before,
    plot_integration_after,
    plot_integration_comparison
)

# Before integration
plot_integration_before.invoke({
    "adata_path": "data_before_integration.h5ad",
    "batch_key": "batch",
    "output_dir": "./integration_plots"
})

# After integration
plot_integration_after.invoke({
    "adata_path": "data_integrated.h5ad",
    "batch_key": "batch",
    "celltype_key": "cell_type",
    "embedding_key": "X_scVI",
    "output_dir": "./integration_plots"
})

# Comparison
plot_integration_comparison.invoke({
    "adata_before_path": "data_before_integration.h5ad",
    "adata_after_path": "data_integrated.h5ad",
    "batch_key": "batch",
    "output_dir": "./integration_plots"
})
```

#### Orchestrated Integration Visualization

```python
from scrna_agent.tools import generate_integration_visualizations

# Generate all integration plots
result = generate_integration_visualizations.invoke({
    "adata_before_path": "data_before_integration.h5ad",
    "adata_after_path": "data_integrated.h5ad",
    "batch_key": "batch",
    "celltype_key": "cell_type",
    "embedding_key": "X_scVI",
    "output_dir": "./integration_plots"
})
```

### Output Files

- `integration_before.png` - Pre-integration analysis
- `integration_after.png` - Post-integration analysis
- `integration_before_after_comparison.png` - Comparison

---

## 4. Condition-wise Analysis

### Overview

Analyze how cell type composition changes across experimental conditions.

### Use Cases

- Young vs Aged comparisons
- Treatment vs Control
- Disease vs Healthy
- Time-course experiments
- Multiple conditions simultaneously

### Usage

```python
from scrna_agent.tools import (
    plot_condition_celltype_distribution,
    generate_condition_analysis
)

# Single condition analysis
plot_condition_celltype_distribution.invoke({
    "adata_path": "data_final.h5ad",
    "condition_key": "age",  # e.g., 'young' vs 'aged'
    "celltype_key": "cell_type",
    "output_dir": "./condition_analysis"
})

# Multiple conditions
generate_condition_analysis.invoke({
    "adata_path": "data_final.h5ad",
    "condition_keys": ["age", "treatment", "timepoint"],
    "celltype_key": "cell_type",
    "output_dir": "./condition_analysis"
})
```

### Output Features

- **Stacked bar plots**: Cell type proportions per condition
- **Heatmaps**: Absolute counts and percentages
- **Fold change plots**: Enrichment/depletion (for 2 conditions)
- **Statistical tables**: Detailed counts and proportions

### Output Files

- `condition_age_celltype_distribution.png`
- `condition_treatment_celltype_distribution.png`
- etc.

---

## 5. Comprehensive Reporting

### Overview

Generate detailed markdown reports with all analysis results and figures.

### Features

- **Executive summary**: High-level overview
- **QC analysis**: Before/after comparison with statistics
- **Annotation results**: Cell type distribution and expert interpretations
- **Integration results**: Batch mixing and quality metrics
- **Condition analysis**: Comparative cell type distribution
- **All figures**: Embedded with captions
- **Methods section**: Documentation of analysis steps

### Usage

```python
from scrna_agent.tools import (
    initialize_report,
    add_qc_metrics,
    add_annotation_results,
    add_integration_results,
    add_figure,
    generate_markdown_report
)

# Initialize report
initialize_report.invoke({
    "project_name": "PBMC_Young_vs_Aged",
    "output_dir": "./reports"
})

# Add QC metrics
add_qc_metrics.invoke({
    "n_cells_before": 50000,
    "n_cells_after": 45000,
    "median_genes_before": 2000,
    "median_genes_after": 2500,
    "median_counts_before": 5000,
    "median_counts_after": 6000,
    "median_mt_before": 8.5,
    "median_mt_after": 5.2
})

# Add annotation results
add_annotation_results.invoke({
    "n_cell_types": 15,
    "annotation_method": "CellTypist + ScType + Literature RAG (Consensus)",
    "consensus_rate": 92.5,
    "cell_type_counts": {
        "CD4 T cells": 12000,
        "CD8 T cells": 8000,
        "B cells": 5000,
        # ... more cell types
    }
})

# Add integration results
add_integration_results.invoke({
    "n_batches": 4,
    "integration_method": "scVI",
    "batch_distribution": {"Batch1": 11000, "Batch2": 11500, ...}
})

# Add figures
add_figure.invoke({
    "figure_name": "QC Comparison",
    "figure_path": "./qc_plots/qc_before_after_comparison.png",
    "figure_caption": "Quality control metrics before and after filtering",
    "figure_category": "qc"
})

# Generate final report
generate_markdown_report.invoke({
    "output_file": "./reports/analysis_report.md"
})
```

### Output Files

- `analysis_report.md` - Markdown report
- `analysis_report.json` - JSON data export

---

## 6. Complete Visualization Pipeline

### Overview

Master orchestration tool that runs the entire visualization and reporting pipeline.

### What It Does

1. Generates all QC visualizations
2. Generates all integration visualizations
3. Generates condition-wise analysis (if requested)
4. Generates comprehensive report
5. Organizes all outputs in structured directories

### Usage

```python
from scrna_agent.tools import run_complete_visualization_pipeline

result = run_complete_visualization_pipeline.invoke({
    "project_name": "PBMC_Young_vs_Aged",
    "adata_raw_path": "data_raw.h5ad",
    "adata_qc_path": "data_qc.h5ad",
    "adata_before_integration_path": "data_normalized.h5ad",
    "adata_final_path": "data_integrated_annotated.h5ad",
    "batch_key": "batch",
    "celltype_key": "cell_type",
    "condition_keys": ["age", "sex"],  # Optional
    "embedding_key": "X_scVI",
    "base_output_dir": "./analysis_output"
})

print(result)
```

### Output Structure

```
analysis_output/
├── qc_plots/
│   ├── qc_before_filtering.png
│   ├── qc_after_filtering.png
│   └── qc_before_after_comparison.png
├── integration_plots/
│   ├── integration_before.png
│   ├── integration_after.png
│   └── integration_before_after_comparison.png
├── condition_analysis/
│   ├── condition_age_celltype_distribution.png
│   └── condition_sex_celltype_distribution.png
└── reports/
    ├── PBMC_Young_vs_Aged_comprehensive_report.md
    └── PBMC_Young_vs_Aged_comprehensive_report.json
```

---

## Example Workflows

### Workflow 1: Complete Analysis with All Features

```python
# 1. Run main analysis pipeline (existing functionality)
# ... QC, normalization, clustering, annotation, integration ...

# 2. Generate all visualizations and reports
from scrna_agent.tools import run_complete_visualization_pipeline

run_complete_visualization_pipeline.invoke({
    "project_name": "MyProject",
    "adata_raw_path": "raw.h5ad",
    "adata_qc_path": "qc.h5ad",
    "adata_before_integration_path": "normalized.h5ad",
    "adata_final_path": "final_integrated.h5ad",
    "batch_key": "sample",
    "celltype_key": "consensus_annotation",
    "condition_keys": ["age", "treatment"],
    "embedding_key": "X_scVI",
    "base_output_dir": "./outputs"
})
```

### Workflow 2: Custom Visualization Pipeline

```python
from scrna_agent.tools import (
    generate_qc_visualizations,
    generate_integration_visualizations,
    generate_condition_analysis
)

# QC only
generate_qc_visualizations.invoke({
    "adata_before_path": "raw.h5ad",
    "adata_after_path": "qc.h5ad",
    "output_dir": "./qc"
})

# Integration only
generate_integration_visualizations.invoke({
    "adata_before_path": "before_integration.h5ad",
    "adata_after_path": "after_integration.h5ad",
    "batch_key": "batch",
    "celltype_key": "cell_type",
    "embedding_key": "X_scVI",
    "output_dir": "./integration"
})

# Condition analysis
generate_condition_analysis.invoke({
    "adata_path": "final.h5ad",
    "condition_keys": ["condition"],
    "celltype_key": "cell_type",
    "output_dir": "./conditions"
})
```

### Workflow 3: Agent-Based Annotation with Domain Experts

```python
from scrna_agent.agents import create_scrna_pipeline_agent

# Create pipeline agent (includes annotation coordinator with domain experts)
agent = create_scrna_pipeline_agent(model_name="gpt-4o")

# Run full pipeline with automatic domain expert consultation
agent.invoke({
    "messages": [{
        "role": "user",
        "content": """
        Run complete pipeline:
        1. Load data from processed.h5ad
        2. Annotate with consensus (CellTypist + ScType + Literature)
        3. Consult domain experts for ambiguous clusters
        4. Integrate batches with scVI
        5. Generate all visualizations and reports
        """
    }]
})
```

---

## Tips and Best Practices

### For Domain Specialists

1. **Be specific about tissue type**: Helps specialists provide better interpretations
2. **Include experimental context**: Infection, aging, time-course, etc.
3. **Report low-confidence clusters**: Domain experts excel at resolving ambiguity

### For Visualizations

1. **Always generate comparison plots**: Before/after comparisons are critical for QC
2. **Check batch mixing**: Visual inspection is essential for integration quality
3. **Use high DPI**: Plots are saved at 300 DPI for publication quality

### For Condition Analysis

1. **Use meaningful condition keys**: "age", "treatment", "timepoint", etc.
2. **Compare 2-3 conditions**: More conditions = more complex plots
3. **Check absolute counts AND proportions**: Both perspectives are important

### For Reports

1. **Add all figures**: Comprehensive reports should include all generated plots
2. **Include expert interpretations**: Add domain expert insights to annotation results
3. **Export both MD and JSON**: Markdown for reading, JSON for programmatic access

---

## Troubleshooting

### Domain Specialists Not Working

**Issue**: Specialists not being consulted during annotation

**Solution**: Ensure you're using the `create_annotation_coordinator()` agent, not individual annotation tools directly.

### Missing Plots

**Issue**: Some plots not generated

**Solution**: Check that required columns exist (batch_key, celltype_key, etc.) in adata.obs

### Low Quality Plots

**Issue**: Plots look pixelated or low resolution

**Solution**: Plots are saved at 300 DPI by default. For higher quality, modify the `_setup_plot_style()` function.

### Integration Comparison Fails

**Issue**: Error when comparing before/after integration

**Solution**: Ensure both AnnData objects have UMAP embeddings. If not, they will be computed automatically.

### Condition Analysis Empty

**Issue**: No condition plots generated

**Solution**: Verify that `condition_key` exists in adata.obs and has expected values.

---

## API Reference

See individual tool docstrings for detailed parameter descriptions:

```python
# QC Visualization
help(plot_qc_metrics_before)
help(plot_qc_metrics_after)
help(plot_qc_comparison)

# Integration Visualization
help(plot_integration_before)
help(plot_integration_after)
help(plot_integration_comparison)
help(plot_condition_celltype_distribution)

# Comprehensive Reporting
help(initialize_report)
help(add_qc_metrics)
help(add_annotation_results)
help(add_integration_results)
help(generate_markdown_report)

# Orchestration
help(generate_qc_visualizations)
help(generate_integration_visualizations)
help(generate_condition_analysis)
help(run_complete_visualization_pipeline)

# Domain Specialists
from scrna_agent.agents import (
    create_hematopoiesis_specialist,
    create_infection_inflammation_specialist,
    create_aging_specialist,
    create_cellstate_plasticity_specialist,
    create_epithelial_specialist,
    create_organoid_stemness_specialist
)
help(create_hematopoiesis_specialist)
# ... etc
```

---

## Summary

The new features provide:

✅ **6 domain specialist agents** for expert biological interpretation
✅ **Comprehensive QC visualization** with before/after comparison
✅ **Integration visualization** with batch mixing analysis
✅ **Condition-wise analysis** for experimental comparisons
✅ **Automated comprehensive reports** with all figures
✅ **End-to-end pipeline orchestration** for convenience

All features are production-ready and fully integrated with the existing scRNA-agent pipeline.
