#!/usr/bin/env python3
"""
Generate Scientific Visualizations using Claude Skills
Prepares specifications for the scientific-visualization skill
"""

import json
from pathlib import Path
import pandas as pd

def create_visualization_specifications(results_dir, adata_file=None):
    """
    Create detailed specifications for scientific visualizations

    Parameters:
    -----------
    results_dir : Path
        Results directory
    adata_file : Path, optional
        Path to h5ad file for data access

    Returns:
    --------
    specifications : dict
        Visualization specifications
    """

    results_dir = Path(results_dir)
    figures_dir = results_dir / 'figures'

    print("="*80)
    print("CREATING SCIENTIFIC VISUALIZATION SPECIFICATIONS")
    print("="*80)

    specifications = {
        'project': 'scRNA-seq Analysis',
        'output_directory': str(figures_dir),
        'style': {
            'theme': 'scientific',
            'dpi': 300,
            'format': 'png',
            'colorblind_safe': True,
            'font_family': 'Arial',
            'figure_facecolor': 'white'
        },
        'visualizations': []
    }

    # Visualization 1: QC Metrics Multi-panel
    specifications['visualizations'].append({
        'name': 'qc_metrics_comprehensive',
        'type': 'multi_panel',
        'title': 'Quality Control Metrics',
        'description': 'Comprehensive QC metrics including genes, UMIs, and mitochondrial percentage',
        'panels': [
            {
                'type': 'violin',
                'data_source': 'adata.obs',
                'x': 'condition',
                'y': 'n_genes_by_counts',
                'title': 'Genes per Cell',
                'ylabel': 'Number of Genes',
                'color_palette': 'Set2'
            },
            {
                'type': 'violin',
                'data_source': 'adata.obs',
                'x': 'condition',
                'y': 'total_counts',
                'title': 'Total UMI Counts',
                'ylabel': 'UMI Counts',
                'color_palette': 'Set2'
            },
            {
                'type': 'violin',
                'data_source': 'adata.obs',
                'x': 'condition',
                'y': 'pct_counts_mt',
                'title': 'Mitochondrial Percentage',
                'ylabel': '% Mitochondrial',
                'color_palette': 'Set2',
                'hline': {'y': 20, 'color': 'red', 'linestyle': '--', 'label': 'Threshold'}
            }
        ],
        'layout': {
            'nrows': 2,
            'ncols': 3,
            'figsize': (18, 10),
            'tight_layout': True
        },
        'output_file': 'qc_metrics_publication.png'
    })

    # Visualization 2: UMAP with Multiple Annotations
    specifications['visualizations'].append({
        'name': 'umap_comprehensive',
        'type': 'multi_panel',
        'title': 'UMAP Dimensionality Reduction',
        'description': 'UMAP projections colored by various annotations',
        'panels': [
            {
                'type': 'scatter',
                'data_source': 'adata.obsm["X_umap"]',
                'color_by': 'cell_type',
                'title': 'Cell Types',
                'legend_location': 'right',
                'point_size': 1,
                'alpha': 0.8
            },
            {
                'type': 'scatter',
                'data_source': 'adata.obsm["X_umap"]',
                'color_by': 'condition',
                'title': 'Experimental Conditions',
                'legend_location': 'right',
                'point_size': 1,
                'alpha': 0.8
            },
            {
                'type': 'scatter',
                'data_source': 'adata.obsm["X_umap"]',
                'color_by': 'n_genes_by_counts',
                'title': 'Gene Count',
                'colormap': 'viridis',
                'colorbar': True,
                'point_size': 1,
                'alpha': 0.8
            },
            {
                'type': 'scatter',
                'data_source': 'adata.obsm["X_umap"]',
                'color_by': 'pct_counts_mt',
                'title': 'Mitochondrial %',
                'colormap': 'plasma',
                'colorbar': True,
                'point_size': 1,
                'alpha': 0.8
            }
        ],
        'layout': {
            'nrows': 2,
            'ncols': 2,
            'figsize': (16, 14),
            'tight_layout': True
        },
        'output_file': 'umap_comprehensive_publication.png'
    })

    # Visualization 3: Marker Gene Expression Heatmap
    specifications['visualizations'].append({
        'name': 'marker_heatmap',
        'type': 'heatmap',
        'title': 'Top Marker Genes by Cluster',
        'description': 'Heatmap showing expression of top marker genes across clusters',
        'data_source': 'marker_genes',
        'clustering': {
            'rows': True,
            'cols': True,
            'method': 'average'
        },
        'colormap': 'RdBu_r',
        'center': 0,
        'cbar_label': 'Scaled Expression',
        'figsize': (12, 16),
        'output_file': 'marker_genes_heatmap_publication.png'
    })

    # Visualization 4: Cell Type Proportions
    specifications['visualizations'].append({
        'name': 'cell_type_proportions',
        'type': 'stacked_bar',
        'title': 'Cell Type Composition by Condition',
        'description': 'Stacked bar plot showing cell type proportions across conditions',
        'data_source': 'adata.obs',
        'x': 'condition',
        'fill': 'cell_type',
        'normalize': True,
        'color_palette': 'tab20',
        'ylabel': 'Proportion (%)',
        'legend_location': 'right',
        'figsize': (12, 8),
        'output_file': 'cell_type_proportions_publication.png'
    })

    # Visualization 5: Marker Gene Dotplot
    specifications['visualizations'].append({
        'name': 'marker_dotplot',
        'type': 'dotplot',
        'title': 'Marker Gene Expression Pattern',
        'description': 'Dot plot showing marker gene expression across cell types',
        'genes': 'top_markers_per_cluster',
        'groupby': 'cell_type',
        'size_title': '% Expressed',
        'color_title': 'Mean Expression',
        'figsize': (16, 10),
        'output_file': 'marker_dotplot_publication.png'
    })

    # Visualization 6: Volcano Plot for Differential Expression
    specifications['visualizations'].append({
        'name': 'volcano_plot',
        'type': 'volcano',
        'title': 'Differential Expression Analysis',
        'description': 'Volcano plot showing differentially expressed genes',
        'x': 'log2_fold_change',
        'y': 'neg_log10_pvalue',
        'significance_thresholds': {
            'fold_change': 1.0,
            'pvalue': 0.05
        },
        'top_genes_label': 10,
        'figsize': (10, 8),
        'output_file': 'volcano_plot_publication.png'
    })

    # Visualization 7: PCA with Loadings
    specifications['visualizations'].append({
        'name': 'pca_biplot',
        'type': 'pca_biplot',
        'title': 'PCA Biplot with Gene Loadings',
        'description': 'PCA visualization with top contributing genes',
        'components': [0, 1],
        'color_by': 'condition',
        'top_loadings': 10,
        'arrow_scale': 0.3,
        'figsize': (12, 10),
        'output_file': 'pca_biplot_publication.png'
    })

    # Visualization 8: Gene Expression Violin Plots
    specifications['visualizations'].append({
        'name': 'gene_expression_violins',
        'type': 'multi_violin',
        'title': 'Key Gene Expression Across Cell Types',
        'description': 'Violin plots for important marker genes',
        'genes': ['ACE2', 'TMPRSS2', 'DPP4', 'TP63', 'FOXJ1', 'SCGB1A1'],
        'groupby': 'cell_type',
        'layout': {
            'nrows': 2,
            'ncols': 3,
            'figsize': (18, 10)
        },
        'output_file': 'gene_expression_violins_publication.png'
    })

    # Save specifications
    spec_file = results_dir / 'visualization_specifications.json'
    with open(spec_file, 'w') as f:
        json.dump(specifications, f, indent=2)

    print(f"\n✓ Saved specifications: {spec_file.name}")
    print(f"  Generated {len(specifications['visualizations'])} visualization specifications")

    return specifications, spec_file


def create_skill_instructions(spec_file, adata_file, results_dir):
    """
    Create instructions for using the scientific-visualization skill

    Parameters:
    -----------
    spec_file : Path
        Path to visualization specifications
    adata_file : Path
        Path to AnnData file
    results_dir : Path
        Results directory
    """

    results_dir = Path(results_dir)

    instructions = f"""# Scientific Visualization Generation

## Overview
Generate publication-quality scientific visualizations using the scientific-visualization skill.

## Data Source
- **AnnData File:** {adata_file}
- **Results Directory:** {results_dir}
- **Specifications:** {spec_file}

## Visualization Requirements

### Style Guidelines
- **DPI:** 300 (publication quality)
- **Format:** PNG
- **Color Scheme:** Colorblind-safe palettes
- **Font:** Arial, consistent sizing
- **Background:** White
- **Grid:** Subtle or none
- **Legend:** Clear and readable
- **Axis Labels:** Descriptive with units

### Figures to Generate

1. **QC Metrics Panel** (qc_metrics_publication.png)
   - Multi-panel violin plots
   - Genes per cell, UMI counts, MT percentage
   - Grouped by condition
   - QC thresholds indicated

2. **UMAP Comprehensive** (umap_comprehensive_publication.png)
   - 4-panel UMAP visualization
   - Cell types, conditions, QC metrics
   - Consistent point size and transparency
   - Publication-ready layout

3. **Marker Gene Heatmap** (marker_genes_heatmap_publication.png)
   - Hierarchically clustered heatmap
   - Top markers per cluster
   - Scaled expression values
   - Clear row and column labels

4. **Cell Type Proportions** (cell_type_proportions_publication.png)
   - Stacked bar plot
   - Normalized proportions
   - Clear legend
   - Error bars if replicates available

5. **Marker Gene Dotplot** (marker_dotplot_publication.png)
   - Dot size = % cells expressing
   - Color = mean expression level
   - Genes vs cell types
   - Clear dual legend

6. **Volcano Plot** (volcano_plot_publication.png)
   - Log2 fold change vs -log10(p-value)
   - Significance thresholds
   - Top genes labeled
   - Color-coded by significance

7. **PCA Biplot** (pca_biplot_publication.png)
   - PCA scatter with loadings
   - Top contributing genes as arrows
   - Colored by condition
   - Variance explained in axis labels

8. **Gene Expression Violins** (gene_expression_violins_publication.png)
   - Multi-panel violin plots
   - Key marker genes
   - Grouped by cell type
   - Statistical comparisons

## Execution Instructions

### Method 1: Using /scientific-visualization skill in Claude Code

```bash
# In Claude Code CLI
/scientific-visualization

# When prompted, provide:
# 1. Data file: {adata_file}
# 2. Specifications: {spec_file}
# 3. Output directory: {results_dir}/figures/
```

The skill will:
- Load the AnnData object
- Parse visualization specifications
- Generate all figures with publication quality
- Apply consistent styling
- Save to specified locations

### Method 2: Batch generation script

Create a file `generate_all_visualizations.txt`:

```
Please use the scientific-visualization skill to generate publication-quality figures.

Data: {adata_file}
Specifications: {spec_file}
Output: {results_dir}/figures/

Generate all visualizations specified in the JSON file with:
- Publication quality (300 DPI)
- Colorblind-safe palettes
- Consistent styling
- Professional layout
- Clear labels and legends
```

Then run:
```bash
claude < generate_all_visualizations.txt
```

## Quality Checklist

For each figure, verify:
- [ ] Resolution is 300 DPI or higher
- [ ] Colors are colorblind-friendly
- [ ] All text is readable (minimum 8pt)
- [ ] Axes are properly labeled with units
- [ ] Legend is clear and complete
- [ ] Figure fits publication column width
- [ ] No overlapping labels
- [ ] Consistent style across figures
- [ ] White background
- [ ] Vector-compatible elements

## Output Files

All figures will be saved in: `{results_dir}/figures/`

Publication-quality figures:
- qc_metrics_publication.png
- umap_comprehensive_publication.png
- marker_genes_heatmap_publication.png
- cell_type_proportions_publication.png
- marker_dotplot_publication.png
- volcano_plot_publication.png
- pca_biplot_publication.png
- gene_expression_violins_publication.png

## Integration with PowerPoint

After generating visualizations, update the PowerPoint:

```bash
python generate_scientific_ppt.py
```

This will:
- Use the scientific-slides skill
- Incorporate new publication-quality figures
- Create professional presentation
- Apply scientific styling

## Notes

- Visualizations are optimized for both digital and print
- Color schemes follow scientific publication standards
- All figures include source data information
- Statistical tests are indicated where appropriate
- Figure legends are comprehensive

## Timeline

- Setup: 5 minutes
- Generation per figure: 2-3 minutes
- Total: ~20-30 minutes for all visualizations
- PowerPoint generation: 5-10 minutes

## Support

If visualizations need adjustments:
1. Modify visualization_specifications.json
2. Re-run the skill
3. Specific figures can be regenerated individually
"""

    instruction_file = results_dir / 'VISUALIZATION_INSTRUCTIONS.md'
    with open(instruction_file, 'w') as f:
        f.write(instructions)

    print(f"✓ Saved instructions: {instruction_file.name}")

    return instruction_file


def main():
    """Main workflow for scientific visualization generation"""

    print("="*80)
    print("SCIENTIFIC VISUALIZATION WORKFLOW")
    print("="*80)

    # Configuration
    results_dir = Path('/mnt2/kwy/scrna_agent/JMV_ALI_results')
    adata_file = results_dir / 'JMV_ALI_harmony_integrated.h5ad'

    # Check if annotated version exists
    annotated_file = results_dir / 'JMV_ALI_harmony_integrated_annotated.h5ad'
    if annotated_file.exists():
        adata_file = annotated_file
        print(f"\nUsing annotated data: {adata_file.name}")
    else:
        print(f"\nUsing integrated data: {adata_file.name}")

    # Create specifications
    print("\n[1/3] Creating visualization specifications...")
    specifications, spec_file = create_visualization_specifications(
        results_dir=results_dir,
        adata_file=adata_file
    )

    # Create instructions
    print("\n[2/3] Creating skill usage instructions...")
    instruction_file = create_skill_instructions(
        spec_file=spec_file,
        adata_file=adata_file,
        results_dir=results_dir
    )

    # Summary
    print("\n[3/3] Workflow setup complete!")
    print("="*80)
    print("NEXT STEPS")
    print("="*80)
    print("\n1. Generate Visualizations:")
    print("   In Claude Code, run: /scientific-visualization")
    print(f"   Provide specification file: {spec_file}")
    print(f"   Provide data file: {adata_file}")
    print("")
    print("2. Generate PowerPoint:")
    print("   In Claude Code, run: /scientific-slides")
    print("   Or run: python generate_scientific_ppt.py")
    print("")
    print("3. Review and Refine:")
    print(f"   Check figures in: {results_dir}/figures/")
    print(f"   Follow instructions in: {instruction_file.name}")
    print("")
    print("All visualizations will be publication-quality (300 DPI)")
    print("with professional scientific styling!")
    print("="*80)


if __name__ == "__main__":
    main()
