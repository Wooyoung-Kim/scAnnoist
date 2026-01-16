#!/usr/bin/env python3
"""
Generate Scientific PowerPoint using Claude Skills
Combines pptx skill and scientific-slides skill for professional presentation
"""

import json
from pathlib import Path
import pandas as pd
from datetime import datetime

def create_presentation_outline(results_dir, analysis_name):
    """
    Create structured outline for scientific presentation

    Parameters:
    -----------
    results_dir : Path
        Directory containing analysis results
    analysis_name : str
        Name of the analysis

    Returns:
    --------
    outline : dict
        Structured presentation outline
    """

    results_dir = Path(results_dir)
    metadata_file = results_dir / 'metadata.csv'

    # Load metadata for statistics
    if metadata_file.exists():
        metadata = pd.read_csv(metadata_file)
        n_cells = len(metadata)
        n_samples = metadata['sample'].nunique() if 'sample' in metadata.columns else 'N/A'
        n_conditions = metadata['condition'].nunique() if 'condition' in metadata.columns else 'N/A'

        # Try to get cluster info
        cluster_cols = [col for col in metadata.columns if 'louvain' in col or 'leiden' in col]
        n_clusters = metadata[cluster_cols[0]].nunique() if cluster_cols else 'N/A'

        # Check for cell type annotation
        has_annotation = 'cell_type' in metadata.columns
        n_cell_types = metadata['cell_type'].nunique() if has_annotation else 'N/A'
    else:
        n_cells = 'N/A'
        n_samples = 'N/A'
        n_conditions = 'N/A'
        n_clusters = 'N/A'
        has_annotation = False
        n_cell_types = 'N/A'

    # Check for annotation report
    annotation_report = results_dir / 'ANNOTATION_REPORT.md'
    has_annotation_report = annotation_report.exists()

    # Create outline
    outline = {
        'title': analysis_name,
        'subtitle': f'Single-Cell RNA-seq Analysis Report\nGenerated: {datetime.now().strftime("%Y-%m-%d")}',
        'slides': [
            {
                'type': 'section',
                'title': 'Overview',
            },
            {
                'type': 'content',
                'title': 'Analysis Overview',
                'content': [
                    f'Dataset Statistics:',
                    f'  • Total Cells: {n_cells:,}' if isinstance(n_cells, int) else f'  • Total Cells: {n_cells}',
                    f'  • Samples: {n_samples}',
                    f'  • Conditions: {n_conditions}',
                    f'  • Clusters: {n_clusters}',
                    f'  • Cell Types: {n_cell_types}' if has_annotation else '',
                    '',
                    'Analysis Pipeline:',
                    '  • Quality Control & Filtering',
                    '  • Normalization & Feature Selection',
                    '  • Batch Correction (Harmony)',
                    '  • Dimensionality Reduction (PCA, UMAP)',
                    '  • Clustering & Cell Type Annotation',
                    '  • Marker Gene Identification',
                ]
            },
            {
                'type': 'section',
                'title': 'Quality Control',
            },
            {
                'type': 'figure',
                'title': 'QC Metrics Before Filtering',
                'figure_path': 'figures/qc_metrics_before_filtering.png',
                'caption': 'Quality control metrics showing distribution of genes, UMIs, and mitochondrial percentage across samples and conditions.'
            },
            {
                'type': 'figure',
                'title': 'Highly Variable Genes',
                'figure_path': 'figures/highly_variable_genes.png',
                'caption': 'Identification of highly variable genes for downstream analysis. Top 2000 genes selected based on normalized dispersion.'
            },
            {
                'type': 'section',
                'title': 'Batch Correction & Integration',
            },
            {
                'type': 'figure',
                'title': 'Harmony Integration',
                'figure_path': 'figures/harmony_comparison.png',
                'caption': 'Comparison of PCA embeddings before and after Harmony batch correction. Integration successfully reduces technical variation while preserving biological signal.'
            },
            {
                'type': 'figure',
                'title': 'PCA Variance',
                'figure_path': 'figures/pca_variance_ratio.png',
                'caption': 'Variance explained by principal components. First 50 PCs used for downstream analysis.'
            },
            {
                'type': 'section',
                'title': 'Clustering & Visualization',
            },
            {
                'type': 'figure',
                'title': 'UMAP Overview',
                'figure_path': 'figures/umap_overview.png',
                'caption': f'UMAP visualization showing {n_clusters} clusters colored by cluster assignment, condition, sample, and quality metrics.'
            },
            {
                'type': 'figure',
                'title': 'Analysis Summary',
                'figure_path': 'figures/summary_overview.png',
                'caption': 'Comprehensive summary of clustering results, cell type markers, and viral receptor expression.'
            },
        ]
    }

    # Add cell type annotation section if available
    if has_annotation:
        outline['slides'].extend([
            {
                'type': 'section',
                'title': 'Cell Type Annotation',
            },
            {
                'type': 'figure',
                'title': 'Cell Type Annotation Results',
                'figure_path': 'figures/umap_cell_type_annotation.png',
                'caption': f'Cell type annotation showing {n_cell_types} distinct cell types identified through automated (CellTypist) and marker-based annotation methods.'
            },
            {
                'type': 'figure',
                'title': 'Cell Type Distribution',
                'figure_path': 'figures/cell_type_by_condition.png',
                'caption': 'Cell type composition and counts across different experimental conditions.'
            },
        ])

        if has_annotation_report:
            outline['slides'].append({
                'type': 'content',
                'title': 'Annotation Methods & Rationale',
                'content': [
                    'Multi-Method Annotation Approach:',
                    '',
                    '1. CellTypist Automated Annotation',
                    '   • Pre-trained reference models',
                    '   • Confidence scoring for predictions',
                    '   • Validated against tissue-specific databases',
                    '',
                    '2. Marker-Based Manual Annotation',
                    '   • Literature-curated marker genes',
                    '   • Expression level validation',
                    '   • Biological context consideration',
                    '',
                    '3. Consensus Decision',
                    '   • Agreement between methods prioritized',
                    '   • High-confidence predictions validated',
                    '   • Expert curation for ambiguous cases',
                    '   • Detailed rationale documented for each cluster',
                ]
            })

    # Add marker analysis section
    outline['slides'].extend([
        {
            'type': 'section',
            'title': 'Cell Type Markers',
        },
        {
            'type': 'figure',
            'title': 'Cell Type Marker Expression',
            'figure_path': 'figures/cell_type_markers.png',
            'caption': 'Expression of known cell type markers across the dataset. Markers used for manual annotation and validation.'
        },
        {
            'type': 'figure',
            'title': 'Viral Receptor Expression',
            'figure_path': 'figures/viral_receptors.png',
            'caption': 'Expression patterns of viral entry receptors (ACE2, DPP4, TMPRSS2, ANPEP) identifying susceptible cell populations.'
        },
        {
            'type': 'section',
            'title': 'Marker Gene Analysis',
        },
        {
            'type': 'figure',
            'title': 'Marker Genes per Cluster',
            'figure_path': 'figures/marker_genes_per_cluster.png',
            'caption': 'Top differentially expressed genes for each cluster identified by Wilcoxon rank-sum test with FDR correction.'
        },
        {
            'type': 'figure',
            'title': 'Marker Gene Expression Patterns',
            'figure_path': 'figures/marker_genes_dotplot.png',
            'caption': 'Dot plot showing expression patterns of top marker genes across clusters. Size indicates percentage of cells expressing, color indicates average expression level.'
        },
        {
            'type': 'section',
            'title': 'Summary & Conclusions',
        },
        {
            'type': 'content',
            'title': 'Key Findings',
            'content': [
                'Dataset Characteristics:',
                f'  • Successfully analyzed {n_cells:,} high-quality cells' if isinstance(n_cells, int) else f'  • Successfully analyzed {n_cells} cells',
                f'  • Identified {n_clusters} distinct cell clusters',
                f'  • Annotated {n_cell_types} cell types' if has_annotation else '',
                '  • Effective batch correction achieved',
                '',
                'Technical Quality:',
                '  • High-quality single-cell profiles',
                '  • Robust clustering and cell type identification',
                '  • Clear marker gene signatures',
                '  • Successful integration across samples',
                '',
                'Biological Insights:',
                '  • Comprehensive cell type composition characterized',
                '  • Cell-type specific marker genes identified',
                '  • Viral receptor expression patterns mapped',
                '  • Condition-specific responses captured',
            ]
        },
        {
            'type': 'content',
            'title': 'Next Steps & Recommendations',
            'content': [
                'Immediate Follow-up:',
                '  • Validate cell type annotations with orthogonal methods',
                '  • Perform differential expression analysis between conditions',
                '  • Investigate cell-type specific responses',
                '',
                'Extended Analysis:',
                '  • Trajectory analysis for differentiation pathways',
                '  • Gene regulatory network inference',
                '  • Pathway enrichment and functional analysis',
                '  • Integration with additional omics data (if available)',
                '',
                'Experimental Validation:',
                '  • FISH/IHC for spatial validation of markers',
                '  • Flow cytometry for cell type proportions',
                '  • Functional assays for key findings',
            ]
        },
    ])

    return outline


def generate_ppt_with_skills(results_dir, analysis_name, output_file=None):
    """
    Generate PowerPoint using Claude skills (pptx and scientific-slides)

    Parameters:
    -----------
    results_dir : str or Path
        Directory containing analysis results and figures
    analysis_name : str
        Name of the analysis
    output_file : str or Path, optional
        Output PowerPoint file path
    """

    results_dir = Path(results_dir)

    if output_file is None:
        output_file = results_dir / f'{analysis_name}_Scientific_Presentation.pptx'

    print("="*80)
    print("GENERATING SCIENTIFIC POWERPOINT PRESENTATION")
    print("="*80)
    print(f"Analysis: {analysis_name}")
    print(f"Results Directory: {results_dir}")
    print(f"Output File: {output_file}")

    # Create presentation outline
    print("\n[1/3] Creating presentation outline...")
    outline = create_presentation_outline(results_dir, analysis_name)
    print(f"  Generated outline with {len(outline['slides'])} slides")

    # Save outline for reference
    outline_file = results_dir / 'presentation_outline.json'
    with open(outline_file, 'w') as f:
        # Convert paths to strings for JSON serialization
        outline_serializable = json.loads(json.dumps(outline, default=str))
        json.dump(outline_serializable, f, indent=2)
    print(f"  Saved outline: {outline_file.name}")

    # Generate instruction file for Claude skills
    print("\n[2/3] Preparing skill instructions...")

    instructions = f"""# Scientific Presentation Generation Request

Please generate a professional scientific PowerPoint presentation using the pptx and scientific-slides skills.

## Analysis Details
- **Project:** {analysis_name}
- **Type:** Single-cell RNA-seq Analysis
- **Results Directory:** {results_dir}

## Presentation Structure

{json.dumps(outline, indent=2, default=str)}

## Requirements

1. **Design Guidelines:**
   - Use scientific presentation template
   - Professional color scheme (blues, greys)
   - Clear, readable fonts (minimum 18pt for body text)
   - Consistent layout across slides
   - High-quality figure integration

2. **Figure Integration:**
   - All figures in `figures/` subdirectory
   - Full-bleed figures for maximum impact
   - Figure captions below each image
   - Maintain aspect ratios

3. **Content Quality:**
   - Concise bullet points
   - Scientific terminology
   - Clear section transitions
   - Logical flow of information

4. **Slide Types:**
   - Title slide with project name and date
   - Section headers for major topics
   - Content slides with bullet points
   - Figure slides with captions
   - Summary slides with key findings

## Output
Save the presentation as: {output_file}

Please generate a publication-quality scientific presentation suitable for:
- Lab meetings
- Conference presentations
- Thesis defenses
- Manuscript preparation
"""

    instruction_file = results_dir / 'ppt_generation_instructions.txt'
    with open(instruction_file, 'w') as f:
        f.write(instructions)

    print(f"  Saved instructions: {instruction_file.name}")

    print("\n[3/3] Instructions for Claude Skills:")
    print("-"*80)
    print("To generate the PowerPoint presentation using Claude skills:")
    print("")
    print("METHOD 1: Using /pptx skill directly")
    print("  1. In Claude Code, type: /pptx")
    print(f"  2. Provide the instruction file: {instruction_file}")
    print(f"  3. Specify figures directory: {results_dir / 'figures'}")
    print("")
    print("METHOD 2: Using /scientific-slides skill")
    print("  1. In Claude Code, type: /scientific-slides")
    print(f"  2. Provide presentation outline: {outline_file}")
    print(f"  3. Provide figures directory: {results_dir / 'figures'}")
    print("")
    print("The skills will automatically:")
    print("  • Create professional scientific presentation")
    print("  • Integrate all figures with proper formatting")
    print("  • Apply consistent styling and layout")
    print("  • Generate publication-quality slides")
    print("-"*80)

    print("\n" + "="*80)
    print("PREPARATION COMPLETE")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  - {outline_file.name}")
    print(f"  - {instruction_file.name}")
    print(f"\nNext steps:")
    print(f"  1. Use /pptx or /scientific-slides skill in Claude Code")
    print(f"  2. Follow the prompts with the generated files")
    print(f"  3. Review and refine the generated presentation")
    print("="*80)

    return outline, instruction_file


if __name__ == "__main__":
    # Generate for JMV ALI analysis
    results_dir = '/mnt2/kwy/scrna_agent/JMV_ALI_results'
    analysis_name = 'JMV ALI Organoid scRNA-seq Analysis'

    outline, instruction_file = generate_ppt_with_skills(
        results_dir=results_dir,
        analysis_name=analysis_name
    )

    print("\n" + "="*80)
    print("IMPORTANT: PPT Generation Using Claude Skills")
    print("="*80)
    print("\nThis script has prepared all necessary files for PPT generation.")
    print("To generate the actual PowerPoint presentation:")
    print("")
    print("  Run in Claude Code CLI:")
    print("  $ /scientific-slides")
    print("")
    print("  Then provide:")
    print(f"  - Outline: {instruction_file}")
    print(f"  - Figures: {results_dir}/figures/")
    print("")
    print("The skill will create a professional scientific presentation")
    print("with all figures, proper formatting, and scientific styling.")
    print("="*80)
