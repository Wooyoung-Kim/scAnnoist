#!/usr/bin/env python3
"""
Multi-Agent Cell Type Annotation System
Uses multiple specialized agents to discuss and determine cell types
"""

import scanpy as sc
import pandas as pd
import numpy as np
import json
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt

def prepare_cluster_information(adata, cluster_key='louvain_0.8', results_dir=None):
    """
    Extract comprehensive information for each cluster for agent discussion

    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    cluster_key : str
        Column name for clusters
    results_dir : Path
        Directory containing marker gene files

    Returns:
    --------
    cluster_info : dict
        Comprehensive information for each cluster
    """

    print("="*80)
    print("PREPARING CLUSTER INFORMATION FOR MULTI-AGENT DISCUSSION")
    print("="*80)

    cluster_info = {}
    results_dir = Path(results_dir) if results_dir else Path('.')

    for cluster_id in sorted(adata.obs[cluster_key].unique()):
        print(f"\nProcessing Cluster {cluster_id}...")

        cluster_mask = adata.obs[cluster_key] == cluster_id
        n_cells = cluster_mask.sum()

        # Get cluster statistics
        cluster_stats = {
            'cluster_id': str(cluster_id),
            'n_cells': int(n_cells),
            'percentage': float((n_cells / adata.n_obs) * 100),
        }

        # QC metrics for this cluster
        qc_metrics = {
            'mean_genes': float(adata.obs.loc[cluster_mask, 'n_genes_by_counts'].mean()),
            'median_genes': float(adata.obs.loc[cluster_mask, 'n_genes_by_counts'].median()),
            'mean_counts': float(adata.obs.loc[cluster_mask, 'total_counts'].mean()),
            'mean_pct_mt': float(adata.obs.loc[cluster_mask, 'pct_counts_mt'].mean()),
        }

        # Load marker genes if available
        marker_file = results_dir / f'markers_cluster_{cluster_id}.csv'
        top_markers = []

        if marker_file.exists():
            markers_df = pd.read_csv(marker_file)
            # Get top 20 markers
            top_20 = markers_df.head(20)
            for idx, row in top_20.iterrows():
                top_markers.append({
                    'gene': row['names'],
                    'log2fc': float(row['logfoldchanges']),
                    'pval_adj': float(row['pvals_adj']),
                    'score': float(row['scores']) if 'scores' in row else None
                })

        # Get condition distribution
        if 'condition' in adata.obs.columns:
            condition_dist = adata.obs.loc[cluster_mask, 'condition'].value_counts().to_dict()
            condition_dist = {str(k): int(v) for k, v in condition_dist.items()}
        else:
            condition_dist = {}

        # Compile all information
        cluster_info[str(cluster_id)] = {
            'statistics': cluster_stats,
            'qc_metrics': qc_metrics,
            'top_markers': top_markers,
            'condition_distribution': condition_dist,
        }

        print(f"  Cells: {n_cells} ({cluster_stats['percentage']:.1f}%)")
        print(f"  Top markers: {', '.join([m['gene'] for m in top_markers[:5]])}")

    # Save cluster information
    output_file = results_dir / 'cluster_information_for_agents.json'
    with open(output_file, 'w') as f:
        json.dump(cluster_info, f, indent=2)

    print(f"\n✓ Saved cluster information: {output_file.name}")
    print("="*80)

    return cluster_info


def create_agent_discussion_prompts(cluster_info, tissue_type='airway epithelium'):
    """
    Create discussion prompts for different specialized agents

    Parameters:
    -----------
    cluster_info : dict
        Cluster information
    tissue_type : str
        Tissue type for context

    Returns:
    --------
    prompts : dict
        Prompts for each agent type
    """

    prompts = {}

    # Format cluster data for agents
    cluster_summaries = []
    for cluster_id, info in cluster_info.items():
        stats = info['statistics']
        qc = info['qc_metrics']
        markers = info['top_markers'][:10]  # Top 10 markers

        summary = f"""
**Cluster {cluster_id}**
- Cells: {stats['n_cells']:,} ({stats['percentage']:.1f}%)
- Mean genes: {qc['mean_genes']:.0f}
- Mean UMIs: {qc['mean_counts']:.0f}
- % MT: {qc['mean_pct_mt']:.2f}%

Top Marker Genes:
"""
        for m in markers:
            summary += f"  - {m['gene']} (log2FC={m['log2fc']:.2f}, padj={m['pval_adj']:.2e})\n"

        if info['condition_distribution']:
            summary += "\nCondition Distribution:\n"
            for cond, count in info['condition_distribution'].items():
                summary += f"  - {cond}: {count} cells\n"

        cluster_summaries.append(summary)

    all_clusters_text = "\n".join(cluster_summaries)

    # Agent 1: Cell Type Expert
    prompts['cell_type_expert'] = f"""You are a cell biology expert specializing in {tissue_type}.

Your task is to annotate cell clusters based on marker gene expression.

## Dataset Information
- Tissue: {tissue_type}
- Number of clusters: {len(cluster_info)}

## Cluster Data
{all_clusters_text}

## Your Task
For each cluster, provide:
1. **Proposed Cell Type:** Most likely cell type based on markers
2. **Confidence Level:** High/Medium/Low
3. **Key Markers:** Which markers support this annotation
4. **Rationale:** Detailed explanation of your reasoning
5. **Alternative Possibilities:** Other potential cell types if uncertain

Consider:
- Known marker genes for {tissue_type} cell types
- Marker gene expression levels (log2FC)
- Statistical significance (adjusted p-values)
- Cell quality metrics
- Biological plausibility

Provide your analysis in structured JSON format:
{{
  "cluster_id": {{
    "proposed_cell_type": "...",
    "confidence": "High/Medium/Low",
    "key_markers": ["gene1", "gene2", ...],
    "rationale": "detailed explanation",
    "alternatives": ["alternative1", "alternative2"]
  }},
  ...
}}
"""

    # Agent 2: Computational Biologist
    prompts['computational_biologist'] = f"""You are a computational biologist expert in single-cell RNA-seq analysis.

Your task is to evaluate cluster annotations from a quantitative perspective.

## Dataset Information
- Tissue: {tissue_type}
- Number of clusters: {len(cluster_info)}

## Cluster Data
{all_clusters_text}

## Your Task
For each cluster, provide:
1. **Quality Assessment:** Evaluate based on QC metrics
2. **Statistical Support:** Strength of marker gene evidence
3. **Uniqueness:** How distinct is this cluster's signature
4. **Potential Issues:** Doublets, low quality, stressed cells
5. **Recommendations:** Merging clusters, removing outliers, etc.

Consider:
- Marker gene statistics (log2FC, p-values)
- Number of cells in cluster
- QC metrics (genes, UMIs, % MT)
- Cluster size and proportion
- Condition-specific effects

Provide your analysis in structured JSON format:
{{
  "cluster_id": {{
    "quality_assessment": "High/Medium/Low",
    "statistical_support": "Strong/Moderate/Weak",
    "uniqueness_score": "High/Medium/Low",
    "potential_issues": ["issue1", "issue2"],
    "recommendations": "...",
    "confidence_in_annotation": "High/Medium/Low"
  }},
  ...
}}
"""

    # Agent 3: Literature Curator
    prompts['literature_curator'] = f"""You are a scientific literature curator specializing in {tissue_type}.

Your task is to validate proposed annotations against published literature.

## Dataset Information
- Tissue: {tissue_type}
- Number of clusters: {len(cluster_info)}

## Cluster Data
{all_clusters_text}

## Your Task
For each cluster, provide:
1. **Literature Support:** Known markers from published studies
2. **Canonical Markers:** Well-established markers for cell types
3. **Novel Findings:** Any unexpected or interesting patterns
4. **Functional Context:** Biological role of the cell type
5. **Validation Recommendations:** Additional markers to check

Consider:
- Published single-cell atlases of {tissue_type}
- Classic immunohistochemistry markers
- Transcription factors defining cell lineages
- Functional markers and receptors
- Recent discoveries in the field

Provide your analysis in structured JSON format:
{{
  "cluster_id": {{
    "literature_support": "Strong/Moderate/Weak/None",
    "canonical_markers": ["marker1", "marker2"],
    "novel_findings": "...",
    "functional_context": "...",
    "validation_markers": ["marker1", "marker2"],
    "references": ["citation1", "citation2"]
  }},
  ...
}}
"""

    # Agent 4: Critical Evaluator
    prompts['critical_evaluator'] = f"""You are a critical evaluator tasked with synthesizing annotations from multiple experts.

Your task is to review annotations from other agents and make final recommendations.

## Dataset Information
- Tissue: {tissue_type}
- Number of clusters: {len(cluster_info)}

## Cluster Data
{all_clusters_text}

You will receive annotations from:
1. Cell Type Expert
2. Computational Biologist
3. Literature Curator

## Your Task
For each cluster, provide:
1. **Agreement Analysis:** Where do experts agree/disagree
2. **Confidence Assessment:** Overall confidence in annotation
3. **Final Recommendation:** Proposed final cell type
4. **Ambiguity Notes:** Any remaining uncertainties
5. **Priority for Validation:** Which clusters need experimental validation

Consider:
- Consistency across expert opinions
- Strength of evidence from each perspective
- Biological plausibility
- Statistical support
- Literature validation

Provide your synthesis in structured JSON format:
{{
  "cluster_id": {{
    "agreement_level": "Full/Partial/Low",
    "final_cell_type": "...",
    "overall_confidence": "High/Medium/Low",
    "consensus_rationale": "...",
    "remaining_ambiguities": "...",
    "validation_priority": "High/Medium/Low",
    "next_steps": "..."
  }},
  ...
}}
"""

    return prompts


def save_agent_prompts(prompts, output_dir):
    """Save agent prompts to files for Task tool usage"""

    output_dir = Path(output_dir)
    prompts_dir = output_dir / 'agent_prompts'
    prompts_dir.mkdir(exist_ok=True)

    prompt_files = {}

    for agent_name, prompt in prompts.items():
        filename = prompts_dir / f'{agent_name}_prompt.txt'
        with open(filename, 'w') as f:
            f.write(prompt)
        prompt_files[agent_name] = filename
        print(f"  Saved: {filename.name}")

    return prompt_files


def generate_discussion_workflow(cluster_info_file, output_dir, tissue_type='airway epithelium'):
    """
    Generate complete multi-agent annotation workflow

    Parameters:
    -----------
    cluster_info_file : str or Path
        Path to cluster information JSON
    output_dir : str or Path
        Output directory
    tissue_type : str
        Tissue type
    """

    output_dir = Path(output_dir)

    print("="*80)
    print("MULTI-AGENT ANNOTATION WORKFLOW GENERATION")
    print("="*80)

    # Load cluster information
    with open(cluster_info_file, 'r') as f:
        cluster_info = json.load(f)

    print(f"\nLoaded information for {len(cluster_info)} clusters")

    # Create agent prompts
    print("\n[1/3] Creating agent discussion prompts...")
    prompts = create_agent_discussion_prompts(cluster_info, tissue_type)
    print(f"  Generated prompts for {len(prompts)} specialized agents")

    # Save prompts
    print("\n[2/3] Saving agent prompts...")
    prompt_files = save_agent_prompts(prompts, output_dir)

    # Create workflow instructions
    print("\n[3/3] Creating workflow instructions...")

    workflow_instructions = f"""# Multi-Agent Cell Type Annotation Workflow

## Overview
This workflow uses 4 specialized agents to discuss and determine cell type annotations:

1. **Cell Type Expert** - Biological knowledge and marker interpretation
2. **Computational Biologist** - Statistical and quantitative analysis
3. **Literature Curator** - Published evidence and validation
4. **Critical Evaluator** - Synthesis and final recommendations

## Execution Steps

### Step 1: Run Cell Type Expert Agent

In Claude Code, use the Task tool:

```
I need you to run a general-purpose agent to analyze cell clusters and propose cell type annotations.

Task: Read the prompt from {prompt_files['cell_type_expert']} and provide detailed cell type annotations for each cluster based on marker genes.

Please analyze all clusters and provide your expert opinion in the requested JSON format.
```

Save the output as: `agent_outputs/cell_type_expert_output.json`

### Step 2: Run Computational Biologist Agent

```
I need you to run a general-purpose agent to evaluate cluster quality and statistical support.

Task: Read the prompt from {prompt_files['computational_biologist']} and provide quantitative assessment of each cluster's annotation quality.

Please evaluate all clusters from a computational perspective and provide analysis in JSON format.
```

Save the output as: `agent_outputs/computational_biologist_output.json`

### Step 3: Run Literature Curator Agent

```
I need you to run a general-purpose agent to validate annotations against published literature.

Task: Read the prompt from {prompt_files['literature_curator']} and provide literature-based validation for proposed annotations.

Please review all clusters and provide literature support in JSON format.
```

Save the output as: `agent_outputs/literature_curator_output.json`

### Step 4: Run Critical Evaluator Agent

```
I need you to run a general-purpose agent to synthesize all expert opinions.

Task:
1. Read outputs from the three previous agents:
   - {prompt_files['cell_type_expert'].parent}/agent_outputs/cell_type_expert_output.json
   - {prompt_files['cell_type_expert'].parent}/agent_outputs/computational_biologist_output.json
   - {prompt_files['cell_type_expert'].parent}/agent_outputs/literature_curator_output.json

2. Read the prompt from {prompt_files['critical_evaluator']}

3. Synthesize all opinions and provide final recommendations

Please analyze the agreements and disagreements, then provide final cell type annotations with rationale.
```

Save the output as: `agent_outputs/critical_evaluator_output.json`

### Step 5: Generate Final Report

After all agents have completed their analysis, run:

```bash
python compile_agent_annotations.py
```

This will:
- Combine all agent outputs
- Generate consensus annotations
- Create detailed annotation report
- Update the h5ad file with annotations
- Regenerate visualizations and PowerPoint

## Output Structure

```
{output_dir}/
├── agent_prompts/
│   ├── cell_type_expert_prompt.txt
│   ├── computational_biologist_prompt.txt
│   ├── literature_curator_prompt.txt
│   └── critical_evaluator_prompt.txt
├── agent_outputs/
│   ├── cell_type_expert_output.json
│   ├── computational_biologist_output.json
│   ├── literature_curator_output.json
│   └── critical_evaluator_output.json
└── MULTI_AGENT_ANNOTATION_REPORT.md
```

## Benefits of Multi-Agent Approach

1. **Multiple Perspectives:** Different expertise areas provide comprehensive analysis
2. **Cross-Validation:** Agreement between agents increases confidence
3. **Transparency:** Each agent's reasoning is documented
4. **Flexibility:** Can add or modify agents for specific needs
5. **Quality Control:** Critical evaluation step catches inconsistencies
6. **Reproducibility:** All decisions and rationale are recorded

## Timeline

- Each agent: ~5-10 minutes
- Total workflow: ~30-45 minutes
- Final compilation: ~5 minutes

## Notes

- Agents can be run in parallel (Steps 1-3)
- Step 4 must wait for Steps 1-3 to complete
- Review each agent's output before proceeding
- Modify prompts if needed for your specific tissue type
"""

    workflow_file = output_dir / 'MULTI_AGENT_WORKFLOW.md'
    with open(workflow_file, 'w') as f:
        f.write(workflow_instructions)

    print(f"  Saved: {workflow_file.name}")

    # Create agent outputs directory
    (output_dir / 'agent_outputs').mkdir(exist_ok=True)

    print("\n" + "="*80)
    print("WORKFLOW GENERATION COMPLETE")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  - 4 agent prompt files in agent_prompts/")
    print(f"  - {workflow_file.name}")
    print(f"\nNext steps:")
    print(f"  1. Review MULTI_AGENT_WORKFLOW.md")
    print(f"  2. Execute agents sequentially using Task tool")
    print(f"  3. Compile results with compile_agent_annotations.py")
    print("="*80)

    return prompt_files, workflow_file


if __name__ == "__main__":
    # Example: Prepare for JMV ALI analysis
    adata_file = '/mnt2/kwy/scrna_agent/JMV_ALI_results/JMV_ALI_harmony_integrated.h5ad'
    results_dir = '/mnt2/kwy/scrna_agent/JMV_ALI_results'

    print("Loading data...")
    adata = sc.read_h5ad(adata_file)

    # Prepare cluster information
    cluster_info = prepare_cluster_information(
        adata=adata,
        cluster_key='louvain_0.8',
        results_dir=results_dir
    )

    # Generate workflow
    cluster_info_file = Path(results_dir) / 'cluster_information_for_agents.json'

    prompt_files, workflow_file = generate_discussion_workflow(
        cluster_info_file=cluster_info_file,
        output_dir=results_dir,
        tissue_type='Human airway epithelium (ALI organoid)'
    )

    print("\n" + "="*80)
    print("READY FOR MULTI-AGENT ANNOTATION")
    print("="*80)
    print(f"\nFollow instructions in: {workflow_file}")
    print("\nThe multi-agent system will provide:")
    print("  • Comprehensive cell type annotations")
    print("  • Multiple expert perspectives")
    print("  • Detailed rationale for each decision")
    print("  • Confidence assessments")
    print("  • Validation recommendations")
    print("="*80)
