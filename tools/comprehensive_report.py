"""Comprehensive Report Generator

Generates detailed analysis reports including:
- QC summary with figures
- Annotation summary with expert interpretations
- Integration summary with batch mixing analysis
- Condition-wise cell type distribution
- All generated figures organized and documented
"""

from typing import Annotated, Optional, Dict, List, Any
from pathlib import Path
import json
from datetime import datetime
import pandas as pd
from langchain_core.tools import tool


# Global state for report data
_report_data = {
    "qc_metrics": {},
    "annotation_results": {},
    "integration_results": {},
    "condition_analysis": {},
    "figures": {},
    "metadata": {}
}


@tool
def initialize_report(
    project_name: Annotated[str, "Name of the project/analysis"],
    output_dir: Annotated[str, "Output directory for reports"] = "./reports",
) -> Annotated[str, "Confirmation message"]:
    """
    Initialize a new comprehensive report.

    Sets up report structure and metadata.
    """
    global _report_data

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    _report_data = {
        "project_name": project_name,
        "output_dir": str(output_path),
        "created_at": datetime.now().isoformat(),
        "qc_metrics": {},
        "annotation_results": {},
        "integration_results": {},
        "condition_analysis": {},
        "figures": {},
        "metadata": {}
    }

    return f"Report initialized for project '{project_name}' in {output_path}"


@tool
def add_qc_metrics(
    n_cells_before: Annotated[int, "Number of cells before QC"],
    n_cells_after: Annotated[int, "Number of cells after QC"],
    median_genes_before: Annotated[float, "Median genes before QC"],
    median_genes_after: Annotated[float, "Median genes after QC"],
    median_counts_before: Annotated[float, "Median counts before QC"],
    median_counts_after: Annotated[float, "Median counts after QC"],
    median_mt_before: Annotated[float, "Median MT% before QC"],
    median_mt_after: Annotated[float, "Median MT% after QC"],
    qc_thresholds: Annotated[Optional[Dict[str, float]], "QC thresholds used"] = None,
) -> Annotated[str, "Confirmation message"]:
    """
    Add QC metrics to the report.
    """
    global _report_data

    _report_data["qc_metrics"] = {
        "n_cells_before": n_cells_before,
        "n_cells_after": n_cells_after,
        "n_cells_removed": n_cells_before - n_cells_after,
        "pct_cells_retained": (n_cells_after / n_cells_before * 100) if n_cells_before > 0 else 0,
        "median_genes_before": median_genes_before,
        "median_genes_after": median_genes_after,
        "median_counts_before": median_counts_before,
        "median_counts_after": median_counts_after,
        "median_mt_before": median_mt_before,
        "median_mt_after": median_mt_after,
        "qc_thresholds": qc_thresholds or {}
    }

    return f"QC metrics added: {n_cells_before:,} → {n_cells_after:,} cells ({_report_data['qc_metrics']['pct_cells_retained']:.1f}% retained)"


@tool
def add_annotation_results(
    n_cell_types: Annotated[int, "Number of distinct cell types identified"],
    annotation_method: Annotated[str, "Primary annotation method used"],
    consensus_rate: Annotated[Optional[float], "Agreement rate between methods (if multiple)"] = None,
    cell_type_counts: Annotated[Optional[Dict[str, int]], "Cell counts per cell type"] = None,
    expert_interpretations: Annotated[Optional[Dict[str, str]], "Expert interpretations for ambiguous clusters"] = None,
) -> Annotated[str, "Confirmation message"]:
    """
    Add annotation results to the report.
    """
    global _report_data

    _report_data["annotation_results"] = {
        "n_cell_types": n_cell_types,
        "annotation_method": annotation_method,
        "consensus_rate": consensus_rate,
        "cell_type_counts": cell_type_counts or {},
        "expert_interpretations": expert_interpretations or {}
    }

    return f"Annotation results added: {n_cell_types} cell types identified using {annotation_method}"


@tool
def add_integration_results(
    n_batches: Annotated[int, "Number of batches integrated"],
    integration_method: Annotated[str, "Integration method used (scVI, Harmony, etc.)"],
    batch_distribution: Annotated[Optional[Dict[str, int]], "Cell counts per batch"] = None,
    integration_metrics: Annotated[Optional[Dict[str, float]], "Integration quality metrics (SCIB)"] = None,
) -> Annotated[str, "Confirmation message"]:
    """
    Add integration results to the report.
    """
    global _report_data

    _report_data["integration_results"] = {
        "n_batches": n_batches,
        "integration_method": integration_method,
        "batch_distribution": batch_distribution or {},
        "integration_metrics": integration_metrics or {}
    }

    return f"Integration results added: {n_batches} batches integrated using {integration_method}"


@tool
def add_condition_analysis(
    condition_key: Annotated[str, "Condition variable name (e.g., 'age', 'treatment')"],
    conditions: Annotated[List[str], "List of condition values"],
    cell_type_distribution: Annotated[Dict[str, Dict[str, int]], "Nested dict: {condition: {celltype: count}}"],
) -> Annotated[str, "Confirmation message"]:
    """
    Add condition-wise cell type distribution analysis to the report.
    """
    global _report_data

    _report_data["condition_analysis"][condition_key] = {
        "conditions": conditions,
        "cell_type_distribution": cell_type_distribution
    }

    return f"Condition analysis added for '{condition_key}': {len(conditions)} conditions"


@tool
def add_figure(
    figure_name: Annotated[str, "Name/identifier for the figure"],
    figure_path: Annotated[str, "Path to the figure file"],
    figure_caption: Annotated[str, "Caption describing the figure"],
    figure_category: Annotated[str, "Category: qc, annotation, integration, or condition"] = "other",
) -> Annotated[str, "Confirmation message"]:
    """
    Add a figure to the report.
    """
    global _report_data

    if "figures" not in _report_data:
        _report_data["figures"] = {}

    if figure_category not in _report_data["figures"]:
        _report_data["figures"][figure_category] = []

    _report_data["figures"][figure_category].append({
        "name": figure_name,
        "path": figure_path,
        "caption": figure_caption
    })

    return f"Figure added: {figure_name} ({figure_category})"


@tool
def generate_markdown_report(
    output_file: Annotated[str, "Output file path for markdown report"] = "./reports/analysis_report.md",
) -> Annotated[str, "Path to generated markdown report"]:
    """
    Generate a comprehensive markdown report from all collected data.

    Creates a detailed, publication-quality report with:
    - Executive summary
    - QC analysis
    - Annotation results
    - Integration results
    - Condition-wise analysis
    - All figures with captions
    - Methods section
    """
    global _report_data

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build markdown content
    md_lines = []

    # Title and metadata
    md_lines.append(f"# {_report_data.get('project_name', 'scRNA-seq Analysis')} - Comprehensive Report")
    md_lines.append("")
    md_lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md_lines.append("")
    md_lines.append("---")
    md_lines.append("")

    # Table of Contents
    md_lines.append("## Table of Contents")
    md_lines.append("")
    md_lines.append("1. [Executive Summary](#executive-summary)")
    md_lines.append("2. [Quality Control Analysis](#quality-control-analysis)")
    md_lines.append("3. [Cell Type Annotation](#cell-type-annotation)")
    md_lines.append("4. [Batch Integration](#batch-integration)")
    md_lines.append("5. [Condition-wise Analysis](#condition-wise-analysis)")
    md_lines.append("6. [Figures](#figures)")
    md_lines.append("7. [Methods](#methods)")
    md_lines.append("")
    md_lines.append("---")
    md_lines.append("")

    # Executive Summary
    md_lines.append("## Executive Summary")
    md_lines.append("")

    if "qc_metrics" in _report_data and _report_data["qc_metrics"]:
        qc = _report_data["qc_metrics"]
        md_lines.append(f"- **Total cells analyzed:** {qc.get('n_cells_before', 'N/A'):,}")
        md_lines.append(f"- **Cells passing QC:** {qc.get('n_cells_after', 'N/A'):,} ({qc.get('pct_cells_retained', 0):.1f}%)")

    if "annotation_results" in _report_data and _report_data["annotation_results"]:
        annot = _report_data["annotation_results"]
        md_lines.append(f"- **Cell types identified:** {annot.get('n_cell_types', 'N/A')}")
        md_lines.append(f"- **Annotation method:** {annot.get('annotation_method', 'N/A')}")

    if "integration_results" in _report_data and _report_data["integration_results"]:
        integ = _report_data["integration_results"]
        md_lines.append(f"- **Batches integrated:** {integ.get('n_batches', 'N/A')}")
        md_lines.append(f"- **Integration method:** {integ.get('integration_method', 'N/A')}")

    md_lines.append("")
    md_lines.append("---")
    md_lines.append("")

    # QC Analysis
    md_lines.append("## Quality Control Analysis")
    md_lines.append("")

    if "qc_metrics" in _report_data and _report_data["qc_metrics"]:
        qc = _report_data["qc_metrics"]

        md_lines.append("### QC Metrics Summary")
        md_lines.append("")
        md_lines.append("| Metric | Before QC | After QC | Change |")
        md_lines.append("|--------|-----------|----------|--------|")
        md_lines.append(f"| **Total Cells** | {qc.get('n_cells_before', 0):,} | {qc.get('n_cells_after', 0):,} | -{qc.get('n_cells_removed', 0):,} ({100-qc.get('pct_cells_retained', 100):.1f}% removed) |")
        md_lines.append(f"| **Median Genes/Cell** | {qc.get('median_genes_before', 0):.0f} | {qc.get('median_genes_after', 0):.0f} | {qc.get('median_genes_after', 0) - qc.get('median_genes_before', 0):.0f} |")
        md_lines.append(f"| **Median UMI Counts/Cell** | {qc.get('median_counts_before', 0):.0f} | {qc.get('median_counts_after', 0):.0f} | {qc.get('median_counts_after', 0) - qc.get('median_counts_before', 0):.0f} |")
        md_lines.append(f"| **Median MT%** | {qc.get('median_mt_before', 0):.2f}% | {qc.get('median_mt_after', 0):.2f}% | {qc.get('median_mt_after', 0) - qc.get('median_mt_before', 0):.2f}% |")
        md_lines.append("")

        if qc.get("qc_thresholds"):
            md_lines.append("### QC Thresholds Applied")
            md_lines.append("")
            for key, value in qc["qc_thresholds"].items():
                md_lines.append(f"- **{key}:** {value}")
            md_lines.append("")

        md_lines.append("### Interpretation")
        md_lines.append("")
        pct_retained = qc.get('pct_cells_retained', 100)
        if pct_retained >= 80:
            md_lines.append(f"✅ **Excellent data quality**: {pct_retained:.1f}% of cells retained after QC filtering.")
        elif pct_retained >= 60:
            md_lines.append(f"✓ **Good data quality**: {pct_retained:.1f}% of cells retained after QC filtering.")
        else:
            md_lines.append(f"⚠️  **Moderate data quality**: {pct_retained:.1f}% of cells retained. Consider reviewing QC thresholds.")

        md_lines.append("")

    md_lines.append("---")
    md_lines.append("")

    # Annotation
    md_lines.append("## Cell Type Annotation")
    md_lines.append("")

    if "annotation_results" in _report_data and _report_data["annotation_results"]:
        annot = _report_data["annotation_results"]

        md_lines.append(f"### Method: {annot.get('annotation_method', 'N/A')}")
        md_lines.append("")

        if annot.get("consensus_rate") is not None:
            md_lines.append(f"**Consensus rate:** {annot['consensus_rate']:.1f}% (agreement between annotation methods)")
            md_lines.append("")

        if annot.get("cell_type_counts"):
            md_lines.append("### Cell Type Distribution")
            md_lines.append("")
            md_lines.append("| Cell Type | Count | Percentage |")
            md_lines.append("|-----------|-------|------------|")

            total_cells = sum(annot["cell_type_counts"].values())
            for celltype, count in sorted(annot["cell_type_counts"].items(), key=lambda x: -x[1]):
                pct = (count / total_cells * 100) if total_cells > 0 else 0
                md_lines.append(f"| {celltype} | {count:,} | {pct:.2f}% |")

            md_lines.append("")

        if annot.get("expert_interpretations"):
            md_lines.append("### Expert Interpretations")
            md_lines.append("")
            md_lines.append("Domain experts provided the following interpretations for ambiguous or complex cell populations:")
            md_lines.append("")

            for cluster, interpretation in annot["expert_interpretations"].items():
                md_lines.append(f"**{cluster}:**")
                md_lines.append(f"> {interpretation}")
                md_lines.append("")

    md_lines.append("---")
    md_lines.append("")

    # Integration
    md_lines.append("## Batch Integration")
    md_lines.append("")

    if "integration_results" in _report_data and _report_data["integration_results"]:
        integ = _report_data["integration_results"]

        md_lines.append(f"### Method: {integ.get('integration_method', 'N/A')}")
        md_lines.append("")

        if integ.get("batch_distribution"):
            md_lines.append("### Batch Distribution")
            md_lines.append("")
            md_lines.append("| Batch | Cell Count |")
            md_lines.append("|-------|------------|")

            for batch, count in sorted(integ["batch_distribution"].items()):
                md_lines.append(f"| {batch} | {count:,} |")

            md_lines.append("")

        if integ.get("integration_metrics"):
            md_lines.append("### Integration Quality Metrics (SCIB)")
            md_lines.append("")
            md_lines.append("| Metric | Score | Interpretation |")
            md_lines.append("|--------|-------|----------------|")

            for metric, score in integ["integration_metrics"].items():
                if score >= 0.8:
                    interp = "Excellent"
                elif score >= 0.6:
                    interp = "Good"
                else:
                    interp = "Needs improvement"

                md_lines.append(f"| {metric} | {score:.3f} | {interp} |")

            md_lines.append("")

    md_lines.append("---")
    md_lines.append("")

    # Condition-wise analysis
    md_lines.append("## Condition-wise Analysis")
    md_lines.append("")

    if "condition_analysis" in _report_data and _report_data["condition_analysis"]:
        for condition_key, data in _report_data["condition_analysis"].items():
            md_lines.append(f"### Condition: {condition_key}")
            md_lines.append("")

            if data.get("cell_type_distribution"):
                # Create a table
                md_lines.append("| Cell Type | " + " | ".join(data["conditions"]) + " |")
                md_lines.append("|-----------|" + "|".join(["---"] * len(data["conditions"])) + "|")

                # Get all cell types
                all_celltypes = set()
                for cond_dict in data["cell_type_distribution"].values():
                    all_celltypes.update(cond_dict.keys())

                for celltype in sorted(all_celltypes):
                    row = [celltype]
                    for cond in data["conditions"]:
                        count = data["cell_type_distribution"].get(cond, {}).get(celltype, 0)
                        row.append(f"{count:,}")
                    md_lines.append("| " + " | ".join(row) + " |")

                md_lines.append("")

    md_lines.append("---")
    md_lines.append("")

    # Figures
    md_lines.append("## Figures")
    md_lines.append("")

    if "figures" in _report_data and _report_data["figures"]:
        for category, figures in _report_data["figures"].items():
            if figures:
                md_lines.append(f"### {category.title()} Figures")
                md_lines.append("")

                for idx, fig in enumerate(figures, 1):
                    md_lines.append(f"#### Figure {category[0].upper()}{idx}: {fig['name']}")
                    md_lines.append("")
                    md_lines.append(f"![{fig['name']}]({fig['path']})")
                    md_lines.append("")
                    md_lines.append(f"**Caption:** {fig['caption']}")
                    md_lines.append("")

        # Create figure index
        md_lines.append("### Figure Index")
        md_lines.append("")
        md_lines.append("| Category | Figure | Path |")
        md_lines.append("|----------|--------|------|")

        for category, figures in _report_data["figures"].items():
            for fig in figures:
                md_lines.append(f"| {category.title()} | {fig['name']} | `{fig['path']}` |")

        md_lines.append("")

    md_lines.append("---")
    md_lines.append("")

    # Methods
    md_lines.append("## Methods")
    md_lines.append("")

    md_lines.append("### Quality Control")
    md_lines.append("")
    md_lines.append("Cells were filtered based on:")
    if "qc_metrics" in _report_data and _report_data["qc_metrics"].get("qc_thresholds"):
        for key, value in _report_data["qc_metrics"]["qc_thresholds"].items():
            md_lines.append(f"- {key}: {value}")
    else:
        md_lines.append("- Standard Scanpy QC metrics (n_genes, n_counts, pct_counts_mt)")
    md_lines.append("")

    if "annotation_results" in _report_data and _report_data["annotation_results"].get("annotation_method"):
        md_lines.append("### Cell Type Annotation")
        md_lines.append("")
        method = _report_data["annotation_results"]["annotation_method"]
        md_lines.append(f"Cell types were annotated using {method}.")
        md_lines.append("")

    if "integration_results" in _report_data and _report_data["integration_results"].get("integration_method"):
        md_lines.append("### Batch Integration")
        md_lines.append("")
        method = _report_data["integration_results"]["integration_method"]
        md_lines.append(f"Batch effects were corrected using {method}.")
        md_lines.append("")

    md_lines.append("---")
    md_lines.append("")

    # Footer
    md_lines.append("## Report Information")
    md_lines.append("")
    md_lines.append(f"- **Generated by:** scRNA-agent comprehensive report generator")
    md_lines.append(f"- **Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md_lines.append(f"- **Project:** {_report_data.get('project_name', 'N/A')}")
    md_lines.append("")

    # Write to file
    with open(output_path, 'w') as f:
        f.write('\n'.join(md_lines))

    # Also save as JSON
    json_path = output_path.with_suffix('.json')
    with open(json_path, 'w') as f:
        json.dump(_report_data, f, indent=2)

    return f"""
Comprehensive report generated:
- Markdown: {output_path}
- JSON data: {json_path}

Report includes:
- QC analysis: {len(_report_data.get('qc_metrics', {}))} metrics
- Annotation: {_report_data.get('annotation_results', {}).get('n_cell_types', 0)} cell types
- Integration: {_report_data.get('integration_results', {}).get('n_batches', 0)} batches
- Figures: {sum(len(figs) for figs in _report_data.get('figures', {}).values())} total
"""


@tool
def get_report_data() -> Annotated[Dict[str, Any], "Current report data"]:
    """
    Get the current report data structure.

    Returns all collected data for inspection or further processing.
    """
    return _report_data.copy()
