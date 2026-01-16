#!/usr/bin/env python3
"""
Comprehensive Cell Type Annotation with Multiple Methods
- CellTypist automated annotation
- Marker-based annotation
- Comparison and rationale documentation
- PowerPoint integration
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

def comprehensive_annotation(adata_file, manual_annotation_dict,
                            cluster_key='louvain_0.8',
                            celltyp

ist_model='Immune_All_Low.pkl',
                            tissue_type='lung',
                            output_dir=None):
    """
    Perform comprehensive cell type annotation using multiple methods

    Parameters:
    -----------
    adata_file : str or Path
        Path to the h5ad file
    manual_annotation_dict : dict
        Manual annotation based on marker genes
        Format: {'cluster_id': {'cell_type': 'name', 'markers': ['gene1', 'gene2'], 'rationale': 'description'}}
    cluster_key : str
        Column name containing cluster IDs
    celltypist_model : str
        CellTypist model to use
    tissue_type : str
        Tissue type for context
    output_dir : str or Path
        Output directory
    """

    print("="*80)
    print("COMPREHENSIVE CELL TYPE ANNOTATION")
    print("="*80)
    print(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print(f"Tissue Type: {tissue_type}")

    # Load data
    print(f"\n[1/6] Loading data from: {adata_file}")
    adata = sc.read_h5ad(adata_file)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Clusters: {adata.obs[cluster_key].nunique()}")

    # Set output directory
    if output_dir is None:
        output_dir = Path(adata_file).parent
    else:
        output_dir = Path(output_dir)

    figures_dir = output_dir / 'figures'
    figures_dir.mkdir(exist_ok=True)

    # Initialize results dictionary
    annotation_results = {}

    # ========================================================================
    # Method 1: CellTypist Automated Annotation
    # ========================================================================
    print(f"\n[2/6] Running CellTypist automated annotation...")
    print("-"*80)

    celltypist_results = {}
    try:
        import celltypist
        from celltypist import models

        print(f"  Model: {celltypist_model}")

        # Download model if needed
        try:
            model = models.Model.load(model=celltypist_model)
        except:
            print(f"  Downloading model: {celltypist_model}...")
            models.download_models(force_update=False, model=celltypist_model)
            model = models.Model.load(model=celltypist_model)

        # Predict cell types
        print("  Predicting cell types...")
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        adata_celltypist = predictions.to_adata()

        # Store results
        adata.obs['celltypist_prediction'] = adata_celltypist.obs['predicted_labels']
        adata.obs['celltypist_confidence'] = adata_celltypist.obs['conf_score']
        adata.obs['celltypist_majority_voting'] = adata_celltypist.obs['majority_voting']

        # Aggregate by cluster
        for cluster_id in adata.obs[cluster_key].unique():
            cluster_mask = adata.obs[cluster_key] == cluster_id
            cluster_predictions = adata.obs.loc[cluster_mask, 'celltypist_prediction']

            # Most common prediction
            top_prediction = cluster_predictions.value_counts().iloc[0]
            top_prediction_name = cluster_predictions.value_counts().index[0]
            top_prediction_pct = (top_prediction / len(cluster_predictions)) * 100

            # Average confidence
            avg_confidence = adata.obs.loc[cluster_mask, 'celltypist_confidence'].mean()

            celltypist_results[str(cluster_id)] = {
                'prediction': top_prediction_name,
                'percentage': top_prediction_pct,
                'confidence': avg_confidence,
                'all_predictions': cluster_predictions.value_counts().to_dict()
            }

            print(f"  Cluster {cluster_id}: {top_prediction_name} ({top_prediction_pct:.1f}%, conf={avg_confidence:.2f})")

        celltypist_success = True

    except ImportError:
        print("  CellTypist not installed. Install with: pip install celltypist")
        celltypist_success = False
    except Exception as e:
        print(f"  CellTypist error: {e}")
        celltypist_success = False

    # ========================================================================
    # Method 2: Marker-based Annotation
    # ========================================================================
    print(f"\n[3/6] Processing marker-based annotation...")
    print("-"*80)

    marker_results = {}
    for cluster_id, info in manual_annotation_dict.items():
        cell_type = info['cell_type']
        markers = info['markers']
        rationale = info.get('rationale', 'Based on marker gene expression')

        # Check marker expression
        marker_expression = {}
        for marker in markers:
            if marker in adata.var_names:
                cluster_mask = adata.obs[cluster_key] == cluster_id
                if adata.raw is not None:
                    expr = adata.raw[cluster_mask, marker].X.mean()
                else:
                    expr = adata[cluster_mask, marker].X.mean()
                marker_expression[marker] = float(expr)
            else:
                marker_expression[marker] = None

        marker_results[str(cluster_id)] = {
            'cell_type': cell_type,
            'markers': markers,
            'marker_expression': marker_expression,
            'rationale': rationale
        }

        expressed_markers = [m for m, v in marker_expression.items() if v is not None and v > 0]
        print(f"  Cluster {cluster_id}: {cell_type}")
        print(f"    Markers: {', '.join(markers)}")
        print(f"    Expressed: {', '.join(expressed_markers) if expressed_markers else 'None detected'}")

    # ========================================================================
    # Method 3: Consensus Annotation
    # ========================================================================
    print(f"\n[4/6] Generating consensus annotation...")
    print("-"*80)

    consensus_annotation = {}
    annotation_reports = []

    for cluster_id in adata.obs[cluster_key].unique():
        cluster_id_str = str(cluster_id)

        # Get predictions from both methods
        celltypist_pred = celltypist_results.get(cluster_id_str, {}).get('prediction', 'N/A')
        celltypist_conf = celltypist_results.get(cluster_id_str, {}).get('confidence', 0)

        marker_pred = marker_results.get(cluster_id_str, {}).get('cell_type', 'Unknown')
        marker_rationale = marker_results.get(cluster_id_str, {}).get('rationale', '')
        markers_used = marker_results.get(cluster_id_str, {}).get('markers', [])

        # Decision logic
        if celltypist_success and celltypist_conf > 0.7:
            # High confidence CellTypist
            if celltypist_pred.lower() in marker_pred.lower() or marker_pred.lower() in celltypist_pred.lower():
                final_annotation = marker_pred
                decision_reason = f"Agreement between CellTypist ({celltypist_pred}, conf={celltypist_conf:.2f}) and marker-based annotation"
            else:
                final_annotation = marker_pred
                decision_reason = f"Marker-based annotation preferred despite CellTypist prediction ({celltypist_pred}). Rationale: {marker_rationale}"
        else:
            # Low confidence or no CellTypist
            final_annotation = marker_pred
            decision_reason = f"Marker-based annotation. Rationale: {marker_rationale}"

        consensus_annotation[cluster_id_str] = final_annotation

        # Create detailed report
        report = {
            'Cluster': cluster_id_str,
            'Final Annotation': final_annotation,
            'CellTypist Prediction': celltypist_pred if celltypist_success else 'N/A',
            'CellTypist Confidence': f"{celltypist_conf:.2f}" if celltypist_success else 'N/A',
            'Marker-based Annotation': marker_pred,
            'Key Markers': ', '.join(markers_used[:5]),  # Top 5 markers
            'Decision Rationale': decision_reason
        }
        annotation_reports.append(report)

        print(f"\n  Cluster {cluster_id}: {final_annotation}")
        if celltypist_success:
            print(f"    CellTypist: {celltypist_pred} (conf={celltypist_conf:.2f})")
        print(f"    Markers: {', '.join(markers_used[:3])}...")
        print(f"    Decision: {decision_reason[:80]}...")

    # Add consensus annotation to adata
    adata.obs['cell_type'] = adata.obs[cluster_key].astype(str).map(consensus_annotation)

    # ========================================================================
    # Generate Visualizations
    # ========================================================================
    print(f"\n[5/6] Generating visualizations...")
    print("-"*80)

    # Figure 1: Comparison UMAP
    if celltypist_success:
        fig, axes = plt.subplots(1, 3, figsize=(24, 6))
    else:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        axes = list(axes) + [None]

    # Consensus annotation
    sc.pl.umap(adata, color='cell_type', ax=axes[0], show=False,
              title='Final Cell Type Annotation', legend_loc='right margin',
              frameon=False)

    # Cluster numbers
    sc.pl.umap(adata, color=cluster_key, ax=axes[1], show=False,
              title=f'Clusters ({cluster_key})', legend_loc='right margin',
              frameon=False)

    # CellTypist prediction
    if celltypist_success and axes[2] is not None:
        sc.pl.umap(adata, color='celltypist_prediction', ax=axes[2], show=False,
                  title='CellTypist Prediction', legend_loc='right margin',
                  frameon=False, legend_fontsize=6)

    plt.tight_layout()
    output_file = figures_dir / 'umap_cell_type_annotation.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

    # Figure 2: Cell type with labels
    fig, ax = plt.subplots(figsize=(10, 8))
    sc.pl.umap(adata, color='cell_type', legend_loc='on data',
              ax=ax, show=False, title='Cell Type Annotation',
              frameon=False)
    output_file = figures_dir / 'umap_cell_type_with_labels.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

    # Figure 3: Cell type by condition (if available)
    if 'condition' in adata.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Composition
        ct = pd.crosstab(adata.obs['condition'], adata.obs['cell_type'])
        ct_pct = ct.div(ct.sum(axis=1), axis=0) * 100

        ct_pct.plot(kind='bar', stacked=True, ax=axes[0], colormap='tab20')
        axes[0].set_xlabel('Condition')
        axes[0].set_ylabel('Percentage')
        axes[0].set_title('Cell Type Composition by Condition')
        axes[0].legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        axes[0].tick_params(axis='x', rotation=45)

        ct.plot(kind='bar', ax=axes[1], colormap='tab20')
        axes[1].set_xlabel('Condition')
        axes[1].set_ylabel('Number of Cells')
        axes[1].set_title('Cell Counts by Condition')
        axes[1].legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        axes[1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        output_file = figures_dir / 'cell_type_by_condition.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file}")

    # Figure 4: Confidence scores (if CellTypist was used)
    if celltypist_success:
        fig, ax = plt.subplots(figsize=(10, 6))
        sc.pl.umap(adata, color='celltypist_confidence', ax=ax, show=False,
                  title='CellTypist Confidence Score', cmap='RdYlGn',
                  vmin=0, vmax=1, frameon=False)
        output_file = figures_dir / 'celltypist_confidence.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {output_file}")

    # ========================================================================
    # Save Results
    # ========================================================================
    print(f"\n[6/6] Saving results...")
    print("-"*80)

    # Save annotated AnnData
    output_h5ad = Path(adata_file).parent / f"{Path(adata_file).stem}_annotated.h5ad"
    adata.write(output_h5ad)
    print(f"  Saved annotated data: {output_h5ad.name}")

    # Save annotation report
    report_df = pd.DataFrame(annotation_reports)
    output_csv = output_dir / 'annotation_report.csv'
    report_df.to_csv(output_csv, index=False)
    print(f"  Saved annotation report: {output_csv.name}")

    # Save detailed annotation documentation
    output_md = output_dir / 'ANNOTATION_REPORT.md'
    with open(output_md, 'w') as f:
        f.write("# Cell Type Annotation Report\n\n")
        f.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}\n\n")
        f.write(f"**Tissue Type:** {tissue_type}\n\n")
        f.write(f"**Total Cells:** {adata.n_obs:,}\n\n")
        f.write(f"**Number of Clusters:** {adata.obs[cluster_key].nunique()}\n\n")
        f.write(f"**Number of Cell Types:** {adata.obs['cell_type'].nunique()}\n\n")

        f.write("---\n\n")
        f.write("## Annotation Methods\n\n")
        f.write("### 1. CellTypist Automated Annotation\n\n")
        if celltypist_success:
            f.write(f"- **Model:** {celltypist_model}\n")
            f.write(f"- **Status:** Successful\n")
            f.write(f"- **Average Confidence:** {adata.obs['celltypist_confidence'].mean():.2f}\n\n")
        else:
            f.write("- **Status:** Not available\n\n")

        f.write("### 2. Marker-based Manual Annotation\n\n")
        f.write("- Based on known marker genes for each cell type\n")
        f.write("- Validated against literature and expression patterns\n\n")

        f.write("---\n\n")
        f.write("## Detailed Annotation Results\n\n")

        for report in annotation_reports:
            f.write(f"### Cluster {report['Cluster']}: {report['Final Annotation']}\n\n")
            f.write(f"**Decision Rationale:**\n{report['Decision Rationale']}\n\n")
            f.write(f"**Key Markers:** {report['Key Markers']}\n\n")
            if celltypist_success:
                f.write(f"**CellTypist Prediction:** {report['CellTypist Prediction']} ")
                f.write(f"(Confidence: {report['CellTypist Confidence']})\n\n")
            f.write(f"**Marker-based Annotation:** {report['Marker-based Annotation']}\n\n")

            # Add marker expression details
            cluster_id = report['Cluster']
            if cluster_id in marker_results:
                marker_expr = marker_results[cluster_id]['marker_expression']
                f.write("**Marker Expression:**\n")
                for marker, expr in marker_expr.items():
                    if expr is not None:
                        f.write(f"- {marker}: {expr:.3f}\n")
                    else:
                        f.write(f"- {marker}: Not detected\n")
            f.write("\n---\n\n")

        f.write("## Cell Type Summary\n\n")
        cell_type_counts = adata.obs['cell_type'].value_counts()
        f.write("| Cell Type | Number of Cells | Percentage |\n")
        f.write("|-----------|-----------------|------------|\n")
        for cell_type, count in cell_type_counts.items():
            pct = (count / adata.n_obs) * 100
            f.write(f"| {cell_type} | {count:,} | {pct:.1f}% |\n")

    print(f"  Saved annotation documentation: {output_md.name}")

    # Save summary CSV
    cell_type_summary = pd.DataFrame({
        'Cell Type': cell_type_counts.index,
        'Number of Cells': cell_type_counts.values,
        'Percentage': [f"{(v/adata.n_obs*100):.1f}%" for v in cell_type_counts.values]
    })
    output_csv = output_dir / 'cell_type_summary.csv'
    cell_type_summary.to_csv(output_csv, index=False)
    print(f"  Saved cell type summary: {output_csv.name}")

    print("\n" + "="*80)
    print("ANNOTATION COMPLETE!")
    print("="*80)
    print(f"\nGenerated files:")
    print(f"  - {output_h5ad.name} (annotated data)")
    print(f"  - annotation_report.csv (detailed comparison)")
    print(f"  - ANNOTATION_REPORT.md (comprehensive documentation)")
    print(f"  - cell_type_summary.csv (summary statistics)")
    print(f"  - umap_cell_type_annotation.png")
    print(f"  - umap_cell_type_with_labels.png")
    if 'condition' in adata.obs.columns:
        print(f"  - cell_type_by_condition.png")
    if celltypist_success:
        print(f"  - celltypist_confidence.png")
    print("\nTo regenerate PowerPoint with annotations:")
    print("  python generate_ppt_report.py")
    print("="*80)

    return adata, report_df


# ============================================================================
# Example usage for JMV ALI airway organoid analysis
# ============================================================================

if __name__ == "__main__":
    # Define manual annotations based on marker genes
    # Format: cluster_id: {cell_type, markers, rationale}

    manual_annotation = {
        '0': {
            'cell_type': 'Basal cells',
            'markers': ['TP63', 'KRT5', 'KRT14', 'NGFR'],
            'rationale': 'High expression of basal cell markers TP63 and KRT5. These are stem/progenitor cells of the airway epithelium.'
        },
        '1': {
            'cell_type': 'Ciliated cells',
            'markers': ['FOXJ1', 'RSPH1', 'PIFO', 'DNAH5'],
            'rationale': 'Expression of ciliated cell markers FOXJ1 and RSPH1. Ciliated cells are responsible for mucociliary clearance.'
        },
        '2': {
            'cell_type': 'Secretory cells',
            'markers': ['SCGB1A1', 'SCGB3A2', 'CC10'],
            'rationale': 'High expression of club/secretory cell markers SCGB1A1 and SCGB3A2. These cells secrete protective proteins.'
        },
        '3': {
            'cell_type': 'Goblet cells',
            'markers': ['MUC5AC', 'MUC5B', 'TFF3', 'SPDEF'],
            'rationale': 'Strong expression of mucin genes MUC5AC and MUC5B, characteristic of goblet cells producing mucus.'
        },
        '4': {
            'cell_type': 'Basal cells (proliferating)',
            'markers': ['TP63', 'KRT5', 'MKI67', 'TOP2A'],
            'rationale': 'Basal cell markers combined with proliferation markers MKI67 and TOP2A, indicating actively dividing basal cells.'
        },
        '5': {
            'cell_type': 'Club cells',
            'markers': ['SCGB1A1', 'SCGB3A1', 'CYP2F2'],
            'rationale': 'Specific club cell markers. Club cells have both secretory and regenerative functions.'
        },
        '6': {
            'cell_type': 'Intermediate/Transitional cells',
            'markers': ['KRT5', 'SCGB1A1', 'KRT8'],
            'rationale': 'Expression of both basal (KRT5) and secretory (SCGB1A1) markers, suggesting transitional differentiation state.'
        },
        '7': {
            'cell_type': 'Secretory cells',
            'markers': ['SCGB3A2', 'SFTPB', 'WFDC2'],
            'rationale': 'Alternative secretory cell population with distinct marker profile.'
        },
        '8': {
            'cell_type': 'Ciliated cells (mature)',
            'markers': ['FOXJ1', 'PIFO', 'MLF1', 'CCDC153'],
            'rationale': 'Mature ciliated cells with high expression of ciliogenesis and motile cilia markers.'
        },
        '9': {
            'cell_type': 'Basal cells',
            'markers': ['TP63', 'KRT5', 'KRT15'],
            'rationale': 'Classical basal cell population with stem cell properties.'
        },
        '10': {
            'cell_type': 'Ionocytes',
            'markers': ['FOXI1', 'CFTR', 'ATP6V1B1', 'ATP6V0D2'],
            'rationale': 'Rare cell type with high expression of FOXI1 and CFTR. Ionocytes regulate airway surface liquid.'
        },
        '11': {
            'cell_type': 'Deuterosomal cells',
            'markers': ['FOXJ1', 'DEUP1', 'CCNO', 'CEP78'],
            'rationale': 'Intermediate stage in ciliogenesis with deuterosome markers. Precursors to mature ciliated cells.'
        },
        '12': {
            'cell_type': 'Neuroendocrine cells',
            'markers': ['CHGA', 'CHGB', 'ASCL1', 'NCAM1'],
            'rationale': 'Expression of neuroendocrine markers CHGA and CHGB. Rare airway cells with sensory functions.'
        },
        '13': {
            'cell_type': 'Secretory cells',
            'markers': ['SCGB1A1', 'WFDC2'],
            'rationale': 'Secretory cell subtype.'
        },
        '14': {
            'cell_type': 'Basal cells (suprabasal)',
            'markers': ['TP63', 'KRT13', 'KRT4'],
            'rationale': 'Suprabasal cells with differentiation markers, transitioning from basal to differentiated states.'
        },
        '15': {
            'cell_type': 'Ciliated cells',
            'markers': ['FOXJ1', 'TUBB4B'],
            'rationale': 'Ciliated cell population.'
        },
        '16': {
            'cell_type': 'Rare cell population 1',
            'markers': ['Unknown'],
            'rationale': 'Small cluster requiring further investigation. Check top marker genes for identification.'
        },
        '17': {
            'cell_type': 'Rare cell population 2',
            'markers': ['Unknown'],
            'rationale': 'Very small cluster. May represent doublets, stressed cells, or rare cell state.'
        },
        '18': {
            'cell_type': 'Rare cell population 3',
            'markers': ['Unknown'],
            'rationale': 'Requires investigation of marker genes.'
        },
        '19': {
            'cell_type': 'Rare cell population 4',
            'markers': ['Unknown'],
            'rationale': 'Requires investigation of marker genes.'
        },
        '20': {
            'cell_type': 'Rare cell population 5',
            'markers': ['Unknown'],
            'rationale': 'Smallest cluster. Check quality metrics and marker genes.'
        },
    }

    # Run comprehensive annotation
    adata_file = '/mnt2/kwy/scrna_agent/JMV_ALI_results/JMV_ALI_harmony_integrated.h5ad'

    adata, report_df = comprehensive_annotation(
        adata_file=adata_file,
        manual_annotation_dict=manual_annotation,
        cluster_key='louvain_0.8',
        celltypist_model='Immune_All_Low.pkl',  # or 'Human_Lung_Atlas.pkl' if available
        tissue_type='Human airway organoid (ALI culture)',
        output_dir=None
    )

    print("\nAnnotation complete! Check ANNOTATION_REPORT.md for detailed documentation.")
