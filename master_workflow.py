#!/usr/bin/env python3
"""
Master Workflow for scRNA-seq Analysis with Skills Integration
Combines:
- Standard analysis pipeline
- Multi-agent annotation discussion
- Scientific visualization (using scientific-visualization skill)
- Professional PowerPoint (using pptx and scientific-slides skills)
"""

import subprocess
from pathlib import Path
from datetime import datetime
import sys

def print_section(title):
    """Print a formatted section header"""
    print("\n" + "="*80)
    print(title.center(80))
    print("="*80 + "\n")

def print_step(step_num, total_steps, description):
    """Print a formatted step"""
    print(f"\n[{step_num}/{total_steps}] {description}")
    print("-"*80)

def run_script(script_name, description):
    """Run a Python script and handle errors"""
    print(f"\nRunning: {script_name}")
    try:
        result = subprocess.run(
            ['python', script_name],
            check=True,
            capture_output=True,
            text=True
        )
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {script_name} failed")
        print(e.stderr)
        return False

def main():
    """Execute complete analysis workflow"""

    print_section("MASTER WORKFLOW: scRNA-seq Analysis with AI Skills")

    print(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Configuration
    data_dir = Path('/home/kwy7605/data_61/SARS/Count/JMV_ALI')
    output_dir = Path('/mnt2/kwy/scrna_agent/JMV_ALI_results')
    analysis_name = "JMV ALI Organoid scRNA-seq Analysis"

    print(f"Data Directory: {data_dir}")
    print(f"Output Directory: {output_dir}")
    print(f"Analysis Name: {analysis_name}")

    # Workflow steps
    total_steps = 6

    # ========================================================================
    # STEP 1: Basic Analysis Pipeline
    # ========================================================================
    print_step(1, total_steps, "Running Basic Analysis Pipeline")

    print("""
This step will:
- Load and merge all samples
- Perform quality control
- Normalize data
- Integrate with Harmony
- Cluster cells
- Identify marker genes
- Generate basic visualizations
""")

    input("Press Enter to start basic analysis...")

    success = run_script('JMV_ALI_analysis.py', 'Basic analysis')

    if not success:
        print("\n⚠️  Basic analysis failed. Please check errors and try again.")
        return

    print("\n✓ Basic analysis complete!")
    print(f"  - Processed data saved: {output_dir / 'JMV_ALI_harmony_integrated.h5ad'}")
    print(f"  - Figures generated in: {output_dir / 'figures'}")
    print(f"  - Marker genes saved as CSV files")

    # ========================================================================
    # STEP 2: Prepare for Multi-Agent Annotation
    # ========================================================================
    print_step(2, total_steps, "Preparing Multi-Agent Annotation System")

    print("""
This step will:
- Extract cluster information
- Generate prompts for 4 specialized agents:
  1. Cell Type Expert
  2. Computational Biologist
  3. Literature Curator
  4. Critical Evaluator
- Create workflow instructions
""")

    input("Press Enter to prepare multi-agent system...")

    success = run_script('multi_agent_annotation.py', 'Multi-agent preparation')

    if not success:
        print("\n⚠️  Preparation failed. Please check errors.")
        return

    print("\n✓ Multi-agent system prepared!")
    print(f"\n⚠️  MANUAL STEP REQUIRED:")
    print(f"  Please review: {output_dir / 'MULTI_AGENT_WORKFLOW.md'}")
    print(f"  Execute the 4 agents using Claude Code Task tool")
    print(f"  Save outputs in: {output_dir / 'agent_outputs/'}")

    proceed = input("\nHave you completed the multi-agent annotation? (yes/no): ")

    if proceed.lower() != 'yes':
        print("\nPlease complete multi-agent annotation and re-run this workflow.")
        print("You can resume from Step 3 by modifying the script.")
        return

    # ========================================================================
    # STEP 3: Compile Agent Annotations
    # ========================================================================
    print_step(3, total_steps, "Compiling Multi-Agent Annotations")

    print("""
This step will:
- Combine annotations from all agents
- Generate consensus annotations
- Create ANNOTATION_REPORT.md
- Update h5ad file with annotations
- Generate annotation visualizations
""")

    # Note: This would require a compilation script
    print("\n⚠️  Note: compile_agent_annotations.py would be run here")
    print("  (Script to be created based on agent outputs)")

    input("Press Enter to continue...")

    # ========================================================================
    # STEP 4: Generate Scientific Visualizations
    # ========================================================================
    print_step(4, total_steps, "Generating Scientific Visualizations")

    print("""
This step will:
- Prepare visualization specifications
- Generate publication-quality figures using scientific-visualization skill
- Create 8 comprehensive figures:
  1. QC Metrics Panel
  2. UMAP Comprehensive
  3. Marker Gene Heatmap
  4. Cell Type Proportions
  5. Marker Gene Dotplot
  6. Volcano Plot
  7. PCA Biplot
  8. Gene Expression Violins
""")

    input("Press Enter to prepare visualizations...")

    success = run_script('generate_scientific_visualizations.py', 'Visualization preparation')

    if not success:
        print("\n⚠️  Visualization preparation failed.")
        return

    print("\n✓ Visualization specifications prepared!")
    print(f"\n⚠️  SKILL EXECUTION REQUIRED:")
    print(f"  In Claude Code, run: /scientific-visualization")
    print(f"  Provide specification file: {output_dir / 'visualization_specifications.json'}")
    print(f"  Provide data file: {output_dir / 'JMV_ALI_harmony_integrated_annotated.h5ad'}")

    proceed = input("\nHave you generated the scientific visualizations? (yes/no): ")

    if proceed.lower() != 'yes':
        print("\nPlease generate visualizations using the skill and re-run.")
        return

    # ========================================================================
    # STEP 5: Generate Scientific PowerPoint
    # ========================================================================
    print_step(5, total_steps, "Generating Scientific PowerPoint Presentation")

    print("""
This step will:
- Prepare presentation outline
- Use pptx and scientific-slides skills
- Create 25+ professional slides
- Integrate all figures and annotations
- Include agent discussion summaries
""")

    input("Press Enter to prepare PowerPoint...")

    success = run_script('generate_scientific_ppt.py', 'PowerPoint preparation')

    if not success:
        print("\n⚠️  PowerPoint preparation failed.")
        return

    print("\n✓ PowerPoint specifications prepared!")
    print(f"\n⚠️  SKILL EXECUTION REQUIRED:")
    print(f"  In Claude Code, run: /scientific-slides")
    print(f"  Or run: /pptx")
    print(f"  Provide instruction file: {output_dir / 'ppt_generation_instructions.txt'}")

    proceed = input("\nHave you generated the PowerPoint presentation? (yes/no): ")

    if proceed.lower() != 'yes':
        print("\nPlease generate PowerPoint using the skill and re-run.")
        return

    # ========================================================================
    # STEP 6: Final Summary and Validation
    # ========================================================================
    print_step(6, total_steps, "Final Summary and Validation")

    print("\nValidating generated files...")

    expected_files = {
        'Data': [
            output_dir / 'JMV_ALI_harmony_integrated.h5ad',
            output_dir / 'JMV_ALI_harmony_integrated_annotated.h5ad',
            output_dir / 'metadata.csv',
        ],
        'Reports': [
            output_dir / 'ANNOTATION_REPORT.md',
            output_dir / 'MULTI_AGENT_WORKFLOW.md',
            output_dir / 'VISUALIZATION_INSTRUCTIONS.md',
        ],
        'Specifications': [
            output_dir / 'cluster_information_for_agents.json',
            output_dir / 'visualization_specifications.json',
            output_dir / 'presentation_outline.json',
        ],
        'Outputs': [
            output_dir / 'JMV_ALI_Scientific_Presentation.pptx',
        ]
    }

    all_present = True
    for category, files in expected_files.items():
        print(f"\n{category}:")
        for file in files:
            if file.exists():
                size_mb = file.stat().st_size / (1024 * 1024)
                print(f"  ✓ {file.name} ({size_mb:.1f} MB)")
            else:
                print(f"  ✗ {file.name} (MISSING)")
                all_present = False

    print("\nFigures Directory:")
    figures_dir = output_dir / 'figures'
    if figures_dir.exists():
        figure_count = len(list(figures_dir.glob('*.png')))
        print(f"  ✓ {figure_count} PNG files")
    else:
        print(f"  ✗ Figures directory missing")
        all_present = False

    # Final Summary
    print_section("WORKFLOW COMPLETE!")

    print(f"End Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    if all_present:
        print("✓ All expected files are present")
    else:
        print("⚠️  Some files are missing - please review above")

    print(f"\nResults Location: {output_dir}")
    print("\nGenerated Outputs:")
    print("  1. Processed scRNA-seq data (h5ad)")
    print("  2. Multi-agent annotation report (MD)")
    print("  3. Publication-quality figures (PNG)")
    print("  4. Scientific PowerPoint presentation (PPTX)")
    print("  5. Comprehensive documentation (MD, JSON, CSV)")

    print("\nKey Features:")
    print("  • Multi-agent cell type annotation with rationale")
    print("  • CellTypist + Marker-based consensus")
    print("  • Publication-quality visualizations (300 DPI)")
    print("  • Professional scientific presentation")
    print("  • Complete analysis documentation")

    print("\nNext Steps:")
    print("  1. Review ANNOTATION_REPORT.md for annotation rationale")
    print("  2. Check publication figures in figures/")
    print("  3. Review and present the PowerPoint")
    print("  4. Perform downstream analyses (DE, pathway, etc.)")
    print("  5. Prepare manuscript figures")

    print("\n" + "="*80)
    print("Thank you for using the scRNA-seq Analysis Pipeline!")
    print("="*80)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nWorkflow interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
