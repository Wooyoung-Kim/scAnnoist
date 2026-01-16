#!/usr/bin/env python
# scRNA-seq Analysis Agent CLI
# Command-line interface with subcommands

import asyncio
import argparse
import sys
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv


def main():
    """Main CLI entry point."""
    load_dotenv(override=True)
    
    parser = argparse.ArgumentParser(
        prog="scrna-agent",
        description="scRNA-seq Analysis Agent - Multi-agent system for cell type annotation and integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run full pipeline
  scrna-agent pipeline --file data.h5ad --batch-key sample --tissue Spleen
  
  # Annotation only
  scrna-agent annotate --file data.h5ad --tissue "Immune system"
  
  # Integration only
  scrna-agent integrate --file data.h5ad --batch-key sample
  
  # Check progress
  scrna-agent status --dataset my_dataset
  
  # Search for reference datasets
  scrna-agent search "mouse spleen B cells"
"""
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # ============================================================
    # Pipeline command
    # ============================================================
    pipeline_parser = subparsers.add_parser(
        "pipeline",
        help="Run the complete analysis pipeline",
        description="Run: Annotation → Integration → scANVI Fine-tuning"
    )
    pipeline_parser.add_argument("--file", "-f", required=True, help="Path to h5ad file")
    pipeline_parser.add_argument("--batch-key", "-b", required=True, help="Batch/sample column")
    pipeline_parser.add_argument("--tissue", "-t", default="Immune system", help="Tissue type")
    pipeline_parser.add_argument("--species", "-s", default="mouse", help="Species (mouse/human)")
    pipeline_parser.add_argument("--output", "-o", default="./output", help="Output directory")
    pipeline_parser.add_argument("--model", "-m", default=None, help="LLM model name")
    pipeline_parser.add_argument("--resume", action="store_true", help="Resume interrupted run")
    
    # ============================================================
    # Annotate command
    # ============================================================
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Run cell type annotation only",
        description="Collaborative annotation using CellTypist, ScType, and Literature RAG"
    )
    annotate_parser.add_argument("--file", "-f", required=True, help="Path to h5ad file")
    annotate_parser.add_argument("--tissue", "-t", default=None, help="Tissue type (applied to all cells)")
    annotate_parser.add_argument("--tissue-map", type=str, help="CSV file with sample->tissue mapping")
    annotate_parser.add_argument("--sample-type", type=str, help="Sample type (applied to all cells)")
    annotate_parser.add_argument("--sample-type-map", type=str, help="CSV file with sample->sample_type mapping")
    annotate_parser.add_argument("--sample-col", default="sample", help="Column name for sample IDs")
    annotate_parser.add_argument("--allow-missing-tissue", action="store_true",
                                 help="Fill missing tissues with 'unknown'")
    annotate_parser.add_argument("--species", "-s", default="mouse", help="Species")
    annotate_parser.add_argument("--output", "-o", default="./annotated.h5ad", help="Output file")
    annotate_parser.add_argument("--methods", nargs="+",
                                 default=["celltypist", "sctype", "literature"],
                                 help="Methods to use")
    
    # ============================================================
    # Integrate command
    # ============================================================
    integrate_parser = subparsers.add_parser(
        "integrate",
        help="Run batch integration only",
        description="Batch correction using scVI, Harmony, Scanorama"
    )
    integrate_parser.add_argument("--file", "-f", required=True, help="Path to h5ad file")
    integrate_parser.add_argument("--batch-key", "-b", required=True, help="Batch column")
    integrate_parser.add_argument("--output", "-o", default="./integrated.h5ad", help="Output file")
    integrate_parser.add_argument("--methods", nargs="+", 
                                  default=["scvi", "harmony"],
                                  help="Methods to use")
    
    # ============================================================
    # Finetune command
    # ============================================================
    finetune_parser = subparsers.add_parser(
        "finetune",
        help="Run scANVI fine-tuning",
        description="Fine-tune integration using cell type labels"
    )
    finetune_parser.add_argument("--file", "-f", required=True, help="Path to h5ad file")
    finetune_parser.add_argument("--batch-key", "-b", required=True, help="Batch column")
    finetune_parser.add_argument("--label-key", "-l", required=True, help="Label column")
    finetune_parser.add_argument("--output", "-o", default="./finetuned.h5ad", help="Output file")
    finetune_parser.add_argument("--epochs", type=int, default=100, help="Max epochs")
    
    # ============================================================
    # Search command (SRAgent integration)
    # ============================================================
    search_parser = subparsers.add_parser(
        "search",
        help="Search SRA for reference datasets",
        description="Search SRA database for scRNA-seq datasets"
    )
    search_parser.add_argument("query", help="Search query")
    search_parser.add_argument("--species", "-s", default="mouse", help="Species filter")
    search_parser.add_argument("--max-results", "-n", type=int, default=10, help="Max results")
    
    # ============================================================
    # Models command
    # ============================================================
    models_parser = subparsers.add_parser(
        "models",
        help="Manage CellTypist models",
        description="Pull and manage CellTypist models"
    )
    models_parser.add_argument("action", choices=["pull", "list"], help="Action to perform")
    models_parser.add_argument("--tissue", help="Tissue type to pull models for")
    models_parser.add_argument("--all", action="store_true", help="Pull all models in registry")

    # ============================================================
    # Report command
    # ============================================================
    report_parser = subparsers.add_parser(
        "report",
        help="Generate PPT presentation from results",
        description="Generate PowerPoint presentation from analysis results"
    )
    report_parser.add_argument("--output-dir", "-o", required=True, help="Directory with analysis results")
    report_parser.add_argument("--title", "-t", default="scRNA-seq Analysis Results", help="Presentation title")

    # ============================================================
    # Status command
    # ============================================================
    status_parser = subparsers.add_parser(
        "status",
        help="Check pipeline status",
        description="View progress and status of pipeline runs"
    )
    status_parser.add_argument("--dataset", "-d", help="Dataset ID to check")
    status_parser.add_argument("--list", "-l", action="store_true", help="List all datasets")
    
    # ============================================================
    # Parse and execute
    # ============================================================
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(0)
    
    # Execute command
    if args.command == "pipeline":
        asyncio.run(run_pipeline_command(args))
    elif args.command == "annotate":
        asyncio.run(run_annotate_command(args))
    elif args.command == "integrate":
        asyncio.run(run_integrate_command(args))
    elif args.command == "finetune":
        asyncio.run(run_finetune_command(args))
    elif args.command == "search":
        run_search_command(args)
    elif args.command == "models":
        run_models_command(args)
    elif args.command == "report":
        run_report_command(args)
    elif args.command == "status":
        run_status_command(args)
    else:
        parser.print_help()


async def run_pipeline_command(args):
    """Run the full pipeline."""
    from scrna_agent.agents.pipeline_agent import run_pipeline
    
    print(f"🚀 Starting scRNA-seq Analysis Pipeline")
    print(f"   File: {args.file}")
    print(f"   Batch key: {args.batch_key}")
    print(f"   Tissue: {args.tissue}")
    print(f"   Species: {args.species}")
    print()
    
    result = await run_pipeline(
        file_path=args.file,
        batch_key=args.batch_key,
        tissue=args.tissue,
        species=args.species,
        output_dir=args.output,
        model_name=args.model,
    )
    
    print("\n✅ Pipeline completed!")
    print(f"   Outputs saved to: {args.output}")


async def run_annotate_command(args):
    """Run annotation only."""
    from scrna_agent.agents.annotation_coordinator import create_annotation_coordinator
    from scrna_agent.tools import load_adata_for_annotation
    from langchain_core.messages import HumanMessage
    
    print(f"🏷️  Starting Cell Type Annotation")
    print(f"   File: {args.file}")
    print(f"   Tissue: {args.tissue}")
    print(f"   Methods: {', '.join(args.methods)}")
    print()
    
    # Load data
    load_adata_for_annotation(args.file)
    
    # Create coordinator and run
    coordinator = create_annotation_coordinator()
    
    task = f"""
    Annotate cells in this {args.species} {args.tissue} dataset.
    
    Methods to use: {', '.join(args.methods)}
    
    1. Run each method
    2. Compare results
    3. Create consensus annotation
    4. Validate with markers
    """
    
    result = await coordinator.ainvoke(
        {"message": task},
        config={"configurable": {"species": args.species}}
    )
    
    print("\n✅ Annotation completed!")
    print(f"   Output: {args.output}")


async def run_integrate_command(args):
    """Run integration only."""
    from scrna_agent.tools import (
        load_adata_for_integration,
        prepare_for_integration,
        run_scvi,
        run_harmony,
        save_integrated_adata,
    )
    
    print(f"🔗 Starting Batch Integration")
    print(f"   File: {args.file}")
    print(f"   Batch key: {args.batch_key}")
    print(f"   Methods: {', '.join(args.methods)}")
    print()
    
    # Load and prepare
    load_adata_for_integration(args.file)
    prepare_for_integration(batch_key=args.batch_key)
    
    # Run methods
    if "scvi" in args.methods:
        print("   Running scVI...")
        run_scvi(batch_key=args.batch_key)
    
    if "harmony" in args.methods:
        print("   Running Harmony...")
        run_harmony(batch_key=args.batch_key)
    
    # Save
    save_integrated_adata(args.output)
    
    print("\n✅ Integration completed!")
    print(f"   Output: {args.output}")


async def run_finetune_command(args):
    """Run scANVI fine-tuning."""
    from scrna_agent.tools import (
        load_data_for_pipeline,
        setup_scanvi_finetuning,
        run_scanvi_finetuning,
        save_pipeline_results,
    )
    
    print(f"🎯 Starting scANVI Fine-tuning")
    print(f"   File: {args.file}")
    print(f"   Batch key: {args.batch_key}")
    print(f"   Label key: {args.label_key}")
    print()
    
    # Load data
    load_data_for_pipeline(args.file, args.batch_key)
    
    # Setup and run
    setup_scanvi_finetuning(annotation_col=args.label_key)
    run_scanvi_finetuning(max_epochs_scanvi=args.epochs)
    
    # Save
    save_pipeline_results(args.output)
    
    print("\n✅ Fine-tuning completed!")
    print(f"   Output: {args.output}")


def run_search_command(args):
    """Search SRA for datasets."""
    from scrna_agent.tools.sragent_tools import search_sra_for_reference
    
    print(f"🔍 Searching SRA: {args.query}")
    print(f"   Species: {args.species}")
    print()
    
    result = search_sra_for_reference(
        query=args.query,
        species=args.species,
        max_results=args.max_results,
    )
    
    print(result)


def run_status_command(args):
    """Check pipeline status."""
    from scrna_agent.tools.state_tools import get_state_manager

    state_mgr = get_state_manager()

    if args.list:
        datasets = state_mgr.list_datasets()
        if not datasets:
            print("No datasets tracked yet.")
        else:
            print("Tracked Datasets:")
            for ds in datasets:
                print(f"  - {ds['id']}: {ds['status']} ({ds['steps_completed']} steps)")
    elif args.dataset:
        progress = state_mgr.get_progress(args.dataset)
        print(f"Dataset: {args.dataset}")
        print(f"Status: {progress['status']}")
        print(f"Progress: {progress['completed_steps']}/{progress['total_steps']} "
              f"({progress['progress_percent']:.1f}%)")
        if progress['remaining_steps']:
            print(f"Next step: {progress['remaining_steps'][0]}")
    else:
        print("Use --list to see all datasets or --dataset ID to check specific dataset")


def run_models_command(args):
    """Manage CellTypist models."""
    from scrna_agent.tools.model_management import pull_celltypist_models, list_available_models

    if args.action == "list":
        print("Available CellTypist models in registry:")
        list_available_models()
    elif args.action == "pull":
        if args.all:
            print("Pulling all models from registry...")
            pull_celltypist_models(tissue=None, pull_all=True)
        elif args.tissue:
            print(f"Pulling models for tissue: {args.tissue}")
            pull_celltypist_models(tissue=args.tissue, pull_all=False)
        else:
            print("Error: Specify --tissue <name> or --all")
            sys.exit(1)


def run_report_command(args):
    """Generate PPT presentation from results."""
    from scrna_agent.tools.report_generator import generate_report

    print(f"📊 Generating presentation from: {args.output_dir}")
    print(f"   Title: {args.title}")
    print()

    try:
        pdf_path = generate_report(
            output_dir=args.output_dir,
            title=args.title
        )
        print()
        print(f"✅ Presentation generated successfully!")
        print(f"   PDF: {pdf_path}")
    except Exception as e:
        print(f"❌ Error generating presentation: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
