#!/usr/bin/env python
"""
scrna-agent CLI
Command-line interface for scRNA-seq analysis pipeline
"""

import argparse
import logging
import sys
import os
from pathlib import Path

from scrna_agent.pipeline import (
    PipelineConfig,
    run_pipeline,
    create_toy_dataset,
)


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def cmd_run(args):
    """Run the analysis pipeline."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Build configuration
    config = PipelineConfig()

    # Override from config file if provided
    if args.config:
        logger.info(f"Loading config from {args.config}")
        config = PipelineConfig.from_yaml(args.config)

    # Override from CLI arguments
    if args.min_genes is not None:
        config.min_genes = args.min_genes
    if args.min_cells is not None:
        config.min_cells = args.min_cells
    if args.max_mt_pct is not None:
        config.max_mt_pct = args.max_mt_pct
    if args.resolution is not None:
        config.resolution = args.resolution
    if args.seed is not None:
        config.seed = args.seed
    if args.clustering:
        config.clustering_method = args.clustering

    # Tissue and annotation parameters
    if args.tissue:
        config.tissue = args.tissue
    if args.tissue_map:
        config.tissue_map_path = args.tissue_map
    if args.sample_col:
        config.sample_col = args.sample_col
    if args.allow_missing_tissue:
        config.allow_missing_tissue = True
    if args.markers:
        config.markers_path = args.markers
    if args.no_annotation:
        config.run_annotation = False

    logger.info("=" * 60)
    logger.info("scRNA-seq Analysis Pipeline")
    logger.info("=" * 60)
    logger.info(f"Input:  {args.input}")
    logger.info(f"Output: {args.out}")
    logger.info("=" * 60)

    # Detect input format
    input_path = Path(args.input)
    if args.format:
        input_format = args.format
    elif input_path.suffix == ".h5ad":
        input_format = "h5ad"
    elif input_path.is_dir():
        input_format = "10x"
    else:
        input_format = "auto"

    try:
        results = run_pipeline(
            input_path=str(args.input),
            output_dir=args.out,
            config=config,
            input_format=input_format
        )

        logger.info("=" * 60)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 60)
        logger.info("Output files:")
        for f in results.files_created:
            logger.info(f"  - {f}")
        logger.info("=" * 60)

        return 0

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_prerun(args):
    """Run pipeline on CellRanger output directory."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    # Check for CellRanger structure
    input_path = Path(args.input)

    # Look for filtered matrix
    possible_paths = [
        input_path / "outs" / "filtered_feature_bc_matrix",
        input_path / "filtered_feature_bc_matrix",
        input_path,
    ]

    matrix_path = None
    for p in possible_paths:
        if p.exists() and (p / "matrix.mtx.gz").exists():
            matrix_path = p
            break
        elif p.exists() and (p / "matrix.mtx").exists():
            matrix_path = p
            break

    if matrix_path is None:
        logger.error(f"Could not find 10x matrix files in {args.input}")
        logger.error("Expected files: matrix.mtx[.gz], barcodes.tsv[.gz], features.tsv[.gz]")
        return 1

    logger.info("=" * 60)
    logger.info("scRNA-seq Analysis Pipeline (CellRanger Input)")
    logger.info("=" * 60)
    logger.info(f"CellRanger dir: {args.input}")
    logger.info(f"Matrix path:    {matrix_path}")
    logger.info(f"Output:         {args.out}")
    logger.info("=" * 60)

    # Build configuration
    config = PipelineConfig()
    if args.config:
        config = PipelineConfig.from_yaml(args.config)

    try:
        results = run_pipeline(
            input_path=str(matrix_path),
            output_dir=args.out,
            config=config,
            input_format="10x"
        )

        logger.info("=" * 60)
        logger.info("Pipeline completed successfully!")
        logger.info("=" * 60)

        return 0

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


def cmd_generate_toy(args):
    """Generate a toy dataset for testing."""
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)

    output_path = args.output or "./data/toy/toy_data.h5ad"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    logger.info(f"Generating toy dataset: {args.n_cells} cells, {args.n_genes} genes")

    adata = create_toy_dataset(
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        n_clusters=args.n_clusters,
        output_path=output_path,
        seed=args.seed
    )

    logger.info(f"Saved to: {output_path}")
    logger.info(f"Shape: {adata.n_obs} cells x {adata.n_vars} genes")

    return 0


def cmd_config(args):
    """Generate a default config file."""
    config = PipelineConfig()
    config.to_yaml(args.output)
    print(f"Default config saved to: {args.output}")
    return 0


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        prog="scrna-agent",
        description="scRNA-seq Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run analysis on h5ad file
  scrna-agent run --input data.h5ad --out ./results

  # Run on 10x CellRanger output
  scrna-agent prerun --input /path/to/cellranger_output --out ./results

  # Generate toy dataset for testing
  scrna-agent generate-toy --output ./data/toy/toy_data.h5ad

  # Generate default config file
  scrna-agent config --output config.yaml

  # Run with custom parameters
  scrna-agent run --input data.h5ad --out ./results --min-genes 300 --resolution 0.8
"""
    )

    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # ============================================================
    # run command
    # ============================================================
    run_parser = subparsers.add_parser(
        "run",
        help="Run the analysis pipeline",
        description="Run scRNA-seq analysis pipeline on input data"
    )
    run_parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to input file (.h5ad) or directory (10x format)"
    )
    run_parser.add_argument(
        "--out", "-o",
        required=True,
        help="Output directory for results"
    )
    run_parser.add_argument(
        "--config", "-c",
        help="Path to YAML config file"
    )
    run_parser.add_argument(
        "--format", "-f",
        choices=["h5ad", "10x", "10x_h5", "auto"],
        help="Input format (default: auto-detect)"
    )
    # QC parameters
    run_parser.add_argument(
        "--min-genes",
        type=int,
        help="Minimum genes per cell (default: 200)"
    )
    run_parser.add_argument(
        "--min-cells",
        type=int,
        help="Minimum cells per gene (default: 3)"
    )
    run_parser.add_argument(
        "--max-mt-pct",
        type=float,
        help="Maximum mitochondrial percentage (default: 20)"
    )
    # Clustering parameters
    run_parser.add_argument(
        "--resolution",
        type=float,
        help="Clustering resolution (default: 0.5)"
    )
    run_parser.add_argument(
        "--clustering",
        choices=["leiden", "louvain"],
        help="Clustering method (default: leiden)"
    )
    # Tissue and annotation parameters
    run_parser.add_argument(
        "--tissue", "-t",
        help="Tissue type for all cells (e.g., 'lung', 'liver')"
    )
    run_parser.add_argument(
        "--tissue-map",
        help="Path to CSV file mapping sample -> tissue (columns: sample, tissue)"
    )
    run_parser.add_argument(
        "--sample-col",
        default="sample",
        help="Column name for sample IDs (default: sample)"
    )
    run_parser.add_argument(
        "--allow-missing-tissue",
        action="store_true",
        help="Allow samples without tissue assignment (fill with 'unknown')"
    )
    run_parser.add_argument(
        "--markers",
        help="Path to custom markers YAML file"
    )
    run_parser.add_argument(
        "--no-annotation",
        action="store_true",
        help="Skip marker-based cell type annotation"
    )
    # Other parameters
    run_parser.add_argument(
        "--seed",
        type=int,
        help="Random seed for reproducibility (default: 42)"
    )
    run_parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    run_parser.set_defaults(func=cmd_run)

    # ============================================================
    # prerun command
    # ============================================================
    prerun_parser = subparsers.add_parser(
        "prerun",
        help="Run pipeline on CellRanger output",
        description="Run analysis on 10x CellRanger output directory"
    )
    prerun_parser.add_argument(
        "--input", "-i",
        required=True,
        help="Path to CellRanger output directory"
    )
    prerun_parser.add_argument(
        "--out", "-o",
        required=True,
        help="Output directory for results"
    )
    prerun_parser.add_argument(
        "--config", "-c",
        help="Path to YAML config file"
    )
    prerun_parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    prerun_parser.set_defaults(func=cmd_prerun)

    # ============================================================
    # generate-toy command
    # ============================================================
    toy_parser = subparsers.add_parser(
        "generate-toy",
        help="Generate toy dataset for testing",
        description="Create a synthetic scRNA-seq dataset"
    )
    toy_parser.add_argument(
        "--output", "-o",
        help="Output path (default: ./data/toy/toy_data.h5ad)"
    )
    toy_parser.add_argument(
        "--n-cells",
        type=int,
        default=500,
        help="Number of cells (default: 500)"
    )
    toy_parser.add_argument(
        "--n-genes",
        type=int,
        default=1000,
        help="Number of genes (default: 1000)"
    )
    toy_parser.add_argument(
        "--n-clusters",
        type=int,
        default=5,
        help="Number of clusters (default: 5)"
    )
    toy_parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)"
    )
    toy_parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    toy_parser.set_defaults(func=cmd_generate_toy)

    # ============================================================
    # config command
    # ============================================================
    config_parser = subparsers.add_parser(
        "config",
        help="Generate default config file",
        description="Create a YAML config file with default parameters"
    )
    config_parser.add_argument(
        "--output", "-o",
        default="config.yaml",
        help="Output path (default: config.yaml)"
    )
    config_parser.set_defaults(func=cmd_config)

    # ============================================================
    # Parse and execute
    # ============================================================
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    # Execute the command
    sys.exit(args.func(args))


if __name__ == "__main__":
    main()
