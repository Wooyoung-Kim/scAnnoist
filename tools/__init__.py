# Tools subpackage
# All tool implementations for the scRNA-seq analysis pipeline

from scrna_agent.tools.annotation_tools import (
    # Data loading
    load_adata_for_annotation,
    # CellTypist
    celltypist_list_models,
    celltypist_download_model,
    celltypist_annotate,
    celltypist_annotate_tissue_aware,
    celltypist_get_probabilities,
    add_tissue_metadata,
    # scANVI annotation
    scanvi_setup_reference,
    scanvi_train_and_transfer,
    # ScType
    sctype_annotate,
    # Comparison
    compare_annotations,
    create_consensus_annotation,
    validate_with_markers,
    save_annotations,
    plot_annotation_umap,
    # Literature RAG
    search_pubmed_markers,
    get_cellmarker_database,
    annotate_with_literature_markers,
    combine_pubmed_and_database_markers,
    # Annotation Integration
    integrate_multiple_annotations,
    compare_celltypist_with_literature,
    create_optimal_annotation,
)

from scrna_agent.tools.mcp_annotation_tools import (
    # MCP-based literature search
    mcp_search_marker_literature,
    mcp_get_marker_evidence,
    mcp_get_full_text,
    mcp_find_related_markers,
    mcp_validate_cell_type_annotation,
    # Hybrid tools
    hybrid_search_markers,
)

from scrna_agent.tools.integration_tools import (
    # Data loading
    load_adata_for_integration,
    prepare_for_integration,
    # scVI/scANVI
    run_scvi,
    run_scanvi_from_scvi,
    # Harmony/Scanorama
    run_harmony,
    run_scanorama,
    # Benchmarking
    benchmark_integration,
    select_best_integration,
    # Post-processing
    compute_umap_from_embedding,
    plot_integration_comparison,
    save_integration_model,
    save_integrated_adata,
)

from scrna_agent.tools.pipeline_tools import (
    # Pipeline state management
    load_data_for_pipeline,
    get_pipeline_state,
    set_pipeline_state,
    get_pipeline_status,
    # scANVI fine-tuning
    setup_scanvi_finetuning,
    run_scanvi_finetuning,
    # Post-processing
    compute_integrated_umap,
    plot_pipeline_results,
    save_pipeline_results,
    # Comprehensive visualization orchestration
    generate_qc_visualizations,
    generate_integration_visualizations,
    generate_condition_analysis,
    generate_comprehensive_report,
    run_complete_visualization_pipeline,
)

from scrna_agent.tools.qc_visualization import (
    plot_qc_metrics_before,
    plot_qc_metrics_after,
    plot_qc_comparison,
    get_qc_plot_paths,
)

from scrna_agent.tools.integration_visualization import (
    plot_integration_before,
    plot_integration_after,
    plot_condition_celltype_distribution,
    plot_integration_comparison,
    get_integration_plot_paths,
)

from scrna_agent.tools.comprehensive_report import (
    initialize_report,
    add_qc_metrics,
    add_annotation_results,
    add_integration_results,
    add_condition_analysis,
    add_figure,
    generate_markdown_report,
    get_report_data,
)

from scrna_agent.tools.report_generator import (
    generate_report,
    PresentationGenerator,
)

__all__ = [
    # Annotation tools
    "load_adata_for_annotation",
    "celltypist_list_models",
    "celltypist_download_model",
    "celltypist_annotate",
    "celltypist_annotate_tissue_aware",
    "celltypist_get_probabilities",
    "add_tissue_metadata",
    "scanvi_setup_reference",
    "scanvi_train_and_transfer",
    "sctype_annotate",
    "compare_annotations",
    "create_consensus_annotation",
    "validate_with_markers",
    "save_annotations",
    "plot_annotation_umap",
    # Literature RAG tools
    "search_pubmed_markers",
    "get_cellmarker_database",
    "annotate_with_literature_markers",
    "combine_pubmed_and_database_markers",
    # MCP-based annotation tools
    "mcp_search_marker_literature",
    "mcp_get_marker_evidence",
    "mcp_get_full_text",
    "mcp_find_related_markers",
    "mcp_validate_cell_type_annotation",
    "hybrid_search_markers",
    # Annotation Integration tools
    "integrate_multiple_annotations",
    "compare_celltypist_with_literature",
    "create_optimal_annotation",
    # Integration tools
    "load_adata_for_integration",
    "prepare_for_integration",
    "run_scvi",
    "run_scanvi_from_scvi",
    "run_harmony",
    "run_scanorama",
    "benchmark_integration",
    "select_best_integration",
    "compute_umap_from_embedding",
    "plot_integration_comparison",
    "save_integration_model",
    "save_integrated_adata",
    # Pipeline tools
    "load_data_for_pipeline",
    "get_pipeline_state",
    "set_pipeline_state",
    "get_pipeline_status",
    "setup_scanvi_finetuning",
    "run_scanvi_finetuning",
    "compute_integrated_umap",
    "plot_pipeline_results",
    "save_pipeline_results",
    # Visualization orchestration tools
    "generate_qc_visualizations",
    "generate_integration_visualizations",
    "generate_condition_analysis",
    "generate_comprehensive_report",
    "run_complete_visualization_pipeline",
    # QC visualization tools
    "plot_qc_metrics_before",
    "plot_qc_metrics_after",
    "plot_qc_comparison",
    "get_qc_plot_paths",
    # Integration visualization tools
    "plot_integration_before",
    "plot_integration_after",
    "plot_condition_celltype_distribution",
    "plot_integration_comparison",
    "get_integration_plot_paths",
    # Comprehensive report tools
    "initialize_report",
    "add_qc_metrics",
    "add_annotation_results",
    "add_integration_results",
    "add_condition_analysis",
    "add_figure",
    "generate_markdown_report",
    "get_report_data",
    # Report generation
    "generate_report",
    "PresentationGenerator",
]
