# scrna-agent
# Single-cell RNA-seq Analysis Pipeline CLI Tool

__version__ = "0.1.0"
__author__ = "BioAnalysis Team"

# MVP Pipeline (always available)
from scrna_agent.pipeline import (
    PipelineConfig,
    PipelineResults,
    run_pipeline,
    create_toy_dataset,
    load_data,
    run_qc,
    normalize_and_scale,
    run_dimensionality_reduction,
    run_clustering,
    find_markers,
)

# Tissue handling and annotation
from scrna_agent.tissue import (
    TissueConfig,
    AnnotationConfig,
    TissueMapError,
    load_tissue_map,
    assign_tissue,
    load_markers,
    get_tissue_markers,
    annotate_cells,
    generate_tissue_summary,
    generate_cell_metadata,
)

# Convenience alias
load_input = load_data

__all__ = [
    # Config
    "PipelineConfig",
    "PipelineResults",
    "TissueConfig",
    "AnnotationConfig",
    # Pipeline functions
    "run_pipeline",
    "create_toy_dataset",
    "load_data",
    "load_input",  # Alias for load_data
    "run_qc",
    "normalize_and_scale",
    "run_dimensionality_reduction",
    "run_clustering",
    "find_markers",
    # Tissue and annotation
    "TissueMapError",
    "load_tissue_map",
    "assign_tissue",
    "load_markers",
    "get_tissue_markers",
    "annotate_cells",
    "generate_tissue_summary",
    "generate_cell_metadata",
]

# Optional: Multi-agent mode (requires langchain)
def _load_agents():
    """Lazy load agent modules (requires [full] dependencies)."""
    try:
        from scrna_agent.agents.pipeline_agent import (
            create_scrna_pipeline_agent,
            run_pipeline as run_agent_pipeline,
        )
        from scrna_agent.agents.annotation_coordinator import (
            create_annotation_coordinator,
        )
        return {
            "create_scrna_pipeline_agent": create_scrna_pipeline_agent,
            "run_agent_pipeline": run_agent_pipeline,
            "create_annotation_coordinator": create_annotation_coordinator,
        }
    except ImportError:
        return None

agents = _load_agents()
