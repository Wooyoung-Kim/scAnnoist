# Main Pipeline Agent
# Orchestrates the complete scRNA-seq analysis workflow

import asyncio
from typing import Annotated, Optional, Callable, Dict, Any

from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables.config import RunnableConfig
from langchain_core.messages import HumanMessage, AIMessage

from scrna_agent.utils import set_model
from scrna_agent.agents.annotation_coordinator import create_annotation_coordinator
from scrna_agent.tools.pipeline_tools import (
    load_data_for_pipeline,
    setup_scanvi_finetuning,
    run_scanvi_finetuning,
    compute_integrated_umap,
    plot_pipeline_results,
    save_pipeline_results,
    get_pipeline_status,
)
from scrna_agent.tools.integration_tools import (
    prepare_for_integration,
    run_scvi,
    run_harmony,
    benchmark_integration,
    select_best_integration,
)
from scrna_agent.skills.skill_tools import (
    list_available_skills,
    get_skill_details,
    invoke_skill,
    search_skills_for_task,
)


def create_scrna_pipeline_agent(model_name: Optional[str] = None):
    """
    Create the main scRNA-seq analysis pipeline agent.
    
    This orchestrator manages the full workflow:
    1. Data loading and QC check
    2. Collaborative annotation (multiple specialists)
    3. Integration with batch correction
    4. scANVI fine-tuning for refined results
    
    Args:
        model_name: Override model name from settings
        
    Returns:
        Configured pipeline agent
    """
    model = set_model(model_name=model_name, agent_name="pipeline")
    
    # Create annotation coordinator (manages specialist agents)
    annotation_coordinator = create_annotation_coordinator(model_name)
    
    # Main pipeline tools
    tools = [
        # Data loading
        load_data_for_pipeline,
        get_pipeline_status,
        # Annotation coordination (manages specialists internally)
        annotation_coordinator,
        # Integration
        prepare_for_integration,
        run_scvi,
        run_harmony,
        benchmark_integration,
        select_best_integration,
        # scANVI fine-tuning
        setup_scanvi_finetuning,
        run_scanvi_finetuning,
        # Post-processing
        compute_integrated_umap,
        plot_pipeline_results,
        save_pipeline_results,
    ]

    # Add Claude Code skill tools
    skill_tools = [
        list_available_skills,
        get_skill_details,
        invoke_skill,
        search_skills_for_task,
    ]
    tools.extend(skill_tools)
    
    prompt = """# Role
You are the scRNA-seq Analysis Pipeline Orchestrator with access to both traditional tools and Claude Code skills.
Your task is to guide data through a complete analysis workflow:
**Annotation → Integration → scANVI Fine-tuning**

# Your Capabilities

## Traditional Tools
- **Data Loading**: load_data_for_pipeline, get_pipeline_status
- **Annotation Coordinator**: Manages multiple annotation specialists
  - CellTypist Specialist: Deep learning annotation
  - ScType Specialist: Marker-based annotation
  - Literature RAG Specialist: PubMed-based annotation
- **Integration**: prepare_for_integration, run_scvi, run_harmony, benchmark_integration
- **Fine-tuning**: setup_scanvi_finetuning, run_scanvi_finetuning
- **Visualization**: compute_integrated_umap, plot_pipeline_results, save_pipeline_results

## Claude Code Skills
You have access to specialized Claude Code skills for advanced analysis:
- Use `list_available_skills()` to discover available skills
- Use `get_skill_details(skill_name)` to learn about a specific skill
- Use `invoke_skill(skill_name, args)` to use a skill
- Use `search_skills_for_task(description)` to find the right skill for a task

**Your Base Skills** (as pipeline orchestrator):
- **scanpy**: Single-cell analysis workflows
- **scvi-tools**: Probabilistic models for scRNA-seq
- **anndata**: Data structure operations
- **scientific-visualization**: Publication-quality figures
- **scientific-slides**: Presentation generation

**When to use Skills vs Traditional Tools:**
- **Traditional tools**: For core pipeline operations (loading, annotation, integration)
- **Skills**: For specialized analysis (pathway enrichment, advanced visualization, cell-cell communication)

# Recommended Workflow

## Phase 1: Data Initialization
1. Load data with load_data_for_pipeline(file_path, batch_key, species)
2. Check data status with get_pipeline_status()

## Phase 2: Collaborative Annotation
3. Invoke annotation_coordinator with tissue type and expected cell types
   - The coordinator will call all three specialists
   - They will discuss and create consensus annotations
   - Results are stored in adata.obs

## Phase 3: Integration
4. Prepare data with prepare_for_integration()
5. Run integration methods (run_scvi, run_harmony)
6. Benchmark with benchmark_integration()
7. Select best method with select_best_integration()

## Phase 4: scANVI Fine-tuning
8. Setup with setup_scanvi_finetuning() using consensus annotations
9. Run fine-tuning with run_scanvi_finetuning()
   - This refines both integration AND annotations
   - Predicts labels for low-confidence cells

## Phase 5: Advanced Analysis (Skills - Optional)
10. **NEW**: Use skills for specialized analysis when needed:
    - `invoke_skill("pathway-enrichment", "GO term analysis for cluster markers")`
    - `invoke_skill("cellchat", "infer cell-cell communication networks")`
    - `invoke_skill("scientific-visualization", "create publication-quality UMAP")`
    - Use `search_skills_for_task(description)` if unsure which skill to use

## Phase 6: Finalization
11. Compute UMAP with compute_integrated_umap()
12. Generate plots with plot_pipeline_results()
13. Save results with save_pipeline_results()

# Key Principles
- **Tool Selection**: Choose the right tool/skill based on the task
- **Collaboration**: Let specialists discuss and reach consensus
- **Validation**: Always validate annotations with markers
- **Iteration**: scANVI fine-tuning can improve initial annotations
- **Skill Discovery**: When unsure about available capabilities, use search_skills_for_task()
- **Documentation**: Report all steps and rationales

# Expected Outputs
- Annotated and integrated AnnData
- UMAP embeddings
- Consensus cell type annotations
- Trained scVI/scANVI models
- Summary plots
- (Optional) Advanced analysis results from skills

# Error Handling
- If a specialist fails, try alternative methods
- Flag persistent issues for user review
- Never proceed without valid annotations
- If a skill is not available, use traditional tools as fallback
"""
    
    agent = create_react_agent(model=model, tools=tools, prompt=prompt)
    
    return agent


async def run_pipeline(
    file_path: str,
    batch_key: str,
    tissue: str = "Immune system",
    species: str = "mouse",
    output_dir: str = "./pipeline_output",
    model_name: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Run the complete scRNA-seq analysis pipeline.
    
    Args:
        file_path: Path to input h5ad file
        batch_key: Column with batch/sample info
        tissue: Tissue type for annotation
        species: Species (mouse/human)
        output_dir: Directory for outputs
        model_name: LLM model to use
        
    Returns:
        Pipeline results dictionary
        
    Example:
        >>> import asyncio
        >>> result = asyncio.run(run_pipeline(
        ...     file_path="/path/to/data.h5ad",
        ...     batch_key="sample",
        ...     tissue="Spleen",
        ...     species="mouse",
        ... ))
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    agent = create_scrna_pipeline_agent(model_name)
    
    task = f"""
    Run the complete scRNA-seq analysis pipeline:
    
    1. Load data from {file_path}
       - Batch key: {batch_key}
       - Species: {species}
    
    2. Annotate cells (collaborate with specialists)
       - Tissue: {tissue}
       - Use CellTypist, ScType, and Literature RAG
       - Create consensus annotation
    
    3. Integrate with scVI and Harmony
       - Benchmark and select best method
    
    4. Fine-tune with scANVI
       - Use consensus annotations as labels
       - Refine predictions for low-confidence cells
    
    5. Save results
       - H5AD: {output_dir}/annotated_integrated.h5ad
       - Models: {output_dir}/models/
       - Plots: {output_dir}/plots/
    """
    
    config = {"configurable": {"species": species, "tissue": tissue}}
    result = await agent.ainvoke(
        {"messages": [HumanMessage(content=task)]}, 
        config=config
    )
    
    return result


# CLI entry point
if __name__ == "__main__":
    import argparse
    from dotenv import load_dotenv
    
    load_dotenv(override=True)
    
    parser = argparse.ArgumentParser(description="scRNA-seq Analysis Pipeline")
    parser.add_argument("--file", "-f", required=True, help="Path to h5ad file")
    parser.add_argument("--batch-key", "-b", required=True, help="Batch column name")
    parser.add_argument("--tissue", "-t", default="Immune system", help="Tissue type")
    parser.add_argument("--species", "-s", default="mouse", help="Species")
    parser.add_argument("--output", "-o", default="./output", help="Output directory")
    
    args = parser.parse_args()
    
    result = asyncio.run(run_pipeline(
        file_path=args.file,
        batch_key=args.batch_key,
        tissue=args.tissue,
        species=args.species,
        output_dir=args.output,
    ))
    
    print("Pipeline completed!")
    print(result)
