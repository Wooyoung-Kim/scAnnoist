# Example: Basic Pipeline Usage
# Run the complete scRNA-seq analysis pipeline

import asyncio
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Import the pipeline
from scrna_agent import run_pipeline, create_scrna_pipeline_agent
from langchain_core.messages import HumanMessage


async def basic_example():
    """
    Basic example: Run the full pipeline with default settings.
    """
    result = await run_pipeline(
        file_path="/path/to/your/data.h5ad",
        batch_key="sample",  # Column containing batch/sample info
        tissue="Spleen",     # Tissue type for annotation
        species="mouse",     # Species
        output_dir="./output",
    )
    
    print("Pipeline completed!")
    return result


async def custom_example():
    """
    Custom example: Use the agent directly for more control.
    """
    # Create the pipeline agent
    agent = create_scrna_pipeline_agent(model_name="gpt-4o")
    
    # Define a custom task
    task = """
    I have mouse spleen B cell data that needs analysis.
    
    1. Load data from /path/to/b_cells.h5ad
       - Batch key is 'sample'
       - This is mouse data
    
    2. For annotation:
       - Focus on B cell subtypes: ImmB, NB1, NB2, MZB, Plasma cells
       - Use CellTypist with Immune_All_High model
       - Also use Literature RAG to find markers for ferret B cells
       - Create consensus annotation
    
    3. Integrate using scVI (primary) and Harmony (comparison)
    
    4. Fine-tune with scANVI to refine annotations
    
    5. Save everything to ./b_cell_analysis/
    """
    
    # Run the agent
    config = {"configurable": {"species": "mouse", "tissue": "Spleen"}}
    result = await agent.ainvoke(
        {"messages": [HumanMessage(content=task)]},
        config=config,
    )
    
    # Print the final response
    print(result["messages"][-1].content)
    return result


async def annotation_only_example():
    """
    Example: Use only the annotation coordinator for cell type annotation.
    """
    from scrna_agent import create_annotation_coordinator
    from scrna_agent.tools import load_adata_for_annotation
    
    # Load data first
    load_adata_for_annotation("/path/to/data.h5ad")
    
    # Create the annotation coordinator
    coordinator = create_annotation_coordinator()
    
    # Run collaborative annotation
    task = """
    Annotate mouse spleen immune cells.
    
    Expected cell types:
    - B cells (ImmB, Naive B, Memory B, Plasma)
    - T cells (CD4, CD8, Tfh)
    - Myeloid cells (Macrophages, DCs)
    
    1. Ask CellTypist to annotate with Immune_All_High model
    2. Ask ScType to annotate with tissue='Spleen'
    3. Ask Literature RAG to search for markers of these cell types
    4. Compare all three methods
    5. Create consensus annotation
    6. Report any disagreements and how they were resolved
    """
    
    config = {"configurable": {"species": "mouse"}}
    result = await coordinator.ainvoke(
        {"message": task},
        config=config,
    )
    
    return result


if __name__ == "__main__":
    # Run the basic example
    asyncio.run(basic_example())
