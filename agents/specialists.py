# Specialist Agents
# Individual expert agents for different annotation methods

from typing import Annotated, Optional, Callable

from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables.config import RunnableConfig
from langchain_core.messages import HumanMessage, AIMessage

from scrna_agent.utils import set_model
from scrna_agent.tools.annotation_tools import (
    celltypist_list_models,
    celltypist_download_model,
    celltypist_annotate,
    celltypist_get_probabilities,
    sctype_annotate,
    validate_with_markers,
    search_pubmed_markers,
    get_cellmarker_database,
    annotate_with_literature_markers,
    combine_pubmed_and_database_markers,
)
from scrna_agent.skills.skill_tools import (
    list_available_skills,
    get_skill_details,
    invoke_skill,
)
from scrna_agent.skills.skill_config import load_agent_skills


# ============================================================
# CellTypist Specialist
# ============================================================

def create_celltypist_specialist(model_name: Optional[str] = None) -> Callable:
    """Create a CellTypist annotation specialist agent."""
    model = set_model(model_name=model_name, agent_name="celltypist_specialist")

    # Traditional CellTypist tools
    tools = [
        celltypist_list_models,
        celltypist_download_model,
        celltypist_annotate,
        celltypist_get_probabilities,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("celltypist_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    # Prepare skill info for prompt
    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n**Your Base Skills**: {', '.join(base_skills)}" if base_skills else ""

    prompt = f"""You are a CellTypist specialist with access to traditional tools and Claude Code skills.
Your expertise is using pre-trained deep learning models for automatic cell type annotation.

# Your Capabilities

## Traditional Tools
- celltypist_list_models: List available pre-trained models
- celltypist_download_model: Download models for use
- celltypist_annotate: Run CellTypist annotation with majority voting
- celltypist_get_probabilities: Get confidence scores

## Claude Code Skills{skill_info}
- Use skills for complementary analysis when traditional tools need support
- Skills available: scanpy (single-cell analysis), celltype-annotation (alternative methods)

# Key Capabilities
- Select appropriate pre-trained models based on tissue type
- Run CellTypist annotation with majority voting
- Provide confidence scores for predictions
- Recommend models (e.g., Immune_All_High.pkl for immune cells)

# When asked to annotate
1. Check available models with celltypist_list_models
2. Download the appropriate model
3. Run annotation with majority_voting=True
4. Report predictions with confidence statistics

# Always report
- Cell type distribution
- Confidence score statistics
- Cells with low confidence (<0.5)
- Your confidence level in the overall annotation

# When to use skills
- If you need additional analysis beyond CellTypist
- For alternative annotation methods (use celltype-annotation skill)
- For advanced data exploration (use scanpy skill)
"""
    
    agent = create_react_agent(model=model, tools=tools, prompt=prompt)
    
    @tool
    async def invoke_celltypist_specialist(
        message: Annotated[str, "Task description for CellTypist annotation"],
        config: RunnableConfig,
    ) -> Annotated[str, "CellTypist annotation results"]:
        """
        Invoke the CellTypist specialist for deep learning-based annotation.
        Best for: immune cells, pre-defined cell types, high-throughput annotation.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="celltypist_specialist")
            ]
        }
    
    return invoke_celltypist_specialist


# ============================================================
# ScType Specialist
# ============================================================

def create_sctype_specialist(model_name: Optional[str] = None) -> Callable:
    """Create a ScType annotation specialist agent."""
    model = set_model(model_name=model_name, agent_name="sctype_specialist")

    # Traditional ScType tools
    tools = [
        sctype_annotate,
        validate_with_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("sctype_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    # Prepare skill info for prompt
    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n**Your Base Skills**: {', '.join(base_skills)}" if base_skills else ""

    prompt = f"""You are a ScType specialist with access to traditional tools and Claude Code skills.
Your expertise is marker gene-based cell type annotation.

# Your Capabilities

## Traditional Tools
- sctype_annotate: Annotate clusters using curated marker databases
- validate_with_markers: Validate annotations with known marker expression

## Claude Code Skills{skill_info}
- Use skills for complementary analysis when traditional tools need support
- Skills available: scanpy (single-cell analysis), celltype-annotation (alternative methods)

# Key Capabilities
- Annotate clusters using curated marker databases
- Validate annotations with known marker expression
- Quick annotation for well-characterized tissues

# When asked to annotate
1. Determine the appropriate tissue type
2. Run ScType annotation on clusters
3. Validate with key marker genes
4. Report cluster-level annotations with scores

# Always report
- Cluster-to-cell-type mapping
- Marker gene expression validation
- Any clusters with ambiguous or low scores
- Your assessment of annotation reliability

# When to use skills
- If you need additional analysis beyond ScType
- For alternative annotation methods (use celltype-annotation skill)
- For advanced marker analysis (use scanpy skill)
"""
    
    agent = create_react_agent(model=model, tools=tools, prompt=prompt)
    
    @tool
    async def invoke_sctype_specialist(
        message: Annotated[str, "Task description for ScType annotation"],
        config: RunnableConfig,
    ) -> Annotated[str, "ScType annotation results"]:
        """
        Invoke the ScType specialist for marker-based annotation.
        Best for: quick annotation, known markers, cluster-level annotation.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="sctype_specialist")
            ]
        }
    
    return invoke_sctype_specialist


# ============================================================
# Literature RAG Specialist
# ============================================================

def create_literature_rag_specialist(model_name: Optional[str] = None) -> Callable:
    """Create a Literature RAG annotation specialist agent."""
    model = set_model(model_name=model_name, agent_name="literature_rag_specialist")

    # Traditional Literature RAG tools
    tools = [
        search_pubmed_markers,
        get_cellmarker_database,
        annotate_with_literature_markers,
        combine_pubmed_and_database_markers,
        validate_with_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("literature_rag_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    # Prepare skill info for prompt
    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n**Your Base Skills**: {', '.join(base_skills)}" if base_skills else ""

    prompt = f"""You are a Literature RAG Specialist with access to traditional tools and Claude Code skills.
Your expertise is using scientific literature (PubMed) to find and validate marker genes for cell type annotation.

# Your Capabilities

## Traditional Tools
- search_pubmed_markers: Search PubMed for marker genes
- get_cellmarker_database: Access curated marker database
- annotate_with_literature_markers: Annotate using literature markers
- combine_pubmed_and_database_markers: Combine multiple marker sources
- validate_with_markers: Validate annotations with expression data

## Claude Code Skills{skill_info}
- Skills available: biopython (PubMed/NCBI queries), scanpy (single-cell analysis)
- Use biopython skill for advanced PubMed searches and NCBI database access
- Use scanpy skill for additional single-cell analysis

# Key Capabilities
- Search PubMed for marker genes of specific cell types
- Access the CellMarker database for curated markers
- Combine literature and database markers for robust annotation
- Validate annotations with expression data

# When asked about markers
1. First check CellMarker database with get_cellmarker_database()
2. Search PubMed for additional markers with search_pubmed_markers()
3. Combine markers with combine_pubmed_and_database_markers()
4. Annotate clusters with annotate_with_literature_markers()
5. Validate with validate_with_markers()

# Always report
- Source of markers (database vs literature)
- Number of papers found in PubMed
- Marker coverage in the data
- Confidence in the annotation

# Best practices
- Use species-specific marker names (Cd19 for mouse, CD19 for human)
- Cross-validate database markers with literature
- Report when literature markers differ from database
- Flag novel markers for expert review

# When to use skills
- For advanced PubMed/NCBI queries (use biopython skill)
- For additional single-cell analysis (use scanpy skill)
- When traditional tools need complementary analysis
"""
    
    agent = create_react_agent(model=model, tools=tools, prompt=prompt)
    
    @tool
    async def invoke_literature_rag_specialist(
        message: Annotated[str, "Task description for literature-based annotation"],
        config: RunnableConfig,
    ) -> Annotated[str, "Literature RAG annotation results"]:
        """
        Invoke the Literature RAG specialist for PubMed-based marker annotation.
        Best for: novel cell types, custom markers, validating other methods.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="literature_rag_specialist")
            ]
        }
    
    return invoke_literature_rag_specialist
