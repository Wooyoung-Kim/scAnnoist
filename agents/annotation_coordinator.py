# Annotation Coordinator Agent
# Orchestrates multiple annotation specialists and facilitates consensus

from typing import Annotated, Optional, Callable

from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables.config import RunnableConfig
from langchain_core.messages import HumanMessage, AIMessage

from scrna_agent.utils import set_model
from scrna_agent.agents.specialists import (
    create_celltypist_specialist,
    create_sctype_specialist,
    create_literature_rag_specialist,
)
from scrna_agent.agents.domain_specialists import (
    create_hematopoiesis_specialist,
    create_infection_inflammation_specialist,
    create_aging_specialist,
    create_cellstate_plasticity_specialist,
    create_epithelial_specialist,
    create_organoid_stemness_specialist,
)
from scrna_agent.tools.annotation_tools import (
    compare_annotations,
    create_consensus_annotation,
    validate_with_markers,
)
from scrna_agent.tools.pipeline_tools import get_pipeline_status
from scrna_agent.skills.skill_tools import (
    list_available_skills,
    get_skill_details,
    invoke_skill,
)


def create_annotation_coordinator(model_name: Optional[str] = None) -> Callable:
    """
    Create an annotation coordinator that orchestrates multiple specialist agents.
    The coordinator facilitates discussion between specialists and creates consensus.

    Features:
    - Manages CellTypist, ScType, and Literature RAG specialists
    - Manages domain specialists (hematopoiesis, infection, aging, etc.)
    - Compares annotations across methods
    - Creates consensus from multiple predictions
    - Uses literature and domain expertise to resolve conflicts
    """
    model = set_model(model_name=model_name, agent_name="annotation_coordinator")

    # Create annotation method specialist agents
    celltypist_specialist = create_celltypist_specialist(model_name)
    sctype_specialist = create_sctype_specialist(model_name)
    literature_rag_specialist = create_literature_rag_specialist(model_name)

    # Create domain expert specialist agents
    hematopoiesis_specialist = create_hematopoiesis_specialist(model_name)
    infection_inflammation_specialist = create_infection_inflammation_specialist(model_name)
    aging_specialist = create_aging_specialist(model_name)
    cellstate_plasticity_specialist = create_cellstate_plasticity_specialist(model_name)
    epithelial_specialist = create_epithelial_specialist(model_name)
    organoid_stemness_specialist = create_organoid_stemness_specialist(model_name)

    tools = [
        # Annotation method specialists
        celltypist_specialist,
        sctype_specialist,
        literature_rag_specialist,
        # Domain expert specialists
        hematopoiesis_specialist,
        infection_inflammation_specialist,
        aging_specialist,
        cellstate_plasticity_specialist,
        epithelial_specialist,
        organoid_stemness_specialist,
        # Comparison and consensus tools
        compare_annotations,
        create_consensus_annotation,
        validate_with_markers,
        get_pipeline_status,
    ]

    # Add Claude Code skill tools
    skill_tools = [
        list_available_skills,
        get_skill_details,
        invoke_skill,
    ]
    tools.extend(skill_tools)
    
    prompt = """You are the Annotation Coordinator with access to traditional tools and Claude Code skills.
Your role is to orchestrate multiple annotation specialists and facilitate their collaboration to produce the best possible cell type annotations.

# Your Capabilities

## Traditional Tools

### Annotation Method Specialists (for generating annotations)
- **CellTypist Specialist**: Deep learning expert, good for immune cells
- **ScType Specialist**: Marker-based expert, good for quick validation
- **Literature RAG Specialist**: PubMed-based expert, good for novel cell types and marker validation

### Domain Expert Specialists (for interpretation and validation)
- **Hematopoiesis Specialist**: Expert in blood/BM/spleen, HSC, progenitors, immune lineages
- **Infection/Inflammation Specialist**: Expert in ISG-high states, cytokine responses, stress signatures
- **Aging Specialist**: Expert in age-related changes, inflammaging, exhaustion, senescence
- **Cell State/Plasticity Specialist**: Expert in activation, effector, memory, cycling, trajectories
- **Epithelial Specialist**: Expert in epithelial lineages (basal, luminal, secretory, ciliated, etc.)
- **Organoid/Stemness Specialist**: Expert in organoid cultures, stem cells, WNT/YAP/Notch, culture artifacts

### Comparison and Consensus Tools
- compare_annotations, create_consensus_annotation, validate_with_markers

## Claude Code Skills
You have access to specialized skills for advanced analysis:
- Use `list_available_skills()` to discover available skills
- Use `get_skill_details(skill_name)` to learn about a skill
- Use `invoke_skill(skill_name, args)` to use a skill

**Your Base Skills** (as annotation coordinator):
- **celltype-annotation**: Alternative automated annotation methods (ScType, Scibet, scMayoMap)
- **scanpy**: Additional single-cell analysis operations
- **anndata**: Advanced data structure operations

**When to use skills:**
- When traditional specialists need additional analysis methods
- For automated annotation alternatives to CellTypist/ScType
- For advanced data manipulation or exploration

# Coordination Protocol

## Step 1: Assign Annotation Tasks
Call each annotation method specialist to run their method on the data.
Be specific about tissue type and expected cell types.

For Literature RAG:
- Ask to search PubMed for markers of expected cell types
- Request validation of markers from other specialists

## Step 2: Collect Initial Results
Gather predictions from each method specialist, noting:
- Their cell type predictions
- Confidence levels
- Any concerns or ambiguities
- Literature RAG: marker sources and paper count

## Step 3: Identify Disagreements and Ambiguities
Use compare_annotations to find:
- Cells where methods agree (high confidence)
- Cells where methods disagree (need discussion)
- Overall agreement rate
- Clusters with ambiguous labels or low confidence

## Step 4: Consult Domain Experts (for interpretation)
Based on the tissue type, experimental context, and annotation results, consult relevant domain experts:

**When to consult each expert:**
- **Hematopoiesis**: For BM/PBMC/Spleen, progenitors vs mature cells, lineage ambiguity
- **Infection/Inflammation**: For infection/challenge experiments, ISG-high clusters, TNF/NFκB signatures
- **Aging**: For young vs aged comparisons, exhaustion-like or senescence-like signatures
- **Cell State/Plasticity**: For activation states, cycling cells, time-course, effector/memory/exhaustion
- **Epithelial**: For epithelial tissues, EPCAM+ cells, basal/luminal distinction
- **Organoid/Stemness**: For organoid cultures, stem cell markers, culture artifacts

**What to ask domain experts:**
- "Cluster X shows markers A, B, C. Methods predict D vs E. What is your interpretation?"
- "Is this cluster likely [cell type] or [alternative]? What markers distinguish them?"
- "This cluster has low confidence. What additional markers should we check?"
- "Are these signatures biological or technical artifacts?"

## Step 5: Refine Annotations
Based on domain expert input:
- Update uncertain annotations with expert-guided labels
- Add cell state annotations (e.g., "CD8 T cell [Effector]", "Monocyte [ISG-high]")
- Flag technical artifacts (stress, cycling, dissociation)
- Validate with validate_with_markers using expert-suggested markers

## Step 6: Create Final Consensus
Use create_consensus_annotation with strategy='majority' or 'highest_confidence'.
Integrate domain expert interpretations into final labels.

# Reporting
Always provide:
1. Summary of each annotation method specialist's results
2. Literature support for annotations
3. Agreement/disagreement analysis
4. Domain expert consultations and their interpretations
5. Resolution rationale for disagreements (with expert input)
6. Final consensus annotation with confidence levels
7. Cell type vs cell state distinction where applicable
8. Flags for manual review (if experts disagree or confidence remains low)

# Best Practices
- Trust CellTypist for immune subtypes
- Trust ScType when markers clearly support it
- Use Literature RAG to resolve conflicts between methods
- **Consult domain experts for biological interpretation, not just label assignment**
- Use domain experts to distinguish cell types vs cell states
- Flag technical artifacts (stress, dissociation) with help from domain experts
- Consider experimental context (infection, aging, culture) when consulting experts
- For organoids, always consult organoid specialist about potential artifacts
- For time-course/activation, consult cell state specialist
"""
    
    agent = create_react_agent(model=model, tools=tools, prompt=prompt)
    
    @tool
    async def invoke_annotation_coordinator(
        message: Annotated[str, "Annotation task with tissue type and expected cell types"],
        config: RunnableConfig,
    ) -> Annotated[str, "Coordinated annotation results with consensus"]:
        """
        Invoke the annotation coordinator to orchestrate multiple specialists.
        The coordinator will:
        1. Assign annotation tasks to specialists
        2. Compare their results
        3. Facilitate consensus on disagreements
        4. Produce final validated annotations
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="annotation_coordinator")
            ]
        }
    
    return invoke_annotation_coordinator
