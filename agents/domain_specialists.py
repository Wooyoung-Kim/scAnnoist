"""Domain Specialist Agents

Expert agents for specific biological domains in scRNA-seq annotation.
These specialists provide domain-specific knowledge and interpretation.
"""

from typing import Annotated, Optional, Callable

from langchain_core.tools import tool
from langgraph.prebuilt import create_react_agent
from langchain_core.runnables.config import RunnableConfig
from langchain_core.messages import HumanMessage, AIMessage

from scrna_agent.utils import set_model
from scrna_agent.tools.annotation_tools import (
    validate_with_markers,
    search_pubmed_markers,
    get_cellmarker_database,
    annotate_with_literature_markers,
    compare_annotations,
)
from scrna_agent.tools.mcp_annotation_tools import (
    mcp_search_marker_literature,
    mcp_get_marker_evidence,
    mcp_find_related_markers,
    mcp_validate_cell_type_annotation,
    hybrid_search_markers,
)
from scrna_agent.skills.skill_tools import (
    list_available_skills,
    get_skill_details,
    invoke_skill,
)
from scrna_agent.skills.mcp_tools import (
    list_mcp_tools,
    get_mcp_tool_details,
    search_mcp_tools,
    get_mcp_tool_recommendations,
)
from scrna_agent.skills.skill_config import load_agent_skills


# ============================================================
# Hematopoiesis Specialist
# ============================================================

def create_hematopoiesis_specialist(model_name: Optional[str] = None) -> Callable:
    """Create a hematopoiesis and lineage hierarchy specialist agent."""
    model = set_model(model_name=model_name, agent_name="hematopoiesis_specialist")

    tools = [
        validate_with_markers,
        search_pubmed_markers,
        get_cellmarker_database,
        annotate_with_literature_markers,
        # MCP-based tools
        mcp_search_marker_literature,
        mcp_get_marker_evidence,
        mcp_find_related_markers,
        mcp_validate_cell_type_annotation,
        hybrid_search_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("hematopoiesis_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    # Add MCP discovery tools if enabled
    if agent_config.get("mcp_discovery_enabled", True):
        tools.extend([list_mcp_tools, get_mcp_tool_details, search_mcp_tools, get_mcp_tool_recommendations])

    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n\n**Claude Code Skills**: {', '.join(base_skills)}\n- Use pathway-enrichment for GO/KEGG analysis\n- Use cellchat/nichenet for cell communication\n- Use scientific-brainstorming for hypothesis generation" if base_skills else ""

    mcp_info = "\n\n**MCP Tools Available**: PubMed literature search and validation\n- Use mcp_search_marker_literature() to find papers about markers\n- Use mcp_validate_cell_type_annotation() to validate annotations\n- Use hybrid_search_markers() for robust literature search"

    prompt = f"""You are a hematopoiesis / stem-cell lineage expert for scRNA-seq annotation with access to Claude Code skills and MCP tools.{skill_info}{mcp_info}

Core responsibilities:
- Distinguish progenitors (HSC/MPP/CLP/GMP/MEP) vs mature immune lineages
- Identify transitional populations and proliferative cycling artifacts
- Provide marker sets (positive + negative markers) and decision rules, not just labels
- When uncertain, propose 2-3 candidate labels with discriminating markers to test

Key marker knowledge:
- HSC/MPP: Cd34, Kit, Sca1 (mouse), CD34, KIT, CD133 (human)
- CLP: Cd127, Flt3, Rag1
- GMP: Csf1r, Mpo, Elane
- MEP: Gata1, Klf1, Epor
- Mature B: Cd19, Ms4a1, Cd79a
- Mature T: Cd3e, Cd4, Cd8a
- Mature Myeloid: Itgam, Cd14, Fcgr3
- Cycling: Mki67, Top2a, Pcna

Output format:
1) Proposed label(s) + confidence level (high/medium/low)
2) Key positive markers supporting this label
3) Key negative markers (what should NOT be expressed)
4) Alternative explanations if confidence is not high
5) Next checks (specific plots/marker tests) to disambiguate

Special considerations:
- Flag potential cycling/proliferative states separately
- Distinguish emergency myelopoiesis from steady-state
- Identify activation-induced progenitor-like signatures in mature cells
- Always consider lineage hierarchy and developmental timing
"""

    agent = create_react_agent(model=model, tools=tools, prompt=prompt)

    @tool
    async def invoke_hematopoiesis_specialist(
        message: Annotated[str, "Question about hematopoietic cell types or lineages"],
        config: RunnableConfig,
    ) -> Annotated[str, "Expert interpretation of hematopoietic populations"]:
        """
        Invoke the hematopoiesis specialist for blood/BM/spleen annotation.
        Best for: HSC, progenitors, immune lineages, BM/PBMC/Spleen datasets.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="hematopoiesis_specialist")
            ]
        }

    return invoke_hematopoiesis_specialist


# ============================================================
# Infection/Inflammation Specialist
# ============================================================

def create_infection_inflammation_specialist(model_name: Optional[str] = None) -> Callable:
    """Create an infection and inflammation specialist agent."""
    model = set_model(model_name=model_name, agent_name="infection_inflammation_specialist")

    tools = [
        validate_with_markers,
        search_pubmed_markers,
        get_cellmarker_database,
        # MCP-based tools
        mcp_search_marker_literature,
        mcp_get_marker_evidence,
        mcp_find_related_markers,
        mcp_validate_cell_type_annotation,
        hybrid_search_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("infection_inflammation_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    # Add MCP discovery tools if enabled
    if agent_config.get("mcp_discovery_enabled", True):
        tools.extend([list_mcp_tools, get_mcp_tool_details, search_mcp_tools, get_mcp_tool_recommendations])

    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n\n**Claude Code Skills**: {', '.join(base_skills)}\n- Use pathway-enrichment for immune pathway analysis\n- Use cellchat for inflammation signaling\n- Use scientific-brainstorming for hypothesis generation" if base_skills else ""

    mcp_info = "\n\n**MCP Tools Available**: PubMed literature search and validation\n- Use mcp_search_marker_literature() for finding inflammation markers\n- Use mcp_validate_cell_type_annotation() for validation\n- Use hybrid_search_markers() for robust literature search"

    prompt = f"""You are an infection/inflammation specialist for scRNA-seq annotation with access to Claude Code skills and MCP tools.{skill_info}{mcp_info}

Core responsibilities:
- Separate true cell types from infection-induced states (ISG-high, TNF/NFκB, stress)
- Identify pDC vs ISG-high monocytes, activated T/NK vs bystander IFN states
- Flag dissociation stress signatures vs biologic inflammation
- Distinguish acute vs chronic inflammation signatures

Key signature knowledge:
- Type I IFN response: Isg15, Ifit1, Ifit3, Mx1, Oas1
- Type II IFN (IFNγ): Cxcl9, Cxcl10, Gbp2
- TNF/NFκB: Tnf, Il1b, Nfkbia, Cxcl1, Cxcl2
- Stress response: Jun, Fos, Hspa1a, Hspa1b, Dnajb1
- pDC: Siglech, Bst2, Tcf4, Irf7, Irf8
- Inflammatory monocytes: Ly6c2, Ccr2, Sell (mouse), LYC6, CCR2, SELL (human)
- Neutrophil activation: S100a8, S100a9, Cxcr2

Output format:
1) Cell type vs cell state distinction (e.g., "Monocyte [ISG-high state]")
2) Infection-driven modules likely present
3) Suggested gating markers & negative controls
4) Caveats (batch/timepoint/dissociation confounds)
5) Recommendations for validation

Special considerations:
- Dissociation stress (Jun/Fos/Hspa) vs true inflammation
- Batch effects in infection time-course experiments
- Bystander IFN activation vs cell-type-specific response
- Distinguish viral vs bacterial signatures when possible
"""

    agent = create_react_agent(model=model, tools=tools, prompt=prompt)

    @tool
    async def invoke_infection_inflammation_specialist(
        message: Annotated[str, "Question about infection or inflammatory states"],
        config: RunnableConfig,
    ) -> Annotated[str, "Expert interpretation of infection/inflammation signatures"]:
        """
        Invoke the infection/inflammation specialist.
        Best for: virus/bacteria challenge, ISG-high states, cytokine responses.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="infection_inflammation_specialist")
            ]
        }

    return invoke_infection_inflammation_specialist


# ============================================================
# Aging Specialist
# ============================================================

def create_aging_specialist(model_name: Optional[str] = None) -> Callable:
    """Create an aging biology specialist agent."""
    model = set_model(model_name=model_name, agent_name="aging_specialist")

    tools = [
        validate_with_markers,
        search_pubmed_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("aging_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n\n**Claude Code Skills**: {', '.join(base_skills)}\n- Use pathway-enrichment for aging pathway analysis\n- Use scientific-brainstorming for hypothesis generation" if base_skills else ""

    prompt = f"""You are an aging biology specialist for scRNA-seq annotation with access to Claude Code skills.{skill_info}

Core responsibilities:
- Interpret age-associated shifts (myeloid bias, reduced lymphopoiesis)
- Distinguish exhaustion/senescence-like programs from acute activation
- Provide age-aware annotation suggestions and pitfalls
- Flag inflammaging signatures and age-related changes

Key aging signatures:
- Myeloid bias: Increased GMP/monocyte progenitors
- Reduced lymphopoiesis: Decreased CLP, pre-B, pro-B
- Inflammaging: Il6, Il1b, Tnf, Ccl2, Cxcl1
- Senescence markers: Cdkn1a (p21), Cdkn2a (p16), Lmnb1
- T cell exhaustion: Pdcd1, Havcr2, Lag3, Tigit, Tox
- DNA damage response: Gadd45a, Mdm2, Trp53

Output format:
1) Age-associated interpretation for clusters (e.g., "likely enriched in aged samples")
2) Markers/programs to check (and what they imply)
3) Risks of mislabeling (activation vs aging programs)
4) Minimal validation checks (compare young vs aged directly)

Special considerations:
- Activation (acute) vs exhaustion/senescence (chronic/age)
- Batch effects when comparing young vs aged samples
- HSC aging signatures (epigenetic drift, clonal expansion)
- Age-associated B cells (ABCs) and other age-specific subsets
- Tissue-specific aging (e.g., thymic involution)
"""

    agent = create_react_agent(model=model, tools=tools, prompt=prompt)

    @tool
    async def invoke_aging_specialist(
        message: Annotated[str, "Question about age-related changes or populations"],
        config: RunnableConfig,
    ) -> Annotated[str, "Expert interpretation of aging-related signatures"]:
        """
        Invoke the aging specialist.
        Best for: young vs aged comparisons, inflammaging, age-skewed hematopoiesis.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="aging_specialist")
            ]
        }

    return invoke_aging_specialist


# ============================================================
# Cell State & Plasticity Specialist
# ============================================================

def create_cellstate_plasticity_specialist(model_name: Optional[str] = None) -> Callable:
    """Create a cell state and plasticity specialist agent."""
    model = set_model(model_name=model_name, agent_name="cellstate_plasticity_specialist")

    tools = [
        validate_with_markers,
        search_pubmed_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("cellstate_plasticity_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n\n**Claude Code Skills**: {', '.join(base_skills)}\n- Use pathway-enrichment for state-specific pathway analysis\n- Use nichenet for cell state communication\n- Use scientific-brainstorming for hypothesis generation" if base_skills else ""

    prompt = f"""You are a cell-state/plasticity expert for scRNA-seq annotation with access to Claude Code skills.{skill_info}

Core responsibilities:
- Label functional states (cycling, effector, memory, exhausted, stress) separately from cell type
- Recommend trajectory-aware checks (RNA velocity/latent time/module scores)
- Provide state-specific markers and how to separate them from batch effects
- Interpret transitional/intermediate states

Key state signatures:
- Cycling/Proliferation: Mki67, Top2a, Pcna, Cdk1, Ccnb1
- Effector (T/NK): Gzmb, Prf1, Ifng, Tnf
- Memory (T): Il7r, Ccr7, Tcf7, Lef1, Sell
- Exhaustion (T): Pdcd1, Havcr2, Lag3, Tigit, Tox
- Activation (early): Cd69, Cd44, Il2ra
- Quiescence: Ccr7, Tcf7, low ribosomal genes
- Stress: Jun, Fos, Hspa1a, Hspa1b, Dusp1

Output format:
1) Cell type label + state label (two-axis annotation, e.g., "CD8 T cell [Effector state]")
2) State-specific markers supporting
3) Competing explanations (e.g., batch vs biology)
4) Next best plots/tests (velocity, pseudotime, module scores)

Special considerations:
- Cycling can occur in any cell type - annotate separately
- Effector vs exhaustion: kinetics matter (acute vs chronic)
- Activation markers can be dissociation artifacts
- Use trajectory methods (Monocle, PAGA, Velocyto) for validation
- Consider EMT-like, epithelial-mesenchymal plasticity
"""

    agent = create_react_agent(model=model, tools=tools, prompt=prompt)

    @tool
    async def invoke_cellstate_plasticity_specialist(
        message: Annotated[str, "Question about cell states, activation, or transitions"],
        config: RunnableConfig,
    ) -> Annotated[str, "Expert interpretation of cell state and plasticity"]:
        """
        Invoke the cell state/plasticity specialist.
        Best for: activation states, trajectories, time-course, effector/memory/exhaustion.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="cellstate_plasticity_specialist")
            ]
        }

    return invoke_cellstate_plasticity_specialist


# ============================================================
# Epithelial Specialist
# ============================================================

def create_epithelial_specialist(model_name: Optional[str] = None) -> Callable:
    """Create an epithelial cell specialist agent."""
    model = set_model(model_name=model_name, agent_name="epithelial_specialist")

    tools = [
        validate_with_markers,
        search_pubmed_markers,
        get_cellmarker_database,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("epithelial_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n\n**Claude Code Skills**: {', '.join(base_skills)}\n- Use pathway-enrichment for epithelial pathway analysis\n- Use scientific-brainstorming for hypothesis generation" if base_skills else ""

    prompt = f"""You are an epithelial lineage expert for scRNA-seq annotation with access to Claude Code skills.{skill_info}

Core responsibilities:
- Identify epithelial compartments (basal, luminal, secretory, ciliated, etc. where applicable)
- Distinguish epithelial immune-like response (IFN, chemokines) from immune contamination
- Provide tissue-agnostic + tissue-adaptive marker logic
- Recognize tissue-specific epithelial subtypes

Key epithelial markers:
- Pan-epithelial: Epcam, Krt8, Krt18, Krt19
- Basal: Krt5, Krt14, Tp63, Itga6
- Luminal: Krt8, Krt18, Cd24a
- Secretory: Muc1, Muc5ac, Muc5b (tissue-dependent)
- Ciliated: Foxj1, Dynlrb2, Tuba1a
- Goblet: Muc2, Tff3, Agr2
- Club (lung): Scgb1a1, Cyp2f2
- AT1 (lung): Pdpn, Aqp5, Hopx
- AT2 (lung): Sftpc, Sftpa1, Sftpb

Output format:
1) Epithelial subtype label(s) with confidence
2) Positive/negative markers + rationale
3) Contamination vs epithelial stress judgment
4) Disambiguation checklist

Special considerations:
- Epithelial cells can express immune genes (e.g., ISG, chemokines) under stress
- Check EPCAM/KRT positivity to confirm epithelial identity
- Tissue-specific naming (e.g., AT2 for lung, enterocyte for gut)
- Dissociation stress can mask epithelial markers
- EMT-like states may have mixed markers
"""

    agent = create_react_agent(model=model, tools=tools, prompt=prompt)

    @tool
    async def invoke_epithelial_specialist(
        message: Annotated[str, "Question about epithelial cell types or states"],
        config: RunnableConfig,
    ) -> Annotated[str, "Expert interpretation of epithelial populations"]:
        """
        Invoke the epithelial specialist.
        Best for: epithelial lineages, EPCAM+ cells, tissue-specific epithelial subtypes.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="epithelial_specialist")
            ]
        }

    return invoke_epithelial_specialist


# ============================================================
# Organoid/Stemness Specialist
# ============================================================

def create_organoid_stemness_specialist(model_name: Optional[str] = None) -> Callable:
    """Create an organoid and stemness specialist agent."""
    model = set_model(model_name=model_name, agent_name="organoid_stemness_specialist")

    tools = [
        validate_with_markers,
        search_pubmed_markers,
    ]

    # Add skill tools based on configuration
    agent_config = load_agent_skills("organoid_stemness_specialist")
    if agent_config.get("discovery_enabled", False):
        tools.extend([list_available_skills, get_skill_details, invoke_skill])

    base_skills = agent_config.get("base_skills", [])
    skill_info = f"\n\n**Claude Code Skills**: {', '.join(base_skills)}\n- Use pathway-enrichment for stemness pathway analysis\n- Use scientific-brainstorming for hypothesis generation" if base_skills else ""

    prompt = f"""You are an organoid/stemness specialist for scRNA-seq annotation with access to Claude Code skills.{skill_info}

Core responsibilities:
- Interpret stemness vs hyperproliferation vs partial differentiation
- Flag culture-specific artifacts (dissociation stress, hypoxia, ECM dependence)
- Suggest annotation strategy for organoid-specific states and QC
- Recognize in vitro culture-induced changes

Key stemness/organoid signatures:
- Stemness: Lgr5, Ascl2, Olfm4 (intestinal), Prom1, Sox2, Nanog
- WNT activity: Axin2, Lgr5, Ascl2, Myc
- YAP/TAZ: Cyr61, Ctgf, Ankrd1
- Notch: Hes1, Hey1, Dll1
- Proliferation: Mki67, Top2a, Pcna
- Partial differentiation: Mix of stem + mature markers
- Culture stress: Jun, Fos, Hspa1a, Atf3
- Hypoxia: Vegfa, Slc2a1, Ldha, Hif1a

Output format:
1) Organoid cell identity + culture/state notes
2) Stemness/differentiation axis markers
3) Likely artifacts + how to test (compare to in vivo, check culture conditions)
4) Minimal recommended QC and controls

Special considerations:
- Culture duration affects differentiation state
- Compare to in vivo counterparts when possible
- Check for culture-induced stress vs biological stemness
- Organoid maturation stage (early vs late culture)
- 3D vs 2D culture artifacts
- ECM/matrix effects on gene expression
- Batch effects from passage number, media composition
"""

    agent = create_react_agent(model=model, tools=tools, prompt=prompt)

    @tool
    async def invoke_organoid_stemness_specialist(
        message: Annotated[str, "Question about organoid identity, stemness, or culture artifacts"],
        config: RunnableConfig,
    ) -> Annotated[str, "Expert interpretation of organoid/stemness states"]:
        """
        Invoke the organoid/stemness specialist.
        Best for: organoids, stem cells, WNT/YAP/Notch states, culture artifacts.
        """
        result = await agent.ainvoke(
            {"messages": [HumanMessage(content=message)]}, config=config
        )
        return {
            "messages": [
                AIMessage(content=result["messages"][-1].content, name="organoid_stemness_specialist")
            ]
        }

    return invoke_organoid_stemness_specialist
