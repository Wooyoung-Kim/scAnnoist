# MCP Integration Guide

## Overview

This guide explains the integration of MCP (Model Context Protocol) tools into the scRNA-seq analysis pipeline, enabling agents to access external data sources like PubMed for literature-based cell type annotation and validation.

## What is MCP?

MCP (Model Context Protocol) is a standard protocol that allows AI agents to interact with external tools and data sources. In this project, we use MCP to:

- **Search PubMed literature** for cell type markers
- **Validate annotations** with published research
- **Extract evidence** from full-text articles
- **Discover related markers** through literature networks

## Architecture

### Components

1. **MCP Registry** (`skills/mcp_registry.py`)
   - Discovers available MCP servers and tools
   - Categorizes tools by domain (biomedical, visualization, etc.)
   - Provides search and recommendation functions

2. **MCP Tool Wrappers** (`skills/mcp_tools.py`)
   - LangChain @tool wrappers for agent integration
   - Discovery tools (list, search, details)
   - Invocation helpers

3. **MCP Annotation Tools** (`tools/mcp_annotation_tools.py`)
   - High-level annotation tools using MCP
   - Hybrid fallback mechanisms
   - Domain-specific helpers

4. **Agent Integration** (`agents/domain_specialists.py`)
   - Specialists automatically get MCP tools
   - Configurable via `agent_skills.yaml`

## Available MCP Tools

### PubMed MCP Server

The PubMed MCP provides access to biomedical literature:

#### Core Search Tools

```python
# Search for articles
mcp_search_marker_literature(
    cell_type="CD8 T cells",
    marker_genes=["CD8A", "CD8B", "GZMB"],
    species="human",
    tissue="blood"
)

# Get detailed article metadata
mcp_get_marker_evidence(
    pmids=["35486828", "34577062"],
    extract_markers=True
)

# Get full-text articles from PMC
mcp_get_full_text(
    pmids=["35486828"]
)
```

#### Discovery Tools

```python
# Find related markers
mcp_find_related_markers(
    known_markers=["CD8A", "CD8B"],
    cell_type="CD8 T cells",
    species="human"
)

# Validate annotation
mcp_validate_cell_type_annotation(
    cell_type="CD8 T cells",
    observed_markers=["CD8A", "CD8B", "CD3D", "GZMB"],
    species="human",
    confidence_threshold=0.7
)
```

#### Hybrid Tools

```python
# Try MCP first, fallback to Bio.Entrez
hybrid_search_markers(
    cell_type="B cells",
    species="human",
    tissue="blood",
    use_mcp=True
)
```

## Configuration

### Agent Configuration

Configure MCP access in `configs/agent_skills.yaml`:

```yaml
# Global MCP settings
mcp:
  enabled: true
  auto_discover: true
  cache_duration_seconds: 600

  servers:
    pubmed:
      enabled: true
      category: "biomedical_literature"
      tools:
        - search_articles
        - get_article_metadata
        - get_full_text_article

# Agent-specific settings
agents:
  hematopoiesis_specialist:
    base_skills:
      - pathway-enrichment
      - cellchat

    # Enable MCP discovery
    mcp_discovery_enabled: true
    mcp_servers:
      - pubmed
```

### Environment Variables

No special environment variables are required for MCP tools. They use the MCP tools available in the Claude Code environment.

## Usage Examples

### Example 1: Literature-Based Marker Search

```python
from scrna_agent.tools import mcp_search_marker_literature

# Search for T cell markers
result = mcp_search_marker_literature(
    cell_type="CD8 T cells",
    marker_genes=["CD8A", "CD8B", "GZMB", "PRF1"],
    species="human",
    tissue="blood",
    max_results=10
)

print(result)
```

### Example 2: Validating Annotation

```python
from scrna_agent.tools import mcp_validate_cell_type_annotation

# Validate a cluster annotation
validation = mcp_validate_cell_type_annotation(
    cell_type="Regulatory T cells",
    observed_markers=["FOXP3", "IL2RA", "CTLA4", "CD25"],
    species="human",
    confidence_threshold=0.7
)

print(validation)
# Output includes:
# - Confidence score (0-1)
# - Papers found for each marker
# - PMIDs for evidence
# - Validation status (VALIDATED / LOW CONFIDENCE)
```

### Example 3: Finding Related Markers

```python
from scrna_agent.tools import mcp_find_related_markers

# Discover additional markers
related = mcp_find_related_markers(
    known_markers=["CD19", "MS4A1"],  # Known B cell markers
    cell_type="B cells",
    species="human"
)

print(related)
# Returns related papers and suggests additional markers
```

### Example 4: Using MCP Discovery in Agents

Agents can dynamically discover and use MCP tools:

```python
from scrna_agent.skills import list_mcp_tools, search_mcp_tools

# List all available MCP tools
tools = list_mcp_tools(server="pubmed")

# Search for tools by task
relevant_tools = search_mcp_tools("validate cell type markers")

# Get recommendations
recommendations = get_mcp_tool_recommendations(
    "I need to find papers about inflammation markers"
)
```

## Integration with Agents

### Domain Specialists

All domain specialist agents automatically have access to MCP tools:

```python
from scrna_agent.agents.domain_specialists import (
    create_hematopoiesis_specialist
)

# Create specialist with MCP tools
specialist = create_hematopoiesis_specialist()

# The specialist can now use:
# - mcp_search_marker_literature()
# - mcp_validate_cell_type_annotation()
# - mcp_find_related_markers()
# - hybrid_search_markers()
# Plus MCP discovery tools
```

### Agent Prompts

Agents are informed about MCP tools in their system prompts:

```
**MCP Tools Available**: PubMed literature search and validation
- Use mcp_search_marker_literature() to find papers about markers
- Use mcp_validate_cell_type_annotation() to validate annotations
- Use hybrid_search_markers() for robust literature search
```

## Adding New MCP Servers

To add support for new MCP servers (e.g., arXiv, bioRxiv):

### 1. Update MCP Registry

Edit `skills/mcp_registry.py` and add tool definitions:

```python
# In _get_known_mcp_tools()
arxiv_tools = [
    MCPToolInfo(
        name="search_papers",
        full_name="mcp__plugin_arxiv_ArXiv__search_papers",
        server="arxiv",
        category="preprints",
        description="Search arXiv preprints",
        when_to_use="When looking for computational biology preprints",
        example='search_papers(query="single-cell", category="q-bio")',
    ),
    # ... more tools
]
tools.extend(arxiv_tools)
```

### 2. Update Configuration

Edit `configs/agent_skills.yaml`:

```yaml
mcp:
  servers:
    arxiv:
      enabled: true
      category: "preprints"
      description: "arXiv preprint search"
      tools:
        - search_papers
        - get_paper_details
```

### 3. Create Tool Wrappers (Optional)

Create high-level wrappers in `tools/` for domain-specific use:

```python
# tools/arxiv_tools.py
@tool
def search_computational_biology_papers(
    query: str,
    max_results: int = 10
) -> str:
    """Search arXiv for computational biology papers."""
    # Use MCP arxiv tools
    pass
```

### 4. Add to Agents

Update specialist agent configurations:

```yaml
agents:
  computational_biology_specialist:
    mcp_servers:
      - pubmed
      - arxiv  # Add new server
```

## Best Practices

### 1. Use Hybrid Tools for Robustness

```python
# Preferred: Tries MCP first, falls back to Bio.Entrez
result = hybrid_search_markers(
    cell_type="T cells",
    species="human"
)

# Instead of: Direct MCP call (no fallback)
result = mcp_search_marker_literature(...)
```

### 2. Cache Literature Searches

The MCP annotation tools automatically cache results:

```python
# First call: searches PubMed
result1 = mcp_search_marker_literature(
    cell_type="B cells",
    species="human"
)

# Second call: uses cached result
result2 = mcp_search_marker_literature(
    cell_type="B cells",
    species="human"
)
```

### 3. Validate Critical Annotations

Always validate key cell type annotations with literature:

```python
# After automated annotation
for cluster in critical_clusters:
    validation = mcp_validate_cell_type_annotation(
        cell_type=cluster.annotation,
        observed_markers=cluster.top_markers,
        species="human"
    )

    if validation.confidence < 0.7:
        # Flag for manual review
        pass
```

### 4. Combine Multiple Evidence Sources

```python
# Use multiple MCP tools together
search_result = mcp_search_marker_literature(...)
evidence = mcp_get_marker_evidence(pmids=search_result.pmids)
related = mcp_find_related_markers(...)

# Combine with traditional tools
celltypist_result = celltypist_annotate(...)
comparison = compare_annotations(
    celltypist=celltypist_result,
    literature=evidence
)
```

## Troubleshooting

### MCP Tools Not Found

If MCP tools are not available:

1. Check if MCP server is installed in Claude Code environment
2. Verify `mcp.enabled: true` in `configs/agent_skills.yaml`
3. Use `list_mcp_servers()` to see available servers

### No Literature Results

If searches return no results:

1. Try broader search terms
2. Check species specification (human vs mouse)
3. Use `hybrid_search_markers()` to fallback to Bio.Entrez
4. Verify PubMed MCP is working: `list_mcp_tools(server="pubmed")`

### Low Validation Confidence

If annotation validation shows low confidence:

1. Check if marker genes are species-appropriate
2. Search for alternative cell type labels
3. Review additional markers from the cluster
4. Use `mcp_get_marker_evidence()` to read supporting papers

### Permission Errors

MCP tools run in the Claude Code sandbox. No special permissions needed.

## API Reference

### MCP Annotation Tools

#### `mcp_search_marker_literature()`
Search PubMed for literature about cell type markers.

**Parameters:**
- `cell_type` (str): Cell type to search
- `marker_genes` (Optional[List[str]]): Specific markers to validate
- `species` (str): "human" or "mouse"
- `tissue` (Optional[str]): Tissue context
- `max_results` (int): Maximum papers to retrieve

**Returns:** Formatted string with search results and PMIDs

#### `mcp_validate_cell_type_annotation()`
Validate a cell type annotation using literature evidence.

**Parameters:**
- `cell_type` (str): Proposed annotation
- `observed_markers` (List[str]): Markers observed in data
- `species` (str): Species context
- `confidence_threshold` (float): Minimum confidence (0-1)

**Returns:** Validation result with confidence score and evidence

#### `mcp_get_marker_evidence()`
Retrieve detailed article information and extract markers.

**Parameters:**
- `pmids` (List[str]): PubMed IDs
- `extract_markers` (bool): Extract marker genes from text

**Returns:** Article details with extracted marker information

#### `mcp_find_related_markers()`
Find related markers through literature networks.

**Parameters:**
- `known_markers` (List[str]): Known marker genes
- `cell_type` (str): Cell type context
- `species` (str): Species context

**Returns:** Related markers and supporting literature

#### `hybrid_search_markers()`
Hybrid search with MCP fallback to Bio.Entrez.

**Parameters:**
- `cell_type` (str): Cell type to search
- `species` (str): Species context
- `tissue` (Optional[str]): Tissue context
- `use_mcp` (bool): Try MCP first

**Returns:** Marker search results from available method

### MCP Discovery Tools

See `skills/mcp_tools.py` for complete API documentation.

## Future Enhancements

Planned improvements:

1. **Additional MCP Servers**
   - arXiv for computational methods
   - bioRxiv/medRxiv for preprints
   - CellMarker database integration

2. **Advanced Features**
   - Semantic search with embeddings
   - Citation graph analysis
   - Automated systematic reviews

3. **Performance**
   - Intelligent caching strategies
   - Parallel literature searches
   - Batch processing

4. **Validation**
   - Multi-source evidence aggregation
   - Confidence scoring improvements
   - Automated evidence quality assessment

## Contributing

To contribute MCP integrations:

1. Add new MCP server definitions to `mcp_registry.py`
2. Create domain-specific tool wrappers
3. Update agent configurations
4. Add tests and documentation
5. Submit pull request

## Support

For questions or issues:

1. Check this documentation
2. Review examples in `examples/`
3. Check agent logs for MCP tool usage
4. Open an issue on GitHub

---

**Version:** 1.0
**Last Updated:** 2026-01-16
**Authors:** scRNA-seq Agent Development Team
