"""
MCP-based Cell Type Annotation Tools

Enhanced annotation tools that leverage NCBI PubMed MCP for marker gene
validation and literature-based cell type identification.

This module provides tools that use the PubMed MCP server to:
- Search for marker genes in literature
- Validate cell type annotations with evidence
- Extract marker information from full-text articles
- Find related markers and cell types
"""

import logging
from typing import Annotated, Optional, List, Dict, Any
from langchain_core.tools import tool

logger = logging.getLogger(__name__)

# Global cache for literature search results
_mcp_literature_cache: Dict[str, Any] = {}


# ============================================================
# MCP-based Literature Search Tools
# ============================================================

@tool
def mcp_search_marker_literature(
    cell_type: Annotated[str, "Cell type to search (e.g., 'CD8 T cells', 'macrophages')"],
    marker_genes: Annotated[Optional[List[str]], "Optional list of marker genes to validate"] = None,
    species: Annotated[str, "Species: 'human' or 'mouse'"] = "human",
    tissue: Annotated[Optional[str], "Tissue context (e.g., 'lung', 'blood')"] = None,
    max_results: Annotated[int, "Maximum number of papers to retrieve"] = 10,
) -> Annotated[str, "Literature search results with marker information"]:
    """
    Search PubMed for literature about cell type markers using NCBI MCP.

    This tool uses the PubMed MCP search_articles function to find relevant
    papers about specific cell types and their marker genes. It constructs
    optimized queries for single-cell literature.

    Args:
        cell_type: The cell type to search for
        marker_genes: Optional list of specific marker genes to validate
        species: Species context (human or mouse)
        tissue: Optional tissue context to narrow search
        max_results: Maximum number of papers to retrieve

    Returns:
        Formatted string with search results, PMIDs, and extracted information

    Example:
        mcp_search_marker_literature(
            cell_type="CD8 T cells",
            marker_genes=["CD8A", "CD8B", "GZMB"],
            species="human",
            tissue="blood"
        )
    """
    from scrna_agent.tools import mcp_tools

    # Build optimized PubMed query
    query_parts = [f'"{cell_type}"[Title/Abstract]']

    # Add marker genes if provided
    if marker_genes:
        marker_query = " OR ".join([f'"{gene}"[Title/Abstract]' for gene in marker_genes[:5]])
        query_parts.append(f"({marker_query})")
    else:
        query_parts.append("(marker OR markers OR signature)[Title/Abstract]")

    # Add species
    query_parts.append(f"{species}[Title/Abstract]")

    # Add tissue if provided
    if tissue:
        query_parts.append(f'"{tissue}"[Title/Abstract]')

    # Focus on single-cell studies
    query_parts.append("(single cell OR scRNA-seq OR RNA-seq OR flow cytometry)[Title/Abstract]")

    # Filter for recent, relevant papers
    query_parts.append("2015:2025[PDAT]")

    query = " AND ".join(query_parts)

    try:
        # Use MCP PubMed search
        search_result = mcp_tools.call_mcp_tool(
            "mcp__plugin_pubmed_PubMed__search_articles",
            {
                "query": query,
                "max_results": max_results,
                "sort": "relevance",
                "date_from": "2015"
            }
        )

        if "error" in search_result:
            return f"❌ PubMed search failed: {search_result['error']}"

        # Cache results
        cache_key = f"{cell_type}_{species}_{tissue}"
        _mcp_literature_cache[cache_key] = search_result

        # Format results
        articles = search_result.get("articles", [])

        if not articles:
            return f"No papers found for '{cell_type}' in {species}. Try broader search terms."

        result = f"Literature Search Results for '{cell_type}' ({species})\n"
        result += "=" * 70 + "\n\n"
        result += f"Found {len(articles)} relevant papers\n\n"

        for i, article in enumerate(articles[:5], 1):
            result += f"{i}. {article.get('title', 'No title')}\n"
            result += f"   PMID: {article.get('pmid', 'N/A')}\n"
            result += f"   Year: {article.get('pub_date', 'N/A')[:4]}\n"
            result += f"   Journal: {article.get('journal', 'N/A')}\n"

            # Extract abstract snippet
            abstract = article.get('abstract', '')
            if abstract and len(abstract) > 200:
                result += f"   Abstract: {abstract[:200]}...\n"

            result += "\n"

        if len(articles) > 5:
            result += f"... and {len(articles) - 5} more papers.\n\n"

        result += f"💡 Use mcp_get_marker_evidence() to extract detailed marker information from these papers.\n"

        return result

    except Exception as e:
        logger.error(f"Error in MCP marker literature search: {e}", exc_info=True)
        return f"❌ Error searching literature: {str(e)}"


@tool
def mcp_get_marker_evidence(
    pmids: Annotated[List[str], "List of PubMed IDs to retrieve"],
    extract_markers: Annotated[bool, "Extract marker genes from text"] = True,
) -> Annotated[str, "Detailed article information with marker genes"]:
    """
    Retrieve detailed article information and extract marker genes using MCP.

    This tool fetches full metadata and optionally full-text articles from
    PubMed Central to extract marker gene information and validation evidence.

    Args:
        pmids: List of PubMed IDs to retrieve
        extract_markers: Whether to extract marker genes from text

    Returns:
        Formatted string with article details and extracted markers

    Example:
        mcp_get_marker_evidence(
            pmids=["35486828", "34577062"],
            extract_markers=True
        )
    """
    from scrna_agent.tools import mcp_tools
    import re

    try:
        # Get detailed metadata
        metadata_result = mcp_tools.call_mcp_tool(
            "mcp__plugin_pubmed_PubMed__get_article_metadata",
            {"pmids": pmids}
        )

        if "error" in metadata_result:
            return f"❌ Failed to retrieve metadata: {metadata_result['error']}"

        articles = metadata_result.get("articles", [])

        result = "Marker Evidence from Literature\n"
        result += "=" * 70 + "\n\n"

        for article in articles:
            pmid = article.get("pmid", "N/A")
            title = article.get("title", "No title")
            abstract = article.get("abstract", "")

            result += f"PMID: {pmid}\n"
            result += f"Title: {title}\n"
            result += f"DOI: {article.get('doi', 'N/A')}\n"
            result += f"Year: {article.get('pub_date', 'N/A')[:4]}\n"
            result += f"Journal: {article.get('journal', 'N/A')}\n\n"

            # Extract marker genes if requested
            if extract_markers and abstract:
                # Pattern for gene symbols (uppercase letters, numbers, hyphens)
                gene_pattern = r'\b[A-Z][A-Z0-9-]{2,10}\b'

                # Common marker-related keywords
                marker_keywords = [
                    "marker", "expressed", "expression", "positive for",
                    "characterized by", "identified by", "CD[0-9]+",
                    "signature", "specific"
                ]

                # Find sentences with marker keywords
                sentences = abstract.split('.')
                marker_sentences = []

                for sentence in sentences:
                    if any(keyword.lower() in sentence.lower() for keyword in marker_keywords):
                        marker_sentences.append(sentence.strip())

                if marker_sentences:
                    result += "Marker-related findings:\n"
                    for sentence in marker_sentences[:3]:
                        # Extract gene symbols from sentence
                        genes = re.findall(gene_pattern, sentence)
                        genes = [g for g in genes if len(g) > 2 and g not in ['THE', 'AND', 'FOR', 'WAS', 'WERE']]

                        if genes:
                            result += f"  • {sentence[:200]}\n"
                            result += f"    Genes mentioned: {', '.join(set(genes))}\n"
                    result += "\n"

            result += "-" * 70 + "\n\n"

        # Try to get full text from PMC if available
        result += "💡 Tip: Use mcp_get_full_text() to access complete articles from PMC.\n"

        return result

    except Exception as e:
        logger.error(f"Error retrieving marker evidence: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


@tool
def mcp_get_full_text(
    pmids: Annotated[List[str], "List of PubMed IDs"],
) -> Annotated[str, "Full-text articles from PubMed Central"]:
    """
    Retrieve full-text articles from PubMed Central using MCP.

    This tool attempts to retrieve complete article text from PMC, which
    provides more detailed marker information than abstracts alone.

    Args:
        pmids: List of PubMed IDs to retrieve full text for

    Returns:
        Full-text content or indication of availability

    Example:
        mcp_get_full_text(pmids=["35486828"])
    """
    from scrna_agent.tools import mcp_tools

    try:
        # First convert PMIDs to PMCIDs
        convert_result = mcp_tools.call_mcp_tool(
            "mcp__plugin_pubmed_PubMed__convert_article_ids",
            {
                "ids": pmids,
                "id_type": "pmid"
            }
        )

        if "error" in convert_result:
            return f"❌ ID conversion failed: {convert_result['error']}"

        conversions = convert_result.get("conversions", [])
        pmc_ids = []

        for conv in conversions:
            pmcid = conv.get("pmcid")
            if pmcid:
                pmc_ids.append(pmcid)

        if not pmc_ids:
            return "No full-text articles available in PubMed Central for these PMIDs."

        # Retrieve full text
        fulltext_result = mcp_tools.call_mcp_tool(
            "mcp__plugin_pubmed_PubMed__get_full_text_article",
            {"pmc_ids": pmc_ids}
        )

        if "error" in fulltext_result:
            return f"❌ Full-text retrieval failed: {fulltext_result['error']}"

        articles = fulltext_result.get("articles", [])

        result = f"Full-Text Articles Retrieved: {len(articles)}\n"
        result += "=" * 70 + "\n\n"

        for article in articles:
            result += f"PMCID: {article.get('pmcid', 'N/A')}\n"
            result += f"Title: {article.get('title', 'No title')}\n"
            result += f"DOI: {article.get('doi', 'N/A')}\n\n"

            # Show section headers available
            sections = article.get("sections", [])
            if sections:
                result += "Available sections:\n"
                for section in sections[:10]:
                    result += f"  • {section.get('title', 'Untitled')}\n"
                result += "\n"

            # Show methods/results snippets if available
            content = article.get("content", "")
            if len(content) > 500:
                result += f"Content preview:\n{content[:500]}...\n\n"

            result += "-" * 70 + "\n\n"

        return result

    except Exception as e:
        logger.error(f"Error retrieving full text: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


@tool
def mcp_find_related_markers(
    known_markers: Annotated[List[str], "Known marker genes"],
    cell_type: Annotated[str, "Cell type context"],
    species: Annotated[str, "Species: 'human' or 'mouse'"] = "human",
) -> Annotated[str, "Related marker genes and literature"]:
    """
    Find related marker genes by exploring connected literature using MCP.

    This tool searches for papers that mention the known markers and extracts
    additional marker genes that are co-mentioned, helping to discover
    comprehensive marker panels.

    Args:
        known_markers: List of known marker genes
        cell_type: Cell type context for search
        species: Species context

    Returns:
        Related markers and supporting literature

    Example:
        mcp_find_related_markers(
            known_markers=["CD8A", "CD8B"],
            cell_type="CD8 T cells",
            species="human"
        )
    """
    from scrna_agent.tools import mcp_tools

    try:
        # Search for papers mentioning known markers
        marker_query = " OR ".join([f'"{m}"[Title/Abstract]' for m in known_markers[:3]])
        query = f'("{cell_type}"[Title/Abstract]) AND ({marker_query}) AND {species}[Title/Abstract] AND (single cell OR scRNA-seq)[Title/Abstract]'

        search_result = mcp_tools.call_mcp_tool(
            "mcp__plugin_pubmed_PubMed__search_articles",
            {
                "query": query,
                "max_results": 10,
                "sort": "relevance"
            }
        )

        if "error" in search_result:
            return f"❌ Search failed: {search_result['error']}"

        articles = search_result.get("articles", [])

        if not articles:
            return f"No related papers found for markers: {', '.join(known_markers)}"

        # Get PMIDs and use find_related_articles
        pmids = [article.get("pmid") for article in articles[:3] if article.get("pmid")]

        if not pmids:
            return "No PMIDs available for related article search"

        related_result = mcp_tools.call_mcp_tool(
            "mcp__plugin_pubmed_PubMed__find_related_articles",
            {
                "pmids": pmids,
                "link_type": "pubmed_pubmed",
                "max_results": 15
            }
        )

        if "error" in related_result:
            return f"❌ Related search failed: {related_result['error']}"

        # Format results
        result = f"Related Markers for {cell_type}\n"
        result += "=" * 70 + "\n\n"
        result += f"Known markers: {', '.join(known_markers)}\n\n"

        related_pmids = related_result.get("related_pmids", [])
        result += f"Found {len(related_pmids)} related papers\n\n"

        # Get metadata for related papers
        if related_pmids[:5]:
            metadata_result = mcp_tools.call_mcp_tool(
                "mcp__plugin_pubmed_PubMed__get_article_metadata",
                {"pmids": related_pmids[:5]}
            )

            if "articles" in metadata_result:
                result += "Top related papers:\n"
                for i, article in enumerate(metadata_result["articles"], 1):
                    result += f"{i}. {article.get('title', 'No title')}\n"
                    result += f"   PMID: {article.get('pmid', 'N/A')}\n"
                    result += f"   Year: {article.get('pub_date', 'N/A')[:4]}\n\n"

        result += "\n💡 Use mcp_get_marker_evidence() with these PMIDs to extract additional markers.\n"

        return result

    except Exception as e:
        logger.error(f"Error finding related markers: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


@tool
def mcp_validate_cell_type_annotation(
    cell_type: Annotated[str, "Proposed cell type annotation"],
    observed_markers: Annotated[List[str], "Marker genes observed in the cluster"],
    species: Annotated[str, "Species: 'human' or 'mouse'"] = "human",
    confidence_threshold: Annotated[float, "Minimum confidence for validation (0-1)"] = 0.7,
) -> Annotated[str, "Validation result with literature support"]:
    """
    Validate a cell type annotation using literature evidence from PubMed MCP.

    This tool searches literature to confirm that the observed markers are
    consistent with the proposed cell type annotation, providing confidence
    scores and evidence citations.

    Args:
        cell_type: The proposed cell type label
        observed_markers: Marker genes observed in the data
        species: Species context
        confidence_threshold: Minimum confidence for validation

    Returns:
        Validation result with confidence score and evidence

    Example:
        mcp_validate_cell_type_annotation(
            cell_type="CD8 T cells",
            observed_markers=["CD8A", "CD8B", "CD3D", "GZMB"],
            species="human"
        )
    """
    from scrna_agent.tools import mcp_tools

    try:
        # Search for papers that mention both cell type and markers
        validation_results = {}

        for marker in observed_markers[:5]:  # Check top 5 markers
            query = f'("{cell_type}"[Title/Abstract]) AND ("{marker}"[Title/Abstract]) AND {species}[Title/Abstract] AND (marker OR expression)[Title/Abstract]'

            search_result = mcp_tools.call_mcp_tool(
                "mcp__plugin_pubmed_PubMed__search_articles",
                {
                    "query": query,
                    "max_results": 5,
                    "sort": "relevance"
                }
            )

            if "articles" in search_result:
                articles = search_result["articles"]
                validation_results[marker] = {
                    "papers_found": len(articles),
                    "pmids": [a.get("pmid") for a in articles[:3]]
                }

        # Calculate confidence score
        total_markers = len(observed_markers[:5])
        validated_markers = sum(1 for m, r in validation_results.items() if r["papers_found"] > 0)
        confidence = validated_markers / total_markers if total_markers > 0 else 0

        # Format results
        result = f"Cell Type Annotation Validation\n"
        result += "=" * 70 + "\n\n"
        result += f"Proposed: {cell_type} ({species})\n"
        result += f"Markers tested: {', '.join(observed_markers[:5])}\n"
        result += f"Confidence: {confidence:.2%}\n"
        result += f"Status: {'✓ VALIDATED' if confidence >= confidence_threshold else '⚠ LOW CONFIDENCE'}\n\n"

        result += "Marker Validation Details:\n"
        result += "-" * 70 + "\n"

        for marker, info in validation_results.items():
            papers = info["papers_found"]
            status = "✓" if papers > 0 else "✗"
            result += f"{status} {marker}: {papers} supporting papers\n"

            if info["pmids"]:
                result += f"   PMIDs: {', '.join(info['pmids'])}\n"

        result += "\n"

        if confidence < confidence_threshold:
            result += "⚠ WARNING: Low confidence. Consider:\n"
            result += "  • Checking if marker genes are species-appropriate\n"
            result += "  • Searching for alternative cell type labels\n"
            result += "  • Reviewing additional markers from the cluster\n\n"

        result += "💡 Use mcp_get_marker_evidence() to read supporting papers.\n"

        return result

    except Exception as e:
        logger.error(f"Error validating annotation: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


# ============================================================
# Hybrid Tools (MCP + Traditional)
# ============================================================

@tool
def hybrid_search_markers(
    cell_type: Annotated[str, "Cell type to search"],
    species: Annotated[str, "Species: 'human' or 'mouse'"] = "human",
    tissue: Annotated[Optional[str], "Tissue context"] = None,
    use_mcp: Annotated[bool, "Use MCP if available, fallback to Bio.Entrez"] = True,
) -> Annotated[str, "Marker search results"]:
    """
    Hybrid marker search that tries MCP first, then falls back to Bio.Entrez.

    This tool provides robustness by attempting to use the MCP PubMed search
    but falling back to traditional Bio.Entrez if MCP is unavailable.

    Args:
        cell_type: Cell type to search for
        species: Species context
        tissue: Optional tissue context
        use_mcp: Whether to try MCP first

    Returns:
        Marker search results from available method
    """
    if use_mcp:
        try:
            # Try MCP first
            result = mcp_search_marker_literature(
                cell_type=cell_type,
                species=species,
                tissue=tissue
            )

            if not result.startswith("❌"):
                return result
        except Exception as e:
            logger.warning(f"MCP search failed, falling back to Bio.Entrez: {e}")

    # Fallback to traditional method
    try:
        from scrna_agent.tools.annotation_tools import search_pubmed_markers

        fallback_result = search_pubmed_markers(
            cell_type=cell_type,
            species=species,
            tissue=tissue
        )

        return f"[Using Bio.Entrez fallback]\n\n{fallback_result}"

    except Exception as e:
        return f"❌ Both MCP and Bio.Entrez searches failed: {str(e)}"


# ============================================================
# Export all tools
# ============================================================

__all__ = [
    "mcp_search_marker_literature",
    "mcp_get_marker_evidence",
    "mcp_get_full_text",
    "mcp_find_related_markers",
    "mcp_validate_cell_type_annotation",
    "hybrid_search_markers",
]
