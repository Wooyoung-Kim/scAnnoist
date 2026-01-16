"""
MCP Tool Registry

Discovers, registers, and manages MCP (Model Context Protocol) tools
available in the Claude Code environment.

This module provides infrastructure for:
- Automatic discovery of installed MCP servers and their tools
- Categorization of MCP tools by domain (biomedical, database, visualization, etc.)
- Search and recommendation of MCP tools for specific tasks
- Integration with the skill system for unified tool access
"""

import logging
import inspect
from typing import Dict, List, Any, Optional, Set
from pathlib import Path
from dataclasses import dataclass, field

logger = logging.getLogger(__name__)


@dataclass
class MCPToolInfo:
    """Information about an MCP tool."""
    name: str
    full_name: str  # e.g., "mcp__plugin_pubmed_PubMed__search_articles"
    server: str  # e.g., "pubmed"
    category: str  # e.g., "biomedical_literature"
    description: str
    parameters: Dict[str, Any] = field(default_factory=dict)
    when_to_use: str = ""
    example: str = ""
    related_tools: List[str] = field(default_factory=list)


class MCPRegistry:
    """
    Registry for MCP tools with automatic discovery and categorization.
    """

    def __init__(self):
        self._tools: Dict[str, MCPToolInfo] = {}
        self._servers: Dict[str, List[str]] = {}  # server -> list of tools
        self._categories: Dict[str, List[str]] = {}  # category -> list of tools
        self._initialized = False

    def initialize(self):
        """
        Initialize the registry by discovering available MCP tools.
        """
        if self._initialized:
            return

        logger.info("Initializing MCP Registry...")

        # Discover MCP tools from the global function namespace
        self._discover_mcp_tools()

        # Categorize tools
        self._categorize_tools()

        self._initialized = True
        logger.info(f"MCP Registry initialized with {len(self._tools)} tools from {len(self._servers)} servers")

    def _discover_mcp_tools(self):
        """
        Discover MCP tools by inspecting the global namespace.

        MCP tools are exposed as functions with names like:
        mcp__plugin_<server>_<Server>__<tool_name>
        """
        import sys

        # Get all functions from the main module's globals
        # In LangChain/LangGraph context, MCP tools are exposed in the tool registry
        discovered_count = 0

        # Try to discover from known MCP tool patterns
        known_mcp_tools = self._get_known_mcp_tools()

        for tool_info in known_mcp_tools:
            self._register_tool(tool_info)
            discovered_count += 1

        logger.debug(f"Discovered {discovered_count} MCP tools")

    def _get_known_mcp_tools(self) -> List[MCPToolInfo]:
        """
        Get list of known MCP tools.

        This includes tools that are commonly available in Claude Code environments.
        In a production system, this would be replaced with dynamic discovery.
        """
        tools = []

        # PubMed MCP Tools
        pubmed_tools = [
            MCPToolInfo(
                name="search_articles",
                full_name="mcp__plugin_pubmed_PubMed__search_articles",
                server="pubmed",
                category="biomedical_literature",
                description="Search PubMed for biomedical and life sciences research articles",
                when_to_use="When you need to find research papers about cell types, markers, diseases, or biological processes",
                example='search_articles(query="CD8 T cells markers human", max_results=10)',
            ),
            MCPToolInfo(
                name="get_article_metadata",
                full_name="mcp__plugin_pubmed_PubMed__get_article_metadata",
                server="pubmed",
                category="biomedical_literature",
                description="Retrieve detailed metadata for PubMed articles",
                when_to_use="When you have PMIDs and need detailed article information, abstracts, authors, citations",
                example='get_article_metadata(pmids=["35486828", "34577062"])',
            ),
            MCPToolInfo(
                name="get_full_text_article",
                full_name="mcp__plugin_pubmed_PubMed__get_full_text_article",
                server="pubmed",
                category="biomedical_literature",
                description="Retrieve full-text articles from PubMed Central",
                when_to_use="When you need complete article text with methods, results, and detailed marker information",
                example='get_full_text_article(pmc_ids=["PMC9046468"])',
            ),
            MCPToolInfo(
                name="find_related_articles",
                full_name="mcp__plugin_pubmed_PubMed__find_related_articles",
                server="pubmed",
                category="biomedical_literature",
                description="Find related articles, genes, proteins, or sequences",
                when_to_use="When you want to discover additional relevant papers or associated biological data",
                example='find_related_articles(pmids=["35486828"], link_type="pubmed_pubmed")',
            ),
            MCPToolInfo(
                name="convert_article_ids",
                full_name="mcp__plugin_pubmed_PubMed__convert_article_ids",
                server="pubmed",
                category="biomedical_literature",
                description="Convert between PMID, PMCID, and DOI identifiers",
                when_to_use="When you need to convert article IDs to access full text or cross-reference papers",
                example='convert_article_ids(ids=["35486828"], id_type="pmid")',
            ),
            MCPToolInfo(
                name="lookup_article_by_citation",
                full_name="mcp__plugin_pubmed_PubMed__lookup_article_by_citation",
                server="pubmed",
                category="biomedical_literature",
                description="Find articles by citation details (journal, year, author, etc.)",
                when_to_use="When you have a citation but need the PMID for further queries",
                example='lookup_article_by_citation(citations=[{"journal": "Nature", "year": 2020, "author": "Smith"}])',
            ),
            MCPToolInfo(
                name="get_copyright_status",
                full_name="mcp__plugin_pubmed_PubMed__get_copyright_status",
                server="pubmed",
                category="biomedical_literature",
                description="Get copyright and licensing information for articles",
                when_to_use="When you need to check if articles are open access or require permissions",
                example='get_copyright_status(pmids=["35891187"])',
            ),
        ]

        tools.extend(pubmed_tools)

        # Add more MCP servers here as they become available
        # Example: Semantic Scholar, arXiv, bioRxiv, etc.

        return tools

    def _register_tool(self, tool_info: MCPToolInfo):
        """Register a single MCP tool."""
        self._tools[tool_info.full_name] = tool_info

        # Add to server index
        if tool_info.server not in self._servers:
            self._servers[tool_info.server] = []
        self._servers[tool_info.server].append(tool_info.full_name)

        # Add to category index
        if tool_info.category not in self._categories:
            self._categories[tool_info.category] = []
        self._categories[tool_info.category].append(tool_info.full_name)

    def _categorize_tools(self):
        """
        Organize tools into meaningful categories for discovery.
        """
        # Categories are already set during tool registration
        # This method can be extended for more sophisticated categorization
        pass

    def get_tool_info(self, tool_name: str) -> Optional[MCPToolInfo]:
        """
        Get information about a specific MCP tool.

        Args:
            tool_name: Full tool name or short name

        Returns:
            MCPToolInfo if found, None otherwise
        """
        if not self._initialized:
            self.initialize()

        # Try exact match
        if tool_name in self._tools:
            return self._tools[tool_name]

        # Try partial match (short name)
        for full_name, tool_info in self._tools.items():
            if tool_info.name == tool_name or full_name.endswith(f"__{tool_name}"):
                return tool_info

        return None

    def list_tools(self, server: Optional[str] = None, category: Optional[str] = None) -> List[MCPToolInfo]:
        """
        List available MCP tools, optionally filtered by server or category.

        Args:
            server: Filter by server name (e.g., "pubmed")
            category: Filter by category (e.g., "biomedical_literature")

        Returns:
            List of MCPToolInfo objects
        """
        if not self._initialized:
            self.initialize()

        tools = []

        if server:
            tool_names = self._servers.get(server, [])
            tools = [self._tools[name] for name in tool_names]
        elif category:
            tool_names = self._categories.get(category, [])
            tools = [self._tools[name] for name in tool_names]
        else:
            tools = list(self._tools.values())

        return sorted(tools, key=lambda t: t.full_name)

    def list_servers(self) -> List[str]:
        """Get list of available MCP servers."""
        if not self._initialized:
            self.initialize()

        return sorted(self._servers.keys())

    def list_categories(self) -> List[str]:
        """Get list of tool categories."""
        if not self._initialized:
            self.initialize()

        return sorted(self._categories.keys())

    def search_tools(self, query: str) -> List[MCPToolInfo]:
        """
        Search for MCP tools by keyword.

        Args:
            query: Search query (matches name, description, when_to_use)

        Returns:
            List of matching MCPToolInfo objects
        """
        if not self._initialized:
            self.initialize()

        query_lower = query.lower()
        matches = []

        for tool_info in self._tools.values():
            # Search in various fields
            searchable = [
                tool_info.name.lower(),
                tool_info.description.lower(),
                tool_info.when_to_use.lower(),
                tool_info.category.lower(),
                tool_info.server.lower(),
            ]

            if any(query_lower in field for field in searchable):
                matches.append(tool_info)

        return sorted(matches, key=lambda t: t.name)

    def get_tools_for_task(self, task_description: str) -> List[MCPToolInfo]:
        """
        Recommend MCP tools for a specific task.

        Args:
            task_description: Natural language description of the task

        Returns:
            List of recommended MCPToolInfo objects
        """
        if not self._initialized:
            self.initialize()

        # Simple keyword-based matching
        # In production, this could use embeddings or LLM-based matching
        keywords = task_description.lower().split()

        tool_scores = {}

        for tool_info in self._tools.values():
            score = 0

            for keyword in keywords:
                if len(keyword) < 3:
                    continue

                # Score based on keyword presence in different fields
                if keyword in tool_info.name.lower():
                    score += 3
                if keyword in tool_info.description.lower():
                    score += 2
                if keyword in tool_info.when_to_use.lower():
                    score += 2
                if keyword in tool_info.category.lower():
                    score += 1

            if score > 0:
                tool_scores[tool_info.full_name] = score

        # Sort by score
        ranked_tools = sorted(tool_scores.items(), key=lambda x: x[1], reverse=True)

        return [self._tools[name] for name, score in ranked_tools[:10]]

    def get_stats(self) -> Dict[str, Any]:
        """Get registry statistics."""
        if not self._initialized:
            self.initialize()

        return {
            "total_tools": len(self._tools),
            "total_servers": len(self._servers),
            "total_categories": len(self._categories),
            "servers": {server: len(tools) for server, tools in self._servers.items()},
            "categories": {cat: len(tools) for cat, tools in self._categories.items()},
        }


# Global registry instance
_mcp_registry = MCPRegistry()


# ============================================================
# Convenience Functions
# ============================================================

def get_mcp_registry() -> MCPRegistry:
    """Get the global MCP registry instance."""
    if not _mcp_registry._initialized:
        _mcp_registry.initialize()
    return _mcp_registry


def discover_mcp_tools() -> List[MCPToolInfo]:
    """Discover all available MCP tools."""
    registry = get_mcp_registry()
    return registry.list_tools()


def search_mcp_tools(query: str) -> List[MCPToolInfo]:
    """Search for MCP tools by keyword."""
    registry = get_mcp_registry()
    return registry.search_tools(query)


def get_mcp_tool_info(tool_name: str) -> Optional[MCPToolInfo]:
    """Get information about a specific MCP tool."""
    registry = get_mcp_registry()
    return registry.get_tool_info(tool_name)


def get_mcp_servers() -> List[str]:
    """Get list of available MCP servers."""
    registry = get_mcp_registry()
    return registry.list_servers()


def get_mcp_categories() -> List[str]:
    """Get list of MCP tool categories."""
    registry = get_mcp_registry()
    return registry.list_categories()
