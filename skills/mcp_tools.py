"""
MCP Tool Wrappers

LangChain @tool wrappers for MCP (Model Context Protocol) tools,
providing a unified interface for agents to discover and use MCP tools
alongside traditional tools and Claude Code skills.

This module provides:
- Discovery tools for finding available MCP tools
- Invocation wrappers for calling MCP tools
- Search and recommendation functions
- Integration with the skill system
"""

import logging
from typing import Annotated, Optional, List, Dict, Any
from langchain_core.tools import tool
from langchain_core.runnables.config import RunnableConfig

from .mcp_registry import get_mcp_registry, MCPToolInfo

logger = logging.getLogger(__name__)


# ============================================================
# MCP Tool Discovery
# ============================================================

@tool
def list_mcp_tools(
    server: Annotated[Optional[str], "Filter by MCP server name (e.g., 'pubmed')"] = None,
    category: Annotated[Optional[str], "Filter by category (e.g., 'biomedical_literature')"] = None,
) -> Annotated[str, "List of available MCP tools"]:
    """
    List all available MCP (Model Context Protocol) tools in the environment.

    MCP tools provide access to external data sources and services like
    PubMed, databases, APIs, and more. Use this to discover what MCP
    tools are available for your analysis.

    Args:
        server: Optional filter by MCP server name
        category: Optional filter by category

    Returns:
        Formatted string with available MCP tools and descriptions

    Example:
        list_mcp_tools()  # List all MCP tools
        list_mcp_tools(server="pubmed")  # List only PubMed tools
        list_mcp_tools(category="biomedical_literature")  # List biomedical tools
    """
    try:
        registry = get_mcp_registry()
        tools = registry.list_tools(server=server, category=category)

        if server:
            header = f"MCP Tools from '{server}' server ({len(tools)})"
        elif category:
            header = f"MCP Tools in '{category}' category ({len(tools)})"
        else:
            header = f"Available MCP Tools ({len(tools)})"

        result = f"{header}\n"
        result += "=" * 70 + "\n\n"

        if not tools:
            result += "No MCP tools found"
            if server:
                result += f" from server '{server}'"
            if category:
                result += f" in category '{category}'"
            result += ".\n"
            return result

        # Group by server
        tools_by_server = {}
        for tool in tools:
            if tool.server not in tools_by_server:
                tools_by_server[tool.server] = []
            tools_by_server[tool.server].append(tool)

        for server_name, server_tools in sorted(tools_by_server.items()):
            result += f"📦 {server_name.upper()} Server\n"
            result += "-" * 70 + "\n\n"

            for tool in sorted(server_tools, key=lambda t: t.name):
                result += f"  • {tool.name}\n"
                result += f"    {tool.description}\n"

                if tool.when_to_use:
                    when = tool.when_to_use[:150]
                    result += f"    When: {when}...\n"

                result += "\n"

            result += "\n"

        result += "=" * 70 + "\n"
        result += "💡 Use get_mcp_tool_details(tool_name) to learn more.\n"
        result += "💡 Use call_mcp_tool(tool_name, params) to invoke a tool.\n"
        result += "💡 Use search_mcp_tools(query) to find tools for a specific task.\n"

        return result

    except Exception as e:
        logger.error(f"Error listing MCP tools: {e}", exc_info=True)
        return f"❌ Error listing MCP tools: {str(e)}"


@tool
def get_mcp_tool_details(
    tool_name: Annotated[str, "Name of the MCP tool (short name or full name)"],
) -> Annotated[str, "Detailed MCP tool information"]:
    """
    Get detailed information about a specific MCP tool.

    Use this to understand what a tool does, when to use it, what
    parameters it needs, and see usage examples.

    Args:
        tool_name: Name of the tool (e.g., "search_articles" or full name)

    Returns:
        Formatted string with detailed tool information

    Example:
        get_mcp_tool_details("search_articles")
        get_mcp_tool_details("get_article_metadata")
    """
    try:
        registry = get_mcp_registry()
        tool_info = registry.get_tool_info(tool_name)

        if not tool_info:
            result = f"❌ MCP tool '{tool_name}' not found.\n\n"

            # Suggest similar tools
            all_tools = registry.list_tools()
            suggestions = [t for t in all_tools if tool_name.lower() in t.name.lower()]

            if suggestions:
                result += "Did you mean one of these?\n"
                for tool in suggestions[:5]:
                    result += f"  • {tool.name} ({tool.server})\n"
            else:
                result += "Use list_mcp_tools() to see all available tools.\n"

            return result

        # Format detailed information
        result = f"MCP Tool: {tool_info.name}\n"
        result += "=" * 70 + "\n\n"

        result += f"Server: {tool_info.server}\n"
        result += f"Category: {tool_info.category}\n"
        result += f"Full name: {tool_info.full_name}\n\n"

        result += f"Description:\n{tool_info.description}\n\n"

        if tool_info.when_to_use:
            result += f"When to use:\n{tool_info.when_to_use}\n\n"

        if tool_info.parameters:
            result += "Parameters:\n"
            for param_name, param_info in tool_info.parameters.items():
                result += f"  • {param_name}: {param_info}\n"
            result += "\n"

        if tool_info.example:
            result += f"Example:\n{tool_info.example}\n\n"

        if tool_info.related_tools:
            result += f"Related tools: {', '.join(tool_info.related_tools)}\n\n"

        result += "=" * 70 + "\n"
        result += f"To use: call_mcp_tool('{tool_info.name}', {{params}})\n"

        return result

    except Exception as e:
        logger.error(f"Error getting MCP tool details: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


@tool
def search_mcp_tools(
    query: Annotated[str, "Search query or task description"],
) -> Annotated[str, "Matching MCP tools"]:
    """
    Search for MCP tools that match a query or task description.

    Use this when you're not sure which MCP tool to use for a specific task.
    The search looks in tool names, descriptions, categories, and usage guidance.

    Args:
        query: Search query or natural language task description

    Returns:
        Formatted string with matching tools

    Example:
        search_mcp_tools("search literature for markers")
        search_mcp_tools("find related papers")
        search_mcp_tools("convert article IDs")
    """
    try:
        registry = get_mcp_registry()
        matches = registry.search_tools(query)

        result = f"MCP Tools matching: '{query}'\n"
        result += "=" * 70 + "\n\n"

        if not matches:
            result += "❌ No matching MCP tools found.\n\n"
            result += "Try:\n"
            result += "  • list_mcp_tools() to see all available tools\n"
            result += "  • Different keywords in your search\n"
            result += "  • get_mcp_tool_recommendations() for task-based suggestions\n"
            return result

        result += f"Found {len(matches)} matching tools:\n\n"

        for i, tool_info in enumerate(matches[:10], 1):
            result += f"{i}. {tool_info.name} ({tool_info.server})\n"
            result += f"   Category: {tool_info.category}\n"
            result += f"   Description: {tool_info.description}\n"

            if tool_info.when_to_use:
                when = tool_info.when_to_use[:150]
                result += f"   When: {when}...\n"

            result += "\n"

        if len(matches) > 10:
            result += f"... and {len(matches) - 10} more matches.\n\n"

        result += "=" * 70 + "\n"
        result += "Use get_mcp_tool_details(tool_name) to learn more.\n"

        return result

    except Exception as e:
        logger.error(f"Error searching MCP tools: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


@tool
def get_mcp_tool_recommendations(
    task: Annotated[str, "Task description (e.g., 'validate cell type annotation')"],
) -> Annotated[str, "Recommended MCP tools for the task"]:
    """
    Get MCP tool recommendations for a specific task.

    This provides smarter recommendations than simple search by analyzing
    the task description and matching it to tool capabilities.

    Args:
        task: Natural language description of what you want to do

    Returns:
        Formatted string with recommended tools

    Example:
        get_mcp_tool_recommendations("I need to find papers about CD8 T cell markers")
        get_mcp_tool_recommendations("validate annotation with literature evidence")
        get_mcp_tool_recommendations("get full text articles about inflammation")
    """
    try:
        registry = get_mcp_registry()
        recommendations = registry.get_tools_for_task(task)

        result = f"MCP Tool Recommendations for:\n'{task}'\n"
        result += "=" * 70 + "\n\n"

        if not recommendations:
            result += "No specific recommendations found.\n\n"
            result += "Try:\n"
            result += "  • list_mcp_tools() to browse all tools\n"
            result += "  • search_mcp_tools() with keywords from your task\n"
            return result

        result += f"Recommended tools ({len(recommendations)}):\n\n"

        for i, tool_info in enumerate(recommendations, 1):
            result += f"{i}. {tool_info.name} ({tool_info.server})\n"
            result += f"   {tool_info.description}\n"

            if tool_info.when_to_use:
                result += f"   When: {tool_info.when_to_use[:150]}...\n"

            if tool_info.example:
                result += f"   Example: {tool_info.example[:100]}...\n"

            result += "\n"

        result += "=" * 70 + "\n"
        result += "Use call_mcp_tool(tool_name, params) to use a tool.\n"

        return result

    except Exception as e:
        logger.error(f"Error getting recommendations: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


# ============================================================
# MCP Tool Invocation
# ============================================================

def call_mcp_tool(
    tool_name: str,
    params: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Call an MCP tool with parameters.

    This is a helper function (not a @tool) that actually invokes MCP tools.
    It's used by the MCP annotation tools and can be called from other tools.

    Args:
        tool_name: Full MCP tool name or short name
        params: Dictionary of parameters for the tool

    Returns:
        Result dictionary from the MCP tool

    Note:
        This function attempts to call the MCP tool directly. In the Claude Code
        environment, MCP tools are available as callable functions.
    """
    try:
        registry = get_mcp_registry()
        tool_info = registry.get_tool_info(tool_name)

        if not tool_info:
            return {"error": f"MCP tool '{tool_name}' not found"}

        # In the Claude Code environment, MCP tools are available as functions
        # We need to import them dynamically or access them from the global namespace

        # Try to get the tool function
        tool_function = _get_mcp_tool_function(tool_info.full_name)

        if tool_function is None:
            # Fallback: return instructions for manual invocation
            return {
                "status": "not_callable",
                "message": f"MCP tool '{tool_name}' found but not directly callable",
                "full_name": tool_info.full_name,
                "params": params,
                "instructions": f"Use the MCP tool '{tool_info.full_name}' with these parameters in your environment"
            }

        # Call the tool
        result = tool_function(**params)

        return {
            "status": "success",
            "result": result,
            "tool": tool_info.full_name
        }

    except Exception as e:
        logger.error(f"Error calling MCP tool '{tool_name}': {e}", exc_info=True)
        return {
            "error": str(e),
            "tool": tool_name,
            "params": params
        }


def _get_mcp_tool_function(full_name: str):
    """
    Attempt to get the MCP tool function from the environment.

    In Claude Code, MCP tools are exposed as functions in the tool context.
    This function tries to access them.
    """
    import sys

    # Try to import from __main__ or globals
    try:
        # Check if the function exists in the current module's globals
        if full_name in globals():
            return globals()[full_name]

        # Try to get from __main__
        if hasattr(sys.modules.get('__main__'), full_name):
            return getattr(sys.modules['__main__'], full_name)

    except Exception as e:
        logger.debug(f"Could not access MCP tool function '{full_name}': {e}")

    return None


# ============================================================
# MCP Statistics and Info
# ============================================================

@tool
def get_mcp_stats() -> Annotated[str, "MCP registry statistics"]:
    """
    Get statistics about available MCP tools and servers.

    Returns:
        Formatted string with registry statistics

    Example:
        get_mcp_stats()
    """
    try:
        registry = get_mcp_registry()
        stats = registry.get_stats()

        result = "MCP Tool Registry Statistics\n"
        result += "=" * 70 + "\n\n"

        result += f"Total MCP tools: {stats['total_tools']}\n"
        result += f"Total MCP servers: {stats['total_servers']}\n"
        result += f"Total categories: {stats['total_categories']}\n\n"

        if stats['servers']:
            result += "Tools by server:\n"
            for server, count in sorted(stats['servers'].items()):
                result += f"  • {server}: {count} tools\n"
            result += "\n"

        if stats['categories']:
            result += "Tools by category:\n"
            for category, count in sorted(stats['categories'].items()):
                result += f"  • {category}: {count} tools\n"

        return result

    except Exception as e:
        logger.error(f"Error getting MCP stats: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


@tool
def list_mcp_servers() -> Annotated[str, "Available MCP servers"]:
    """
    List all available MCP servers.

    Returns:
        Formatted string with server names

    Example:
        list_mcp_servers()
    """
    try:
        registry = get_mcp_registry()
        servers = registry.list_servers()

        result = f"Available MCP Servers ({len(servers)})\n"
        result += "=" * 70 + "\n\n"

        for server in servers:
            tools = registry.list_tools(server=server)
            result += f"• {server} ({len(tools)} tools)\n"

        result += "\n"
        result += "Use list_mcp_tools(server='name') to see tools for a specific server.\n"

        return result

    except Exception as e:
        logger.error(f"Error listing MCP servers: {e}", exc_info=True)
        return f"❌ Error: {str(e)}"


# ============================================================
# Export all tools
# ============================================================

__all__ = [
    "list_mcp_tools",
    "get_mcp_tool_details",
    "search_mcp_tools",
    "get_mcp_tool_recommendations",
    "call_mcp_tool",
    "get_mcp_stats",
    "list_mcp_servers",
]
