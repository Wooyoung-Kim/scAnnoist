"""
Claude Code Skills Integration Module

This module provides integration between LangGraph agents and Claude Code skills,
enabling agents to discover and invoke specialized capabilities dynamically.
Includes MCP tool discovery and integration.
"""

from .skill_registry import SkillRegistry, SkillMetadata, SkillCategory
from .skill_manager import SkillManager
from .skill_config import load_skill_config, load_agent_skills, get_category_skills

# MCP tool integration
from .mcp_tools import (
    list_mcp_tools,
    get_mcp_tool_details,
    search_mcp_tools,
    get_mcp_tool_recommendations,
    call_mcp_tool,
    get_mcp_stats,
    list_mcp_servers,
)

from .mcp_registry import (
    get_mcp_registry,
    discover_mcp_tools,
    get_mcp_tool_info,
    get_mcp_servers,
    get_mcp_categories,
)

__all__ = [
    "SkillRegistry",
    "SkillMetadata",
    "SkillCategory",
    "SkillManager",
    "load_skill_config",
    "load_agent_skills",
    "get_category_skills",
    # MCP tools
    "list_mcp_tools",
    "get_mcp_tool_details",
    "search_mcp_tools",
    "get_mcp_tool_recommendations",
    "call_mcp_tool",
    "get_mcp_stats",
    "list_mcp_servers",
    "get_mcp_registry",
    "discover_mcp_tools",
    "get_mcp_tool_info",
    "get_mcp_servers",
    "get_mcp_categories",
]
