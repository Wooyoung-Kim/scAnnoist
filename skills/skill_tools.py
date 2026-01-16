"""
Skill Tool Wrappers

LangChain @tool wrappers for Claude Code skills, providing the interface
for agents to discover and invoke skills.
"""

import logging
from typing import Annotated, Optional
from langchain_core.tools import tool
from langchain_core.runnables.config import RunnableConfig

from .skill_manager import SkillManager

logger = logging.getLogger(__name__)

# Initialize global skill manager (shared across all tools)
_skill_manager = SkillManager()


# ============================================================
# Discovery Tools
# ============================================================

@tool
def list_available_skills(
    category: Annotated[Optional[str], "Filter by category: core_analysis, communication, visualization, supporting"] = None
) -> Annotated[str, "List of available Claude Code skills"]:
    """
    List all Claude Code skills available in this environment.

    Use this tool to discover what specialized capabilities are available
    beyond your traditional tools.

    Categories:
    - core_analysis: scanpy, scvi-tools, anndata, celltype-annotation, pathway-enrichment
    - communication: cellchat, nichenet, cellxgene-census
    - visualization: scientific-visualization, matplotlib, seaborn, scientific-slides
    - supporting: biopython, scientific-brainstorming, scientific-critical-thinking

    Args:
        category: Optional category filter to narrow results

    Returns:
        Formatted string listing available skills with descriptions

    Example:
        list_available_skills()  # List all skills
        list_available_skills(category="core_analysis")  # List only core analysis skills
    """
    try:
        if category:
            skills = _skill_manager.discover_skills(category=category)
            category_display = category.replace("_", " ").title()
            header = f"Available Skills in {category_display} ({len(skills)})"
        else:
            skills = _skill_manager.discover_skills()
            header = f"Available Claude Code Skills ({len(skills)})"

        result = f"{header}\n"
        result += "=" * 70 + "\n\n"

        if not skills:
            result += f"No skills found"
            if category:
                result += f" in category '{category}'"
            result += ".\n"
            return result

        # Get detailed info for each skill
        for skill_name in sorted(skills):
            info = _skill_manager.get_skill_info(skill_name)

            if "error" not in info:
                result += f"{skill_name}\n"
                result += f"  Category: {info['category']}\n"
                result += f"  Description: {info['description']}\n"

                # Show when to use (truncated)
                when_to_use = info.get('when_to_use', '')[:150]
                if when_to_use:
                    result += f"  When to use: {when_to_use}...\n"

                result += "\n"

        result += f"\nUse get_skill_details(skill_name) to learn more about a specific skill.\n"
        result += f"Use invoke_skill(skill_name, args) to use a skill.\n"

        return result

    except Exception as e:
        logger.error(f"Error listing skills: {e}", exc_info=True)
        return f"Error listing skills: {str(e)}"


@tool
def get_skill_details(
    skill_name: Annotated[str, "Name of the skill to get details about"]
) -> Annotated[str, "Detailed skill information"]:
    """
    Get detailed information about a specific Claude Code skill.

    Use this to understand when and how to use a particular skill,
    including usage guidance and examples.

    Args:
        skill_name: Name of the skill (e.g., "scanpy", "cellchat", "scientific-visualization")

    Returns:
        Formatted string with detailed skill information

    Example:
        get_skill_details("scanpy")
        get_skill_details("pathway-enrichment")
    """
    try:
        info = _skill_manager.get_skill_info(skill_name)

        if "error" in info:
            result = f"Error: {info['error']}\n\n"

            if "available_skills" in info:
                result += "Did you mean one of these?\n"
                for skill in info["available_skills"]:
                    result += f"  - {skill}\n"

            result += "\nUse list_available_skills() to see all available skills."
            return result

        # Format detailed information
        result = f"Skill: {info['name']}\n"
        result += "=" * 70 + "\n\n"

        result += f"Description:\n{info['description']}\n\n"

        result += f"Category: {info['category']}\n\n"

        result += f"When to use:\n{info['when_to_use']}\n\n"

        if info.get('requires_python'):
            result += "Requirements: Python environment\n"
        if info.get('requires_r'):
            result += "Requirements: R environment\n"

        if info.get('relevant_agents'):
            result += f"\nRelevant for agents: {', '.join(info['relevant_agents'])}\n"

        if info.get('example'):
            result += f"\nExample invocation:\n{info['example']}\n"

        result += "\n" + "-" * 70 + "\n"
        result += "To use this skill, call: invoke_skill(skill_name, args)\n"

        return result

    except Exception as e:
        logger.error(f"Error getting skill details: {e}", exc_info=True)
        return f"Error getting skill details: {str(e)}"


# ============================================================
# Invocation Tool
# ============================================================

@tool
async def invoke_skill(
    skill_name: Annotated[str, "Name of the Claude Code skill to invoke"],
    args: Annotated[Optional[str], "Optional arguments for the skill"] = None,
    config: RunnableConfig = None,
) -> Annotated[str, "Skill invocation result or instructions"]:
    """
    Invoke a Claude Code skill dynamically.

    This tool provides a bridge to Claude Code skills, allowing agents to
    access specialized capabilities beyond their base tools.

    IMPORTANT: This function prepares skill invocation. The actual skill
    execution happens through Claude Code's Skill tool. This function will
    provide instructions on how to invoke the skill properly.

    Args:
        skill_name: Name of the skill (e.g., "scanpy", "pathway-enrichment", "cellchat")
        args: Optional arguments string to pass to the skill
        config: Optional configuration (automatically provided by LangChain)

    Returns:
        Invocation instructions or prepared invocation result

    Examples:
        invoke_skill("scanpy", "analyze quality control metrics")
        invoke_skill("cellchat", "infer cell-cell communication networks")
        invoke_skill("pathway-enrichment", "GO term enrichment for cluster markers")
        invoke_skill("scientific-visualization", "create publication-quality UMAP plot")
    """
    try:
        # Extract context from config if available
        context = {}
        if config and "configurable" in config:
            context = config.get("configurable", {})

        # Prepare invocation
        result = _skill_manager.invoke_skill(skill_name, args, context)

        if result["status"] == "error":
            output = f"❌ Error: {result['message']}\n\n"

            if "suggestion" in result:
                output += f"💡 Suggestion: {result['suggestion']}\n"

            return output

        # Format invocation instructions
        output = f"✓ Skill Invocation Prepared: {skill_name}\n"
        output += "=" * 70 + "\n\n"

        output += f"📋 Instructions:\n{result['instructions']}\n\n"

        # Include skill metadata for context
        metadata = result.get("metadata", {})
        if metadata and "description" not in result.get("error", ""):
            output += f"📖 About this skill:\n"
            output += f"  Description: {metadata.get('description', 'N/A')}\n"
            output += f"  Category: {metadata.get('category', 'N/A')}\n\n"

        # Add context if provided
        if args:
            output += f"📝 Arguments: {args}\n\n"

        if context:
            output += f"🔧 Context: {', '.join(f'{k}={v}' for k, v in context.items())}\n\n"

        output += "=" * 70 + "\n"
        output += "Note: The skill will be executed in the Claude Code environment.\n"
        output += "Results will be available after skill execution completes.\n"

        return output

    except Exception as e:
        logger.error(f"Error invoking skill: {e}", exc_info=True)
        return f"❌ Error invoking skill '{skill_name}': {str(e)}"


# ============================================================
# Search Tool
# ============================================================

@tool
def search_skills_for_task(
    task_description: Annotated[str, "Description of the task you want to accomplish"]
) -> Annotated[str, "Recommended skills for this task"]:
    """
    Search for Claude Code skills that match a specific task or analysis need.

    Use this when you're not sure which skill to use for a particular task.
    The tool will search skill descriptions and usage guidance to find matches.

    Args:
        task_description: Natural language description of what you want to do

    Returns:
        Formatted string with matching skills and recommendations

    Examples:
        search_skills_for_task("I need to perform pathway enrichment analysis")
        search_skills_for_task("visualize cell-cell communication networks")
        search_skills_for_task("create publication-quality scientific figures")
        search_skills_for_task("annotate cell types automatically")
    """
    try:
        # Extract keywords for searching
        keywords = task_description.lower().split()

        # Search for matching skills
        all_matches = []

        for keyword in keywords:
            if len(keyword) < 3:  # Skip very short keywords
                continue

            matches = _skill_manager.search_skills(keyword)
            all_matches.extend(matches)

        # Remove duplicates and score by frequency
        skill_scores = {}
        for match in all_matches:
            name = match['name']
            skill_scores[name] = skill_scores.get(name, 0) + 1

        # Sort by score
        ranked_skills = sorted(skill_scores.items(), key=lambda x: x[1], reverse=True)

        # Format results
        result = f"Skills matching: '{task_description}'\n"
        result += "=" * 70 + "\n\n"

        if not ranked_skills:
            result += "❌ No matching skills found.\n\n"
            result += "Try:\n"
            result += "  - list_available_skills() to see all available skills\n"
            result += "  - More specific keywords in your task description\n"
            result += "  - Different phrasing of your task\n"
            return result

        result += f"Found {len(ranked_skills)} matching skills:\n\n"

        # Show top 5 matches
        for i, (skill_name, score) in enumerate(ranked_skills[:5], 1):
            info = _skill_manager.get_skill_info(skill_name)

            if "error" not in info:
                result += f"{i}. {skill_name} (relevance: {score})\n"
                result += f"   Category: {info['category']}\n"
                result += f"   Description: {info['description']}\n"

                # Truncate when_to_use for brevity
                when_to_use = info.get('when_to_use', '')[:150]
                if when_to_use:
                    result += f"   When to use: {when_to_use}...\n"

                result += "\n"

        if len(ranked_skills) > 5:
            result += f"... and {len(ranked_skills) - 5} more matches.\n\n"

        result += "=" * 70 + "\n"
        result += "To learn more: get_skill_details(skill_name)\n"
        result += "To use a skill: invoke_skill(skill_name, args)\n"

        return result

    except Exception as e:
        logger.error(f"Error searching skills: {e}", exc_info=True)
        return f"❌ Error searching skills: {str(e)}"


# ============================================================
# Utility Tools
# ============================================================

@tool
def get_agent_skills(
    agent_name: Annotated[str, "Name of the agent (e.g., 'pipeline', 'celltypist_specialist')"]
) -> Annotated[str, "List of skills assigned to this agent"]:
    """
    Get all skills assigned to a specific agent.

    Use this to see what skills are configured for a particular agent type.

    Args:
        agent_name: Name of the agent

    Returns:
        Formatted string with agent's skills

    Example:
        get_agent_skills("pipeline")
        get_agent_skills("hematopoiesis_specialist")
    """
    try:
        skills = _skill_manager.get_skills_for_agent(agent_name)

        result = f"Skills for agent '{agent_name}':\n"
        result += "=" * 70 + "\n\n"

        if not skills:
            result += f"No skills assigned to agent '{agent_name}'.\n"
            result += "This agent will use only its traditional tools.\n"
            return result

        result += f"Total: {len(skills)} skills\n\n"

        for skill_name in sorted(skills):
            info = _skill_manager.get_skill_info(skill_name)
            if "error" not in info:
                result += f"  • {skill_name}\n"
                result += f"    {info['description']}\n\n"

        return result

    except Exception as e:
        logger.error(f"Error getting agent skills: {e}", exc_info=True)
        return f"❌ Error getting agent skills: {str(e)}"


@tool
def get_skill_stats() -> Annotated[str, "Statistics about available skills"]:
    """
    Get statistics about the skill registry.

    Returns:
        Formatted string with statistics

    Example:
        get_skill_stats()
    """
    try:
        stats = _skill_manager.get_stats()

        result = "Skill Registry Statistics\n"
        result += "=" * 70 + "\n\n"

        result += f"Total skills registered: {stats['total_skills']}\n"
        result += f"Total agents configured: {stats['total_agents']}\n"
        result += f"Skills directory: {stats['skills_dir']}\n\n"

        result += "Skills by category:\n"
        for category, count in stats['skills_by_category'].items():
            result += f"  • {category}: {count}\n"

        return result

    except Exception as e:
        logger.error(f"Error getting skill stats: {e}", exc_info=True)
        return f"❌ Error getting skill stats: {str(e)}"


# ============================================================
# Export all tools
# ============================================================

__all__ = [
    "list_available_skills",
    "get_skill_details",
    "invoke_skill",
    "search_skills_for_task",
    "get_agent_skills",
    "get_skill_stats",
]
