# Utility functions for model setup and configuration

import os
from typing import Optional
from langchain_openai import ChatOpenAI
from langchain_anthropic import ChatAnthropic


# Default model configurations
DEFAULT_MODELS = {
    "pipeline": "gpt-4o",
    "annotation_coordinator": "gpt-4o",
    "celltypist_specialist": "gpt-4o-mini",
    "sctype_specialist": "gpt-4o-mini",
    "literature_rag_specialist": "gpt-4o-mini",
    "integration": "gpt-4o-mini",
}


def set_model(
    model_name: Optional[str] = None,
    agent_name: str = "default",
    temperature: float = 0.1,
):
    """
    Create a language model instance for an agent.
    
    Args:
        model_name: Override model name. If None, uses default for agent.
        agent_name: Name of the agent (for default model lookup)
        temperature: Model temperature
        
    Returns:
        Configured LLM instance
    """
    # Determine model name
    if model_name is None:
        model_name = DEFAULT_MODELS.get(agent_name, "gpt-4o-mini")
    
    # Check for API keys and select provider
    if "claude" in model_name.lower() or "anthropic" in model_name.lower():
        return ChatAnthropic(
            model=model_name,
            temperature=temperature,
            max_tokens=4096,
        )
    else:
        # Default to OpenAI
        return ChatOpenAI(
            model=model_name,
            temperature=temperature,
        )


def get_env_or_raise(key: str) -> str:
    """Get environment variable or raise error."""
    value = os.getenv(key)
    if not value:
        raise ValueError(f"Environment variable '{key}' not set")
    return value


# ============================================================
# Skill Integration Utilities
# ============================================================

def get_agent_skills(agent_name: str) -> list:
    """
    Get list of base skills for an agent.

    Args:
        agent_name: Name of the agent

    Returns:
        List of skill names assigned to the agent

    Example:
        >>> skills = get_agent_skills("pipeline")
        >>> print(skills)
        ['scanpy', 'scvi-tools', 'anndata', 'scientific-visualization', 'scientific-slides']
    """
    try:
        from scrna_agent.skills.skill_config import load_agent_skills
        config = load_agent_skills(agent_name)
        return config.get("base_skills", [])
    except Exception as e:
        import logging
        logging.warning(f"Failed to load agent skills: {e}")
        return []


def is_skill_available(skill_name: str) -> bool:
    """
    Check if a skill is available in the environment.

    Args:
        skill_name: Name of the skill to check

    Returns:
        True if skill is available, False otherwise

    Example:
        >>> if is_skill_available("scanpy"):
        ...     print("Scanpy skill is available")
    """
    try:
        from scrna_agent.skills.skill_manager import SkillManager
        manager = SkillManager()
        skills = manager.discover_skills()
        return skill_name in skills
    except Exception as e:
        import logging
        logging.warning(f"Failed to check skill availability: {e}")
        return False


def get_available_skills(category: Optional[str] = None) -> list:
    """
    Get list of all available skills, optionally filtered by category.

    Args:
        category: Optional category filter (core_analysis, communication, visualization, supporting)

    Returns:
        List of available skill names

    Example:
        >>> all_skills = get_available_skills()
        >>> core_skills = get_available_skills("core_analysis")
    """
    try:
        from scrna_agent.skills.skill_manager import SkillManager
        manager = SkillManager()
        return manager.discover_skills(category=category)
    except Exception as e:
        import logging
        logging.warning(f"Failed to get available skills: {e}")
        return []


def validate_skill_environment() -> dict:
    """
    Validate that the skill integration environment is properly configured.

    Returns:
        Dictionary with validation results:
        - skills_available: Boolean indicating if skills are available
        - skills_dir_exists: Boolean indicating if skills directory exists
        - total_skills: Number of available skills
        - issues: List of any issues found

    Example:
        >>> validation = validate_skill_environment()
        >>> if not validation['skills_available']:
        ...     print("Skills not available:", validation['issues'])
    """
    issues = []
    skills_available = False
    skills_dir_exists = False
    total_skills = 0

    try:
        from scrna_agent.skills.skill_manager import SkillManager
        from scrna_agent.skills.skill_config import SkillConfig
        import os

        # Check if skills directory exists
        skills_dir_exists = os.path.exists(SkillConfig.SKILLS_DIR)
        if not skills_dir_exists:
            issues.append(f"Skills directory not found: {SkillConfig.SKILLS_DIR}")

        # Try to create manager and discover skills
        manager = SkillManager()
        skills = manager.discover_skills()
        total_skills = len(skills)

        if total_skills == 0:
            issues.append("No skills found in skills directory")
        else:
            skills_available = True

    except Exception as e:
        issues.append(f"Failed to validate skill environment: {str(e)}")

    return {
        "skills_available": skills_available,
        "skills_dir_exists": skills_dir_exists,
        "total_skills": total_skills,
        "issues": issues
    }
