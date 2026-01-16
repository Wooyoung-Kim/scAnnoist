"""
Skill Configuration Management

Handles loading and parsing of agent skills configuration from YAML files.
"""

import os
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List
import logging

logger = logging.getLogger(__name__)


class SkillConfig:
    """Global skill configuration constants"""

    # Paths
    SKILLS_DIR = os.getenv(
        "SCRNA_AGENT_SKILLS_DIR",
        os.path.join(os.getcwd(), ".claude", "skills")
    )

    CONFIG_PATH = os.getenv(
        "SCRNA_AGENT_SKILL_CONFIG",
        os.path.join(os.getcwd(), "configs", "agent_skills.yaml")
    )

    # Discovery settings
    CACHE_DURATION = int(os.getenv("SKILL_CACHE_DURATION", "300"))
    DISCOVERY_ENABLED = os.getenv("SKILL_DISCOVERY_ENABLED", "true").lower() == "true"

    # Permissions
    AUTO_APPROVE_SKILLS = os.getenv("SKILL_AUTO_APPROVE", "false").lower() == "true"
    STRICT_PERMISSIONS = os.getenv("SKILL_STRICT_PERMISSIONS", "true").lower() == "true"


def load_skill_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load agent skills configuration from YAML file.

    Args:
        config_path: Optional path to configuration file. If None, uses default path.

    Returns:
        Configuration dictionary

    Raises:
        FileNotFoundError: If configuration file doesn't exist
        yaml.YAMLError: If configuration file is invalid YAML
    """
    if config_path is None:
        config_path = SkillConfig.CONFIG_PATH

    config_path = Path(config_path)

    if not config_path.exists():
        logger.warning(
            f"Skill configuration file not found: {config_path}. "
            f"Using default configuration."
        )
        return _get_default_config()

    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)

        logger.info(f"Loaded skill configuration from {config_path}")
        return config

    except yaml.YAMLError as e:
        logger.error(f"Failed to parse skill configuration: {e}")
        raise
    except Exception as e:
        logger.error(f"Failed to load skill configuration: {e}")
        raise


def load_agent_skills(agent_name: str, config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load skills configuration for a specific agent.

    Args:
        agent_name: Name of the agent (e.g., "pipeline_agent", "celltypist_specialist")
        config_path: Optional path to configuration file

    Returns:
        Agent-specific configuration dictionary with keys:
        - base_skills: List of skill names
        - discovery_enabled: Boolean
        - discovery_categories: List of category names
        - description: Agent description
    """
    config = load_skill_config(config_path)

    agent_config = config.get("agents", {}).get(agent_name, {})

    if not agent_config:
        logger.warning(
            f"Agent '{agent_name}' not configured in agent_skills.yaml. "
            f"Using empty configuration."
        )
        return {
            "base_skills": [],
            "discovery_enabled": False,
            "discovery_categories": [],
            "description": f"Agent: {agent_name}",
        }

    # Ensure all required keys exist
    agent_config.setdefault("base_skills", [])
    agent_config.setdefault("discovery_enabled", False)
    agent_config.setdefault("discovery_categories", [])
    agent_config.setdefault("description", f"Agent: {agent_name}")

    logger.debug(
        f"Loaded configuration for agent '{agent_name}': "
        f"{len(agent_config['base_skills'])} base skills, "
        f"discovery={'enabled' if agent_config['discovery_enabled'] else 'disabled'}"
    )

    return agent_config


def get_category_skills(category: str, config_path: Optional[str] = None) -> List[str]:
    """
    Get all skills in a specific category.

    Args:
        category: Category name (e.g., "tier1_core", "tier2_communication")
        config_path: Optional path to configuration file

    Returns:
        List of skill names in the category
    """
    config = load_skill_config(config_path)

    categories = config.get("categories", {})
    skills = categories.get(category, [])

    if not skills:
        logger.warning(f"Category '{category}' not found or empty in configuration")

    return skills


def get_all_configured_skills(config_path: Optional[str] = None) -> List[str]:
    """
    Get all skills configured across all categories.

    Args:
        config_path: Optional path to configuration file

    Returns:
        List of all unique skill names
    """
    config = load_skill_config(config_path)

    all_skills = set()
    categories = config.get("categories", {})

    for category_skills in categories.values():
        all_skills.update(category_skills)

    return sorted(list(all_skills))


def get_agent_list(config_path: Optional[str] = None) -> List[str]:
    """
    Get list of all configured agents.

    Args:
        config_path: Optional path to configuration file

    Returns:
        List of agent names
    """
    config = load_skill_config(config_path)
    return list(config.get("agents", {}).keys())


def _get_default_config() -> Dict[str, Any]:
    """
    Get default configuration when config file is not found.

    Returns:
        Default configuration dictionary
    """
    return {
        "version": "1.0",
        "discovery": {
            "enabled": True,
            "cache_duration_seconds": 300,
            "skill_directory": ".claude/skills",
        },
        "categories": {
            "tier1_core": [
                "scanpy", "scvi-tools", "anndata", "celltype-annotation"
            ],
            "tier2_communication": [
                "cellchat", "nichenet"
            ],
            "tier3_visualization": [
                "scientific-visualization", "matplotlib", "seaborn"
            ],
            "tier4_supporting": [
                "biopython", "scientific-brainstorming"
            ],
        },
        "agents": {},
        "permissions": {
            "required_in_settings": [],
            "auto_approve": False,
        },
    }


# Convenience functions for common operations
def is_skill_enabled_for_agent(agent_name: str, skill_name: str, config_path: Optional[str] = None) -> bool:
    """
    Check if a skill is enabled for an agent (either in base_skills or via discovery categories).

    Args:
        agent_name: Name of the agent
        skill_name: Name of the skill
        config_path: Optional path to configuration file

    Returns:
        True if skill is enabled for agent, False otherwise
    """
    agent_config = load_agent_skills(agent_name, config_path)

    # Check base skills
    if skill_name in agent_config.get("base_skills", []):
        return True

    # Check discovery categories if enabled
    if agent_config.get("discovery_enabled", False):
        config = load_skill_config(config_path)
        categories = config.get("categories", {})

        for category_name in agent_config.get("discovery_categories", []):
            category_skills = categories.get(category_name, [])
            if skill_name in category_skills:
                return True

    return False


def get_all_skills_for_agent(agent_name: str, config_path: Optional[str] = None) -> List[str]:
    """
    Get all skills available to an agent (base skills + discovery category skills).

    Args:
        agent_name: Name of the agent
        config_path: Optional path to configuration file

    Returns:
        List of all available skill names
    """
    agent_config = load_agent_skills(agent_name, config_path)
    config = load_skill_config(config_path)

    all_skills = set(agent_config.get("base_skills", []))

    # Add skills from discovery categories if enabled
    if agent_config.get("discovery_enabled", False):
        categories = config.get("categories", {})
        for category_name in agent_config.get("discovery_categories", []):
            category_skills = categories.get(category_name, [])
            all_skills.update(category_skills)

    return sorted(list(all_skills))
