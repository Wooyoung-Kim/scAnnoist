"""
Skill Manager

Core logic for Claude Code skill discovery, validation, and invocation.
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional, List

from .skill_registry import SkillRegistry, SkillCategory, SkillMetadata
from .skill_config import SkillConfig

logger = logging.getLogger(__name__)


class SkillManager:
    """
    Manages Claude Code skill discovery and invocation.

    This class provides the main interface for agents to interact with skills:
    - Discover available skills
    - Get skill metadata and information
    - Prepare skill invocations
    - Validate permissions
    """

    def __init__(self, skills_dir: Optional[str] = None, config_override: Optional[Dict[str, Any]] = None):
        """
        Initialize the Skill Manager.

        Args:
            skills_dir: Optional path to skills directory. If None, uses default from config.
            config_override: Optional configuration dictionary to override defaults
        """
        self.skills_dir = skills_dir or SkillConfig.SKILLS_DIR
        self.config_override = config_override or {}
        self.registry = SkillRegistry(skills_dir=self.skills_dir)

        logger.info(f"SkillManager initialized with skills_dir: {self.skills_dir}")
        self._validate_environment()

    def _validate_environment(self):
        """
        Ensure Claude Code skills are available.

        Raises:
            RuntimeError: If skills directory is not found or empty
        """
        skills_path = Path(self.skills_dir)

        if not skills_path.exists():
            raise RuntimeError(
                f"Claude Code skills directory not found: {self.skills_dir}\n"
                f"Expected structure: .claude/skills/\n"
                f"Make sure you're running from the project root directory."
            )

        # Check for SKILL.md files
        skill_files = list(skills_path.glob("*/SKILL.md"))
        if len(skill_files) == 0:
            logger.warning(
                f"No skills found in {self.skills_dir}. "
                f"Expected to find SKILL.md files in subdirectories."
            )

        logger.info(f"Found {len(skill_files)} skill files in {self.skills_dir}")

    def discover_skills(self, category: Optional[str] = None) -> List[str]:
        """
        Discover all available Claude Code skills.

        Args:
            category: Optional category filter (core_analysis, communication, visualization, supporting)

        Returns:
            List of skill names

        Example:
            >>> manager = SkillManager()
            >>> skills = manager.discover_skills(category="core_analysis")
            >>> print(skills)
            ['scanpy', 'scvi-tools', 'anndata', 'celltype-annotation', ...]
        """
        if category:
            try:
                skill_cat = SkillCategory(category)
                skills = self.registry.get_skills_by_category(skill_cat)
                return [s.name for s in skills]
            except ValueError:
                logger.warning(f"Invalid category: {category}")
                return []
        else:
            return self.registry.list_all_skills()

    def get_skill_info(self, skill_name: str) -> Dict[str, Any]:
        """
        Get detailed information about a skill.

        Args:
            skill_name: Name of the skill

        Returns:
            Dictionary with skill information:
            - name: Skill name
            - description: Short description
            - category: Skill category
            - when_to_use: Usage guidance
            - example: Example invocation
            - error: Error message if skill not found

        Example:
            >>> manager = SkillManager()
            >>> info = manager.get_skill_info("scanpy")
            >>> print(info['description'])
            'Single-cell RNA-seq analysis toolkit'
        """
        metadata = self.registry.get_skill_by_name(skill_name)

        if not metadata:
            return {
                "error": f"Skill '{skill_name}' not found",
                "available_skills": self.discover_skills()[:10],  # Show first 10
                "suggestion": "Use discover_skills() to see all available skills"
            }

        return {
            "name": metadata.name,
            "description": metadata.description,
            "category": metadata.category.value,
            "when_to_use": metadata.when_to_use,
            "example": metadata.example_invocation,
            "requires_python": metadata.requires_python,
            "requires_r": metadata.requires_r,
            "relevant_agents": metadata.relevant_agents,
        }

    def invoke_skill(
        self,
        skill_name: str,
        args: Optional[str] = None,
        context: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Prepare skill invocation.

        Note: This function PREPARES skill invocation rather than actually invoking it.
        Actual skill execution happens through the Skill tool provided by Claude Code CLI.

        Args:
            skill_name: Name of the skill to invoke
            args: Optional arguments string
            context: Optional context dictionary (e.g., species, tissue type)

        Returns:
            Dictionary with invocation metadata:
            - status: "ready" or "error"
            - skill: Skill name
            - args: Arguments
            - instructions: How to invoke the skill
            - metadata: Skill information
            - message: Error message if status is "error"

        Example:
            >>> manager = SkillManager()
            >>> result = manager.invoke_skill("scanpy", args="quality control analysis")
            >>> print(result['instructions'])
            To invoke this skill, use: Skill(skill='scanpy', args='quality control analysis')
        """
        # Check if skill exists
        metadata = self.registry.get_skill_by_name(skill_name)
        if not metadata:
            return {
                "status": "error",
                "message": (
                    f"Skill '{skill_name}' not found in registry. "
                    f"Available skills: {', '.join(self.discover_skills()[:5])}..."
                ),
                "suggestion": "Use list_available_skills() to see all available skills"
            }

        # Check permissions (if strict mode enabled)
        if SkillConfig.STRICT_PERMISSIONS:
            try:
                if not self.check_skill_permissions(skill_name):
                    return {
                        "status": "error",
                        "message": (
                            f"Skill '{skill_name}' not permitted. "
                            f"Add 'Skill({skill_name})' to .claude/settings.local.json permissions"
                        )
                    }
            except Exception as e:
                logger.warning(f"Permission check failed: {e}")

        # Prepare invocation
        return {
            "status": "ready",
            "skill": skill_name,
            "args": args or "",
            "context": context or {},
            "instructions": (
                f"To invoke this skill, the agent should use the Skill tool "
                f"provided by Claude Code with:\n"
                f"  Skill(skill='{skill_name}', args='{args or ''}')\n\n"
                f"The skill will be executed in the Claude Code environment."
            ),
            "metadata": self.get_skill_info(skill_name)
        }

    def check_skill_permissions(self, skill_name: str) -> bool:
        """
        Check if a skill is permitted in .claude/settings.local.json.

        Args:
            skill_name: Name of the skill

        Returns:
            True if skill is permitted, False otherwise

        Raises:
            PermissionError: If strict permissions enabled and skill not permitted
        """
        settings_path = Path(self.skills_dir).parent / "settings.local.json"

        if not settings_path.exists():
            # No settings file, assume all skills allowed
            logger.debug(f"No settings file found at {settings_path}, assuming all skills allowed")
            return True

        try:
            with open(settings_path, 'r') as f:
                settings = json.load(f)

            allowed = settings.get("permissions", {}).get("allow", [])
            skill_permission = f"Skill({skill_name})"

            is_allowed = skill_permission in allowed

            if not is_allowed and SkillConfig.STRICT_PERMISSIONS:
                raise PermissionError(
                    f"Skill '{skill_name}' not permitted. "
                    f"Add 'Skill({skill_name})' to .claude/settings.local.json permissions array"
                )

            return is_allowed

        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse settings file: {e}")
            return True  # Fail open if can't parse settings
        except Exception as e:
            logger.error(f"Error checking permissions: {e}")
            return True  # Fail open on errors

    def list_skills_by_category(self, category: str) -> List[str]:
        """
        List all skills in a specific category.

        Args:
            category: Category name (core_analysis, communication, visualization, supporting)

        Returns:
            List of skill names in the category

        Example:
            >>> manager = SkillManager()
            >>> core_skills = manager.list_skills_by_category("core_analysis")
            >>> print(core_skills)
            ['scanpy', 'scvi-tools', 'anndata', 'celltype-annotation', ...]
        """
        try:
            skill_cat = SkillCategory(category)
            skills = self.registry.get_skills_by_category(skill_cat)
            return sorted([s.name for s in skills])
        except ValueError:
            logger.warning(f"Invalid category: {category}")
            valid_categories = self.registry.list_categories()
            logger.info(f"Valid categories: {valid_categories}")
            return []

    def search_skills(self, query: str, category: Optional[str] = None) -> List[Dict[str, str]]:
        """
        Search for skills matching a query.

        Args:
            query: Search query string
            category: Optional category filter

        Returns:
            List of dictionaries with skill name and description

        Example:
            >>> manager = SkillManager()
            >>> results = manager.search_skills("pathway analysis")
            >>> for skill in results:
            ...     print(f"{skill['name']}: {skill['description']}")
        """
        category_enum = None
        if category:
            try:
                category_enum = SkillCategory(category)
            except ValueError:
                logger.warning(f"Invalid category: {category}")

        matches = self.registry.search_skills(query, category_enum)

        return [
            {
                "name": skill.name,
                "description": skill.description,
                "category": skill.category.value,
                "when_to_use": skill.when_to_use[:200] + "..."  # Truncate for brevity
            }
            for skill in matches
        ]

    def get_skills_for_agent(self, agent_name: str) -> List[str]:
        """
        Get all skills assigned to a specific agent.

        Args:
            agent_name: Name of the agent (e.g., "pipeline", "celltypist_specialist")

        Returns:
            List of skill names assigned to the agent

        Example:
            >>> manager = SkillManager()
            >>> skills = manager.get_skills_for_agent("pipeline")
            >>> print(skills)
            ['scanpy', 'scvi-tools', 'anndata', 'scientific-visualization', ...]
        """
        skills = self.registry.get_skills_for_agent(agent_name)
        return [s.name for s in skills]

    def get_stats(self) -> Dict[str, Any]:
        """
        Get statistics about the skill registry.

        Returns:
            Dictionary with statistics

        Example:
            >>> manager = SkillManager()
            >>> stats = manager.get_stats()
            >>> print(f"Total skills: {stats['total_skills']}")
        """
        return self.registry.get_stats()

    def refresh(self):
        """Force refresh of the skill registry"""
        logger.info("Refreshing skill registry...")
        self.registry.refresh()

    def validate_skill_environment(self, skill_name: str) -> Dict[str, Any]:
        """
        Validate that a skill's environment requirements are met.

        Args:
            skill_name: Name of the skill

        Returns:
            Dictionary with validation results:
            - valid: Boolean indicating if environment is valid
            - issues: List of issues found
            - recommendations: List of recommendations

        Example:
            >>> manager = SkillManager()
            >>> validation = manager.validate_skill_environment("scanpy")
            >>> if not validation['valid']:
            ...     print(validation['issues'])
        """
        metadata = self.registry.get_skill_by_name(skill_name)

        if not metadata:
            return {
                "valid": False,
                "issues": [f"Skill '{skill_name}' not found"],
                "recommendations": ["Check skill name spelling"]
            }

        issues = []
        recommendations = []

        # Check if skill directory exists
        skill_path = Path(metadata.skill_path)
        if not skill_path.exists():
            issues.append(f"Skill directory not found: {skill_path}")
            recommendations.append("Ensure Claude Code skills are properly installed")

        # Check Python/R requirements
        if metadata.requires_python:
            # Could add actual Python environment checks here
            pass

        if metadata.requires_r:
            # Could add actual R environment checks here
            pass

        # Check permissions
        try:
            if not self.check_skill_permissions(skill_name):
                issues.append(f"Skill '{skill_name}' not permitted")
                recommendations.append(f"Add 'Skill({skill_name})' to .claude/settings.local.json")
        except Exception as e:
            issues.append(f"Permission check failed: {e}")

        return {
            "valid": len(issues) == 0,
            "issues": issues,
            "recommendations": recommendations if issues else ["Skill environment is valid"]
        }
