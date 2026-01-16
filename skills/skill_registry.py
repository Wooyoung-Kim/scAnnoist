"""
Skill Registry

Centralized registry for Claude Code skills metadata, discovery, and agent mappings.
"""

import os
import re
import time
import logging
from pathlib import Path
from typing import Dict, List, Optional, Set
from dataclasses import dataclass
from enum import Enum
from functools import lru_cache

logger = logging.getLogger(__name__)


class SkillCategory(Enum):
    """Skill categories matching the tier system"""
    CORE_ANALYSIS = "core_analysis"          # Tier 1: Core scRNA-seq analysis
    COMMUNICATION = "communication"           # Tier 2: Cell-cell communication & integration
    VISUALIZATION = "visualization"           # Tier 3: Visualization & publication
    SUPPORTING = "supporting"                 # Tier 4: Supporting tools


@dataclass
class SkillMetadata:
    """Metadata for a single Claude Code skill"""
    name: str
    description: str
    category: SkillCategory
    when_to_use: str
    relevant_agents: List[str]
    skill_path: str = ""
    requires_python: bool = False
    requires_r: bool = False
    example_invocation: str = ""


class SkillRegistry:
    """
    Central registry for Claude Code skills.

    Provides skill discovery, metadata management, and agent-to-skill mappings.
    """

    # Agent-to-skill mapping based on user requirements
    AGENT_SKILL_MAPPING = {
        "pipeline": [
            "scanpy", "scvi-tools", "anndata",
            "scientific-visualization", "scientific-slides"
        ],
        "annotation_coordinator": [
            "celltype-annotation", "scanpy", "anndata"
        ],
        "celltypist_specialist": [
            "scanpy", "celltype-annotation"
        ],
        "sctype_specialist": [
            "scanpy", "celltype-annotation"
        ],
        "literature_rag_specialist": [
            "biopython", "scanpy"
        ],
        "hematopoiesis_specialist": [
            "pathway-enrichment", "cellchat", "nichenet", "scientific-brainstorming"
        ],
        "infection_inflammation_specialist": [
            "pathway-enrichment", "cellchat", "scientific-brainstorming"
        ],
        "aging_specialist": [
            "pathway-enrichment", "scientific-brainstorming"
        ],
        "cellstate_plasticity_specialist": [
            "pathway-enrichment", "nichenet", "scientific-brainstorming"
        ],
        "epithelial_specialist": [
            "pathway-enrichment", "scientific-brainstorming"
        ],
        "organoid_stemness_specialist": [
            "pathway-enrichment", "scientific-brainstorming"
        ],
    }

    def __init__(self, skills_dir: Optional[str] = None, cache_duration: int = 300):
        """
        Initialize the skill registry.

        Args:
            skills_dir: Path to Claude Code skills directory. If None, uses .claude/skills/
            cache_duration: Cache duration in seconds (default: 300)
        """
        self.skills_dir = skills_dir or os.path.join(os.getcwd(), ".claude", "skills")
        self.cache_duration = cache_duration
        self._skills: Dict[str, SkillMetadata] = {}
        self._agent_skills: Dict[str, Set[str]] = {}
        self._last_scan: float = 0

        logger.info(f"Initializing SkillRegistry with skills_dir: {self.skills_dir}")
        self._initialize_registry()

    def _initialize_registry(self):
        """
        Initialize skill registry by scanning skills directory.

        Scans .claude/skills/ for SKILL.md files and parses their metadata.
        """
        current_time = time.time()

        # Use cached data if recent enough
        if self._skills and (current_time - self._last_scan < self.cache_duration):
            logger.debug("Using cached skill registry")
            return

        logger.info("Scanning skills directory for metadata...")
        self._scan_skills()
        self._build_agent_mappings()
        self._last_scan = current_time

    def _scan_skills(self):
        """Scan skills directory and parse SKILL.md files"""
        skills_path = Path(self.skills_dir)

        if not skills_path.exists():
            logger.warning(
                f"Skills directory not found: {self.skills_dir}. "
                f"Skill discovery will be limited."
            )
            return

        skill_files = list(skills_path.glob("*/SKILL.md"))
        logger.info(f"Found {len(skill_files)} skill files")

        for skill_file in skill_files:
            try:
                metadata = self._parse_skill_metadata(skill_file)
                if metadata:
                    self._skills[metadata.name] = metadata
                    logger.debug(f"Registered skill: {metadata.name}")
            except Exception as e:
                logger.warning(f"Failed to parse {skill_file}: {e}")

        logger.info(f"Successfully registered {len(self._skills)} skills")

    def _parse_skill_metadata(self, skill_file: Path) -> Optional[SkillMetadata]:
        """
        Parse SKILL.md file to extract metadata.

        Args:
            skill_file: Path to SKILL.md file

        Returns:
            SkillMetadata if successfully parsed, None otherwise
        """
        try:
            with open(skill_file, 'r', encoding='utf-8') as f:
                content = f.read()

            # Extract frontmatter
            frontmatter_match = re.search(r'^---\s*\n(.*?)\n---\s*\n', content, re.DOTALL | re.MULTILINE)
            if not frontmatter_match:
                logger.warning(f"No frontmatter found in {skill_file}")
                return None

            frontmatter = frontmatter_match.group(1)

            # Parse name and description from frontmatter
            name_match = re.search(r'^name:\s*(.+)$', frontmatter, re.MULTILINE)
            desc_match = re.search(r'^description:\s*(.+)$', frontmatter, re.MULTILINE)

            if not name_match or not desc_match:
                logger.warning(f"Missing name or description in {skill_file}")
                return None

            name = name_match.group(1).strip().strip('"\'')
            description = desc_match.group(1).strip().strip('"\'')

            # Extract body (when to use)
            body = content[frontmatter_match.end():].strip()
            when_to_use = body[:500] if body else "Use this skill for specialized tasks."

            # Determine category based on name
            category = self._categorize_skill(name)

            # Determine which agents should have this skill
            relevant_agents = self._get_relevant_agents(name)

            # Check if requires Python or R
            requires_python = "python" in body.lower() or ".py" in body.lower()
            requires_r = " r " in body.lower() or ".r " in body.lower()

            return SkillMetadata(
                name=name,
                description=description,
                category=category,
                when_to_use=when_to_use,
                relevant_agents=relevant_agents,
                skill_path=str(skill_file.parent),
                requires_python=requires_python,
                requires_r=requires_r,
                example_invocation=f"invoke_skill('{name}', args='...')"
            )

        except Exception as e:
            logger.error(f"Error parsing {skill_file}: {e}")
            return None

    def _categorize_skill(self, skill_name: str) -> SkillCategory:
        """Categorize a skill based on its name"""
        skill_lower = skill_name.lower()

        # Tier 1: Core analysis
        tier1_keywords = [
            "scanpy", "seurat", "scvi", "anndata", "scanvi",
            "deseq", "celltype", "annotation", "pathway", "enrichment"
        ]
        if any(kw in skill_lower for kw in tier1_keywords):
            return SkillCategory.CORE_ANALYSIS

        # Tier 2: Communication
        tier2_keywords = ["cellchat", "nichenet", "census", "communication"]
        if any(kw in skill_lower for kw in tier2_keywords):
            return SkillCategory.COMMUNICATION

        # Tier 3: Visualization
        tier3_keywords = [
            "visualization", "matplotlib", "seaborn", "slides", "schematics",
            "writing", "pptx"
        ]
        if any(kw in skill_lower for kw in tier3_keywords):
            return SkillCategory.VISUALIZATION

        # Tier 4: Supporting
        return SkillCategory.SUPPORTING

    def _get_relevant_agents(self, skill_name: str) -> List[str]:
        """Determine which agents should have access to this skill"""
        relevant = []

        for agent, skills in self.AGENT_SKILL_MAPPING.items():
            if skill_name in skills:
                relevant.append(agent)

        return relevant

    def _build_agent_mappings(self):
        """Build agent-to-skills mapping"""
        self._agent_skills.clear()

        for agent, skill_names in self.AGENT_SKILL_MAPPING.items():
            self._agent_skills[agent] = set(skill_names)

    def get_skills_for_agent(self, agent_name: str) -> List[SkillMetadata]:
        """
        Get all skills assigned to a specific agent.

        Args:
            agent_name: Name of the agent

        Returns:
            List of SkillMetadata objects for the agent
        """
        self._initialize_registry()  # Refresh if needed

        if agent_name not in self._agent_skills:
            logger.warning(f"Agent '{agent_name}' not found in registry")
            return []

        skill_names = self._agent_skills[agent_name]
        skills = []

        for name in skill_names:
            if name in self._skills:
                skills.append(self._skills[name])
            else:
                logger.debug(f"Skill '{name}' not found in registry for agent '{agent_name}'")

        return skills

    @lru_cache(maxsize=128)
    def get_skill_by_name(self, skill_name: str) -> Optional[SkillMetadata]:
        """
        Get metadata for a specific skill (cached).

        Args:
            skill_name: Name of the skill

        Returns:
            SkillMetadata if found, None otherwise
        """
        self._initialize_registry()
        return self._skills.get(skill_name)

    def search_skills(
        self,
        query: str,
        category: Optional[SkillCategory] = None
    ) -> List[SkillMetadata]:
        """
        Search for skills by keyword or category.

        Args:
            query: Search query string
            category: Optional category filter

        Returns:
            List of matching SkillMetadata objects
        """
        self._initialize_registry()

        query_lower = query.lower()
        matches = []

        for skill in self._skills.values():
            # Category filter
            if category and skill.category != category:
                continue

            # Keyword matching
            searchable_text = (
                skill.name + " " +
                skill.description + " " +
                skill.when_to_use
            ).lower()

            if query_lower in searchable_text:
                matches.append(skill)

        return matches

    def list_categories(self) -> List[str]:
        """
        List all skill categories.

        Returns:
            List of category names
        """
        return [cat.value for cat in SkillCategory]

    def list_all_skills(self) -> List[str]:
        """
        List all registered skill names.

        Returns:
            Sorted list of skill names
        """
        self._initialize_registry()
        return sorted(self._skills.keys())

    def get_skills_by_category(self, category: SkillCategory) -> List[SkillMetadata]:
        """
        Get all skills in a specific category.

        Args:
            category: Skill category

        Returns:
            List of SkillMetadata objects in the category
        """
        self._initialize_registry()
        return [skill for skill in self._skills.values() if skill.category == category]

    def refresh(self):
        """Force refresh of the skill registry"""
        self._last_scan = 0
        self._initialize_registry()

    def get_stats(self) -> Dict[str, any]:
        """
        Get registry statistics.

        Returns:
            Dictionary with statistics
        """
        self._initialize_registry()

        category_counts = {cat.value: 0 for cat in SkillCategory}
        for skill in self._skills.values():
            category_counts[skill.category.value] += 1

        return {
            "total_skills": len(self._skills),
            "total_agents": len(self._agent_skills),
            "skills_by_category": category_counts,
            "last_scan": self._last_scan,
            "skills_dir": self.skills_dir,
        }
