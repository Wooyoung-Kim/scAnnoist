"""
Biological context validation for scRNA-seq analysis.
Validates cell type annotations against sample type constraints.

Key insight from JMV_ALI analysis: organoid samples were incorrectly annotated
as 80.6% neutrophils because S100A8/S100A9 stress markers were misinterpreted
as immune cell markers. This module prevents such errors by enforcing
biological constraints based on sample type.
"""

import yaml
from pathlib import Path
import logging
from typing import Dict, List, Optional, Any

logger = logging.getLogger(__name__)


class BiologicalValidator:
    """
    Validates cell type annotations against biological context.

    Sample types:
        - organoid: Only epithelial/stem cells (no immune/stromal)
        - primary_tissue: All cell types allowed
        - cell_line: Homogeneous population expected
        - pbmc: Only immune cells
    """

    def __init__(self, config_path: Optional[Path] = None):
        """
        Initialize the validator with configuration.

        Args:
            config_path: Path to sample_types.yaml. If None, uses default location.
        """
        if config_path is None:
            config_path = Path(__file__).parent.parent / "data" / "sample_types.yaml"
        self.config_path = config_path
        self.config = self._load_config(config_path)

    def _load_config(self, path: Path) -> Dict:
        """Load configuration from YAML file."""
        if not path.exists():
            logger.warning(f"Config file not found: {path}. Using defaults.")
            return self._get_default_config()

        with open(path) as f:
            return yaml.safe_load(f)

    def _get_default_config(self) -> Dict:
        """Return default configuration if file not found."""
        return {
            'sample_types': {
                'organoid': {
                    'forbidden_cell_types': ['immune', 'stromal', 'endothelial'],
                    'validation_warning': 'Organoids should NOT contain immune/stromal cells'
                },
                'primary_tissue': {
                    'forbidden_cell_types': [],
                    'validation_warning': None
                }
            },
            'immune_keywords': ['neutrophil', 'monocyte', 'macrophage', 't_cell', 'b_cell', 'nk'],
            'epithelial_keywords': ['basal', 'ciliated', 'secretory', 'club', 'goblet', 'epithelial']
        }

    def validate_annotation(self, adata, sample_type: str, celltype_col: str = 'celltype') -> Dict:
        """
        Validate cell type annotations against sample type constraints.

        Args:
            adata: AnnData object with cell type annotations
            sample_type: Type of sample ('organoid', 'primary_tissue', 'cell_line', 'pbmc')
            celltype_col: Column name containing cell type annotations

        Returns:
            Dict with keys:
                - valid: bool, whether annotation passes validation
                - issues: List of critical issues (>5% forbidden cell types)
                - warnings: List of minor warnings (<5% forbidden cell types)
                - sample_type: The sample type used for validation
                - summary: Human-readable summary
        """
        if sample_type not in self.config.get('sample_types', {}):
            return {
                'valid': True,
                'issues': [],
                'warnings': [f"Unknown sample type: {sample_type}. Skipping validation."],
                'sample_type': sample_type,
                'summary': f"Unknown sample type '{sample_type}' - validation skipped"
            }

        constraints = self.config['sample_types'][sample_type]
        issues = []
        warnings = []

        # Check for forbidden cell types
        if celltype_col in adata.obs.columns:
            celltype_counts = adata.obs[celltype_col].value_counts()
            forbidden_types = constraints.get('forbidden_cell_types', [])

            for forbidden in forbidden_types:
                # Find cell types matching forbidden keyword
                matching = [ct for ct in celltype_counts.index
                           if forbidden.lower() in ct.lower()]

                for ct in matching:
                    count = celltype_counts[ct]
                    pct = count / len(adata) * 100

                    if pct > 5:  # >5% is a critical issue
                        issues.append({
                            'type': 'forbidden_celltype',
                            'celltype': ct,
                            'count': int(count),
                            'percentage': round(pct, 1),
                            'forbidden_category': forbidden,
                            'message': f"{ct}: {count} cells ({pct:.1f}%) - UNEXPECTED in {sample_type}!"
                        })
                        logger.warning(f"Biological validation issue: {ct} ({pct:.1f}%) unexpected in {sample_type}")
                    elif pct > 0:  # <5% is a warning
                        warnings.append({
                            'celltype': ct,
                            'count': int(count),
                            'percentage': round(pct, 1),
                            'message': f"{ct}: {pct:.1f}% - unusual for {sample_type}, possible contamination"
                        })

            # Check using immune/epithelial keywords
            if sample_type == 'organoid':
                immune_keywords = self.config.get('immune_keywords', [])
                for ct in celltype_counts.index:
                    ct_lower = ct.lower()
                    for kw in immune_keywords:
                        if kw in ct_lower:
                            pct = celltype_counts[ct] / len(adata) * 100
                            if pct > 5 and ct not in [i['celltype'] for i in issues]:
                                issues.append({
                                    'type': 'immune_in_organoid',
                                    'celltype': ct,
                                    'percentage': round(pct, 1),
                                    'message': f"{ct}: {pct:.1f}% - Immune cells UNEXPECTED in organoid!"
                                })

        is_valid = len(issues) == 0

        # Generate summary
        if is_valid:
            summary = f"Annotation validated for {sample_type}: No biological constraint violations"
        else:
            issue_summary = ", ".join([f"{i['celltype']} ({i['percentage']}%)" for i in issues])
            summary = f"Validation FAILED for {sample_type}: Unexpected cell types - {issue_summary}"

        return {
            'valid': is_valid,
            'issues': issues,
            'warnings': [w['message'] if isinstance(w, dict) else w for w in warnings],
            'sample_type': sample_type,
            'summary': summary
        }

    def get_recommended_celltypist_model(self, tissue_type: str) -> str:
        """
        Get recommended CellTypist model for tissue type.

        Args:
            tissue_type: Type of tissue ('lung_airway', 'pbmc', 'gut', etc.)

        Returns:
            Model filename (e.g., 'Cells_Lung_Airway.pkl')
        """
        models = self.config.get('celltypist_models', {})
        if tissue_type in models:
            model = models[tissue_type]['model']
            logger.info(f"Recommended CellTypist model for {tissue_type}: {model}")
            return model

        # Default fallback
        default = "Immune_All_Low.pkl"
        logger.info(f"No specific model for {tissue_type}, using default: {default}")
        return default

    def check_ambiguous_markers(self, marker: str, sample_type: str) -> Dict:
        """
        Check if a marker has context-dependent interpretation.

        Key example: S100A8/S100A9 are commonly interpreted as neutrophil markers,
        but in organoid samples they indicate epithelial stress response.

        Args:
            marker: Gene symbol (e.g., 'S100A8')
            sample_type: Type of sample for context

        Returns:
            Dict with keys:
                - is_ambiguous: bool
                - common_interpretation: Default interpretation
                - context_interpretation: Context-specific interpretation
                - recommendation: Action to take
        """
        ambiguous = self.config.get('ambiguous_markers', {})

        if marker in ambiguous:
            info = ambiguous[marker]
            context_interp = info.get('context_dependent', {}).get(sample_type)

            return {
                'is_ambiguous': True,
                'marker': marker,
                'common_interpretation': info.get('common_interpretation'),
                'context_interpretation': context_interp,
                'recommendation': info.get('recommendation'),
                'warning': f"Marker {marker} is ambiguous in {sample_type} context!"
            }

        return {'is_ambiguous': False, 'marker': marker}

    def suggest_reannotation(self, adata, sample_type: str, celltype_col: str = 'celltype') -> Dict:
        """
        Suggest whether reannotation is needed based on validation results.

        Args:
            adata: AnnData object
            sample_type: Type of sample
            celltype_col: Cell type column name

        Returns:
            Dict with recommendation and details
        """
        validation = self.validate_annotation(adata, sample_type, celltype_col)

        if validation['valid']:
            return {
                'needs_reannotation': False,
                'reason': 'Annotation passes biological validation',
                'validation': validation
            }

        # Calculate severity
        total_bad_pct = sum(i['percentage'] for i in validation['issues'])

        if total_bad_pct > 50:
            severity = 'critical'
            action = 'Strongly recommended to reannotate using marker-based method for organoids'
        elif total_bad_pct > 20:
            severity = 'high'
            action = 'Consider reannotation or manual review of suspicious clusters'
        else:
            severity = 'moderate'
            action = 'Review suspicious clusters, may be contamination or doublets'

        return {
            'needs_reannotation': True,
            'severity': severity,
            'bad_cell_percentage': round(total_bad_pct, 1),
            'action': action,
            'validation': validation
        }

    def get_sample_type_info(self, sample_type: str) -> Dict:
        """Get information about a sample type and its constraints."""
        if sample_type not in self.config.get('sample_types', {}):
            return {'error': f'Unknown sample type: {sample_type}'}

        info = self.config['sample_types'][sample_type]
        return {
            'sample_type': sample_type,
            'description': info.get('description', ''),
            'allowed_cell_types': info.get('allowed_cell_types', []),
            'forbidden_cell_types': info.get('forbidden_cell_types', []),
            'validation_warning': info.get('validation_warning')
        }

    def list_sample_types(self) -> List[str]:
        """List all available sample types."""
        return list(self.config.get('sample_types', {}).keys())
