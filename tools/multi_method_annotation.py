"""
Multi-method cell type annotation with consensus decision.
Combines CellTypist (ML-based) + marker-based annotation with biological validation.

This approach was developed during JMV_ALI analysis where single-method annotation
(using only marker-based scoring) incorrectly identified 80.6% of cells as neutrophils.
The consensus approach with biological context awareness correctly identified the
cells as epithelial (Basal, Ciliated, Secretory/Club).
"""

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

logger = logging.getLogger(__name__)

# Default airway epithelial markers
DEFAULT_AIRWAY_MARKERS = {
    'Basal': ['KRT5', 'TP63', 'KRT15', 'KRT14', 'NGFR', 'KRT17'],
    'Ciliated': ['FOXJ1', 'DNAH5', 'TPPP3', 'CAPS', 'SNTN', 'PIFO'],
    'Secretory_Club': ['SCGB1A1', 'SCGB3A2', 'CYP2F1', 'LYPD2'],
    'Goblet': ['MUC5AC', 'MUC5B', 'TFF3', 'SPDEF', 'FCGBP'],
    'AT2_like': ['SFTPC', 'SFTPA1', 'SFTPB', 'LAMP3', 'ABCA3'],
    'Ionocyte': ['FOXI1', 'CFTR', 'ATP6V1C2', 'ATP6V0D2'],
    'General_Epithelial': ['EPCAM', 'KRT18', 'KRT8', 'KRT19', 'CDH1']
}


class MultiMethodAnnotator:
    """
    Multi-method cell type annotation with consensus decision logic.

    Combines:
    1. CellTypist (ML-based): Uses pre-trained models for cell type prediction
    2. Marker-based: Scores cells based on canonical marker expression
    3. Consensus: Combines both with biological context awareness

    For organoid samples, immune cell annotations from CellTypist are overridden
    with marker-based epithelial annotations.
    """

    def __init__(self, sample_type: str = 'primary_tissue'):
        """
        Initialize the multi-method annotator.

        Args:
            sample_type: Type of sample ('organoid', 'primary_tissue', 'cell_line', 'pbmc')
                        This affects the consensus decision logic.
        """
        self.sample_type = sample_type

        # Lazy import of BiologicalValidator to avoid circular imports
        try:
            from .biological_validator import BiologicalValidator
            self.validator = BiologicalValidator()
        except ImportError:
            self.validator = None
            logger.warning("BiologicalValidator not available")

        # Track results
        self.celltypist_results = None
        self.marker_results = None
        self.consensus_results = None

    def run_celltypist(self, adata, model_name: Optional[str] = None) -> sc.AnnData:
        """
        Run CellTypist annotation using specified or auto-selected model.

        Args:
            adata: AnnData object (should be log-normalized)
            model_name: CellTypist model name. If None, auto-selects based on tissue.

        Returns:
            AnnData with CellTypist annotations in obs
        """
        try:
            import celltypist
            from celltypist import models
        except ImportError:
            logger.error("celltypist not installed. Run: pip install celltypist")
            raise

        # Auto-select model based on sample type
        if model_name is None:
            if self.validator:
                model_name = self.validator.get_recommended_celltypist_model(
                    self._infer_tissue_type(adata)
                )
            else:
                model_name = 'Cells_Lung_Airway.pkl'

        logger.info(f"Running CellTypist with model: {model_name}")

        # Load model
        try:
            model = models.Model.load(model=model_name)
        except Exception as e:
            logger.info(f"Downloading model {model_name}...")
            model_base = model_name.replace('.pkl', '')
            models.download_models(model=[model_base])
            model = models.Model.load(model=model_name)

        logger.info(f"Model loaded: {len(model.cell_types)} cell types available")

        # Run prediction with majority voting for robustness
        predictions = celltypist.annotate(
            adata,
            model=model,
            majority_voting=True
        )

        adata_ct = predictions.to_adata()

        # Log summary
        ct_counts = adata_ct.obs['majority_voting'].value_counts()
        logger.info(f"CellTypist results ({len(ct_counts)} cell types):")
        for ct, count in ct_counts.head(5).items():
            pct = count / len(adata_ct) * 100
            logger.info(f"  {ct}: {count} cells ({pct:.1f}%)")

        self.celltypist_results = adata_ct
        return adata_ct

    def run_marker_based(self, adata, markers_dict: Optional[Dict] = None) -> Tuple[pd.Series, pd.Series]:
        """
        Run marker-based annotation by scoring each cell type.

        Args:
            adata: AnnData object
            markers_dict: Dictionary of cell type -> marker genes list

        Returns:
            Tuple of (cell_types, confidence_scores) as pandas Series
        """
        if markers_dict is None:
            markers_dict = DEFAULT_AIRWAY_MARKERS

        logger.info("Running marker-based annotation")

        scores = {}
        markers_found = {}

        for celltype, markers in markers_dict.items():
            present = [m for m in markers if m in adata.var_names]
            markers_found[celltype] = len(present)

            if not present:
                scores[celltype] = np.zeros(adata.n_obs)
                logger.debug(f"  {celltype}: 0/{len(markers)} markers found")
                continue

            logger.debug(f"  {celltype}: {len(present)}/{len(markers)} markers found")

            # Get expression matrix
            expr = adata[:, present].X
            if issparse(expr):
                expr = expr.toarray()

            # Mean expression across markers
            scores[celltype] = expr.mean(axis=1)

        # Create DataFrame and assign cell types
        scores_df = pd.DataFrame(scores, index=adata.obs_names)
        cell_types = scores_df.idxmax(axis=1)
        confidence = scores_df.max(axis=1)

        # Log summary
        type_counts = cell_types.value_counts()
        logger.info(f"Marker-based results ({len(type_counts)} cell types):")
        for ct, count in type_counts.items():
            pct = count / len(adata) * 100
            avg_conf = confidence[cell_types == ct].mean()
            logger.info(f"  {ct}: {count} cells ({pct:.1f}%), avg score={avg_conf:.3f}")

        self.marker_results = (cell_types, confidence, scores_df)
        return cell_types, confidence

    def consensus_annotation(
        self,
        adata,
        celltypist_results: Optional[sc.AnnData] = None,
        marker_types: Optional[pd.Series] = None,
        marker_conf: Optional[pd.Series] = None
    ) -> Tuple[List[str], List[str], List[str]]:
        """
        Create consensus annotation combining CellTypist and marker-based results.

        The decision logic depends on sample_type:
        - organoid: Prioritize marker-based, reject immune annotations
        - primary_tissue: Standard consensus favoring high-confidence results
        - pbmc: Prioritize CellTypist for immune cells

        Args:
            adata: Original AnnData
            celltypist_results: AnnData from run_celltypist()
            marker_types: Series from run_marker_based()
            marker_conf: Series from run_marker_based()

        Returns:
            Tuple of (consensus_types, confidence_levels, agreement_types)
        """
        # Use stored results if not provided
        if celltypist_results is None:
            celltypist_results = self.celltypist_results
        if marker_types is None and self.marker_results:
            marker_types, marker_conf, _ = self.marker_results

        if celltypist_results is None or marker_types is None:
            raise ValueError("Must run both CellTypist and marker-based annotation first")

        logger.info(f"Creating consensus annotation for sample_type={self.sample_type}")

        consensus_types = []
        confidence_levels = []
        agreement_types = []

        ct_types = celltypist_results.obs['majority_voting']
        ct_conf = celltypist_results.obs['conf_score']

        # Immune keywords for organoid filtering
        immune_keywords = ['neutrophil', 'monocyte', 'macrophage', 't_cell', 'b_cell',
                          'nk', 'dendritic', 'mast', 'plasma', 'lymphocyte']

        for idx in range(len(adata)):
            ct_type = ct_types.iloc[idx]
            ct_c = ct_conf.iloc[idx]
            mb_type = marker_types.iloc[idx]
            mb_c = marker_conf.iloc[idx]

            # Apply decision logic based on sample type
            if self.sample_type == 'organoid':
                final_type, final_conf, agree = self._organoid_decision(
                    ct_type, ct_c, mb_type, mb_c, immune_keywords
                )
            elif self.sample_type == 'pbmc':
                final_type, final_conf, agree = self._pbmc_decision(
                    ct_type, ct_c, mb_type, mb_c
                )
            else:
                final_type, final_conf, agree = self._standard_decision(
                    ct_type, ct_c, mb_type, mb_c
                )

            consensus_types.append(final_type)
            confidence_levels.append(final_conf)
            agreement_types.append(agree)

        # Log summary
        logger.info("Consensus decision summary:")
        agree_counts = pd.Series(agreement_types).value_counts()
        for decision, count in agree_counts.items():
            pct = count / len(adata) * 100
            logger.info(f"  {decision}: {count} cells ({pct:.1f}%)")

        logger.info("Final consensus cell types:")
        consensus_counts = pd.Series(consensus_types).value_counts()
        for ct, count in consensus_counts.head(10).items():
            pct = count / len(adata) * 100
            logger.info(f"  {ct}: {count} cells ({pct:.1f}%)")

        self.consensus_results = (consensus_types, confidence_levels, agreement_types)
        return consensus_types, confidence_levels, agreement_types

    def _organoid_decision(self, ct_type, ct_conf, mb_type, mb_conf, immune_keywords):
        """
        Decision logic for organoid samples.

        Key principle: Organoids should NOT have immune cells. If CellTypist
        predicts immune cell, override with marker-based epithelial annotation.
        """
        # Check if CellTypist predicted immune cell
        ct_is_immune = any(kw in ct_type.lower() for kw in immune_keywords)

        if ct_is_immune:
            # CellTypist says immune, but this is organoid - use marker-based
            return mb_type, "Medium", "Marker_Priority_Organoid"

        # Both methods agree
        if mb_type.lower() in ct_type.lower() or ct_type.lower() in mb_type.lower():
            return ct_type, "High", "Agree"

        # Marker-based has high confidence
        if mb_conf > 0.5:
            return mb_type, "Medium" if ct_conf < 0.5 else "High", "Marker_Priority"

        # CellTypist has high confidence and predicts epithelial
        epithelial_keywords = ['basal', 'ciliated', 'secretory', 'club', 'goblet', 'epithelial']
        ct_is_epithelial = any(kw in ct_type.lower() for kw in epithelial_keywords)
        if ct_conf > 0.7 and ct_is_epithelial:
            return ct_type, "High", "CellTypist_Priority"

        # Default to marker-based for organoids
        return mb_type, "Low", "Default_Marker"

    def _pbmc_decision(self, ct_type, ct_conf, mb_type, mb_conf):
        """Decision logic for PBMC samples - prioritize CellTypist for immune."""
        # For PBMCs, CellTypist is usually more reliable
        if ct_conf > 0.5:
            return ct_type, "High" if ct_conf > 0.7 else "Medium", "CellTypist_Priority"

        if mb_type.lower() in ct_type.lower() or ct_type.lower() in mb_type.lower():
            return ct_type, "High", "Agree"

        return ct_type, "Low", "Uncertain"

    def _standard_decision(self, ct_type, ct_conf, mb_type, mb_conf):
        """Standard decision logic for primary tissue."""
        # Methods agree
        if mb_type.lower() in ct_type.lower() or ct_type.lower() in mb_type.lower():
            return ct_type, "High", "Agree"

        # CellTypist has high confidence
        if ct_conf > 0.7:
            return ct_type, "High", "CellTypist_Priority"

        # Marker-based has high confidence
        if mb_conf > 0.5:
            return mb_type, "Medium", "Marker_Priority"

        # Default to CellTypist
        return ct_type, "Low", "Uncertain"

    def _infer_tissue_type(self, adata) -> str:
        """Infer tissue type from adata metadata."""
        # Check common metadata fields
        for col in ['tissue', 'tissue_type', 'organ']:
            if col in adata.obs.columns:
                tissue = str(adata.obs[col].iloc[0]).lower()
                if 'lung' in tissue or 'airway' in tissue:
                    return 'lung_airway'
                elif 'blood' in tissue or 'pbmc' in tissue:
                    return 'pbmc'
                elif 'gut' in tissue or 'intestin' in tissue:
                    return 'gut'

        # Default based on sample_type
        if self.sample_type == 'pbmc':
            return 'pbmc'

        return 'lung_airway'

    def annotate(
        self,
        adata,
        markers_dict: Optional[Dict] = None,
        celltypist_model: Optional[str] = None,
        add_to_adata: bool = True
    ) -> Dict:
        """
        Run full multi-method annotation pipeline.

        Args:
            adata: AnnData object (should be log-normalized)
            markers_dict: Custom markers dictionary
            celltypist_model: CellTypist model name
            add_to_adata: If True, add results to adata.obs

        Returns:
            Dict with all annotation results and validation
        """
        logger.info("=" * 60)
        logger.info(f"MULTI-METHOD ANNOTATION (sample_type={self.sample_type})")
        logger.info("=" * 60)

        # Run CellTypist
        adata_ct = self.run_celltypist(adata, model_name=celltypist_model)

        # Run marker-based
        marker_types, marker_conf = self.run_marker_based(adata, markers_dict)

        # Create consensus
        consensus, conf, agreement = self.consensus_annotation(
            adata, adata_ct, marker_types, marker_conf
        )

        # Add to adata if requested
        if add_to_adata:
            adata.obs['celltypist'] = adata_ct.obs['majority_voting'].values
            adata.obs['celltypist_conf'] = adata_ct.obs['conf_score'].values
            adata.obs['marker_based'] = marker_types.values
            adata.obs['marker_conf'] = marker_conf.values
            adata.obs['consensus_celltype'] = consensus
            adata.obs['annotation_confidence'] = conf
            adata.obs['annotation_agreement'] = agreement

        # Validate results
        validation = None
        if self.validator:
            adata.obs['celltype'] = consensus  # Temporarily set for validation
            validation = self.validator.validate_annotation(adata, self.sample_type)
            if not validation['valid']:
                logger.warning("Biological validation issues detected:")
                for issue in validation['issues']:
                    logger.warning(f"  - {issue['message']}")

        return {
            'celltypist_adata': adata_ct,
            'marker_types': marker_types,
            'marker_confidence': marker_conf,
            'consensus_types': consensus,
            'consensus_confidence': conf,
            'agreement_types': agreement,
            'validation': validation
        }


def run_multi_method_annotation(
    adata,
    sample_type: str = 'primary_tissue',
    markers_dict: Optional[Dict] = None,
    celltypist_model: Optional[str] = None,
    output_dir: Optional[Path] = None
) -> sc.AnnData:
    """
    Convenience function to run multi-method annotation.

    Args:
        adata: AnnData object
        sample_type: 'organoid', 'primary_tissue', 'cell_line', or 'pbmc'
        markers_dict: Custom marker dictionary
        celltypist_model: CellTypist model name
        output_dir: Optional output directory for saving results

    Returns:
        AnnData with annotations added to obs
    """
    annotator = MultiMethodAnnotator(sample_type=sample_type)
    results = annotator.annotate(
        adata,
        markers_dict=markers_dict,
        celltypist_model=celltypist_model,
        add_to_adata=True
    )

    # Save comparison CSV if output_dir provided
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)

        comparison_df = pd.DataFrame({
            'cell_id': adata.obs_names,
            'celltypist': adata.obs['celltypist'],
            'celltypist_conf': adata.obs['celltypist_conf'],
            'marker_based': adata.obs['marker_based'],
            'marker_conf': adata.obs['marker_conf'],
            'consensus': adata.obs['consensus_celltype'],
            'confidence': adata.obs['annotation_confidence'],
            'agreement': adata.obs['annotation_agreement']
        })
        comparison_df.to_csv(output_dir / 'multi_method_comparison.csv', index=False)
        logger.info(f"Saved comparison to {output_dir / 'multi_method_comparison.csv'}")

    return adata
