"""
Tests for statistical methods in the scRNA-seq pipeline.
"""

import os
import tempfile
import pytest
import numpy as np
import pandas as pd
import scanpy as sc

from scrna_agent.pipeline import create_toy_dataset, PipelineConfig, run_pipeline
from scrna_agent.stats import (
    compute_mad,
    mad_to_std,
    compute_adaptive_qc_thresholds,
    apply_adaptive_qc,
    select_pcs_elbow,
    select_pcs_variance,
    compute_pc_selection,
    compute_qc_distributions,
    run_full_deg_analysis,
    QCThresholds,
    PCSelectionResult,
)


class TestMADComputation:
    """Test MAD (Median Absolute Deviation) computation."""

    def test_mad_normal_distribution(self):
        """Test MAD on normally distributed data."""
        np.random.seed(42)
        # For normal distribution, MAD ≈ 0.6745 * σ
        data = np.random.normal(0, 1, 10000)
        mad = compute_mad(data)
        # MAD should be approximately 0.67 for standard normal
        assert 0.6 < mad < 0.75

    def test_mad_with_outliers(self):
        """Test that MAD is robust to outliers."""
        data = np.array([1, 2, 3, 4, 5, 100])  # 100 is outlier
        mad = compute_mad(data)
        # MAD should be close to 1.5 (median of |x - 3.5| = |[2.5,1.5,0.5,0.5,1.5,96.5]| = 1.5)
        assert 0.5 <= mad <= 2.0

    def test_mad_to_std_conversion(self):
        """Test MAD to standard deviation conversion."""
        # For normal distribution: σ ≈ 1.4826 * MAD
        mad = 0.6745
        std = mad_to_std(mad)
        assert 0.9 < std < 1.1  # Should be close to 1


class TestAdaptiveQCThresholds:
    """Test adaptive QC threshold computation."""

    @pytest.fixture
    def toy_adata_with_qc(self):
        """Create toy data with QC metrics computed."""
        adata = create_toy_dataset(n_cells=500, n_genes=1000, seed=42)
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        return adata

    def test_compute_adaptive_thresholds_mad(self, toy_adata_with_qc):
        """Test MAD-based threshold computation."""
        thresholds = compute_adaptive_qc_thresholds(
            toy_adata_with_qc,
            nmads=3.0,
            method='mad'
        )

        assert isinstance(thresholds, QCThresholds)
        assert thresholds.min_genes > 0
        assert thresholds.max_genes > thresholds.min_genes
        assert thresholds.max_mt_pct > 0
        assert thresholds.method == 'mad'
        assert thresholds.nmads == 3.0

    def test_compute_adaptive_thresholds_quantile(self, toy_adata_with_qc):
        """Test quantile-based threshold computation."""
        thresholds = compute_adaptive_qc_thresholds(
            toy_adata_with_qc,
            method='quantile'
        )

        assert isinstance(thresholds, QCThresholds)
        assert thresholds.method == 'quantile'

    def test_apply_adaptive_qc(self, toy_adata_with_qc):
        """Test applying adaptive QC filters."""
        n_cells_before = toy_adata_with_qc.n_obs

        thresholds = compute_adaptive_qc_thresholds(toy_adata_with_qc)
        stats = apply_adaptive_qc(toy_adata_with_qc, thresholds, min_cells=3)

        assert stats['cells_before'] == n_cells_before
        assert stats['cells_after'] <= n_cells_before
        assert toy_adata_with_qc.n_obs == stats['cells_after']


class TestPCSelection:
    """Test PC selection methods."""

    @pytest.fixture
    def variance_ratio(self):
        """Create mock variance ratio array."""
        # Simulate typical PCA variance explained
        # First few PCs explain most variance, then tapers off
        n_pcs = 50
        ratios = np.exp(-np.linspace(0, 5, n_pcs))
        ratios = ratios / ratios.sum()
        return ratios

    def test_select_pcs_elbow(self, variance_ratio):
        """Test elbow method for PC selection."""
        n_pcs = select_pcs_elbow(variance_ratio, min_pcs=5, max_pcs=50)

        assert 5 <= n_pcs <= 50
        assert isinstance(n_pcs, (int, np.integer))

    def test_select_pcs_variance(self, variance_ratio):
        """Test variance threshold method."""
        # Select PCs explaining 80% variance
        n_pcs = select_pcs_variance(
            variance_ratio,
            variance_threshold=0.80,
            min_pcs=5,
            max_pcs=50
        )

        # Verify we reach threshold
        cum_var = np.cumsum(variance_ratio)
        assert cum_var[n_pcs - 1] >= 0.80 or n_pcs == 50

    def test_pc_selection_result(self):
        """Test PCSelectionResult dataclass."""
        result = PCSelectionResult(
            n_pcs_selected=20,
            method='elbow',
            variance_explained=[0.1, 0.08, 0.06],
            cumulative_variance=[0.1, 0.18, 0.24],
            elbow_point=15
        )

        assert result.n_pcs_selected == 20
        result_dict = result.to_dict()
        assert 'n_pcs_selected' in result_dict
        assert 'method' in result_dict


class TestQCDistributions:
    """Test QC distribution statistics."""

    @pytest.fixture
    def toy_adata_with_qc(self):
        """Create toy data with QC metrics."""
        adata = create_toy_dataset(n_cells=500, n_genes=1000, seed=42)
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        return adata

    def test_compute_distributions(self, toy_adata_with_qc):
        """Test distribution statistics computation."""
        distributions = compute_qc_distributions(toy_adata_with_qc)

        assert 'n_genes_by_counts' in distributions
        assert 'total_counts' in distributions
        assert 'pct_counts_mt' in distributions

        for metric, stats in distributions.items():
            assert 'mean' in stats
            assert 'median' in stats
            assert 'std' in stats
            assert 'mad' in stats
            assert 'q05' in stats
            assert 'q95' in stats
            assert 'skewness' in stats
            assert 'kurtosis' in stats


class TestDEGAnalysis:
    """Test differential expression analysis."""

    @pytest.fixture
    def clustered_adata(self):
        """Create toy data with clusters."""
        adata = create_toy_dataset(n_cells=300, n_genes=500, n_clusters=3, seed=42)
        # Normalize
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        # Add cluster column
        adata.obs['cluster'] = adata.obs['true_cluster']
        return adata

    def test_full_deg_analysis(self, clustered_adata):
        """Test full DEG analysis with proper statistics."""
        deg_df = run_full_deg_analysis(
            clustered_adata,
            groupby='cluster',
            method='wilcoxon',
            min_fold_change=0.25
        )

        assert isinstance(deg_df, pd.DataFrame)
        assert 'cluster' in deg_df.columns
        assert 'gene' in deg_df.columns
        assert 'log2fc' in deg_df.columns
        assert 'pval' in deg_df.columns
        assert 'pval_adj' in deg_df.columns
        assert 'significant' in deg_df.columns
        assert 'pct_in' in deg_df.columns
        assert 'pct_out' in deg_df.columns

    def test_deg_multiple_testing_correction(self, clustered_adata):
        """Test that multiple testing correction is applied."""
        deg_df = run_full_deg_analysis(
            clustered_adata,
            groupby='cluster',
            method='wilcoxon'
        )

        # Adjusted p-values should be >= raw p-values
        valid_mask = (deg_df['pval'] > 0) & (deg_df['pval_adj'] > 0)
        assert (deg_df.loc[valid_mask, 'pval_adj'] >= deg_df.loc[valid_mask, 'pval']).all()


class TestEndToEndStatistical:
    """End-to-end tests with statistical features."""

    def test_pipeline_with_adaptive_qc(self):
        """Test full pipeline with adaptive QC."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            create_toy_dataset(
                n_cells=300,
                n_genes=500,
                n_clusters=3,
                output_path=input_path
            )

            config = PipelineConfig(
                adaptive_qc=True,
                qc_method='mad',
                qc_nmads=3.0,
                min_genes=50,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=200,
                n_pcs=30,
                pc_selection_method='elbow',
                resolution=0.3,
                seed=42
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check outputs exist
            assert os.path.exists(os.path.join(output_dir, "run_metadata.json"))
            assert os.path.exists(os.path.join(output_dir, "hvg_list.txt"))
            assert os.path.exists(os.path.join(output_dir, "markers_full.csv"))
            assert os.path.exists(os.path.join(output_dir, "pca_variance.png"))

            # Check run_metadata content
            import json
            with open(os.path.join(output_dir, "run_metadata.json")) as f:
                metadata = json.load(f)

            assert 'summary' in metadata
            assert 'qc_thresholds' in metadata or 'cells_after_qc' in metadata

    def test_pipeline_variance_ratio_pc_selection(self):
        """Test pipeline with variance ratio PC selection."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            create_toy_dataset(
                n_cells=200,
                n_genes=400,
                output_path=input_path
            )

            config = PipelineConfig(
                adaptive_qc=False,  # Use fixed thresholds
                min_genes=20,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=150,
                n_pcs=30,
                pc_selection_method='variance_ratio',
                variance_threshold=0.85,
                resolution=0.3,
                seed=42
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check PC selection was recorded
            assert 'pc_selection' in results.adata.uns
            assert results.adata.uns['pc_selection']['method'] == 'variance_ratio'

    def test_full_deg_table_saved(self):
        """Test that full DEG table is saved with all statistics."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            create_toy_dataset(
                n_cells=200,
                n_genes=300,
                n_clusters=3,
                output_path=input_path
            )

            config = PipelineConfig(
                min_genes=20,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=100,
                n_pcs=20,
                resolution=0.3,
                deg_method='wilcoxon',
                seed=42
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check full DEG table
            markers_full = pd.read_csv(os.path.join(output_dir, "markers_full.csv"))
            assert 'pct_in' in markers_full.columns
            assert 'pct_out' in markers_full.columns
            assert 'significant' in markers_full.columns

            # Should have many more genes than top markers
            markers_top = pd.read_csv(os.path.join(output_dir, "markers.csv"))
            assert len(markers_full) > len(markers_top)


class TestStatisticalRationale:
    """Test that statistical rationale is properly generated."""

    def test_qc_rationale_in_report(self):
        """Test that QC rationale appears in report."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            create_toy_dataset(
                n_cells=200,
                n_genes=300,
                output_path=input_path
            )

            config = PipelineConfig(
                adaptive_qc=True,
                min_genes=20,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=100,
                n_pcs=20,
                seed=42
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check report contains statistical methods section
            with open(os.path.join(output_dir, "report.md")) as f:
                report = f.read()

            assert "Statistical Methods" in report or "QC Threshold Rationale" in report
            assert "MAD" in report or "Median Absolute Deviation" in report


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
