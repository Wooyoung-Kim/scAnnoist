"""
Tests for scrna-agent pipeline
"""

import os
import tempfile
import pytest
import numpy as np
import pandas as pd
import scanpy as sc

from scrna_agent.pipeline import (
    PipelineConfig,
    create_toy_dataset,
    load_data,
    run_qc,
    normalize_and_scale,
    run_dimensionality_reduction,
    run_clustering,
    find_markers,
    run_pipeline,
)


@pytest.fixture
def toy_adata():
    """Create a toy AnnData object for testing."""
    return create_toy_dataset(n_cells=200, n_genes=500, n_clusters=3, seed=42)


@pytest.fixture
def config():
    """Create a default config for testing."""
    return PipelineConfig(
        min_genes=50,  # Lower threshold for toy data
        min_cells=1,
        max_mt_pct=50.0,  # Higher threshold for synthetic data
        n_top_genes=200,
        n_pcs=20,
        resolution=0.3,
        seed=42
    )


class TestPipelineConfig:
    """Test PipelineConfig class."""

    def test_default_config(self):
        """Test default configuration values."""
        config = PipelineConfig()
        assert config.min_genes == 200
        assert config.min_cells == 3
        assert config.max_mt_pct == 20.0
        assert config.clustering_method == "leiden"
        assert config.seed == 42

    def test_config_yaml_roundtrip(self):
        """Test saving and loading config from YAML."""
        with tempfile.NamedTemporaryFile(suffix='.yaml', delete=False) as f:
            yaml_path = f.name

        try:
            config = PipelineConfig(min_genes=300, resolution=0.8)
            config.to_yaml(yaml_path)

            loaded_config = PipelineConfig.from_yaml(yaml_path)
            assert loaded_config.min_genes == 300
            assert loaded_config.resolution == 0.8
        finally:
            os.unlink(yaml_path)


class TestToyDataset:
    """Test toy dataset generation."""

    def test_create_toy_dataset_shape(self):
        """Test that toy dataset has correct shape."""
        adata = create_toy_dataset(n_cells=100, n_genes=200)
        assert adata.n_obs == 100
        assert adata.n_vars == 200

    def test_create_toy_dataset_clusters(self):
        """Test that toy dataset has cluster information."""
        adata = create_toy_dataset(n_cells=100, n_genes=200, n_clusters=4)
        assert 'true_cluster' in adata.obs.columns
        assert adata.obs['true_cluster'].nunique() == 4

    def test_create_toy_dataset_has_mt_genes(self):
        """Test that toy dataset has MT genes."""
        adata = create_toy_dataset(n_cells=100, n_genes=200)
        mt_genes = [g for g in adata.var_names if g.startswith('MT-')]
        assert len(mt_genes) > 0

    def test_create_toy_dataset_save(self):
        """Test saving toy dataset to file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "test.h5ad")
            adata = create_toy_dataset(n_cells=50, n_genes=100, output_path=output_path)
            assert os.path.exists(output_path)

            # Load and verify
            loaded = sc.read_h5ad(output_path)
            assert loaded.n_obs == 50
            assert loaded.n_vars == 100


class TestQC:
    """Test QC functions."""

    def test_run_qc_filters_cells(self, toy_adata, config):
        """Test that QC filtering removes cells."""
        n_before = toy_adata.n_obs
        qc_summary = run_qc(toy_adata, config)

        assert isinstance(qc_summary, pd.DataFrame)
        assert 'metric' in qc_summary.columns
        assert 'value' in qc_summary.columns

    def test_run_qc_adds_metrics(self, toy_adata, config):
        """Test that QC adds metrics to adata."""
        run_qc(toy_adata, config)
        assert 'n_genes_by_counts' in toy_adata.obs.columns
        assert 'total_counts' in toy_adata.obs.columns
        assert 'pct_counts_mt' in toy_adata.obs.columns


class TestNormalization:
    """Test normalization functions."""

    def test_normalize_and_scale(self, toy_adata, config):
        """Test normalization adds expected layers."""
        run_qc(toy_adata, config)
        normalize_and_scale(toy_adata, config)

        assert 'counts' in toy_adata.layers
        assert 'highly_variable' in toy_adata.var.columns
        assert toy_adata.raw is not None


class TestDimensionalityReduction:
    """Test dimensionality reduction."""

    def test_pca_umap(self, toy_adata, config):
        """Test PCA and UMAP computation."""
        run_qc(toy_adata, config)
        normalize_and_scale(toy_adata, config)
        run_dimensionality_reduction(toy_adata, config)

        assert 'X_pca' in toy_adata.obsm
        assert 'X_umap' in toy_adata.obsm
        assert toy_adata.obsm['X_pca'].shape[1] == config.n_pcs
        assert toy_adata.obsm['X_umap'].shape[1] == 2


class TestClustering:
    """Test clustering functions."""

    def test_leiden_clustering(self, toy_adata, config):
        """Test Leiden clustering."""
        run_qc(toy_adata, config)
        normalize_and_scale(toy_adata, config)
        run_dimensionality_reduction(toy_adata, config)

        clusters = run_clustering(toy_adata, config)

        assert isinstance(clusters, pd.DataFrame)
        assert 'cell_id' in clusters.columns
        assert 'cluster' in clusters.columns
        assert 'leiden' in toy_adata.obs.columns

    def test_louvain_clustering(self, toy_adata, config):
        """Test Louvain clustering."""
        config.clustering_method = "louvain"
        run_qc(toy_adata, config)
        normalize_and_scale(toy_adata, config)
        run_dimensionality_reduction(toy_adata, config)

        clusters = run_clustering(toy_adata, config)

        assert 'louvain' in toy_adata.obs.columns


class TestMarkers:
    """Test marker gene detection."""

    def test_find_markers(self, toy_adata, config):
        """Test marker gene detection."""
        run_qc(toy_adata, config)
        normalize_and_scale(toy_adata, config)
        run_dimensionality_reduction(toy_adata, config)
        run_clustering(toy_adata, config)

        markers, markers_full = find_markers(toy_adata, config)

        assert isinstance(markers, pd.DataFrame)
        assert 'cluster' in markers.columns
        assert 'gene' in markers.columns
        assert 'log2fc' in markers.columns
        assert 'pval_adj' in markers.columns

        # Check full markers table has additional columns
        assert isinstance(markers_full, pd.DataFrame)
        assert 'pct_in' in markers_full.columns
        assert 'pct_out' in markers_full.columns
        assert 'significant' in markers_full.columns


class TestEndToEnd:
    """End-to-end pipeline tests."""

    def test_full_pipeline(self):
        """Test running the complete pipeline on toy data."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create toy data
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            adata = create_toy_dataset(
                n_cells=200,
                n_genes=500,
                n_clusters=3,
                output_path=input_path
            )

            # Run pipeline with relaxed config for toy data
            config = PipelineConfig(
                min_genes=50,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=200,
                n_pcs=20,
                resolution=0.3,
                seed=42
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check all required outputs exist
            assert os.path.exists(os.path.join(output_dir, "qc_summary.csv"))
            assert os.path.exists(os.path.join(output_dir, "umap.png"))
            assert os.path.exists(os.path.join(output_dir, "clusters.csv"))
            assert os.path.exists(os.path.join(output_dir, "markers.csv"))
            assert os.path.exists(os.path.join(output_dir, "report.md"))
            assert os.path.exists(os.path.join(output_dir, "processed.h5ad"))

            # Check results object
            assert results.adata is not None
            assert results.qc_summary is not None
            assert results.clusters is not None
            assert results.markers is not None
            assert len(results.files_created) >= 6

    def test_pipeline_reproducibility(self):
        """Test that pipeline produces reproducible results with same seed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir1 = os.path.join(tmpdir, "output1")
            output_dir2 = os.path.join(tmpdir, "output2")

            create_toy_dataset(n_cells=100, n_genes=200, output_path=input_path, seed=42)

            config = PipelineConfig(
                min_genes=20,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=100,
                n_pcs=10,
                seed=42
            )

            results1 = run_pipeline(input_path, output_dir1, config)
            results2 = run_pipeline(input_path, output_dir2, config)

            # Check cluster assignments are identical
            clusters1 = pd.read_csv(os.path.join(output_dir1, "clusters.csv"))
            clusters2 = pd.read_csv(os.path.join(output_dir2, "clusters.csv"))

            assert clusters1['cluster'].equals(clusters2['cluster'])


class TestDataLoading:
    """Test data loading functions."""

    def test_load_h5ad(self):
        """Test loading h5ad file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test.h5ad")
            adata = create_toy_dataset(n_cells=50, n_genes=100, output_path=path)

            loaded = load_data(path, input_format="h5ad")
            assert loaded.n_obs == 50
            assert loaded.n_vars == 100

    def test_load_auto_detect(self):
        """Test auto-detection of h5ad format."""
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "test.h5ad")
            create_toy_dataset(n_cells=50, n_genes=100, output_path=path)

            loaded = load_data(path, input_format="auto")
            assert loaded.n_obs == 50


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
