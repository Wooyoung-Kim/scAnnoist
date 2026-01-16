"""
Tests for tissue handling and marker-based annotation.
"""

import os
import tempfile
import pytest
import numpy as np
import pandas as pd
import scanpy as sc

from scrna_agent.pipeline import create_toy_dataset, PipelineConfig, run_pipeline
from scrna_agent.tissue import (
    TissueMapError,
    load_tissue_map,
    assign_tissue,
    load_markers,
    get_tissue_markers,
    score_cell_types,
    annotate_cells,
    generate_tissue_summary,
    generate_cell_metadata,
)


@pytest.fixture
def toy_adata():
    """Create a toy AnnData object with sample column."""
    adata = create_toy_dataset(n_cells=200, n_genes=500, n_clusters=3, seed=42)
    # Add sample column
    adata.obs['sample'] = np.random.choice(['sample_A', 'sample_B'], size=adata.n_obs)
    return adata


@pytest.fixture
def tissue_map_csv(tmp_path):
    """Create a temporary tissue map CSV file."""
    csv_path = tmp_path / "tissue_map.csv"
    df = pd.DataFrame({
        'sample': ['sample_A', 'sample_B'],
        'tissue': ['lung', 'liver']
    })
    df.to_csv(csv_path, index=False)
    return str(csv_path)


class TestLoadTissueMap:
    """Test tissue map loading."""

    def test_load_valid_tissue_map(self, tissue_map_csv):
        """Test loading a valid tissue map CSV."""
        tissue_map = load_tissue_map(tissue_map_csv)
        assert tissue_map == {'sample_A': 'lung', 'sample_B': 'liver'}

    def test_load_tissue_map_custom_sample_col(self, tmp_path):
        """Test loading tissue map with custom sample column."""
        csv_path = tmp_path / "tissue_map.csv"
        df = pd.DataFrame({
            'sample_id': ['S1', 'S2'],
            'tissue': ['brain', 'heart']
        })
        df.to_csv(csv_path, index=False)

        tissue_map = load_tissue_map(str(csv_path), sample_col='sample_id')
        assert tissue_map == {'S1': 'brain', 'S2': 'heart'}

    def test_load_tissue_map_missing_sample_col(self, tmp_path):
        """Test error when sample column is missing."""
        csv_path = tmp_path / "tissue_map.csv"
        df = pd.DataFrame({
            'wrong_col': ['S1', 'S2'],
            'tissue': ['brain', 'heart']
        })
        df.to_csv(csv_path, index=False)

        with pytest.raises(TissueMapError, match="must have 'sample' column"):
            load_tissue_map(str(csv_path))

    def test_load_tissue_map_missing_tissue_col(self, tmp_path):
        """Test error when tissue column is missing."""
        csv_path = tmp_path / "tissue_map.csv"
        df = pd.DataFrame({
            'sample': ['S1', 'S2'],
            'wrong_col': ['brain', 'heart']
        })
        df.to_csv(csv_path, index=False)

        with pytest.raises(TissueMapError, match="must have 'tissue' column"):
            load_tissue_map(str(csv_path))


class TestAssignTissue:
    """Test tissue assignment to AnnData."""

    def test_assign_single_tissue(self, toy_adata):
        """Test assigning single tissue to all cells."""
        adata = assign_tissue(toy_adata, tissue='lung')
        assert 'tissue' in adata.obs.columns
        assert (adata.obs['tissue'] == 'lung').all()

    def test_assign_tissue_from_map(self, toy_adata, tissue_map_csv):
        """Test assigning tissue from map."""
        tissue_map = load_tissue_map(tissue_map_csv)
        adata = assign_tissue(toy_adata, tissue_map=tissue_map, sample_col='sample')

        assert 'tissue' in adata.obs.columns
        # Check that sample_A cells have 'lung' tissue
        sample_a_mask = adata.obs['sample'] == 'sample_A'
        assert (adata.obs.loc[sample_a_mask, 'tissue'] == 'lung').all()
        # Check that sample_B cells have 'liver' tissue
        sample_b_mask = adata.obs['sample'] == 'sample_B'
        assert (adata.obs.loc[sample_b_mask, 'tissue'] == 'liver').all()

    def test_assign_tissue_missing_sample_col(self, toy_adata):
        """Test error when sample column is not in adata."""
        tissue_map = {'sample_A': 'lung'}
        adata = toy_adata.copy()
        del adata.obs['sample']

        with pytest.raises(TissueMapError, match="Sample column 'sample' not found"):
            assign_tissue(adata, tissue_map=tissue_map)

    def test_assign_tissue_missing_sample_in_map(self, toy_adata):
        """Test error when sample is missing from tissue map."""
        tissue_map = {'sample_A': 'lung'}  # Missing sample_B

        with pytest.raises(TissueMapError, match="Samples missing from tissue map"):
            assign_tissue(toy_adata, tissue_map=tissue_map)

    def test_assign_tissue_allow_missing(self, toy_adata):
        """Test allowing missing samples with unknown tissue."""
        tissue_map = {'sample_A': 'lung'}  # Missing sample_B

        adata = assign_tissue(
            toy_adata,
            tissue_map=tissue_map,
            allow_missing=True
        )

        # sample_B should be 'unknown'
        sample_b_mask = adata.obs['sample'] == 'sample_B'
        assert (adata.obs.loc[sample_b_mask, 'tissue'] == 'unknown').all()

    def test_assign_no_tissue_default_unknown(self, toy_adata):
        """Test that no tissue assignment defaults to 'unknown'."""
        adata = assign_tissue(toy_adata)
        assert (adata.obs['tissue'] == 'unknown').all()


class TestLoadMarkers:
    """Test marker loading."""

    def test_load_default_markers(self):
        """Test loading default markers file."""
        markers_data = load_markers()

        assert 'markers' in markers_data
        assert 'settings' in markers_data
        assert 'default' in markers_data['markers']

    def test_load_custom_markers(self, tmp_path):
        """Test loading custom markers file."""
        yaml_path = tmp_path / "custom_markers.yaml"
        yaml_content = """
markers:
  default:
    MyCell: [GENE1, GENE2, GENE3]
settings:
  min_markers_expressed: 1
"""
        yaml_path.write_text(yaml_content)

        markers_data = load_markers(str(yaml_path))
        assert 'MyCell' in markers_data['markers']['default']
        assert markers_data['settings']['min_markers_expressed'] == 1


class TestGetTissueMarkers:
    """Test tissue-specific marker retrieval."""

    def test_get_default_markers(self):
        """Test getting default markers for unknown tissue."""
        markers_data = load_markers()
        tissue_markers = get_tissue_markers(markers_data, 'unknown_tissue')

        # Should return default markers
        assert 'T_cell' in tissue_markers
        assert 'B_cell' in tissue_markers

    def test_get_lung_markers(self):
        """Test getting lung-specific markers."""
        markers_data = load_markers()
        tissue_markers = get_tissue_markers(markers_data, 'lung')

        # Should have lung-specific markers
        assert 'AT1' in tissue_markers
        assert 'AT2' in tissue_markers
        # Should also have default markers
        assert 'T_cell' in tissue_markers

    def test_get_liver_markers(self):
        """Test getting liver-specific markers."""
        markers_data = load_markers()
        tissue_markers = get_tissue_markers(markers_data, 'liver')

        assert 'Hepatocyte' in tissue_markers
        assert 'Kupffer_cell' in tissue_markers


class TestAnnotateCells:
    """Test cell type annotation."""

    def test_annotate_cells_adds_columns(self, toy_adata):
        """Test that annotation adds required columns."""
        # First normalize the data (required for annotation)
        sc.pp.normalize_total(toy_adata, target_sum=1e4)
        sc.pp.log1p(toy_adata)
        toy_adata.obs['tissue'] = 'unknown'

        annotate_cells(toy_adata)

        assert 'celltype' in toy_adata.obs.columns
        assert 'celltype_score' in toy_adata.obs.columns

    def test_annotate_cells_no_tissue_warning(self, toy_adata):
        """Test that annotation works without tissue column."""
        sc.pp.normalize_total(toy_adata, target_sum=1e4)
        sc.pp.log1p(toy_adata)

        # Should not raise error, just use 'unknown'
        annotate_cells(toy_adata)
        assert 'tissue' in toy_adata.obs.columns
        assert (toy_adata.obs['tissue'] == 'unknown').all()


class TestGenerateSummary:
    """Test summary generation functions."""

    def test_generate_tissue_summary(self, toy_adata):
        """Test tissue summary generation."""
        toy_adata.obs['tissue'] = np.random.choice(['lung', 'liver'], size=toy_adata.n_obs)
        toy_adata.obs['celltype'] = np.random.choice(['T_cell', 'B_cell'], size=toy_adata.n_obs)

        summary = generate_tissue_summary(toy_adata)

        assert isinstance(summary, pd.DataFrame)
        assert 'Total' in summary.index
        assert 'Total' in summary.columns

    def test_generate_cell_metadata(self, toy_adata):
        """Test cell metadata generation."""
        toy_adata.obs['tissue'] = 'lung'
        toy_adata.obs['celltype'] = 'T_cell'
        toy_adata.obs['celltype_score'] = 1.5
        toy_adata.obs['leiden'] = '0'

        metadata = generate_cell_metadata(toy_adata)

        assert isinstance(metadata, pd.DataFrame)
        assert 'tissue' in metadata.columns
        assert 'celltype' in metadata.columns
        assert 'celltype_score' in metadata.columns
        assert 'leiden' in metadata.columns


class TestEndToEndWithTissue:
    """End-to-end tests with tissue features."""

    def test_pipeline_with_single_tissue(self):
        """Test full pipeline with single tissue annotation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            create_toy_dataset(
                n_cells=200,
                n_genes=500,
                n_clusters=3,
                output_path=input_path
            )

            config = PipelineConfig(
                min_genes=50,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=200,
                n_pcs=20,
                resolution=0.3,
                seed=42,
                tissue='lung',
                run_annotation=True
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check tissue was assigned
            assert 'tissue' in results.adata.obs.columns
            assert (results.adata.obs['tissue'] == 'lung').all()

            # Check annotation was performed
            assert 'celltype' in results.adata.obs.columns
            assert 'celltype_score' in results.adata.obs.columns

            # Check outputs exist
            assert os.path.exists(os.path.join(output_dir, "cell_metadata.csv"))
            assert os.path.exists(os.path.join(output_dir, "tissue_summary.csv"))

    def test_pipeline_without_annotation(self):
        """Test pipeline with annotation disabled."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")

            create_toy_dataset(
                n_cells=100,
                n_genes=200,
                output_path=input_path
            )

            config = PipelineConfig(
                min_genes=20,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=100,
                n_pcs=10,
                run_annotation=False
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Annotation columns should not exist (or be default)
            # Pipeline should complete without errors
            assert results.adata is not None

    def test_pipeline_with_tissue_map(self):
        """Test pipeline with tissue map file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.h5ad")
            output_dir = os.path.join(tmpdir, "output")
            tissue_map_path = os.path.join(tmpdir, "tissue_map.csv")

            # Create toy data with sample column
            adata = create_toy_dataset(
                n_cells=200,
                n_genes=500,
                n_clusters=3,
                output_path=input_path
            )
            # Add sample column
            adata.obs['sample'] = np.random.choice(['S1', 'S2'], size=adata.n_obs)
            adata.write(input_path)

            # Create tissue map
            pd.DataFrame({
                'sample': ['S1', 'S2'],
                'tissue': ['lung', 'liver']
            }).to_csv(tissue_map_path, index=False)

            config = PipelineConfig(
                min_genes=50,
                min_cells=1,
                max_mt_pct=50.0,
                n_top_genes=200,
                n_pcs=20,
                resolution=0.3,
                seed=42,
                tissue_map_path=tissue_map_path,
                sample_col='sample',
                run_annotation=True
            )

            results = run_pipeline(
                input_path=input_path,
                output_dir=output_dir,
                config=config
            )

            # Check tissue was assigned correctly
            assert 'tissue' in results.adata.obs.columns
            tissues = results.adata.obs['tissue'].unique()
            assert set(tissues) == {'lung', 'liver'}


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
