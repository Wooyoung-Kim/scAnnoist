"""
Tests for tissue-aware CellTypist annotation functionality
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import os
import sys
import importlib.util

# Helper to import module without going through __init__.py
def import_module_from_path(module_name, file_path):
    """Import a module directly from a file path."""
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

# Import model_management directly
package_dir = Path(__file__).parent.parent
model_mgmt = import_module_from_path(
    "model_management",
    package_dir / "tools" / "model_management.py"
)
annotation_tools_path = package_dir / "tools" / "annotation_tools.py"


# ============================================================
# Unit Tests: Model Selection Priority
# ============================================================

def test_model_selection_priority():
    """Test deterministic model selection priority."""
    registry = model_mgmt.load_model_registry()

    # Test Priority 1: Override (tissue|sample_type)
    model, rule = model_mgmt.select_model_for_tissue("Blood", "PBMC", registry)
    assert rule.startswith("override"), f"Expected override, got {rule}"
    assert model == "Healthy_COVID19_PBMC.pkl"

    # Test Priority 2: Tissue-specific
    model, rule = select_model_for_tissue("Spleen", None, registry)
    assert rule.startswith("tissue"), f"Expected tissue, got {rule}"
    assert model == "Immune_All_High.pkl"

    # Test Priority 3: Sample type-specific
    model, rule = select_model_for_tissue("Unknown_Tissue", "PBMC", registry)
    assert rule.startswith("sample_type"), f"Expected sample_type, got {rule}"
    assert model == "Healthy_COVID19_PBMC.pkl"

    # Test Priority 4: Default
    model, rule = select_model_for_tissue("Unknown_Tissue", None, registry)
    assert rule == "default", f"Expected default, got {rule}"
    assert model == "Immune_All_High.pkl"


def test_override_beats_tissue():
    """Test that override priority beats tissue priority."""
    from tools.model_management import select_model_for_tissue

    # Blood|PBMC should use override, not just Blood tissue
    model, rule = select_model_for_tissue("Blood", "PBMC")
    assert rule.startswith("override"), "Override should beat tissue priority"
    assert "Blood|PBMC" in rule


def test_tissue_beats_sample_type():
    """Test that tissue priority beats sample_type priority."""
    from tools.model_management import select_model_for_tissue

    # Spleen tissue should use tissue-specific model, even if sample_type matches
    model, rule = select_model_for_tissue("Spleen", "tumor")
    assert rule.startswith("tissue"), "Tissue should beat sample_type priority"


def test_sample_type_beats_default():
    """Test that sample_type priority beats default."""
    from tools.model_management import select_model_for_tissue

    # Unknown tissue with PBMC sample_type should use PBMC model, not default
    model, rule = select_model_for_tissue("Unknown", "PBMC")
    assert rule.startswith("sample_type"), "Sample type should beat default"
    assert model == "Healthy_COVID19_PBMC.pkl"


# ============================================================
# Unit Tests: Tissue Metadata Validation
# ============================================================

def create_mock_adata(n_cells=1000, n_genes=500, samples=None):
    """Create mock AnnData for testing."""
    import scanpy as sc

    if samples is None:
        samples = ["sample1", "sample2", "sample3"]

    # Create random data
    X = np.random.randn(n_cells, n_genes)
    obs = pd.DataFrame({
        "sample": np.random.choice(samples, n_cells),
    }, index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    adata = sc.AnnData(X=X, obs=obs, var=var)
    return adata


def test_tissue_map_parsing():
    """Test parsing of tissue mapping CSV."""
    from tools.model_management import validate_and_add_tissue_metadata

    adata = create_mock_adata(samples=["S1", "S2", "S3"])

    # Create tissue map
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("sample,tissue\n")
        f.write("S1,Spleen\n")
        f.write("S2,Blood\n")
        f.write("S3,Liver\n")
        tissue_map_path = f.name

    try:
        validate_and_add_tissue_metadata(
            adata,
            tissue_map=tissue_map_path,
            sample_col="sample"
        )

        assert "tissue" in adata.obs.columns
        assert set(adata.obs["tissue"].unique()) == {"Spleen", "Blood", "Liver"}
    finally:
        os.unlink(tissue_map_path)


def test_unknown_sample_error():
    """Test that unknown sample IDs raise error without allow_missing_tissue."""
    from tools.model_management import validate_and_add_tissue_metadata

    adata = create_mock_adata(samples=["S1", "S2", "S3", "S4"])

    # Create tissue map missing S4
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("sample,tissue\n")
        f.write("S1,Spleen\n")
        f.write("S2,Blood\n")
        f.write("S3,Liver\n")
        tissue_map_path = f.name

    try:
        with pytest.raises(ValueError, match="Unknown sample IDs"):
            validate_and_add_tissue_metadata(
                adata,
                tissue_map=tissue_map_path,
                sample_col="sample",
                allow_missing_tissue=False
            )
    finally:
        os.unlink(tissue_map_path)


def test_allow_missing_tissue():
    """Test that allow_missing_tissue fills unknowns with 'unknown'."""
    from tools.model_management import validate_and_add_tissue_metadata

    adata = create_mock_adata(samples=["S1", "S2", "S3", "S4"])

    # Create tissue map missing S4
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("sample,tissue\n")
        f.write("S1,Spleen\n")
        f.write("S2,Blood\n")
        f.write("S3,Liver\n")
        tissue_map_path = f.name

    try:
        validate_and_add_tissue_metadata(
            adata,
            tissue_map=tissue_map_path,
            sample_col="sample",
            allow_missing_tissue=True
        )

        assert "tissue" in adata.obs.columns
        # S4 cells should be marked as "unknown"
        s4_cells = adata.obs[adata.obs["sample"] == "S4"]
        assert (s4_cells["tissue"] == "unknown").all()
    finally:
        os.unlink(tissue_map_path)


def test_apply_tissue_to_all():
    """Test applying single tissue to all cells."""
    from tools.model_management import validate_and_add_tissue_metadata

    adata = create_mock_adata()

    validate_and_add_tissue_metadata(
        adata,
        tissue="Spleen",
        sample_col="sample"
    )

    assert "tissue" in adata.obs.columns
    assert (adata.obs["tissue"] == "Spleen").all()


def test_sample_type_map():
    """Test sample_type mapping."""
    from tools.model_management import validate_and_add_tissue_metadata

    adata = create_mock_adata(samples=["S1", "S2", "S3"])

    # Create sample_type map
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("sample,sample_type\n")
        f.write("S1,PBMC\n")
        f.write("S2,tumor\n")
        f.write("S3,organoid\n")
        sample_type_map_path = f.name

    try:
        validate_and_add_tissue_metadata(
            adata,
            tissue="Blood",
            sample_type_map=sample_type_map_path,
            sample_col="sample"
        )

        assert "sample_type" in adata.obs.columns
        assert set(adata.obs["sample_type"].dropna().unique()) == {"PBMC", "tumor", "organoid"}
    finally:
        os.unlink(sample_type_map_path)


# ============================================================
# Integration Tests: Per-Tissue Annotation Merge
# ============================================================

def test_per_tissue_subset_annotation_mock():
    """Test per-tissue subset annotation with mock predictions."""
    from tools.model_management import get_models_for_adata

    adata = create_mock_adata(n_cells=1000, samples=["S1", "S2", "S3"])

    # Add tissue metadata
    tissue_mapping = {"S1": "Spleen", "S2": "Blood", "S3": "Liver"}
    adata.obs["tissue"] = adata.obs["sample"].map(tissue_mapping)
    adata.obs["sample_type"] = None

    # Get models for each tissue
    tissue_info = get_models_for_adata(adata, sample_col="sample")

    assert len(tissue_info) == 3
    assert "Spleen" in tissue_info
    assert "Blood" in tissue_info
    assert "Liver" in tissue_info

    # Check that each tissue has correct info
    for tissue, info in tissue_info.items():
        assert "model" in info
        assert "rule" in info
        assert "n_cells" in info
        assert info["n_cells"] > 0


def test_consistent_output_columns():
    """Test that output columns are consistent across tissues."""
    adata = create_mock_adata(n_cells=500, samples=["S1", "S2"])

    # Add tissue metadata
    adata.obs["tissue"] = adata.obs["sample"].map({"S1": "Spleen", "S2": "Blood"})
    adata.obs["sample_type"] = None

    # Mock annotation outputs
    adata.obs["celltypist_label"] = np.random.choice(["B cells", "T cells"], len(adata))
    adata.obs["celltypist_confidence"] = np.random.rand(len(adata))
    adata.obs["celltypist_model_used"] = "Immune_All_High.pkl"

    # Check columns exist
    assert "celltypist_label" in adata.obs.columns
    assert "celltypist_confidence" in adata.obs.columns
    assert "celltypist_model_used" in adata.obs.columns


def test_merge_preserves_cell_order():
    """Test that merging per-tissue results preserves cell order."""
    adata = create_mock_adata(n_cells=500, samples=["S1", "S2", "S3"])
    original_index = adata.obs.index.copy()

    # Add tissue metadata
    tissue_mapping = {"S1": "Spleen", "S2": "Blood", "S3": "Liver"}
    adata.obs["tissue"] = adata.obs["sample"].map(tissue_mapping)

    # Mock per-tissue annotation
    for tissue in ["Spleen", "Blood", "Liver"]:
        tissue_mask = adata.obs["tissue"] == tissue
        adata.obs.loc[tissue_mask, "celltypist_label"] = f"{tissue}_celltype"

    # Check that index is preserved
    assert (adata.obs.index == original_index).all()


# ============================================================
# End-to-End Tests
# ============================================================

def test_end_to_end_tissue_metadata_workflow():
    """Test complete workflow: metadata -> model selection -> outputs."""
    from tools.model_management import (
        validate_and_add_tissue_metadata,
        get_models_for_adata,
        save_run_metadata
    )

    adata = create_mock_adata(n_cells=1000, samples=["S1", "S2", "S3"])

    # Step 1: Add tissue metadata
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
        f.write("sample,tissue\n")
        f.write("S1,Spleen\n")
        f.write("S2,Blood\n")
        f.write("S3,Liver\n")
        tissue_map_path = f.name

    try:
        validate_and_add_tissue_metadata(
            adata,
            tissue_map=tissue_map_path,
            sample_col="sample"
        )

        # Step 2: Get model assignments
        tissue_info = get_models_for_adata(adata, sample_col="sample")

        assert len(tissue_info) == 3
        assert all("model" in info for info in tissue_info.values())

        # Step 3: Save metadata
        with tempfile.TemporaryDirectory() as tmpdir:
            save_run_metadata(tmpdir, tissue_info)

            metadata_file = Path(tmpdir) / "run_metadata.json"
            assert metadata_file.exists()

            import json
            with open(metadata_file) as f:
                metadata = json.load(f)

            assert "tissue_model_mapping" in metadata
            assert len(metadata["tissue_model_mapping"]) == 3

    finally:
        os.unlink(tissue_map_path)


def test_output_csv_generation():
    """Test generation of prediction and tissue usage CSVs."""
    # Import directly from annotation_tools module, not through __init__
    from tools.annotation_tools import save_celltypist_outputs
    from tools.model_management import get_models_for_adata

    adata = create_mock_adata(n_cells=500, samples=["S1", "S2"])

    # Add metadata
    adata.obs["tissue"] = adata.obs["sample"].map({"S1": "Spleen", "S2": "Blood"})
    adata.obs["sample_type"] = None

    # Mock annotation results
    adata.obs["celltypist_label"] = np.random.choice(["B cells", "T cells"], len(adata))
    adata.obs["celltypist_confidence"] = np.random.rand(len(adata))
    adata.obs["celltypist_model_used"] = "Immune_All_High.pkl"

    # Get tissue info
    tissue_info = get_models_for_adata(adata, sample_col="sample")

    # Save outputs
    with tempfile.TemporaryDirectory() as tmpdir:
        save_celltypist_outputs(adata, tissue_info, tmpdir, sample_col="sample")

        # Check files exist
        predictions_file = Path(tmpdir) / "celltypist_predictions.csv"
        tissue_usage_file = Path(tmpdir) / "tissue_model_usage.csv"

        assert predictions_file.exists()
        assert tissue_usage_file.exists()

        # Check predictions CSV structure
        pred_df = pd.read_csv(predictions_file)
        assert "cell_id" in pred_df.columns
        assert "sample" in pred_df.columns
        assert "tissue" in pred_df.columns
        assert "label" in pred_df.columns
        assert "confidence" in pred_df.columns
        assert "model_used" in pred_df.columns

        # Check tissue usage CSV structure
        tissue_df = pd.read_csv(tissue_usage_file)
        assert "tissue" in tissue_df.columns
        assert "model_used" in tissue_df.columns
        assert "selection_rule" in tissue_df.columns
        assert "n_cells" in tissue_df.columns
        assert len(tissue_df) == 2  # S1 and S2


def test_model_registry_loading():
    """Test that model registry loads correctly."""
    from tools.model_management import load_model_registry

    registry = load_model_registry()

    assert "default" in registry
    assert "tissue" in registry
    assert "sample_type" in registry
    assert "overrides" in registry

    # Check default model exists
    assert "model" in registry["default"]

    # Check some expected tissues
    assert "Spleen" in registry["tissue"]
    assert "Blood" in registry["tissue"]


# ============================================================
# Run Tests
# ============================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
