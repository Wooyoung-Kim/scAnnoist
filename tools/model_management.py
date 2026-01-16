"""
CellTypist Model Management and Metadata Handling
Handles tissue/sample-aware model selection and metadata validation
"""

import os
from pathlib import Path
from typing import Dict, Optional, Any, Tuple
import yaml
import pandas as pd
import json


# Cache directory for CellTypist models
def get_cache_dir() -> Path:
    """Get cache directory for CellTypist models."""
    cache_dir = os.getenv("SCRNA_AGENT_CACHE_DIR")
    if cache_dir:
        cache_path = Path(cache_dir) / "celltypist"
    else:
        cache_path = Path.home() / ".cache" / "scrna-agent" / "celltypist"

    cache_path.mkdir(parents=True, exist_ok=True)
    return cache_path


def load_model_registry() -> Dict[str, Any]:
    """Load the CellTypist model registry from YAML file."""
    # Find registry file
    package_dir = Path(__file__).parent.parent
    registry_path = package_dir / "resources" / "celltypist_models.yaml"

    if not registry_path.exists():
        raise FileNotFoundError(
            f"Model registry not found at {registry_path}. "
            "Please ensure celltypist_models.yaml exists in resources/"
        )

    with open(registry_path, 'r') as f:
        registry = yaml.safe_load(f)

    return registry


def select_model_for_tissue(
    tissue: str,
    sample_type: Optional[str] = None,
    registry: Optional[Dict[str, Any]] = None
) -> Tuple[str, str]:
    """
    Select appropriate CellTypist model based on tissue and sample_type.

    Selection priority:
        1. overrides(tissue|sample_type) - exact match
        2. tissue - tissue-specific model
        3. sample_type - sample type-specific model
        4. default - fallback model

    Args:
        tissue: Tissue type
        sample_type: Optional sample type
        registry: Model registry (loads if not provided)

    Returns:
        Tuple of (model_name, selection_rule)
    """
    if registry is None:
        registry = load_model_registry()

    # Priority 1: Check overrides for exact tissue|sample_type match
    if sample_type and "overrides" in registry:
        override_key = f"{tissue}|{sample_type}"
        if override_key in registry["overrides"]:
            model = registry["overrides"][override_key]["model"]
            return model, f"override({override_key})"

    # Priority 2: Check tissue-specific model
    if "tissue" in registry and tissue in registry["tissue"]:
        model = registry["tissue"][tissue]["model"]
        return model, f"tissue({tissue})"

    # Priority 3: Check sample_type-specific model
    if sample_type and "sample_type" in registry:
        if sample_type in registry["sample_type"]:
            model = registry["sample_type"][sample_type]["model"]
            return model, f"sample_type({sample_type})"

    # Priority 4: Default model
    if "default" in registry:
        model = registry["default"]["model"]
        return model, "default"

    raise ValueError(
        f"No model found for tissue='{tissue}', sample_type='{sample_type}' "
        "and no default model configured in registry."
    )


def validate_and_add_tissue_metadata(
    adata,
    tissue: Optional[str] = None,
    tissue_map: Optional[str] = None,
    sample_type: Optional[str] = None,
    sample_type_map: Optional[str] = None,
    sample_col: str = "sample",
    allow_missing_tissue: bool = False
) -> None:
    """
    Validate and add tissue/sample_type metadata to AnnData object.

    Args:
        adata: AnnData object
        tissue: Tissue type to apply to all cells
        tissue_map: Path to CSV file with sample->tissue mapping
        sample_type: Sample type to apply to all cells
        sample_type_map: Path to CSV file with sample->sample_type mapping
        sample_col: Column name containing sample IDs
        allow_missing_tissue: If True, fill missing tissues with "unknown"

    Raises:
        ValueError: If validation fails
    """
    # Validate sample_col exists
    if sample_col not in adata.obs.columns:
        raise ValueError(
            f"Sample column '{sample_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Get unique samples
    samples = adata.obs[sample_col].unique()

    # Handle tissue metadata
    if tissue is not None:
        # Apply tissue to all cells
        adata.obs["tissue"] = tissue
    elif tissue_map is not None:
        # Load tissue mapping
        if not os.path.exists(tissue_map):
            raise FileNotFoundError(f"Tissue map file not found: {tissue_map}")

        tissue_df = pd.read_csv(tissue_map)
        if len(tissue_df.columns) < 2:
            raise ValueError(
                f"Tissue map must have at least 2 columns (sample, tissue). "
                f"Found: {tissue_df.columns.tolist()}"
            )

        # Use first two columns as sample, tissue
        tissue_df.columns = ["sample", "tissue"] + list(tissue_df.columns[2:])
        tissue_mapping = dict(zip(tissue_df["sample"], tissue_df["tissue"]))

        # Validate all samples are in mapping
        unknown_samples = set(samples) - set(tissue_mapping.keys())
        if unknown_samples:
            if not allow_missing_tissue:
                raise ValueError(
                    f"Unknown sample IDs in tissue mapping: {unknown_samples}. "
                    f"Expected samples: {set(tissue_mapping.keys())}. "
                    f"Use --allow-missing-tissue to fill with 'unknown'."
                )
            else:
                # Fill missing with "unknown"
                for sample in unknown_samples:
                    tissue_mapping[sample] = "unknown"

        # Map tissue to cells
        adata.obs["tissue"] = adata.obs[sample_col].map(tissue_mapping)
    else:
        if not allow_missing_tissue:
            raise ValueError(
                "Must provide either --tissue or --tissue-map. "
                "Use --allow-missing-tissue to skip tissue annotation."
            )
        else:
            adata.obs["tissue"] = "unknown"

    # Handle sample_type metadata
    if sample_type is not None:
        # Apply sample_type to all cells
        adata.obs["sample_type"] = sample_type
    elif sample_type_map is not None:
        # Load sample_type mapping
        if not os.path.exists(sample_type_map):
            raise FileNotFoundError(f"Sample type map file not found: {sample_type_map}")

        sample_type_df = pd.read_csv(sample_type_map)
        if len(sample_type_df.columns) < 2:
            raise ValueError(
                f"Sample type map must have at least 2 columns (sample, sample_type). "
                f"Found: {sample_type_df.columns.tolist()}"
            )

        # Use first two columns
        sample_type_df.columns = ["sample", "sample_type"] + list(sample_type_df.columns[2:])
        sample_type_mapping = dict(zip(sample_type_df["sample"], sample_type_df["sample_type"]))

        # Validate all samples are in mapping
        unknown_samples = set(samples) - set(sample_type_mapping.keys())
        if unknown_samples:
            # For sample_type, we allow missing values (optional metadata)
            for sample in unknown_samples:
                sample_type_mapping[sample] = None

        # Map sample_type to cells
        adata.obs["sample_type"] = adata.obs[sample_col].map(sample_type_mapping)
    else:
        # sample_type is optional
        adata.obs["sample_type"] = None


def get_models_for_adata(
    adata,
    sample_col: str = "sample"
) -> Dict[str, Dict[str, Any]]:
    """
    Determine which models to use for each tissue in the dataset.

    Args:
        adata: AnnData object with tissue and sample_type columns
        sample_col: Column name for samples

    Returns:
        Dictionary mapping tissue -> {
            "model": model_name,
            "rule": selection_rule,
            "n_cells": count,
            "samples": list of samples
        }
    """
    if "tissue" not in adata.obs.columns:
        raise ValueError("Tissue metadata not found. Run validate_and_add_tissue_metadata first.")

    registry = load_model_registry()
    tissue_info = {}

    # Group by tissue
    for tissue in adata.obs["tissue"].unique():
        tissue_mask = adata.obs["tissue"] == tissue
        tissue_cells = adata.obs[tissue_mask]

        # Get representative sample_type (most common in this tissue)
        sample_types = tissue_cells["sample_type"].dropna()
        if len(sample_types) > 0:
            sample_type = sample_types.mode()[0] if len(sample_types.mode()) > 0 else None
        else:
            sample_type = None

        # Select model
        model, rule = select_model_for_tissue(tissue, sample_type, registry)

        tissue_info[tissue] = {
            "model": model,
            "rule": rule,
            "n_cells": tissue_mask.sum(),
            "samples": tissue_cells[sample_col].unique().tolist(),
            "sample_type": sample_type
        }

    return tissue_info


def save_run_metadata(
    output_dir: str,
    tissue_info: Dict[str, Dict[str, Any]],
    additional_info: Optional[Dict[str, Any]] = None
) -> None:
    """
    Save metadata about model selection and run info.

    Args:
        output_dir: Output directory
        tissue_info: Tissue information from get_models_for_adata
        additional_info: Additional metadata to save
    """
    os.makedirs(output_dir, exist_ok=True)

    metadata = {
        "tissue_model_mapping": tissue_info,
        "total_tissues": len(tissue_info),
        "total_cells": sum(info["n_cells"] for info in tissue_info.values()),
    }

    if additional_info:
        metadata.update(additional_info)

    metadata_path = os.path.join(output_dir, "run_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"Saved run metadata to {metadata_path}")


def pull_celltypist_models(
    tissue: Optional[str] = None,
    pull_all: bool = False
) -> None:
    """
    Pull CellTypist models from repository.

    Args:
        tissue: Specific tissue to pull model for
        pull_all: Pull all models in registry
    """
    try:
        from celltypist import models
    except ImportError:
        raise ImportError(
            "CellTypist not installed. "
            "Install with: pip install celltypist"
        )

    registry = load_model_registry()
    cache_dir = get_cache_dir()

    models_to_pull = set()

    if pull_all:
        # Collect all unique models
        if "default" in registry:
            models_to_pull.add(registry["default"]["model"])
        if "tissue" in registry:
            for tissue_data in registry["tissue"].values():
                models_to_pull.add(tissue_data["model"])
        if "sample_type" in registry:
            for st_data in registry["sample_type"].values():
                models_to_pull.add(st_data["model"])
        if "overrides" in registry:
            for override_data in registry["overrides"].values():
                models_to_pull.add(override_data["model"])
    elif tissue:
        # Pull model for specific tissue
        model, rule = select_model_for_tissue(tissue, registry=registry)
        models_to_pull.add(model)
    else:
        raise ValueError("Must specify --tissue or --all")

    # Download models
    print(f"Pulling {len(models_to_pull)} model(s)...")
    for model_name in models_to_pull:
        # Skip if it's a path (custom model)
        if "/" in model_name or model_name.startswith("."):
            print(f"  Skipping {model_name} (custom path)")
            continue

        try:
            print(f"  Downloading {model_name}...")
            models.download_models(model=model_name, force_update=False)
            print(f"    ✓ {model_name}")
        except Exception as e:
            print(f"    ✗ Failed to download {model_name}: {e}")

    print("Done!")


def list_available_models() -> None:
    """List all models in the registry."""
    registry = load_model_registry()

    print("\nDefault Model:")
    print(f"  {registry['default']['model']} - {registry['default']['description']}")

    if "tissue" in registry:
        print("\nTissue-Specific Models:")
        for tissue, data in registry["tissue"].items():
            print(f"  {tissue}: {data['model']}")

    if "sample_type" in registry:
        print("\nSample Type-Specific Models:")
        for st, data in registry["sample_type"].items():
            print(f"  {st}: {data['model']}")

    if "overrides" in registry:
        print("\nOverrides:")
        for override, data in registry["overrides"].items():
            print(f"  {override}: {data['model']}")
