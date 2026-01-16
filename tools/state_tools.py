# State Management and Tracking
# Persistent state storage for tracking pipeline progress

import os
import json
from datetime import datetime
from typing import Dict, Any, Optional, List
from pathlib import Path

# Default state file location
DEFAULT_STATE_DIR = Path.home() / ".scrna_agent"
DEFAULT_STATE_FILE = DEFAULT_STATE_DIR / "pipeline_state.json"


class PipelineState:
    """
    Persistent state manager for the scRNA-seq analysis pipeline.
    Tracks completed steps, saves results, and enables resumption.
    """
    
    def __init__(self, state_file: Optional[Path] = None):
        self.state_file = state_file or DEFAULT_STATE_FILE
        self.state_file.parent.mkdir(parents=True, exist_ok=True)
        self._state = self._load_state()
    
    def _load_state(self) -> Dict[str, Any]:
        """Load state from file or create new."""
        if self.state_file.exists():
            try:
                with open(self.state_file, "r") as f:
                    return json.load(f)
            except json.JSONDecodeError:
                return self._default_state()
        return self._default_state()
    
    def _default_state(self) -> Dict[str, Any]:
        """Create default state structure."""
        return {
            "version": "0.1.0",
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat(),
            "datasets": {},  # dataset_id -> dataset state
            "current_dataset": None,
        }
    
    def save(self):
        """Save state to file."""
        self._state["updated_at"] = datetime.now().isoformat()
        with open(self.state_file, "w") as f:
            json.dump(self._state, f, indent=2, default=str)
    
    def get_dataset(self, dataset_id: str) -> Dict[str, Any]:
        """Get state for a specific dataset."""
        if dataset_id not in self._state["datasets"]:
            self._state["datasets"][dataset_id] = {
                "id": dataset_id,
                "created_at": datetime.now().isoformat(),
                "status": "initialized",
                "steps_completed": [],
                "current_step": None,
                "file_path": None,
                "output_dir": None,
                "results": {},
                "errors": [],
            }
        return self._state["datasets"][dataset_id]
    
    def update_dataset(self, dataset_id: str, **kwargs):
        """Update dataset state."""
        dataset = self.get_dataset(dataset_id)
        dataset.update(kwargs)
        dataset["updated_at"] = datetime.now().isoformat()
        self.save()
    
    def mark_step_complete(self, dataset_id: str, step_name: str, result: Any = None):
        """Mark a step as complete for a dataset."""
        dataset = self.get_dataset(dataset_id)
        
        if step_name not in dataset["steps_completed"]:
            dataset["steps_completed"].append(step_name)
        
        if result:
            dataset["results"][step_name] = {
                "completed_at": datetime.now().isoformat(),
                "summary": str(result)[:500],  # Truncate for storage
            }
        
        dataset["current_step"] = None
        self.save()
    
    def set_current_step(self, dataset_id: str, step_name: str):
        """Set the current step being executed."""
        dataset = self.get_dataset(dataset_id)
        dataset["current_step"] = step_name
        dataset["status"] = "running"
        self.save()
    
    def log_error(self, dataset_id: str, error_message: str, step_name: Optional[str] = None):
        """Log an error for a dataset."""
        dataset = self.get_dataset(dataset_id)
        dataset["errors"].append({
            "timestamp": datetime.now().isoformat(),
            "step": step_name or dataset.get("current_step"),
            "message": error_message,
        })
        dataset["status"] = "error"
        self.save()
    
    def get_progress(self, dataset_id: str) -> Dict[str, Any]:
        """Get progress summary for a dataset."""
        dataset = self.get_dataset(dataset_id)
        
        all_steps = [
            "load_data",
            "annotation_celltypist",
            "annotation_sctype",
            "annotation_literature",
            "consensus_annotation",
            "prepare_integration",
            "run_scvi",
            "run_harmony",
            "benchmark_integration",
            "select_integration",
            "setup_scanvi",
            "run_scanvi_finetuning",
            "compute_umap",
            "plot_results",
            "save_results",
        ]
        
        completed = set(dataset["steps_completed"])
        
        return {
            "dataset_id": dataset_id,
            "status": dataset["status"],
            "current_step": dataset.get("current_step"),
            "completed_steps": len(completed),
            "total_steps": len(all_steps),
            "progress_percent": len(completed) / len(all_steps) * 100,
            "remaining_steps": [s for s in all_steps if s not in completed],
            "errors": dataset.get("errors", []),
        }
    
    def list_datasets(self) -> List[Dict[str, Any]]:
        """List all tracked datasets."""
        return [
            {
                "id": did,
                "status": data.get("status"),
                "file_path": data.get("file_path"),
                "steps_completed": len(data.get("steps_completed", [])),
                "created_at": data.get("created_at"),
            }
            for did, data in self._state["datasets"].items()
        ]
    
    def clear_dataset(self, dataset_id: str):
        """Clear state for a specific dataset."""
        if dataset_id in self._state["datasets"]:
            del self._state["datasets"][dataset_id]
            self.save()


# Global state instance
_pipeline_state_manager: Optional[PipelineState] = None


def get_state_manager() -> PipelineState:
    """Get or create the global state manager."""
    global _pipeline_state_manager
    if _pipeline_state_manager is None:
        _pipeline_state_manager = PipelineState()
    return _pipeline_state_manager


# ============================================================
# Tools for state tracking
# ============================================================

from langchain_core.tools import tool
from typing import Annotated


@tool
def track_pipeline_progress(
    dataset_id: Annotated[str, "Unique identifier for the dataset"],
) -> Annotated[str, "Progress report"]:
    """
    Get the current progress of a pipeline run for a dataset.
    Shows completed steps, current step, and remaining steps.
    """
    state_mgr = get_state_manager()
    progress = state_mgr.get_progress(dataset_id)
    
    result = f"""
Pipeline Progress: {dataset_id}
================================
Status: {progress['status']}
Current Step: {progress['current_step'] or 'None'}
Progress: {progress['completed_steps']}/{progress['total_steps']} ({progress['progress_percent']:.1f}%)

Remaining Steps:
"""
    for step in progress['remaining_steps'][:5]:
        result += f"  - {step}\n"
    
    if progress['errors']:
        result += f"\nErrors ({len(progress['errors'])}):\n"
        for err in progress['errors'][-3:]:  # Last 3 errors
            result += f"  - [{err['step']}] {err['message'][:100]}\n"
    
    return result


@tool
def list_tracked_datasets() -> Annotated[str, "List of tracked datasets"]:
    """
    List all datasets that have been processed or are being processed.
    """
    state_mgr = get_state_manager()
    datasets = state_mgr.list_datasets()
    
    if not datasets:
        return "No datasets tracked yet."
    
    result = "Tracked Datasets\n================\n"
    for ds in datasets:
        result += f"""
- {ds['id']}
  Status: {ds['status']}
  Steps: {ds['steps_completed']}
  File: {ds['file_path'] or 'Not set'}
"""
    
    return result


@tool
def resume_pipeline(
    dataset_id: Annotated[str, "Dataset ID to resume"],
) -> Annotated[str, "Resume instructions"]:
    """
    Get instructions for resuming a paused or failed pipeline.
    Shows the last completed step and what to run next.
    """
    state_mgr = get_state_manager()
    progress = state_mgr.get_progress(dataset_id)
    
    if not progress['remaining_steps']:
        return f"Pipeline for {dataset_id} is already complete!"
    
    next_step = progress['remaining_steps'][0]
    
    result = f"""
Resume Pipeline: {dataset_id}
============================
Last Status: {progress['status']}
Completed: {progress['completed_steps']} steps

Next Step: {next_step}

To resume, run the pipeline with:
  - Skip completed steps
  - Start from: {next_step}

Recent Errors:
"""
    if progress['errors']:
        for err in progress['errors'][-2:]:
            result += f"  - {err['message'][:150]}\n"
    else:
        result += "  None\n"
    
    return result
