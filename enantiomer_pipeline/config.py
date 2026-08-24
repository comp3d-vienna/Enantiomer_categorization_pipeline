"""Pipeline configuration and path resolution."""

from __future__ import annotations

import importlib.util
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


def _default_project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _resolve_user_path(project_root: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = project_root / path
    return path.resolve()


def resolve_data_dir(project_root: Path, mmcif_source: Optional[str] = None) -> Path:
    """Return Data/ or Data/test/ (when mmCIF input is Data/test/mmCIF)."""
    project_root = Path(project_root).resolve()
    production = project_root / "Data"
    trial = production / "test"
    if mmcif_source:
        src = _resolve_user_path(project_root, mmcif_source)
        if src == (trial / "mmCIF").resolve():
            return trial
    override = os.environ.get("PIPELINE_DATA_DIR")
    if override:
        return _resolve_user_path(project_root, override)
    return production


def _load_preprocess_paths(project_root: Path, data_dir: Path, mmcif_source: Optional[str]):
    """Load Data_preprocess.paths with the selected data root."""
    prev_data = os.environ.get("PIPELINE_DATA_DIR")
    prev_mmcif = os.environ.get("MMCIF_SOURCE")
    os.environ["PIPELINE_DATA_DIR"] = str(data_dir)
    if mmcif_source:
        os.environ["MMCIF_SOURCE"] = str(_resolve_user_path(project_root, mmcif_source))
    try:
        paths_file = project_root / "Data_preprocess" / "paths.py"
        spec = importlib.util.spec_from_file_location("preprocess_paths", paths_file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if prev_data is None:
            os.environ.pop("PIPELINE_DATA_DIR", None)
        else:
            os.environ["PIPELINE_DATA_DIR"] = prev_data
        if prev_mmcif is None:
            os.environ.pop("MMCIF_SOURCE", None)
        else:
            os.environ["MMCIF_SOURCE"] = prev_mmcif


@dataclass
class PipelineConfig:
    """Root configuration for pipeline runs and summaries."""

    project_root: Path = field(default_factory=_default_project_root)
    mmcif_source: Optional[str] = None
    python_cmd: Optional[str] = None
    schrodinger: Optional[str] = None

    def __post_init__(self) -> None:
        self.project_root = Path(self.project_root).resolve()

    @property
    def data_dir(self) -> Path:
        return resolve_data_dir(self.project_root, self.mmcif_source)

    @property
    def paths(self):
        return _load_preprocess_paths(self.project_root, self.data_dir, self.mmcif_source)

    def run_script(self, stage: str) -> Path:
        """Return the bash driver script for a pipeline stage."""
        scripts = {
            "preprocess": "run_data_preprocess.sh",
            "alignment": "run_alignment.sh",
            "pharmacophore": "run_pharmacophore.sh",
            "categorization": "run_categorization.sh",
        }
        return self.project_root / scripts[stage]

    def subprocess_env(self) -> dict[str, str]:
        env = os.environ.copy()
        env.pop("PIPELINE_TEST", None)
        env["PIPELINE_DATA_DIR"] = str(self.data_dir)
        if self.mmcif_source:
            env["MMCIF_SOURCE"] = str(_resolve_user_path(self.project_root, self.mmcif_source))
        if self.python_cmd:
            env["PYTHON_CMD"] = self.python_cmd
        if self.schrodinger:
            env["SCHRODINGER"] = self.schrodinger
        return env

    def stage_dir(self, relative: str) -> Path:
        return self.data_dir / relative
