"""Pipeline configuration and path resolution."""

from __future__ import annotations

import importlib.util
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Mapping, Optional


def _default_project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def _load_preprocess_paths(project_root: Path, test_mode: bool):
    """Load Data_preprocess.paths with PIPELINE_TEST set as needed."""
    prev = os.environ.get("PIPELINE_TEST")
    if test_mode:
        os.environ["PIPELINE_TEST"] = "1"
    else:
        os.environ.pop("PIPELINE_TEST", None)
    try:
        paths_file = project_root / "Data_preprocess" / "paths.py"
        spec = importlib.util.spec_from_file_location("preprocess_paths", paths_file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        if prev is None:
            os.environ.pop("PIPELINE_TEST", None)
        else:
            os.environ["PIPELINE_TEST"] = prev


@dataclass
class PipelineConfig:
    """Root configuration for pipeline runs and summaries."""

    project_root: Path = field(default_factory=_default_project_root)
    test_mode: bool = False
    mmcif_source: Optional[str] = None
    python_cmd: Optional[str] = None
    schrodinger: Optional[str] = None

    def __post_init__(self) -> None:
        self.project_root = Path(self.project_root).resolve()

    @property
    def data_dir(self) -> Path:
        return self.test_dir if self.test_mode else self.project_root / "Data"

    @property
    def test_dir(self) -> Path:
        return self.project_root / "test"

    @property
    def paths(self):
        return _load_preprocess_paths(self.project_root, self.test_mode)

    def run_script(self, stage: str) -> Path:
        """Return the bash driver script for a pipeline stage."""
        scripts = {
            "preprocess": (
                "run_data_preprocess_test.sh" if self.test_mode else "run_data_preprocess.sh"
            ),
            "alignment": (
                "run_alignment_test.sh" if self.test_mode else "run_alignment.sh"
            ),
            "pharmacophore": (
                "run_pharmacophore_test.sh" if self.test_mode else "run_pharmacophore.sh"
            ),
            "categorization": (
                "run_categorization_test.sh" if self.test_mode else "run_categorization.sh"
            ),
        }
        name = scripts[stage]
        base = self.test_dir if self.test_mode else self.project_root
        return base / name

    def subprocess_env(self) -> dict[str, str]:
        env = os.environ.copy()
        if self.test_mode:
            env["PIPELINE_TEST"] = "1"
        else:
            env.pop("PIPELINE_TEST", None)
        if self.mmcif_source:
            env["MMCIF_SOURCE"] = self.mmcif_source
        if self.python_cmd:
            env["PYTHON_CMD"] = self.python_cmd
        if self.schrodinger:
            env["SCHRODINGER"] = self.schrodinger
        return env

    def stage_dir(self, relative: str) -> Path:
        return self.data_dir / relative
