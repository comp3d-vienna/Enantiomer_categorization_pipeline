"""Run pipeline stages by invoking existing shell drivers."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from typing import Any, Optional, Sequence

from ._results import STAGE_NAMES, collect_stage_result
from .config import PipelineConfig

STAGES = ("preprocess", "alignment", "pharmacophore", "categorization")


@dataclass
class StageResult:
    stage: str
    returncode: int
    script: str
    result: dict[str, Any]

    @property
    def ok(self) -> bool:
        return self.returncode == 0


def run_stage(
    config: PipelineConfig,
    stage: str,
    extra_args: Optional[Sequence[str]] = None,
) -> StageResult:
    if stage not in STAGES:
        raise ValueError(f"Unknown stage: {stage}. Choose from: {', '.join(STAGES)}")

    script = config.run_script(stage)
    if not script.is_file():
        raise FileNotFoundError(f"Pipeline script not found: {script}")

    cmd = ["bash", str(script)]
    if extra_args:
        cmd.extend(extra_args)

    proc = subprocess.run(
        cmd,
        cwd=str(config.project_root),
        env=config.subprocess_env(),
        check=False,
    )

    return StageResult(
        stage=stage,
        returncode=proc.returncode,
        script=str(script),
        result=collect_stage_result(config, stage),
    )


def run_all(
    config: PipelineConfig,
    stages: Optional[Sequence[str]] = None,
    stage_args: Optional[dict[str, list[str]]] = None,
    *,
    stop_on_error: bool = True,
) -> list[StageResult]:
    """Run multiple stages in order."""
    to_run = list(stages) if stages else list(STAGES)
    stage_args = stage_args or {}
    results: list[StageResult] = []

    for stage in to_run:
        extra = stage_args.get(stage, [])
        result = run_stage(config, stage, extra)
        results.append(result)
        if stop_on_error and not result.ok:
            break
    return results


def format_stage_result(result: StageResult) -> str:
    lines = [
        f"=== {STAGE_NAMES.get(result.stage, result.stage)} ===",
        f"Script: {result.script}",
        f"Exit code: {result.returncode}",
        "",
        result.result.get("text", ""),
    ]
    return "\n".join(lines)
