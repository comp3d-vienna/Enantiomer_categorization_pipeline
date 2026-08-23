"""Command-line interface for the enantiomer pipeline package."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Optional, Sequence

from . import __version__
from .config import PipelineConfig
from .runner import STAGES, format_stage_result, run_all, run_stage


def _build_config(args: argparse.Namespace) -> PipelineConfig:
    return PipelineConfig(
        project_root=args.project_root,
        test_mode=args.test,
        mmcif_source=args.mmcif_source,
        python_cmd=args.python_cmd,
        schrodinger=args.schrodinger,
    )


def _print_result(stage_result, *, as_json: bool = False) -> None:
    if as_json:
        payload = {
            "stage": stage_result.stage,
            "returncode": stage_result.returncode,
            "script": stage_result.script,
            **stage_result.result,
        }
        print(json.dumps(payload, indent=2, default=str))
    else:
        print(format_stage_result(stage_result))


def _cmd_run(args: argparse.Namespace) -> int:
    config = _build_config(args)
    extra_args: list[str] = list(args.extra_args or [])

    if args.stage == "all":
        stage_args = {s: extra_args for s in STAGES} if extra_args else None
        results = run_all(config, stage_args=stage_args, stop_on_error=not args.continue_on_error)
        if args.json:
            payload = [
                {
                    "stage": r.stage,
                    "returncode": r.returncode,
                    "script": r.script,
                    **r.result,
                }
                for r in results
            ]
            print(json.dumps(payload, indent=2, default=str))
        else:
            for result in results:
                print(format_stage_result(result))
        failed = [r for r in results if not r.ok]
        if failed:
            if not args.json:
                print(f"\nPipeline stopped with {len(failed)} failed stage(s).", file=sys.stderr)
            return failed[0].returncode or 1
        if not args.json:
            print("\nFull pipeline completed successfully.")
        return 0

    result = run_stage(config, args.stage, extra_args)
    _print_result(result, as_json=args.json)
    return result.returncode


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="enantiomer-pipeline",
        description=(
            "Run the enantiomer binding-mode pipeline "
            "(4 stages: preprocess, alignment, pharmacophore, categorization). "
            "Each run prints a result summary when the stage finishes."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--project-root",
        type=str,
        default=None,
        help="Path to Organize_script repository root (default: auto-detect).",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Use test/ outputs and test/run_*_test.sh drivers (sets PIPELINE_TEST=1).",
    )
    parser.add_argument(
        "--mmcif-source",
        dest="mmcif_source",
        default=None,
        help="PDB mmCIF archive root for preprocess (MMCIF_SOURCE; default: Data/mmCIF).",
    )
    parser.add_argument(
        "--python-cmd",
        dest="python_cmd",
        default=None,
        help="Python executable for pipeline scripts (PYTHON_CMD).",
    )
    parser.add_argument(
        "--schrodinger",
        default=None,
        help="Schrödinger installation root (SCHRODINGER).",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Print stage result as JSON (includes metrics and output paths).",
    )
    parser.add_argument(
        "stage",
        choices=[*STAGES, "all"],
        help="Pipeline stage to run, or 'all' for the full workflow.",
    )
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="When running all stages, continue after a failed stage.",
    )
    parser.add_argument(
        "extra_args",
        nargs=argparse.REMAINDER,
        help="Extra flags passed to the underlying bash driver (e.g. -- --skip-uniprot --skip-prepwizard).",
    )
    parser.set_defaults(func=_cmd_run)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.project_root is None:
        from .config import _default_project_root

        args.project_root = _default_project_root()

    if getattr(args, "extra_args", None):
        if args.extra_args and args.extra_args[0] == "--":
            args.extra_args = args.extra_args[1:]

    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
