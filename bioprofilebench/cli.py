"""BioProfileBench command-line interface."""

from __future__ import annotations

import argparse
import sys

from .batch import run_batch, run_single
from .io import load_config
from .legacy import format_existing_bk_results


def add_common_metadata_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--profile-id", default="profile_001", help="Prediction profile ID for single-run mode.")
    parser.add_argument("--dataset", default="Dataset1", help="Dataset label for single-run mode.")
    parser.add_argument("--method", default="Prediction", help="Method label for single-run mode.")
    parser.add_argument("--kind", default=None, help="Profile kind label, e.g. ARGs, Bac, Bacteria, VFs.")
    parser.add_argument("--gene-kind", default=None, help="Deprecated alias for --kind; kept for backward compatibility.")
    parser.add_argument("--parameter-tag", default="default", help="Parameter tag for single-run mode.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="BioProfileBench",
        description="Unified benchmarking for host-resolved microbiome abundance profiles.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run = subparsers.add_parser("run", help="Benchmark one prediction profile against one truth profile.")
    run.add_argument("--truth", required=True, help="Truth abundance TSV.")
    run.add_argument("--pred", required=True, help="Prediction abundance TSV.")
    run.add_argument("--config", required=True, help="YAML/JSON config file.")
    run.add_argument("--out", required=True, help="Output directory or benchmark_long TSV path.")
    add_common_metadata_args(run)

    batch = subparsers.add_parser("batch", help="Benchmark profiles listed in a manifest TSV.")
    batch.add_argument("--manifest", required=True, help="Batch manifest TSV.")
    batch.add_argument("--config", required=True, help="YAML/JSON config file.")
    batch.add_argument("--out", required=True, help="Output directory.")
    batch.add_argument("--threads", type=int, default=1, help="Number of parallel worker threads.")

    optimize = subparsers.add_parser("optimize", help="Alias for batch mode with optimized-filter reporting.")
    optimize.add_argument("--manifest", required=True, help="Batch manifest TSV.")
    optimize.add_argument("--config", required=True, help="YAML/JSON config file.")
    optimize.add_argument("--out", required=True, help="Output directory.")
    optimize.add_argument("--threads", type=int, default=1, help="Number of parallel worker threads.")

    legacy = subparsers.add_parser("format-bk", help="Convert existing .BK result files to an R-compatible table.")
    legacy.add_argument("--input", nargs="+", required=True, help="Existing .BK file(s) or directories.")
    legacy.add_argument("--out", required=True, help="Output AllBK_filter-compatible TSV.")
    legacy.add_argument("--config", help="Optional config for digits and sample group-size metadata.")

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        config = load_config(args.config) if getattr(args, "config", None) else {}
        if args.command == "run":
            paths = run_single(args, config)
        elif args.command in {"batch", "optimize"}:
            paths = run_batch(args, config)
        elif args.command == "format-bk":
            out = format_existing_bk_results(args, config)
            print("BioProfileBench completed.", file=sys.stderr)
            print(f"allbk_filter_compatible: {out}", file=sys.stderr)
            return 0
        else:
            parser.error(f"Unsupported command: {args.command}")
            return 2
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print("BioProfileBench completed.", file=sys.stderr)
    if "benchmark" in paths:
        print(f"benchmark: {paths['benchmark']}", file=sys.stderr)
    print(f"filter_trace: {paths['filter_trace']}", file=sys.stderr)
    if "pair_filter_ranking" in paths:
        print(f"pair_filter_ranking: {paths['pair_filter_ranking']}", file=sys.stderr)
    print(f"benchmark_long: {paths['benchmark_long']}", file=sys.stderr)
    return 0
