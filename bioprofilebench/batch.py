"""Single and batch runners for BioProfileBench."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import pandas as pd

from .benchmark import (
    benchmark_profile_pair,
    build_method_comparison,
    select_best_filters,
    summarize_benchmark,
)
from .bk_metadata import build_r_compatible_table
from .io import build_profile, load_manifest, write_table
from .report import (
    build_benchmark_summary,
    build_filter_trace_table,
    build_pair_filter_level_ranking,
    build_pair_filter_ranking,
    build_simple_benchmark_table,
    prepare_public_table,
    split_summary_by_level,
)


def output_paths(out: str | Path) -> dict[str, Path]:
    """Resolve standard output paths from an output directory or main TSV path."""
    out = Path(out)
    if out.suffix.lower() == ".tsv":
        base = out.with_suffix("")
        return {
            "benchmark": out,
            "benchmark_long": base.with_suffix(".benchmark_long.tsv"),
            "benchmark_summary": base.with_suffix(".summary.tsv"),
            "filter_trace": base.with_suffix(".filter_trace.tsv"),
            "filter_diagnostics": base.with_suffix(".filter_diagnostics.tsv"),
            "filter_config_resolved": base.with_suffix(".filter_config_resolved.tsv"),
            "pair_filter_ranking": base.with_suffix(".pair_filter_ranking.tsv"),
            "pair_filter_level_ranking": base.with_suffix(".pair_filter_level_ranking.tsv"),
            "best_filters": base.with_suffix(".best_filters.tsv"),
            "method_comparison": base.with_suffix(".method_comparison.tsv"),
            "failed_or_skipped_runs": base.with_suffix(".failed_or_skipped_runs.tsv"),
            "run_metadata": base.with_suffix(".run_metadata.tsv"),
            "fp_fn_decomposition": base.with_suffix(".fp_fn_decomposition.tsv"),
            "allbk_filter_compatible": base.with_suffix(".allbk_filter_compatible.tsv"),
        }
    return {
        "benchmark": out / "benchmark.tsv",
        "benchmark_long": out / "benchmark_long.tsv",
        "benchmark_summary": out / "benchmark_summary.tsv",
        "filter_trace": out / "filter_trace.tsv",
        "filter_diagnostics": out / "filter_diagnostics.tsv",
        "filter_config_resolved": out / "filter_config_resolved.tsv",
        "pair_filter_ranking": out / "pair_filter_ranking.tsv",
        "pair_filter_level_ranking": out / "pair_filter_level_ranking.tsv",
        "best_filters": out / "best_filters.tsv",
        "method_comparison": out / "method_comparison.tsv",
        "failed_or_skipped_runs": out / "failed_or_skipped_runs.tsv",
        "run_metadata": out / "run_metadata.tsv",
        "fp_fn_decomposition": out / "fp_fn_decomposition.tsv",
        "allbk_filter_compatible": out / "allbk_filter_compatible.tsv",
    }




def summary_level_path(summary_path: Path, level: str) -> Path:
    """Return benchmark_summary.<level>.tsv beside the main summary file."""
    if summary_path.name == "benchmark_summary.tsv":
        return summary_path.with_name(f"benchmark_summary.{level}.tsv")
    return summary_path.with_suffix(f".{level}.tsv")


def build_single_profiles(args: Any) -> tuple[Any, Any]:
    truth = build_profile(
        args.truth,
        profile_id=f"{args.profile_id}_truth",
        dataset=args.dataset,
        method="Truth",
        gene_kind=args.gene_kind,
        parameter_tag=args.parameter_tag,
    )
    prediction = build_profile(
        args.pred,
        profile_id=args.profile_id,
        dataset=args.dataset,
        method=args.method,
        gene_kind=args.gene_kind,
        parameter_tag=args.parameter_tag,
        extra_metadata={"prediction_path": Path(args.pred).name, "truth_path": Path(args.truth).name},
    )
    return truth, prediction


def write_result_bundle(
    out: str | Path,
    config: dict[str, Any],
    benchmark: pd.DataFrame,
    trace: pd.DataFrame,
    filter_configs: pd.DataFrame,
    skipped: pd.DataFrame,
    decomposition: pd.DataFrame,
    run_metadata: pd.DataFrame,
) -> dict[str, Path]:
    """Write all standard output files."""
    paths = output_paths(out)
    digits = int(config.get("benchmark", {}).get("digits", 6))
    internal_summary = summarize_benchmark(benchmark)
    best = select_best_filters(internal_summary, config)
    method_comparison = build_method_comparison(internal_summary)
    simple_benchmark = build_simple_benchmark_table(benchmark)
    filter_trace = build_filter_trace_table(filter_configs)
    benchmark_summary = build_benchmark_summary(simple_benchmark, filter_trace)
    pair_ranking = build_pair_filter_ranking(simple_benchmark, filter_trace)
    pair_level_ranking = build_pair_filter_level_ranking(simple_benchmark, filter_trace)
    report_config = config.get("report", {})
    sample_config = config.get("sample_metadata", {})
    group_size_map = report_config.get("group_size_map", sample_config.get("group_size_map", {}))
    r_compatible = build_r_compatible_table(benchmark, group_size_map)
    failed = pd.concat(
        [
            skipped,
            benchmark.loc[benchmark["Status"].isin(["SKIPPED", "EMPTY_PREDICTION", "EMPTY_BOTH"])]
            if not benchmark.empty and "Status" in benchmark.columns
            else pd.DataFrame(),
        ],
        ignore_index=True,
    )

    write_table(prepare_public_table(simple_benchmark), paths["benchmark"], digits)
    write_table(prepare_public_table(benchmark), paths["benchmark_long"], digits)
    write_table(prepare_public_table(benchmark_summary), paths["benchmark_summary"], digits)
    for level, level_summary in split_summary_by_level(benchmark_summary).items():
        write_table(prepare_public_table(level_summary), summary_level_path(paths["benchmark_summary"], level), digits)
    write_table(prepare_public_table(filter_trace), paths["filter_trace"], digits)
    write_table(prepare_public_table(trace), paths["filter_diagnostics"], digits)
    write_table(prepare_public_table(filter_configs), paths["filter_config_resolved"], digits)
    write_table(prepare_public_table(pair_ranking), paths["pair_filter_ranking"], digits)
    write_table(prepare_public_table(pair_level_ranking), paths["pair_filter_level_ranking"], digits)
    write_table(prepare_public_table(best), paths["best_filters"], digits)
    write_table(prepare_public_table(method_comparison), paths["method_comparison"], digits)
    write_table(prepare_public_table(failed), paths["failed_or_skipped_runs"], digits)
    write_table(prepare_public_table(run_metadata), paths["run_metadata"], digits)
    write_table(prepare_public_table(decomposition), paths["fp_fn_decomposition"], digits)
    write_table(prepare_public_table(r_compatible), paths["allbk_filter_compatible"], digits)
    return paths


def run_single(args: Any, config: dict[str, Any]) -> dict[str, Path]:
    truth, prediction = build_single_profiles(args)
    benchmark, trace, filter_configs, skipped, decomposition = benchmark_profile_pair(truth, prediction, config)
    run_metadata = pd.DataFrame([{**prediction.metadata, "truth_path": Path(args.truth).name, "prediction_path": Path(args.pred).name}])
    return write_result_bundle(args.out, config, benchmark, trace, filter_configs, skipped, decomposition, run_metadata)


def row_to_prediction_profile(row: pd.Series):
    known = {"prediction_path", "truth_path"}
    extra_metadata = {key: value for key, value in row.to_dict().items() if key not in known}
    extra_metadata["prediction_path"] = Path(row["prediction_path"]).name
    extra_metadata["truth_path"] = Path(row["truth_path"]).name
    return build_profile(
        row["prediction_path"],
        profile_id=row["profile_id"],
        dataset=row.get("dataset", ""),
        method=row.get("method", ""),
        gene_kind=row.get("gene_kind", ""),
        parameter_tag=row.get("parameter_tag", ""),
        extra_metadata=extra_metadata,
    )


def row_to_truth_profile(row: pd.Series, cache: dict[str, Any]):
    truth_path = row["truth_path"]
    if truth_path not in cache:
        cache[truth_path] = build_profile(
            truth_path,
            profile_id=f"truth_{Path(truth_path).stem}",
            dataset=row.get("dataset", ""),
            method="Truth",
            gene_kind=row.get("gene_kind", ""),
            parameter_tag="",
        )
    return cache[truth_path]


def benchmark_manifest_row(row: pd.Series, truth_profile: Any, config: dict[str, Any]):
    prediction = row_to_prediction_profile(row)
    return benchmark_profile_pair(truth_profile, prediction, config)


def run_batch(args: Any, config: dict[str, Any]) -> dict[str, Path]:
    manifest = load_manifest(args.manifest)
    truth_cache: dict[str, Any] = {}
    for _, row in manifest.iterrows():
        row_to_truth_profile(row, truth_cache)

    results = []
    threads = max(1, int(getattr(args, "threads", 1)))
    if threads == 1:
        for _, row in manifest.iterrows():
            truth = row_to_truth_profile(row, truth_cache)
            results.append(benchmark_manifest_row(row, truth, config))
    else:
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = []
            for _, row in manifest.iterrows():
                truth = row_to_truth_profile(row, truth_cache)
                futures.append(executor.submit(benchmark_manifest_row, row.copy(), truth, config))
            for future in as_completed(futures):
                results.append(future.result())

    benchmark = pd.concat([item[0] for item in results], ignore_index=True) if results else pd.DataFrame()
    trace = pd.concat([item[1] for item in results], ignore_index=True) if results else pd.DataFrame()
    filter_configs = (
        pd.concat([item[2] for item in results], ignore_index=True).drop_duplicates()
        if results
        else pd.DataFrame()
    )
    skipped = pd.concat([item[3] for item in results], ignore_index=True) if results else pd.DataFrame()
    decomposition = pd.concat([item[4] for item in results], ignore_index=True) if results else pd.DataFrame()
    run_metadata = manifest.copy()
    return write_result_bundle(args.out, config, benchmark, trace, filter_configs, skipped, decomposition, run_metadata)
