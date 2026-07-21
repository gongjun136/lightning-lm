#!/usr/bin/env python3
"""Validate completeness, provenance and share-readiness of the formal suite."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ANALYSES = {
    "m3dgr_frontend": 48,
    # The backend evaluator reuses the 12 Lightning frontend runs so that
    # frontend/legacy/new/full trajectories share one accuracy support.
    # These are analysis rows, not 12 additional algorithm executions.
    "m3dgr_backend": 48,
    "m3dgr_localization": 24,
    "sany_20260701_mapping": 27,
    "sany_20260716_relocalization": 12,
}
EXPECTED_FIGURES = 10
FRONTEND_COMPUTE_RUNS = 48
EXPECTED_COMPUTE_FIGURES = 4
STATISTICAL_RUN_COUNT = 195


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def parse_metadata(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    return values


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    formal = repo / "runs" / "formal_report_20260717"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--formal-root", type=Path, default=formal)
    parser.add_argument(
        "--figure-dir", type=Path, default=repo / "doc" / "assets" / "formal_report_20260717",
    )
    parser.add_argument(
        "--report", type=Path, default=repo / "doc" / "SLAM建图定位算法报告.md",
    )
    parser.add_argument("--require-final-report", action="store_true")
    args = parser.parse_args()

    failures: list[str] = []
    warnings: list[str] = []
    evidence: dict[str, object] = {}
    analysis_root = args.formal_root / "analysis"
    for name, expected_runs in ANALYSES.items():
        directory = analysis_root / name
        required = [
            directory / "validation.json",
            directory / "run_metrics.csv",
            directory / "summary_metrics.csv",
            directory / "aligned_samples.csv",
            directory / "metrics.json",
            directory / "analysis_manifest.json",
        ]
        if name == "m3dgr_backend":
            required.extend([
                directory / "backend_activity.csv",
                directory / "backend_activity_summary.csv",
            ])
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            failures.extend(f"missing analysis artifact: {path}" for path in missing)
            continue
        validation = json.loads(required[0].read_text(encoding="utf-8"))
        run_rows = read_csv(required[1])
        summary_rows = read_csv(required[2])
        if validation.get("status") not in {"passed", "passed_with_warnings"}:
            failures.append(f"{name}: validation status={validation.get('status')}")
        if validation.get("required_run_count") != expected_runs:
            failures.append(f"{name}: required_run_count={validation.get('required_run_count')} expected={expected_runs}")
        if validation.get("evaluated_run_count") != expected_runs or len(run_rows) != expected_runs:
            failures.append(
                f"{name}: evaluated={validation.get('evaluated_run_count')} csv_rows={len(run_rows)} expected={expected_runs}"
            )
        bad_repeats = [row for row in summary_rows if row.get("successful_repeats") != "3"]
        if bad_repeats:
            failures.append(f"{name}: {len(bad_repeats)} summary cells do not contain 3 successful repeats")
        for index, row in enumerate(run_rows, 2):
            for field in ("ate_rmse_m", "realtime_factor", "peak_rss_mb"):
                try:
                    value = float(row[field])
                except (KeyError, TypeError, ValueError):
                    failures.append(f"{name}: run_metrics.csv:{index} invalid {field}")
                    continue
                if not math.isfinite(value) or value < 0:
                    failures.append(f"{name}: run_metrics.csv:{index} non-finite/negative {field}")
        if required[3].stat().st_size <= 100:
            failures.append(f"{name}: aligned_samples.csv is empty")
        analysis_warnings = validation.get("warnings", [])
        warnings.extend(f"{name}: {item}" for item in analysis_warnings)
        evidence[name] = {
            "validation_status": validation.get("status"),
            "run_count": len(run_rows),
            "summary_cell_count": len(summary_rows),
            "warning_count": len(analysis_warnings),
            "artifact_sha256": {path.name: sha256(path) for path in required},
        }

    compute_name = "m3dgr_frontend_compute"
    compute_directory = analysis_root / compute_name
    compute_required = [
        compute_directory / "validation.json",
        compute_directory / "run_timing_summary.csv",
        compute_directory / "sequence_method_summary.csv",
        compute_directory / "frame_timing.csv",
        compute_directory / "metrics.json",
        compute_directory / "analysis_manifest.json",
    ]
    compute_missing = [str(path) for path in compute_required if not path.is_file()]
    if compute_missing:
        failures.extend(f"missing analysis artifact: {path}" for path in compute_missing)
    else:
        validation = json.loads(compute_required[0].read_text(encoding="utf-8"))
        run_rows = read_csv(compute_required[1])
        summary_rows = read_csv(compute_required[2])
        if validation.get("status") not in {"passed", "passed_with_warnings"}:
            failures.append(f"{compute_name}: validation status={validation.get('status')}")
        if validation.get("expected_run_count") != FRONTEND_COMPUTE_RUNS:
            failures.append(
                f"{compute_name}: expected_run_count={validation.get('expected_run_count')} "
                f"expected={FRONTEND_COMPUTE_RUNS}"
            )
        if validation.get("evaluated_run_count") != FRONTEND_COMPUTE_RUNS or len(run_rows) != FRONTEND_COMPUTE_RUNS:
            failures.append(
                f"{compute_name}: evaluated={validation.get('evaluated_run_count')} "
                f"csv_rows={len(run_rows)} expected={FRONTEND_COMPUTE_RUNS}"
            )
        bad_repeats = [row for row in summary_rows if row.get("successful_repeats") != "3"]
        if bad_repeats:
            failures.append(
                f"{compute_name}: {len(bad_repeats)} summary cells do not contain 3 successful repeats"
            )
        for index, row in enumerate(run_rows, 2):
            for field in ("total_ms_mean", "total_ms_p95", "timing_coverage", "repeat_frame_coverage"):
                try:
                    value = float(row[field])
                except (KeyError, TypeError, ValueError):
                    failures.append(f"{compute_name}: run_timing_summary.csv:{index} invalid {field}")
                    continue
                if not math.isfinite(value) or value < 0:
                    failures.append(f"{compute_name}: run_timing_summary.csv:{index} non-finite/negative {field}")
        if compute_required[3].stat().st_size <= 100:
            failures.append(f"{compute_name}: frame_timing.csv is empty")
        compute_warnings = validation.get("warnings", [])
        warnings.extend(f"{compute_name}: {item}" for item in compute_warnings)
        evidence[compute_name] = {
            "validation_status": validation.get("status"),
            "run_count": len(run_rows),
            "summary_cell_count": len(summary_rows),
            "frame_row_count": validation.get("frame_row_count"),
            "included_frame_count": validation.get("included_frame_count"),
            "warning_count": len(compute_warnings),
            "artifact_sha256": {path.name: sha256(path) for path in compute_required},
        }

    references: dict[str, object] = {}
    reference_root = args.formal_root / "sany_voxel114_reference"
    for dataset in ("data_20260701", "data1", "data2"):
        metadata_path = reference_root / dataset / "run_metadata.txt"
        trajectory_path = reference_root / dataset / "results" / "trajectory_voxel_opt.tum"
        resources_path = reference_root / dataset / "resource_summary.json"
        if not all(path.is_file() for path in (metadata_path, trajectory_path, resources_path)):
            failures.append(f"missing SANY Voxel114 reference artifacts: {dataset}")
            continue
        metadata = parse_metadata(metadata_path)
        if metadata.get("completion") != "reached_final_lidar":
            failures.append(f"SANY Voxel114 reference incomplete: {dataset}")
        try:
            ratio = float(metadata["output_ratio"])
        except (KeyError, ValueError):
            failures.append(f"SANY Voxel114 reference invalid output ratio: {dataset}")
            ratio = math.nan
        if not (0.80 <= ratio <= 1.01):
            failures.append(f"SANY Voxel114 reference output ratio out of contract: {dataset}={ratio}")
        references[dataset] = {
            "completion": metadata.get("completion"),
            "output_ratio": ratio,
            "trajectory_sha256": sha256(trajectory_path),
            "metadata_sha256": sha256(metadata_path),
        }
    evidence["sany_voxel114_references"] = references

    environment_path = args.formal_root / "_environment" / "environment.json"
    if not environment_path.is_file():
        failures.append(f"missing environment manifest: {environment_path}")
    else:
        environment = json.loads(environment_path.read_text(encoding="utf-8"))
        for key in ("cpu", "memory", "gpu"):
            if environment.get("host", {}).get(key, {}).get("returncode") != 0:
                failures.append(f"environment host.{key} capture failed")
        for key in ("ubuntu_22_04", "ubuntu_20_04", "repository"):
            if environment.get(key, {}).get("returncode") != 0:
                failures.append(f"environment {key} capture failed")
        evidence["environment_sha256"] = sha256(environment_path)

    figure_manifest = args.figure_dir / "figure_manifest.json"
    if not figure_manifest.is_file():
        failures.append(f"missing figure manifest: {figure_manifest}")
    else:
        figures = json.loads(figure_manifest.read_text(encoding="utf-8"))
        output_hashes = figures.get("outputs", {})
        if len(output_hashes) != EXPECTED_FIGURES:
            failures.append(f"figure count={len(output_hashes)} expected={EXPECTED_FIGURES}")
        for filename, expected_hash in output_hashes.items():
            path = args.figure_dir / filename
            if not path.is_file() or sha256(path) != expected_hash:
                failures.append(f"missing or changed figure: {path}")
        evidence["figure_manifest_sha256"] = sha256(figure_manifest)

    compute_figure_manifest = args.figure_dir / "frontend_compute_figure_manifest.json"
    if not compute_figure_manifest.is_file():
        failures.append(f"missing compute figure manifest: {compute_figure_manifest}")
    else:
        compute_figures = json.loads(compute_figure_manifest.read_text(encoding="utf-8"))
        outputs = compute_figures.get("figures", {})
        if len(outputs) != EXPECTED_COMPUTE_FIGURES:
            failures.append(
                f"compute figure count={len(outputs)} expected={EXPECTED_COMPUTE_FIGURES}"
            )
        for name, artifact in outputs.items():
            try:
                path = Path(artifact["path"])
                expected_hash = artifact["sha256"]
            except (KeyError, TypeError):
                failures.append(f"invalid compute figure manifest entry: {name}")
                continue
            if not path.is_file() or sha256(path) != expected_hash:
                failures.append(f"missing or changed compute figure: {path}")
        evidence["frontend_compute_figure_manifest_sha256"] = sha256(compute_figure_manifest)

    if args.require_final_report:
        if not args.report.is_file():
            failures.append(f"missing final report: {args.report}")
        else:
            report_text = args.report.read_text(encoding="utf-8")
            if "TODO" in report_text:
                failures.append("final report still contains TODO")
            for filename in (
                "m3dgr_frontend_summary.png", "m3dgr_backend_summary.png",
                "m3dgr_localization_summary.png", "sany_mapping_summary.png",
                "sany_relocalization_summary.png", "m3dgr_frontend_compute_summary.png",
                "lightning_lio_timing_stage_summary.png", "lightning_lio_timing_distribution.png",
                "lightning_lio_timing_timeseries.png",
            ):
                if filename not in report_text:
                    failures.append(f"final report does not reference {filename}")
            evidence["report_sha256"] = sha256(args.report)

    status = "failed" if failures else ("passed_with_warnings" if warnings else "passed")
    fallacy_scan = {
        "coverage": "11/11",
        "items": [
            {"fallacy": "Simpson's paradox", "severity": "NOTE",
             "finding": "Primary results remain sequence-stratified; any aggregate is shown only with per-sequence values."},
            {"fallacy": "Ecological fallacy", "severity": "NOTE",
             "finding": "No inference from sequence-level metrics to individual frames, vehicles or users."},
            {"fallacy": "Berkson's paradox", "severity": "CAUTION",
             "finding": "The four challenge sequences were purposefully selected and do not represent the full deployment distribution."},
            {"fallacy": "Collider bias", "severity": "NOTE",
             "finding": "No regression adjustment or conditioning on a shared outcome variable is used."},
            {"fallacy": "Base-rate neglect", "severity": "NOTE",
             "finding": "The report does not estimate diagnostic sensitivity, specificity, PPV or NPV."},
            {"fallacy": "Regression to the mean", "severity": "NOTE",
             "finding": "All predeclared repeats are retained; runs are not selected because of an extreme first result."},
            {"fallacy": "Survivorship bias", "severity": "CAUTION",
             "finding": "Accuracy uses common temporal support, so output deficits, coverage and continuity warnings must be read alongside ATE/RPE."},
            {"fallacy": "Look-elsewhere effect", "severity": "NOTE",
             "finding": "Sequences, methods and metrics were frozen before this formal rerun; no p-value significance claims are made."},
            {"fallacy": "Garden of forking paths", "severity": "CAUTION",
             "finding": "The evaluated datasets were seen during earlier development; this is a frozen regression/acceptance rerun, not an unseen holdout."},
            {"fallacy": "Correlation is not causation", "severity": "NOTE",
             "finding": "Comparative claims are limited to observed behavior under the documented hardware, data and configuration."},
            {"fallacy": "Reverse causality", "severity": "NOTE",
             "finding": "No directional observational association is interpreted causally."},
        ],
    }
    result = {
        "schema_version": 2,
        "status": status,
        "required_formal_run_count": 198,
        "statistical_run_count": STATISTICAL_RUN_COUNT,
        "proxy_reference_run_count": 3,
        "failures": failures,
        "warnings": warnings,
        "methodology_fallacy_scan": fallacy_scan,
        "evidence": evidence,
    }
    output = args.formal_root / "suite_validation.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": status, "failures": len(failures), "warnings": len(warnings)}, ensure_ascii=False))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
