#!/usr/bin/env python3
"""Run the frozen map-only cross-source QC on all 21 validated SANY cells."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from statistics import median
from typing import Any


EXPECTED_IDS_BY_VARIANT = {
    "C0_single114_noise": (0,),
    "C1_four_noise": (0, 1, 2, 3),
    "C2_four_no_noise": (0, 1, 2, 3),
    "C3_drop114_lidar": (1, 2, 3),
    "C4_drop127": (0, 2, 3),
    "C5_drop187": (0, 1, 3),
    "C6_drop195": (0, 1, 2),
}


def one_run(analyzer: Path, map_path: Path, output: Path, expected_ids: tuple[int, ...]) -> dict[str, Any]:
    result = subprocess.run(
        [
            sys.executable, str(analyzer), str(map_path), "--output", str(output),
            "--expected-ids", ",".join(str(value) for value in expected_ids),
        ],
        text=True, capture_output=True,
    )
    if not output.is_file():
        raise RuntimeError(f"map QC produced no result for {map_path}: {result.stderr[-2000:]}")
    data = json.loads(output.read_text(encoding="utf-8"))
    if result.returncode not in (0, 2):
        raise RuntimeError(f"map QC crashed for {map_path}: {result.stderr[-2000:]}")
    if bool(data["contract"]["passed"]) != (result.returncode == 0):
        raise RuntimeError(f"map QC exit status disagrees with contract for {map_path}")
    return data


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--experiment-root", type=Path, required=True)
    parser.add_argument(
        "--analyzer", type=Path,
        default=Path(__file__).resolve().parent / "analyze_sany_map_cross_source.py",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4:
        raise SystemExit("workers must be within [1,4]")
    jobs = []
    for cell_dir in sorted((args.experiment_root / "runs").glob("*/repeat_[0-9][0-9]")):
        selection_path = cell_dir / "selected_attempt.json"
        if not selection_path.is_file():
            continue
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
        attempt = Path(selection["selected_attempt"])
        manifest = json.loads((attempt / "completion_manifest.json").read_text(encoding="utf-8"))
        if not manifest["output_validation"]["passed"]:
            continue
        variant = manifest["variant"]
        if variant not in EXPECTED_IDS_BY_VARIANT:
            raise SystemExit(f"unknown variant has no lidar_id contract: {variant}")
        repeat = int(manifest["repeat"])
        output = args.output_dir / "per_run" / f"{variant}_repeat_{repeat:02d}.json"
        jobs.append((variant, repeat, attempt / "results" / "map_lio.pcd", output))
    if len(jobs) != 21:
        raise SystemExit(f"expected 21 validated maps, found {len(jobs)}")
    (args.output_dir / "per_run").mkdir(parents=True, exist_ok=True)
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                one_run, args.analyzer, map_path, output, EXPECTED_IDS_BY_VARIANT[variant]
            ): (variant, repeat, map_path, output)
            for variant, repeat, map_path, output in jobs
        }
        for future in as_completed(futures):
            variant, repeat, map_path, output = futures[future]
            data = future.result()
            aggregate = data["aggregate_symmetric_point_to_plane_m"]
            coverage = data["coverage"]
            rows.append({
                "variant": variant, "repeat": repeat, "map": str(map_path), "result": str(output),
                "contract_passed": data["contract"]["passed"],
                "contract_errors": data["contract"]["errors"],
                "expected_ids": data["lidar_id_contract"]["expected_ids"],
                "missing_expected_ids": data["lidar_id_contract"]["missing_expected_ids"],
                "unexpected_ids": data["lidar_id_contract"]["unexpected_ids"],
                "raw_point_count": data["raw_point_count"],
                "nonfinite_xyz_point_count": data["nonfinite_xyz_point_count"],
                "invalid_lidar_id_point_count": data["invalid_lidar_id_point_count"],
                "cross_source_metric_status": data["cross_source_metric"]["status"],
                "residual_count": aggregate["count"],
                "residual_median_m": aggregate["median"],
                "residual_p95_m": aggregate["p95"],
                "coverage_at_least_2": coverage["fraction_observed_by_at_least_2_lidars"],
                "coverage_at_least_3": coverage["fraction_observed_by_at_least_3_lidars"],
                "coverage_all_4": coverage["fraction_observed_by_all_4_lidars"],
            })
    rows.sort(key=lambda row: (row["variant"], row["repeat"]))
    variants = {}
    for variant in sorted({row["variant"] for row in rows}):
        selected = [row for row in rows if row["variant"] == variant]
        passed = [row for row in selected if row["contract_passed"]]
        valid_residual = [row for row in passed if row["residual_median_m"] is not None]
        variants[variant] = {
            "runs": len(selected),
            "contract_passed_runs": sum(row["contract_passed"] for row in selected),
            "raw_point_count_median": median(row["raw_point_count"] for row in selected),
            "residual_median_m_median": median(row["residual_median_m"] for row in valid_residual) if valid_residual else None,
            "residual_p95_m_median": median(row["residual_p95_m"] for row in valid_residual) if valid_residual else None,
            "coverage_at_least_2_median": median(row["coverage_at_least_2"] for row in passed) if passed else None,
            "coverage_at_least_3_median": median(row["coverage_at_least_3"] for row in passed) if passed else None,
            "coverage_all_4_median": median(row["coverage_all_4"] for row in passed) if passed else None,
        }
    summary = {
        "runs": rows, "variants": variants,
        "contract_passed_runs": sum(row["contract_passed"] for row in rows),
        "interpretation_limit": "跨源平面一致性与覆盖率 QC，不是绝对地图精度。",
    }
    (args.output_dir / "sany_map_qc_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"runs": len(rows), "variants": variants}, ensure_ascii=False))
    return 0 if summary["contract_passed_runs"] == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
