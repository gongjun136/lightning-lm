#!/usr/bin/env python3
"""Freeze the selected passing SANY online/fault evidence and retained failures."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def identity(path: Path) -> dict:
    return {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--clock-reanalysis", type=Path)
    args = parser.parse_args()
    selected = {
        "missing_imu114_all": args.root / "fault_runs_final" / "missing_imu114_all",
        "interrupt_imu114_40_50": args.root / "fault_runs_final" / "interrupt_imu114_40_50",
        "interrupt_lidar127_40_60": args.root / "fault_runs_final" / "interrupt_lidar127_40_60",
        "interrupt_lidar114_40_60": args.root / "fault_runs_final_retry_clock" / "interrupt_lidar114_40_60",
    }
    scenarios = {}
    for name, directory in selected.items():
        metadata_path = directory / "run_metadata.json"
        contract_path = directory / "fault_recovery_contract.json"
        metadata = load(metadata_path)
        contract = load(contract_path)
        if metadata.get("status") != "completed" or contract.get("contract_passed") is not True:
            raise SystemExit(f"selected evidence is not passing: {name}")
        scenarios[name] = {
            "directory": str(directory),
            "experiment_fingerprint": metadata["experiment_fingerprint"],
            "exit_codes": metadata["exit_codes"],
            "contract_passed": True,
            "realtime_contract_passed": contract.get("realtime_contract_passed"),
            "message_counts": contract["message_counts"],
            "receive_rates_hz": contract["receive_rates_hz"],
            "details": contract["details"],
            "files": [identity(metadata_path), identity(contract_path)],
        }
    retained_failure = args.root / "fault_runs_final" / "interrupt_lidar114_40_60"
    failed_contract = retained_failure / "fault_recovery_contract.json"
    failure_payload = load(failed_contract)
    clock_reanalysis = load(args.clock_reanalysis) if args.clock_reanalysis else None
    result = {
        "selected_fault_scenarios": scenarios,
        "retained_nonselected_failure": {
            "directory": str(retained_failure),
            "contract_failures": failure_payload["contract_failures"],
            "maximum_clock_backward_jump_s": (
                clock_reanalysis["details"].get("maximum_clock_backward_jump_s")
                if clock_reanalysis else failure_payload["details"].get("maximum_clock_backward_jump_s")
            ),
            "file": identity(failed_contract),
            "diagnostic_reanalysis": identity(args.clock_reanalysis) if args.clock_reanalysis else None,
        },
        "normal_online": {
            "directory": str(args.root / "online_contract_normal_attempt4"),
            "contract": identity(args.root / "online_contract_normal_attempt4" / "topic_contract.json"),
            "metadata": identity(args.root / "online_contract_normal_attempt4" / "run_metadata.json"),
        },
        "offline": {
            "posthoc_audit": identity(args.root / "results" / "posthoc_strict_audit.json"),
            "trajectory_proxy": identity(args.root / "results" / "trajectory_proxy" / "sany_proxy_trajectory_evaluation.json"),
            "map_qc": identity(args.root / "results" / "map_qc" / "sany_map_qc_summary.json"),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"scenarios": len(scenarios), "passed": True}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
