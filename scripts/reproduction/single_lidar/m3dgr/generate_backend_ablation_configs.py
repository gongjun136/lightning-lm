#!/usr/bin/env python3
"""Generate self-contained M3DGR backend ablation YAML files."""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path

import yaml


VARIANTS = {
    "frontend_only": {"mode": "disabled"},
    "frontend_velocity_propagation": {
        "mode": "disabled",
        "fasterlio": {
            "propagate_velocity": True,
            "adaptive_velocity_propagation": False,
            "max_update_velocity_step": 0.5,
        },
    },
    "frontend_velocity_propagation_relaxed": {
        "mode": "disabled",
        "fasterlio": {
            "propagate_velocity": True,
            "adaptive_velocity_propagation": False,
            "max_update_velocity_step": 2.0,
        },
    },
    "frontend_innovation_adaptive": {
        "mode": "disabled",
        "fasterlio": {
            "propagate_velocity": False,
            "max_update_velocity_step": 2.0,
            "adaptive_velocity_propagation": True,
            "velocity_innovation_ema_alpha": 0.05,
            "velocity_innovation_enable_threshold": 0.16,
            "velocity_innovation_disable_threshold": 0.10,
            "velocity_propagation_max_active_updates": 200,
            "velocity_propagation_cooldown_updates": 100,
        },
    },
    "frontend_innovation_adaptive_wide": {
        "mode": "disabled",
        "fasterlio": {
            "propagate_velocity": False,
            "max_update_velocity_step": 2.0,
            "adaptive_velocity_propagation": True,
            "velocity_innovation_ema_alpha": 0.05,
            "velocity_innovation_enable_threshold": 0.14,
            "velocity_innovation_disable_threshold": 0.08,
            "velocity_propagation_max_active_updates": 200,
            "velocity_propagation_cooldown_updates": 100,
        },
    },
    "legacy": {
        "mode": "legacy",
        "fasterlio": {
            "propagate_velocity": False,
            "adaptive_velocity_propagation": False,
        },
    },
    "ba_btc_hba": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": True,
    },
    "ba_btc_hba_velocity_propagation": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": True,
        "fasterlio": {
            "propagate_velocity": True,
            "adaptive_velocity_propagation": False,
            "max_update_velocity_step": 0.5,
        },
    },
    "btc_context20": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": True,
        "btc_overrides": {
            "descriptor_submap_size": 20,
        },
    },
    "btc_safe_drift": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": True,
        "btc_overrides": {
            "max_drift_ratio": 0.005,
        },
    },
    "no_local_ba": {
        "mode": "ba_btc_hba",
        "local_ba": False,
        "btc": True,
        "hba": True,
    },
    "no_btc": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": False,
        "hba": True,
    },
    "no_hba": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": False,
    },
    "btc_high_recall": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": True,
    },
    "btc_high_recall_no_hba": {
        "mode": "ba_btc_hba",
        "local_ba": True,
        "btc": True,
        "hba": False,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--variants",
        nargs="+",
        choices=tuple(VARIANTS),
        default=tuple(VARIANTS),
    )
    return parser.parse_args()


def set_variant(base: dict, name: str) -> dict:
    config = deepcopy(base)
    backend = config.setdefault("backend", {})
    settings = VARIANTS[name]
    backend["mode"] = settings["mode"]
    for module in ("local_ba", "btc", "hba"):
        if module in settings:
            backend.setdefault(module, {})["enabled"] = settings[module]
    for key, value in settings.get("btc_overrides", {}).items():
        backend.setdefault("btc", {})[key] = value
    for key, value in settings.get("fasterlio", {}).items():
        config.setdefault("fasterlio", {})[key] = value
    if name in {"btc_high_recall", "btc_high_recall_no_hba"}:
        btc = backend.setdefault("btc", {})
        btc["min_loop_score"] = 0.29
        btc["downsample_leaf_size"] = 0.0
        btc["plane_icp_min_observability"] = 3.0
        btc["allow_degenerate_plane_icp"] = True
        btc["degenerate_min_loop_score"] = 0.29
        btc["degenerate_min_matches"] = 60
        btc["confirmation_count"] = 2
        btc["confirmation_max_current_gap"] = 10
        btc["confirmation_max_history_gap"] = 10
        btc["loop_cooldown_descriptors"] = 10
        btc["max_odom_revisit_distance"] = 0.0
        btc["min_optimization_translation"] = 0.20
        btc["min_optimization_rotation_deg"] = 2.0
        descriptor = btc.setdefault("descriptor", {})
        descriptor["useful_corner_num"] = 200
        descriptor["candidate_count"] = 50
        descriptor["rough_distance_threshold"] = 0.03
        descriptor["similarity_threshold"] = 0.50
    config.setdefault("system", {})["with_loop_closing"] = settings["mode"] != "disabled"
    config.setdefault("common", {})["benchmark_backend_variant"] = name
    return config


def main() -> None:
    args = parse_args()
    with args.base.open("r", encoding="utf-8") as stream:
        base = yaml.safe_load(stream)
    if not isinstance(base, dict) or "backend" not in base:
        raise SystemExit(f"base YAML does not contain a backend mapping: {args.base}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for name in args.variants:
        output = args.output_dir / f"lightning_m3dgr_{name}.yaml"
        with output.open("w", encoding="utf-8") as stream:
            stream.write(
                "# Generated by generate_backend_ablation_configs.py; "
                f"variant={name}\n"
            )
            yaml.safe_dump(
                set_variant(base, name),
                stream,
                allow_unicode=True,
                sort_keys=False,
                default_flow_style=False,
            )
        print(output.resolve())


if __name__ == "__main__":
    main()
