#!/usr/bin/env python3
"""Map-only cross-source QC using reciprocal planar correspondences and coverage."""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def lzf_decompress(data: bytes, expected_size: int) -> bytes:
    output = bytearray()
    index = 0
    while index < len(data):
        control = data[index]
        index += 1
        if control < 32:
            length = control + 1
            output.extend(data[index:index + length])
            index += length
            continue
        length = control >> 5
        reference = len(output) - ((control & 0x1F) << 8) - 1
        if length == 7:
            length += data[index]
            index += 1
        reference -= data[index]
        index += 1
        length += 2
        if reference < 0:
            raise ValueError("invalid LZF back-reference")
        for _ in range(length):
            output.append(output[reference])
            reference += 1
    if len(output) != expected_size:
        raise ValueError(f"LZF size mismatch: {len(output)} != {expected_size}")
    return bytes(output)


def load_binary_compressed_pcd(path: Path) -> dict[str, np.ndarray]:
    header: dict[str, list[str]] = {}
    with path.open("rb") as stream:
        while True:
            line = stream.readline()
            if not line:
                raise ValueError("PCD DATA header missing")
            text = line.decode("ascii").strip()
            if not text or text.startswith("#"):
                continue
            tokens = text.split()
            header[tokens[0].upper()] = tokens[1:]
            if tokens[0].upper() == "DATA":
                break
        if header["DATA"] != ["binary_compressed"]:
            raise ValueError("only binary_compressed PCD is supported")
        compressed_size, uncompressed_size = struct.unpack("<II", stream.read(8))
        compressed = stream.read(compressed_size)
    raw = lzf_decompress(compressed, uncompressed_size)
    fields = header["FIELDS"]
    sizes = [int(value) for value in header["SIZE"]]
    types = header["TYPE"]
    counts = [int(value) for value in header.get("COUNT", ["1"] * len(fields))]
    points = int(header["POINTS"][0])
    type_map = {
        ("F", 4): "<f4", ("F", 8): "<f8", ("U", 1): "u1", ("U", 2): "<u2",
        ("U", 4): "<u4", ("I", 1): "i1", ("I", 2): "<i2", ("I", 4): "<i4",
    }
    arrays: dict[str, np.ndarray] = {}
    offset = 0
    for name, size, kind, count in zip(fields, sizes, types, counts):
        length = points * size * count
        dtype = np.dtype(type_map[(kind, size)])
        array = np.frombuffer(raw, dtype=dtype, count=points * count, offset=offset)
        arrays[name] = array.reshape(points, count)[:, 0] if count > 1 else array
        offset += length
    if offset != len(raw):
        raise ValueError("PCD field blocks do not consume the uncompressed payload")
    return arrays


def voxel_sample(points: np.ndarray, resolution: float, max_points: int, seed: int) -> np.ndarray:
    keys = np.floor(points / resolution).astype(np.int64)
    _, indices = np.unique(keys, axis=0, return_index=True)
    sampled = points[np.sort(indices)]
    if len(sampled) > max_points:
        rng = np.random.default_rng(seed)
        sampled = sampled[np.sort(rng.choice(len(sampled), max_points, replace=False))]
    return sampled


def distribution(values: np.ndarray) -> dict[str, float | int | None]:
    if not len(values):
        return {"count": 0, "median": None, "p90": None, "p95": None, "max": None}
    return {
        "count": int(len(values)), "median": float(np.median(values)),
        "p90": float(np.percentile(values, 90)), "p95": float(np.percentile(values, 95)),
        "max": float(np.max(values)),
    }


def parse_expected_ids(value: str) -> tuple[int, ...]:
    try:
        expected = tuple(sorted({int(item.strip()) for item in value.split(",") if item.strip()}))
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected IDs must be comma-separated integers") from error
    if not expected or any(lidar_id not in range(4) for lidar_id in expected):
        raise argparse.ArgumentTypeError("expected IDs must be a non-empty subset of 0,1,2,3")
    return expected


def directional_planar_residual(
    source: np.ndarray, target: np.ndarray, max_distance: float, max_pairs: int, seed: int,
) -> tuple[np.ndarray, int, float]:
    if len(source) < 20 or len(target) < 20:
        return np.empty(0), 0, 0.0
    target_tree = cKDTree(target)
    source_tree = cKDTree(source)
    distances, target_indices = target_tree.query(source, k=1, workers=-1)
    candidates = np.flatnonzero(distances <= max_distance)
    if len(candidates) > max_pairs:
        rng = np.random.default_rng(seed)
        candidates = np.sort(rng.choice(candidates, max_pairs, replace=False))
    matched_target = target_indices[candidates]
    _, reverse_indices = source_tree.query(target[matched_target], k=1, workers=-1)
    reciprocal = reverse_indices == candidates
    source_indices = candidates[reciprocal]
    target_indices = matched_target[reciprocal]
    if not len(source_indices):
        return np.empty(0), 0, float(np.mean(distances <= max_distance))
    _, neighbor_indices = target_tree.query(target[target_indices], k=12, workers=-1)
    residuals = []
    planar_count = 0
    for source_index, target_index, neighbors in zip(source_indices, target_indices, neighbor_indices):
        neighborhood = target[neighbors]
        centered = neighborhood - np.mean(neighborhood, axis=0)
        eigenvalues, eigenvectors = np.linalg.eigh(centered.T @ centered / len(centered))
        if eigenvalues[-1] <= 0.0 or eigenvalues[0] / eigenvalues[-1] > 0.08:
            continue
        planar_count += 1
        normal = eigenvectors[:, 0]
        residuals.append(abs(float(np.dot(source[source_index] - target[target_index], normal))))
    return np.asarray(residuals), planar_count, float(np.mean(distances <= max_distance))


def coverage(points_by_id: dict[int, np.ndarray], resolution: float) -> dict[str, Any]:
    voxel_sets = {
        lidar_id: {tuple(row) for row in np.unique(np.floor(points / resolution).astype(np.int64), axis=0)}
        for lidar_id, points in points_by_id.items()
    }
    union = set().union(*voxel_sets.values()) if voxel_sets else set()
    multiplicity = np.asarray([sum(voxel in voxels for voxels in voxel_sets.values()) for voxel in union], dtype=np.int64)
    return {
        "resolution_m": resolution,
        "union_voxel_count": len(union),
        "per_lidar_voxel_count": {str(key): len(value) for key, value in voxel_sets.items()},
        "fraction_observed_by_at_least_2_lidars": float(np.mean(multiplicity >= 2)) if len(multiplicity) else 0.0,
        "fraction_observed_by_at_least_3_lidars": float(np.mean(multiplicity >= 3)) if len(multiplicity) else 0.0,
        "fraction_observed_by_all_4_lidars": float(np.mean(multiplicity >= 4)) if len(multiplicity) else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("map", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-ids", type=parse_expected_ids, required=True,
                        help="Comma-separated lidar_id contract, for example 0,1,3")
    parser.add_argument("--voxel-resolution", type=float, default=0.30)
    parser.add_argument("--max-correspondence-distance", type=float, default=0.50)
    parser.add_argument("--max-points-per-lidar", type=int, default=150000)
    parser.add_argument("--max-pairs", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260711)
    args = parser.parse_args()
    arrays = load_binary_compressed_pcd(args.map)
    required = {"x", "y", "z", "lidar_id"}
    if not required.issubset(arrays):
        raise SystemExit(f"map lacks fields: {sorted(required - set(arrays))}")
    xyz_original = np.column_stack((arrays["x"], arrays["y"], arrays["z"])).astype(np.float64)
    ids_original = np.asarray(arrays["lidar_id"], dtype=np.float64)
    original_point_count = len(xyz_original)
    finite_xyz = np.all(np.isfinite(xyz_original), axis=1)
    finite_id = np.isfinite(ids_original)
    integral_id = np.zeros(original_point_count, dtype=bool)
    integral_id[finite_id] = ids_original[finite_id] == np.floor(ids_original[finite_id])
    ids_integral = ids_original[integral_id].astype(np.int64)
    expected_ids = set(args.expected_ids)
    accepted = finite_xyz & integral_id & np.isin(ids_original, args.expected_ids)
    xyz = xyz_original[accepted]
    ids = ids_original[accepted].astype(np.int64)
    observed_integral_ids = set(int(value) for value in np.unique(ids_integral))
    unexpected_ids = sorted(observed_integral_ids - expected_ids)
    missing_ids = sorted(expected_ids - set(int(value) for value in np.unique(ids)))
    invalid_id_mask = (~integral_id) | (integral_id & ~np.isin(ids_original, args.expected_ids))
    points_by_id = {
        lidar_id: voxel_sample(xyz[ids == lidar_id], args.voxel_resolution, args.max_points_per_lidar, args.seed + lidar_id)
        for lidar_id in args.expected_ids
    }
    pairs: dict[str, Any] = {}
    all_residuals = []
    for left, right in combinations(args.expected_ids, 2):
        lr, lr_planar, lr_overlap = directional_planar_residual(
            points_by_id[left], points_by_id[right], args.max_correspondence_distance, args.max_pairs, args.seed + left * 10 + right,
        )
        rl, rl_planar, rl_overlap = directional_planar_residual(
            points_by_id[right], points_by_id[left], args.max_correspondence_distance, args.max_pairs, args.seed + right * 10 + left,
        )
        residuals = np.concatenate((lr, rl))
        all_residuals.append(residuals)
        pairs[f"{left}-{right}"] = {
            "symmetric_point_to_plane_m": distribution(residuals),
            "planar_reciprocal_count": int(lr_planar + rl_planar),
            "left_to_right_distance_overlap_ratio": lr_overlap,
            "right_to_left_distance_overlap_ratio": rl_overlap,
        }
    aggregate = np.concatenate(all_residuals) if all_residuals else np.empty(0)
    contract_errors = []
    if np.count_nonzero(~finite_xyz):
        contract_errors.append("nonfinite_xyz_points")
    if np.count_nonzero(~finite_id):
        contract_errors.append("nonfinite_lidar_id_points")
    if np.count_nonzero(invalid_id_mask):
        contract_errors.append("invalid_or_unexpected_lidar_id_points")
    if missing_ids:
        contract_errors.append("expected_lidar_ids_missing")
    cross_source_applicable = len(args.expected_ids) >= 2
    if cross_source_applicable and not len(aggregate):
        contract_errors.append("zero_valid_cross_source_correspondences")
    result = {
        "map": str(args.map), "map_sha256": sha256(args.map),
        "analyzer_sha256": sha256(Path(__file__).resolve()),
        "raw_point_count": int(original_point_count),
        "original_point_count": int(original_point_count),
        "nonfinite_point_count": int(np.count_nonzero(~(finite_xyz & finite_id))),
        "nonfinite_xyz_point_count": int(np.count_nonzero(~finite_xyz)),
        "nonfinite_lidar_id_point_count": int(np.count_nonzero(~finite_id)),
        "nonintegral_lidar_id_point_count": int(np.count_nonzero(finite_id & ~integral_id)),
        "invalid_lidar_id_point_count": int(np.count_nonzero(invalid_id_mask)),
        "invalid_or_unexpected_lidar_id_values": unexpected_ids,
        "finite_point_count": int(np.count_nonzero(finite_xyz & finite_id)),
        "accepted_point_count": int(len(xyz)),
        "original_points_by_integral_lidar_id": {
            str(lidar_id): int(np.count_nonzero(ids_integral == lidar_id))
            for lidar_id in sorted(observed_integral_ids)
        },
        "raw_points_by_lidar": {
            str(lidar_id): int(np.count_nonzero(ids_original[integral_id] == lidar_id))
            for lidar_id in args.expected_ids
        },
        "sampled_points_by_lidar": {str(key): len(value) for key, value in points_by_id.items()},
        "lidar_id_contract": {
            "expected_ids": list(args.expected_ids),
            "observed_integral_ids": sorted(observed_integral_ids),
            "missing_expected_ids": missing_ids,
            "unexpected_ids": unexpected_ids,
        },
        "parameters": {
            "voxel_resolution_m": args.voxel_resolution,
            "max_correspondence_distance_m": args.max_correspondence_distance,
            "max_points_per_lidar": args.max_points_per_lidar,
            "max_pairs_per_direction": args.max_pairs,
            "seed": args.seed,
        },
        "pair_metrics": pairs,
        "aggregate_symmetric_point_to_plane_m": distribution(aggregate),
        "cross_source_metric": {
            "applicable": cross_source_applicable,
            "status": "valid" if cross_source_applicable and len(aggregate) else (
                "not_applicable" if not cross_source_applicable else "invalid"
            ),
            "reason": None if cross_source_applicable else "fewer_than_two_expected_lidar_sources",
        },
        "coverage": coverage(points_by_id, 1.0),
        "contract": {"passed": not contract_errors, "errors": contract_errors},
        "interpretation_limit": (
            "This is a robust map-only cross-source consistency/coverage QC metric, not absolute map accuracy. "
            "Reciprocal locally planar correspondences reduce but do not eliminate dynamic-object and self-inclusion effects."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "contract": result["contract"],
        "aggregate": result["aggregate_symmetric_point_to_plane_m"],
        "output": str(args.output),
    }))
    return 0 if result["contract"]["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
