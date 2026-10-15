#!/usr/bin/env python3
"""Regression tests for output-point translation and frozen sweep contracts."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

import accuracy_reference as ar

ae = ar.ae
FIELDS = ["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw", "reference_valid", "reason", "epoch", "time_basis", "position_sigma", "raw_quality"]


def rotation(roll, pitch, yaw):
    a, b, c = np.radians([roll, pitch, yaw])
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()

    def document(self, name, value):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, ensure_ascii=False), encoding="utf-8")
        return path

    def source(self, mode="fixed_point", invalid=False, multi_epoch=False):
        q = ae.matrix_quat(rotation(27, -18, 76))
        row = dict(zip(FIELDS, [1_000_000_123, 4, -2, 7, *q, 1, "", "e0", "receive_proxy_ns", ".0123", "FIXED,heading_valid"]))
        rows = [row]
        if invalid:
            rows += [{**row, "timestamp_ns": 1_500_000_000, "reference_valid": 0, "reason": "RTK_FLOAT", "x": "nan", "qw": ""}]
        intervals = [dict(start_ns=0, end_ns=2_000_000_000, epoch="e0")]
        if multi_epoch:
            rows += [{**row, "timestamp_ns": 123, "epoch": "restart"}]
            intervals += [dict(start_ns=0, end_ns=2_000_000_000, epoch="restart")]
        path = self.root / "source.csv"
        self.csv(path, rows)
        meta = dict(schema_version=1, source_recording_id="synthetic_C", source_recording_sha256="recording-hash",
                    world_frame_id="enu:fixed_origin", body_frame_id="cgi_output_flu" if mode == "fixed_point" else "rear_axle_flu",
                    time_basis="receive_proxy_ns", recorded_intervals=intervals,
                    conditions=["SOURCE_MEASUREMENT_EPOCH_UNVERIFIED"], evidence={"georef": "survey", "clock": "receive proxy"})
        if mode == "fixed_point":
            meta["reference_point"] = dict(position="fixed_point", point_id="cgi_origin", attitude_axes="FLU",
                                           device_lever_arm_already_applied=False, confirmed=True,
                                           rigid_to_rear_axle_confirmed=True, evidence="synthetic rigid output contract")
        config = dict(schema_version=1, mode=mode, source_point_id="cgi_origin", confirmed=True,
                      evidence="synthetic known rear-body axes and original output point", rigid_to_rear_axle_confirmed=True,
                      lever_arm_m=[1.2, -.7, .4])
        if mode == "rear_axle_delta":
            config.update(previous_device_lever_arm_m=[.5, .1, -.3],
                          previous_compensation_is_rigid_translation_confirmed=True,
                          previous_compensation_evidence="synthetic old settings: only rigid output translation")
        if mode == "already_rear":
            config["lever_arm_m"] = [0, 0, 0]
        return path, self.document("source.metadata.json", meta), self.document("normalization.yaml", config), rows

    def csv(self, path, rows):
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, FIELDS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    def read_rows(self, path):
        with path.open(newline="", encoding="utf-8") as stream:
            return list(csv.DictReader(stream))

    def test_nontrivial_body_rotation_and_covariance_evidence(self):
        source, meta, config, rows = self.source()
        result = ar.normalize(source, meta, config, self.root / "rear.csv")
        actual = self.read_rows(self.root / "rear.csv")[0]
        expected = np.array([4, -2, 7]) + rotation(27, -18, 76) @ np.array([1.2, -.7, .4])
        np.testing.assert_allclose([float(actual[k]) for k in ("x", "y", "z")], expected, atol=1e-12)
        original = self.read_rows(source)[0]
        self.assertEqual({k: actual[k] for k in FIELDS if k not in ("x", "y", "z")},
                         {k: original[k] for k in FIELDS if k not in ("x", "y", "z")})
        self.assertIn(ar.COVARIANCE, result["conditions"])
        self.assertEqual(result["reference_normalization"]["target_covariance_status"], "not_propagated")
        self.assertEqual(result["reference_normalization"]["quality_fields_reference_point"]["point_id"], "cgi_origin")
        self.assertEqual(result["reference_normalization"]["output_csv_sha256"], ae.sha256(self.root / "rear.csv"))
        self.assertEqual(len(ae.read_reference(self.root / "rear.csv", ae.read_metadata(self.root / "rear.metadata.json", True))), 1)

    def test_old_rear_output_delta_equals_fixed_point_translation(self):
        source, meta, config, _ = self.source()
        ar.normalize(source, meta, config, self.root / "fixed_rear.csv")
        original = self.read_rows(source)
        old = np.array([.5, .1, -.3])
        old_position = np.array([float(original[0][k]) for k in ("x", "y", "z")]) + rotation(27, -18, 76) @ old
        original[0].update(zip(("x", "y", "z"), old_position))
        self.csv(self.root / "old_output.csv", original)
        old_meta = ae.read_document(meta)
        old_meta["body_frame_id"] = "rear_axle_flu"
        del old_meta["reference_point"]  # legacy product compatibility
        old_config = ae.read_document(config)
        old_config.update(mode="rear_axle_delta", previous_device_lever_arm_m=old.tolist(),
                          previous_compensation_is_rigid_translation_confirmed=True,
                          previous_compensation_evidence="known device output translation with matching epoch/attitude")
        ar.normalize(self.root / "old_output.csv", self.document("old.metadata.json", old_meta),
                     self.document("delta.yaml", old_config), self.root / "delta_rear.csv")
        a, b = self.read_rows(self.root / "fixed_rear.csv")[0], self.read_rows(self.root / "delta_rear.csv")[0]
        np.testing.assert_allclose([float(a[k]) for k in ("x", "y", "z")], [float(b[k]) for k in ("x", "y", "z")], atol=1e-12)
        self.assertEqual([a[k] for k in ("qx", "qy", "qz", "qw")], [b[k] for k in ("qx", "qy", "qz", "qw")])

    def test_legacy_identity_and_exclusive_output(self):
        source, meta, config, _ = self.source("already_rear")
        result = ar.normalize(source, meta, config, self.root / "rear.csv")
        self.assertEqual(self.read_rows(source), self.read_rows(self.root / "rear.csv"))
        self.assertNotIn(ar.COVARIANCE, result["conditions"])
        with self.assertRaisesRegex(ValueError, "already exists"):
            ar.normalize(source, meta, config, self.root / "rear.csv")
        self.document("blocked.metadata.json", {"do_not_touch": True})
        with self.assertRaisesRegex(ValueError, "already exists"):
            ar.normalize(source, meta, config, self.root / "blocked.csv")
        self.assertFalse((self.root / "blocked.csv").exists())
        with self.assertRaises(ValueError):
            ar.normalize(source, meta, config, source)

    def test_invalid_rows_epochs_and_provenance_preserved(self):
        source, meta, config, _ = self.source(invalid=True, multi_epoch=True)
        before = ae.read_document(meta)
        result = ar.normalize(source, meta, config, self.root / "rear.csv")
        rows = self.read_rows(self.root / "rear.csv")
        self.assertEqual(rows[1], self.read_rows(source)[1])
        for field in ("recorded_intervals", "source_recording_id", "source_recording_sha256", "world_frame_id", "time_basis", "evidence"):
            self.assertEqual(result[field], before[field])
        self.assertIn(before["conditions"][0], result["conditions"])
        chain = result["reference_normalization"]["source"]
        self.assertEqual(chain["sha256"], ae.sha256(source))
        self.assertEqual(chain["metadata_sha256"], ae.sha256(meta))
        self.assertEqual(chain["metadata"], before)
        self.assertEqual([p.epoch for p in ae.read_reference(self.root / "rear.csv", result)], ["e0", "e0", "restart"])

    def test_current_exporter_rear_axle_center_metadata_compatible(self):
        source, metadata, config, _ = self.source("already_rear")
        meta = ae.read_document(metadata)
        meta["reference_point"] = dict(position="rear_axle_center", attitude_axes="FLU", confirmed=True,
                                       device_lever_arm_already_applied=True, evidence="current exporter rear-axle contract")
        self.document("source.metadata.json", meta)
        result = ar.normalize(source, metadata, config, self.root / "current_rear.csv")
        self.assertEqual(result["reference_point"]["position"], "rear_axle_center")
        self.assertEqual(self.read_rows(source), self.read_rows(self.root / "current_rear.csv"))
        ae.read_reference(self.root / "current_rear.csv", ae.read_metadata(self.root / "current_rear.metadata.json", True))

    def test_wrong_point_axes_frame_and_double_compensation_rejected(self):
        source, meta, config, _ = self.source()
        original = ae.read_document(meta)
        variants = [dict(body_frame_id="rear_axle_flu")]
        for field, value in (("point_id", "other"), ("attitude_axes", "FRD"),
                             ("device_lever_arm_already_applied", True), ("rigid_to_rear_axle_confirmed", False)):
            variants.append({"reference_point": {**original["reference_point"], field: value}})
        for i, change in enumerate(variants):
            with self.subTest(change=change), self.assertRaises(ValueError):
                ar.normalize(source, self.document(f"bad{i}.json", {**original, **change}), config, self.root / f"bad{i}.csv")
        ar.normalize(source, meta, config, self.root / "rear.csv")
        already = ae.read_document(config)
        already.update(mode="already_rear", lever_arm_m=[0, 0, 0])
        with self.assertRaisesRegex(ValueError, "already normalized"):
            ar.normalize(self.root / "rear.csv", self.root / "rear.metadata.json", self.document("already.yaml", already), self.root / "twice.csv")

    def test_invalid_lever_unknown_old_compensation_and_articulation_rejected(self):
        source, meta, config, _ = self.source("rear_axle_delta")
        original = ae.read_document(config)
        changes = [{"lever_arm_m": [float("nan"), 0, 0]}, {"lever_arm_m": [0, 0, float("inf")]},
                   {"lever_arm_m": [True, 0, 0]}, {"previous_device_lever_arm_m": None},
                   {"previous_compensation_is_rigid_translation_confirmed": False},
                   {"previous_compensation_evidence": ""}, {"rigid_to_rear_axle_confirmed": False},
                   {"mode": "already_rear"}]
        for i, change in enumerate(changes):
            with self.subTest(change=change), self.assertRaises(ValueError):
                ar.normalize(source, meta, self.document(f"bad{i}.yaml", {**original, **change}), self.root / f"bad{i}.csv")

    def test_invalid_valid_pose_and_time_contract_rejected(self):
        source, meta, config, rows = self.source()
        for field, value in (("x", "nan"), ("qw", 2), ("time_basis", "another_clock"), ("timestamp_ns", 4_000_000_000)):
            self.csv(source, [{**rows[0], field: value}])
            with self.subTest(field=field), self.assertRaises(ValueError):
                ar.normalize(source, meta, config, self.root / f"bad_{field}.csv")

    def sweep_fixture(self, empty_estimate=False, all_invalid=False):
        # B and C are independent source recordings, with arbitrary ENU/map SE3
        # and a known nonzero clock correction. No test-set alignment is fitted.
        r_map, t_map = rotation(13, -8, 71), np.array([32., -18., 6.])
        lever, offset = np.array([1.2, -.7, .4]), 35_000_000
        (self.root / "data").mkdir()
        (self.root / "configs").mkdir()
        for recording, start in (("B", 100), ("C", 200)):
            reference, tum = [], []
            for i in range(9):
                timestamp = (start * ae.NS) + i * ae.NS // 2
                p = np.array([i, np.sin(i*.8)*3, .1*i*i])
                r = rotation(i*2, -i, 15*i)
                rear, q = r_map @ p + t_map, ae.matrix_quat(r_map @ r)
                source = rear if recording == "B" else rear - (r_map @ r) @ lever
                valid = recording == "B" or (not all_invalid and i != 3)
                reference.append(dict(zip(FIELDS, [timestamp, *source, *q, int(valid), "" if valid else "SYNTHETIC_BAD_REFERENCE",
                                                   "e0", "synthetic_ns", .01, "SYNTHETIC"])))
                est_time = timestamp - offset
                tum.append(f"{est_time//ae.NS}.{est_time%ae.NS:09d} " + " ".join(format(v, ".17g") for v in np.r_[p, ae.matrix_quat(r)]))
            self.csv(self.root / "data" / f"{recording}.csv", reference)
            (self.root / "data" / f"{recording}.tum").write_text("" if recording == "C" and empty_estimate else "\n".join(tum)+"\n", encoding="utf-8")
            interval = dict(start_ns=start*ae.NS, end_ns=(start+5)*ae.NS, epoch="e0")
            meta = dict(schema_version=1, source_recording_id=f"SYNTHETIC_{recording}", source_recording_sha256=f"source-{recording}",
                        world_frame_id="enu:site", body_frame_id="rear_axle_flu" if recording == "B" else "cgi_output_flu",
                        time_basis="synthetic_ns", recorded_intervals=[interval], conditions=["SYNTHETIC_NOT_FIELD_DATA"])
            if recording == "C":
                meta["reference_point"] = dict(position="fixed_point", point_id="cgi_origin", attitude_axes="FLU",
                                               device_lever_arm_already_applied=False, confirmed=True, rigid_to_rear_axle_confirmed=True,
                                               evidence="synthetic rear-body pose at fixed output origin")
            self.document(f"data/{recording}.reference.json", meta)
            emeta = {**meta, "body_frame_id": "rear_axle_flu", "world_frame_id": "map:asset", "trajectory_kind": "localization",
                     "clock_epoch": "e0", "recorded_intervals": ae.shifted([interval], -offset)}
            emeta.pop("reference_point", None)
            self.document(f"data/{recording}.estimate.json", emeta)
            cfg = dict(schema_version=1, role="calibration" if recording == "B" else "evaluation",
                       reference=dict(csv=f"../data/{recording}.csv", metadata=f"../data/{recording}.reference.json"),
                       estimates=[dict(name="localization", tum=f"../data/{recording}.tum", metadata=f"../data/{recording}.estimate.json", alignment="../frozen.json")],
                       analysis_intervals=[interval], time=dict(offset_status="known", estimate_to_reference_offset_ns=offset,
                                                               reference_epoch="e0", evidence="synthetic known offset"),
                       pairing=dict(max_reference_gap_sec=.6, availability_hold_sec=.6), conditions=["SYNTHETIC_CONFIG"],
                       thresholds={"position_3d_rmse_m": 100}, fit=dict(method="position_attitude", min_pairs=5,
                       min_secondary_span_m=.1, attitude_weight_m=1, max_attitude_dispersion_deg=1))
            self.document(f"configs/{recording}.yaml", cfg)
        ae.fit(ae.load_config(self.root / "configs/B.yaml"), "localization", self.root / "frozen.json")
        normal = dict(schema_version=1, mode="fixed_point", source_point_id="cgi_origin", confirmed=True,
                      evidence="synthetic fixed output to rear axle", rigid_to_rear_axle_confirmed=True, lever_arm_m=lever.tolist())
        self.document("configs/normal.yaml", normal)
        sweep = dict(schema_version=1, normalization_config="normal.yaml", evaluation_config="C.yaml",
                     candidates=[dict(name="measured", lever_arm_m=lever.tolist()), dict(name="perturbed", lever_arm_m=(lever + [.2, -.1, .3]).tolist())])
        self.document("configs/sweep.yaml", sweep)
        return self.root / "data/C.csv", self.root / "data/C.reference.json", self.root / "configs/sweep.yaml"

    def test_sweep_frozen_alignment_same_samples_attitude_and_no_acceptance(self):
        args = self.sweep_fixture()
        before = ae.sha256(self.root / "frozen.json")
        # A sweep must never call fit, including with calibration options present.
        with patch.object(ae, "fit", side_effect=AssertionError("sweep called fit")):
            result = ar.sweep(*args, self.root / "sweep")
        self.assertEqual(result["report_type"], "lever_arm_sensitivity")
        self.assertFalse(result["automatic_selection"])
        self.assertTrue(all(result["invariants"].values()))
        self.assertEqual(before, ae.sha256(self.root / "frozen.json"))
        measured, perturbed = [c["estimates"][0] for c in result["candidates"]]
        self.assertLess(measured["metrics"]["position_3d_m"]["rmse"], 1e-12)
        self.assertAlmostEqual(perturbed["metrics"]["position_3d_m"]["rmse"], np.linalg.norm([.2, -.1, .3]), places=10)
        for candidate in result["candidates"]:
            estimate = candidate["estimates"][0]
            self.assertEqual(estimate["acceptance"]["status"], "NOT_ASSESSED_NO_THRESHOLDS")
            self.assertIsNone(estimate["acceptance"]["pass"])
            self.assertIn(ar.SENSITIVITY, estimate["conditions"])
            self.assertIn(ar.COVARIANCE, estimate["conditions"])
            self.assertEqual(estimate["coverage"]["recorded_duration_sec"], 5)
            self.assertEqual(estimate["metrics"]["attitude_geodesic_deg"], measured["metrics"]["attitude_geodesic_deg"])
        self.assertEqual(measured["sample_counts"], perturbed["sample_counts"])
        self.assertEqual(measured["sample_counts"]["unpaired"], 1)
        for name in ("sweep.json", "sweep.md", "summary.csv"):
            self.assertTrue((self.root / "sweep" / name).is_file())
        written_config = ae.read_document(self.root / "sweep/measured/evaluation.yaml")
        self.assertIsNone(written_config["thresholds"])
        self.assertEqual(written_config["time"]["estimate_to_reference_offset_ns"], 35_000_000)
        self.assertEqual(Path(written_config["estimates"][0]["tum"]), self.root / "data/C.tum")

    def test_sweep_complete_estimate_outage_keeps_denominator(self):
        args = self.sweep_fixture(empty_estimate=True)
        result = ar.sweep(*args, self.root / "sweep")
        self.assertEqual(result["status"], "NO_PAIRED_DATA")
        for candidate in result["candidates"]:
            estimate = candidate["estimates"][0]
            self.assertEqual(estimate["acceptance"]["status"], "NO_PAIRED_DATA")
            self.assertEqual(estimate["coverage"]["recorded_duration_sec"], 5)
            self.assertEqual(estimate["coverage"]["estimate_availability_fraction"], 0)
            self.assertIsNone(estimate["metrics"]["position_3d_m"]["rmse"])

    def test_sweep_all_invalid_reference_keeps_failure_population(self):
        args = self.sweep_fixture(all_invalid=True)
        result = ar.sweep(*args, self.root / "sweep")
        self.assertEqual(result["status"], "NO_PAIRED_DATA")
        self.assertEqual(result["candidates"][0]["estimates"][0]["sample_counts"]["unpaired"], 9)
        self.assertEqual(result["candidates"][0]["estimates"][0]["coverage"]["reference_valid_fraction"], 0)

    def test_sweep_rejects_calibration_identity_bad_names_and_reused_output(self):
        args = self.sweep_fixture()
        original = ae.read_document(args[2])
        for i, name in enumerate(("../escape", "CON", "", "with space")):
            wrong = {**original, "candidates": [dict(name=name, lever_arm_m=[0, 0, 0])]}
            with self.subTest(name=name), self.assertRaises(ValueError):
                ar.sweep(args[0], args[1], self.document(f"configs/wrong{i}.yaml", wrong), self.root / f"bad{i}")
        with self.assertRaisesRegex(ValueError, "role: evaluation"):
            ar.sweep(args[0], args[1], self.document("configs/calibration_sweep.yaml", {**original, "evaluation_config": "B.yaml"}), self.root / "bad_role")
        (self.root / "existing").mkdir()
        with self.assertRaisesRegex(ValueError, "already exists"):
            ar.sweep(*args, self.root / "existing")
        normal = ae.read_document(self.root / "configs/normal.yaml")
        self.document("configs/normal.yaml", {**normal, "mode": "already_rear", "lever_arm_m": [0, 0, 0]})
        with self.assertRaisesRegex(ValueError, "not a lever sweep mode"):
            ar.sweep(*args, self.root / "identity")

    def test_sweep_rejects_changed_frozen_input(self):
        args = self.sweep_fixture()
        real_evaluate = ae.evaluate
        def mutate_after_evaluation(config, output):
            report = real_evaluate(config, output)
            with (self.root / "frozen.json").open("a", encoding="utf-8") as stream:
                stream.write(" ")
            return report
        with patch.object(ae, "evaluate", side_effect=mutate_after_evaluation), self.assertRaisesRegex(ValueError, "changed during processing"):
            ar.sweep(*args, self.root / "changed")
        self.assertFalse((self.root / "changed/sweep.json").exists())

    def test_sweep_rejects_changed_sample_status(self):
        args = self.sweep_fixture()
        real_signature, calls = ar.sample_signature, []
        def bad_signature(report, directory):
            value = real_signature(report, directory)
            if calls:
                value["localization"]["samples"][0][3] = "unexpected_different_pairing"
            calls.append(True)
            return value
        with patch.object(ar, "sample_signature", side_effect=bad_signature), self.assertRaisesRegex(ValueError, "changed pairing/status"):
            ar.sweep(*args, self.root / "changed")

    def test_normalize_cli(self):
        source, meta, config, _ = self.source()
        command = [sys.executable, "-B", str(Path(ar.__file__)), "normalize", "--source", str(source),
                   "--metadata", str(meta), "--config", str(config), "--output", str(self.root / "cli.csv")]
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(result.stdout)["status"], "NORMALIZED")


if __name__ == "__main__":
    unittest.main()
