"""Meaningful regression cases for frame, time, quality and denominator contracts."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import accuracy_eval as ae


def rotation(roll=0, pitch=0, yaw=0):
    a, b, c = np.radians([roll, pitch, yaw])
    rx = np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
    ry = np.array([[np.cos(b), 0, np.sin(b)], [0, 1, 0], [-np.sin(b), 0, np.cos(b)]])
    rz = np.array([[np.cos(c), -np.sin(c), 0], [np.sin(c), np.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


def pose(t, p=(0, 0, 0), r=None, valid=True, reason="", epoch="e0"):
    return ae.Pose(int(t*ae.NS), np.array(p, dtype=float), ae.matrix_quat(np.eye(3) if r is None else r), epoch, valid, reason)


class GeometryTests(unittest.TestCase):
    def test_quaternion_roundtrip_and_slerp_near_pi(self):
        for r in (rotation(35, -25, 179.99), rotation(180, 0, 0), rotation(0, 90, 0)):
            np.testing.assert_allclose(ae.quat_matrix(ae.matrix_quat(r)), r, atol=1e-12)
        q = ae.slerp(ae.matrix_quat(rotation(yaw=179)), ae.matrix_quat(rotation(yaw=-179)), .5)
        self.assertAlmostEqual(abs(ae.euler_xyz(ae.quat_matrix(q))[2]), 180)

    def test_wrapped_euler_and_so3(self):
        error = ae.disagreement(pose(0, r=rotation(yaw=-179)), pose(0, r=rotation(yaw=179)), np.eye(3), np.zeros(3))
        self.assertAlmostEqual(error["yaw_deg"], 2)
        self.assertAlmostEqual(error["attitude_geodesic_deg"], 2)
        coupled = ae.disagreement(pose(0, r=rotation(40, 30, 20)), pose(0), np.eye(3), np.zeros(3))
        self.assertNotAlmostEqual(coupled["attitude_geodesic_deg"], np.linalg.norm([40, 30, 20]))

    def test_body_error_is_reference_flu(self):
        error = ae.disagreement(pose(0, (2, 3, 4)), pose(0, r=rotation(yaw=90)), np.eye(3), np.zeros(3))
        self.assertAlmostEqual(error["longitudinal_m"], 3)
        self.assertAlmostEqual(error["lateral_m"], -2)
        self.assertAlmostEqual(error["vertical_m"], 4)

    def test_exported_lever_arm_is_not_applied_twice(self):
        # Exporter contract: T_WB = T_WI * T_IL * inverse(T_BL).
        r_wi, p_wi = rotation(8, -4, 60), np.array([10., 20., 3.])
        r_il, p_il = rotation(2, 3, -8), np.array([.4, -.2, .7])
        r_bl, p_bl = rotation(-1, 5, 12), np.array([1.7, .3, 2.1])
        r_wb = r_wi @ r_il @ r_bl.T
        p_wb = p_wi + r_wi @ p_il - r_wb @ p_bl
        canonical = pose(0, p_wb, r_wb)
        error = ae.disagreement(canonical, canonical, np.eye(3), np.zeros(3))
        self.assertEqual(error["position_3d_m"], 0)
        self.assertAlmostEqual(error["attitude_geodesic_deg"], 0)
        self.assertGreater(np.linalg.norm(p_wb-p_wi), 1)

    def test_invalid_reference_gap_is_never_bridged(self):
        rows = [pose(0), pose(1, valid=False, reason="RTK_FLOAT"), pose(2, (2, 0, 0)), pose(3, (3, 0, 0))]
        reference = ae.Reference(rows, [{"start_ns": 0, "end_ns": 4*ae.NS, "epoch": "e0"}], 3*ae.NS)
        self.assertIsNone(reference.at("e0", int(.5*ae.NS))[0])
        self.assertIsNone(reference.at("e0", int(1.5*ae.NS))[0])
        self.assertIsNotNone(reference.at("e0", int(2.5*ae.NS))[0])
        self.assertEqual(reference.valid_intervals(), [{"start_ns": 2*ae.NS, "end_ns": 3*ae.NS, "epoch": "e0"}])

    def test_reference_max_gap_epoch_and_recording_interruptions(self):
        reference = ae.Reference([pose(0), pose(2)], [{"start_ns": 0, "end_ns": 3*ae.NS, "epoch": "e0"}], ae.NS)
        self.assertEqual(reference.at("e0", ae.NS)[1], "reference_gap_exceeded")
        self.assertIsNone(reference.at("other", 0)[0])
        intervals = [{"start_ns": 0, "end_ns": ae.NS, "epoch": "e0"}, {"start_ns": 2*ae.NS, "end_ns": 3*ae.NS, "epoch": "e0"}]
        reference = ae.Reference([pose(.5), pose(2.5)], intervals, 3*ae.NS)
        self.assertEqual(reference.at("e0", ae.NS)[1], "reference_recording_interruption")

    def test_invalid_duplicate_timestamp_rejects_exact_match(self):
        reference = ae.Reference([pose(0), pose(1), pose(1, valid=False, reason="duplicate"), pose(2)],
                                 [{"start_ns": 0, "end_ns": 3*ae.NS, "epoch": "e0"}], 2*ae.NS)
        self.assertIsNone(reference.at("e0", ae.NS)[0])
        self.assertIsNone(reference.at("e0", ae.NS//2)[0])
        self.assertIsNone(reference.at("e0", ae.NS*3//2)[0])
        self.assertEqual(reference.valid_intervals(), [])

    def test_outages_denominator_is_full_recorded_duration(self):
        analysis = [{"start_ns": 0, "end_ns": 10*ae.NS, "epoch": "e0"}]
        available = ae.estimate_availability([pose(2), pose(3), pose(7)], analysis, ae.NS, 0, "e0")
        reference_valid = [{"start_ns": ae.NS, "end_ns": 9*ae.NS, "epoch": "e0"}]
        result, segments = ae.coverage(analysis, reference_valid, available)
        self.assertEqual(result["recorded_duration_sec"], 10)
        self.assertEqual(result["reference_valid_fraction"], .8)
        self.assertEqual(result["estimate_availability_fraction"], .3)
        self.assertEqual(result["estimate_unavailable_duration_sec"], 7)
        self.assertEqual(result["estimate_availability_given_valid_reference"], 3/8)
        self.assertEqual(sum(s["duration_sec"] for s in segments), 10)

    def test_straight_calibration_path_is_rejected_even_with_attitude(self):
        pairs = [(pose(t, (t, 0, 0)), pose(t, (t+3, 0, 0))) for t in range(10)]
        with self.assertRaisesRegex(ValueError, "near-collinear"):
            ae.fit_rigid(pairs, {"min_pairs": 3, "min_secondary_span_m": .1,
                                 "method": "position_attitude", "attitude_weight_m": 1, "max_attitude_dispersion_deg": 10})


class FileAndIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.r = rotation(7, -11, 123)
        self.translation = np.array([125., -60., 8.])

    def write_json(self, name, value):
        path = self.root/name
        path.write_text(json.dumps(value), encoding="utf-8")
        return path

    def dataset(self, name, start, offset_ns=0, bias=(0, 0, 0), recording_id=None, empty=False):
        points, reference_rows, estimate_rows = [], [], []
        for i in range(21):
            phase = i/20
            timestamp = int((start + i*.5)*ae.NS)
            p = np.array([12*phase, 4*np.sin(phase*4), phase**2])
            attitude = rotation(3*phase, -5*phase, 30*phase)
            ref_p = self.r @ p + self.translation
            ref_q = ae.matrix_quat(self.r @ attitude)
            reference_rows.append([timestamp, *ref_p, *ref_q, 1, "", "e0", "unix_ns"])
            est_p = p + self.r.T @ np.asarray(bias)
            est_q = ae.matrix_quat(attitude)
            estimate_rows.append(" ".join([format((timestamp-offset_ns)/ae.NS, ".9f"), *[format(v, ".16g") for v in np.r_[est_p, est_q]]]))
        csv_path = self.root/(name+".csv")
        with csv_path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw", "reference_valid", "reason", "epoch", "time_basis"])
            writer.writerows(reference_rows)
        tum = self.root/(name+".tum")
        tum.write_text("" if empty else "\n".join(estimate_rows)+"\n", encoding="utf-8")
        interval = {"start_ns": int(start*ae.NS), "end_ns": int((start+10.5)*ae.NS), "epoch": "e0"}
        meta = {"schema_version": 1, "source_recording_id": recording_id or name, "world_frame_id": "enu:site1",
                "body_frame_id": ae.BODY_FRAME, "time_basis": "unix_ns", "recorded_intervals": [interval]}
        ref_meta = self.write_json(name+".reference.json", meta)
        est_meta = {**meta, "world_frame_id": "map:asset1", "trajectory_kind": "localization", "clock_epoch": "e0",
                    "recorded_intervals": ae.shifted([interval], -offset_ns)}
        estimate_meta = self.write_json(name+".estimate.json", est_meta)
        return {"csv": csv_path.name, "metadata": ref_meta.name}, {"name": name, "tum": tum.name, "metadata": estimate_meta.name}, interval

    def config(self, name, role, reference, estimate, interval, offset_ns=0, **extra):
        raw = {"schema_version": 1, "role": role, "reference": reference, "estimates": [estimate], "analysis_intervals": [interval],
               "time": {"offset_status": "known", "estimate_to_reference_offset_ns": offset_ns, "evidence": "Synthetic known clock relation", "reference_epoch": "e0"},
               "pairing": {"max_reference_gap_sec": .6, "availability_hold_sec": .6},
               "fit": {"method": "position_attitude", "min_pairs": 10, "min_secondary_span_m": .1, "attitude_weight_m": 2, "max_attitude_dispersion_deg": 1}, **extra}
        path = self.root/(name+".yaml")
        path.write_text(yaml.safe_dump(raw), encoding="utf-8")
        return ae.load_config(path)

    def calibrate(self, offset_ns=0):
        ref, est, interval = self.dataset("calibration", 10, offset_ns)
        config = self.config("calibration", "calibration", ref, est, interval, offset_ns)
        path = self.root/"frozen.json"
        result = ae.fit(config, "calibration", path)
        return path, result, config

    def evaluation(self, frozen, offset_ns=0, **data_options):
        ref, est, interval = self.dataset("evaluation", 30, offset_ns, **data_options)
        est["alignment"] = frozen.name
        return self.config("evaluation", "evaluation", ref, est, interval, offset_ns)

    def test_arbitrary_map_fixed_alignment_and_full_report(self):
        path, fitted, _ = self.calibrate()
        np.testing.assert_allclose(ae.quat_matrix(np.array(fitted["transform"]["quaternion_xyzw"])), self.r, atol=1e-12)
        np.testing.assert_allclose(fitted["transform"]["translation_m"], self.translation, atol=1e-11)
        config = self.evaluation(path, bias=(.1, -.2, .3))
        report = ae.evaluate(config, self.root/"report")
        result = report["estimates"][0]
        self.assertAlmostEqual(result["metrics"]["position_3d_m"]["rmse"], np.linalg.norm([.1, -.2, .3]))
        self.assertAlmostEqual(result["metrics"]["x_m"]["bias"], .1)
        self.assertLess(result["metrics"]["attitude_geodesic_deg"]["max"], 1e-10)
        self.assertEqual(result["acceptance"]["status"], "NOT_ASSESSED_NO_THRESHOLDS")
        self.assertIsNone(result["acceptance"]["pass"])
        self.assertTrue((self.root/"report"/"report.md").is_file())
        json.loads((self.root/"report"/"report.json").read_text(encoding="utf-8"))

    def test_nonzero_time_offset_is_applied_with_declared_sign(self):
        offset = 250_000_000
        path, _, _ = self.calibrate(offset)
        config = self.evaluation(path, offset)
        good = ae.evaluate(config, self.root/"good")["estimates"][0]
        self.assertLess(good["metrics"]["position_3d_m"]["rmse"], 1e-10)
        config["time"]["estimate_to_reference_offset_ns"] = 0
        wrong = ae.evaluate(config, self.root/"wrong")["estimates"][0]
        self.assertGreater(wrong["metrics"]["position_3d_m"]["rmse"], .2)

    def test_identical_calibration_file_is_rejected(self):
        path, _, config = self.calibrate()
        config["role"] = "evaluation"
        config["estimates"][0]["alignment"] = path.name
        with self.assertRaisesRegex(ValueError, "file hashes overlap"):
            ae.evaluate(config, self.root/"bad")
        self.assertFalse((self.root/"bad").exists())

    def test_canonical_hash_rejects_reformatted_reused_data(self):
        path, _, config = self.calibrate()
        original = self.root/config["reference"]["csv"]
        copied = self.root/"reformatted.csv"
        copied.write_bytes(original.read_bytes().replace(b"\r\n", b"\n"))
        self.assertNotEqual(ae.sha256(copied), ae.sha256(original))
        config["reference"]["csv"] = copied.name
        config["role"] = "evaluation"
        config["estimates"][0]["alignment"] = path.name
        with self.assertRaisesRegex(ValueError, "canonical pose hashes overlap"):
            ae.evaluate(config, self.root/"bad")

    def test_same_recording_overlapping_ranges_are_rejected(self):
        path, _, _ = self.calibrate()
        ref, est, interval = self.dataset("changed_copy", 12, recording_id="calibration", bias=(1, 0, 0))
        est["alignment"] = path.name
        config = self.config("bad", "evaluation", ref, est, interval)
        with self.assertRaisesRegex(ValueError, "recording intervals overlap"):
            ae.evaluate(config, self.root/"bad")

    def test_disjoint_reserved_segment_same_run_is_allowed(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path, recording_id="calibration")
        self.assertEqual(ae.evaluate(config, self.root/"disjoint")["estimates"][0]["sample_counts"]["paired"], 21)

    def test_different_map_frame_is_rejected(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path)
        meta_path = self.root/config["estimates"][0]["metadata"]
        meta = json.loads(meta_path.read_text())
        meta["world_frame_id"] = "map:other"
        meta_path.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, "another map/run frame"):
            ae.evaluate(config, self.root/"bad")

    def test_canonical_body_frame_required(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path)
        meta_path = self.root/config["estimates"][0]["metadata"]
        meta = json.loads(meta_path.read_text())
        meta["body_frame_id"] = "imu"
        meta_path.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, "already be rear_axle_flu"):
            ae.evaluate(config, self.root/"bad")

    def test_tum_preserves_nanoseconds_at_unix_epoch_magnitude(self):
        timestamp = 1_800_000_000_123_456_789
        tum = self.root/"precise.tum"
        tum.write_text("1800000000.123456789 0 0 0 0 0 0 1\n")
        meta = {"clock_epoch": "e0", "recorded_intervals": [{"start_ns": timestamp-1, "end_ns": timestamp+1, "epoch": "e0"}]}
        self.assertEqual(ae.read_tum(tum, meta)[0].timestamp, timestamp)

    def test_complete_estimate_failure_is_reported(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path, empty=True)
        result = ae.evaluate(config, self.root/"empty")["estimates"][0]
        self.assertEqual(result["coverage"]["recorded_duration_sec"], 10.5)
        self.assertEqual(result["coverage"]["estimate_availability_fraction"], 0)
        self.assertEqual(result["coverage"]["estimate_unavailable_duration_sec"], 10.5)
        self.assertEqual(result["acceptance"]["status"], "NO_PAIRED_DATA")
        self.assertIsNone(result["metrics"]["position_3d_m"]["rmse"])

    def test_unknown_offset_fails_closed_and_assumed_is_conditional(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path)
        config["time"]["offset_status"] = "unknown"
        with self.assertRaisesRegex(ValueError, "unknown time offset"):
            ae.evaluate(config, self.root/"unknown")
        config["time"].update(offset_status="assumed", allow_conditional=True)
        config["thresholds"] = {"position_3d_rmse_m": 1}
        result = ae.evaluate(config, self.root/"conditional")["estimates"][0]
        self.assertEqual(result["acceptance"]["status"], "CONDITIONAL")
        self.assertIsNone(result["acceptance"]["pass"])

    def test_reference_and_frozen_conditions_survive_known_clock_override(self):
        path, fitted, _ = self.calibrate()
        fitted["conditions"] = ["calibration_field_measurement_epoch_unverified"]
        path.write_text(json.dumps(fitted), encoding="utf-8")
        config = self.evaluation(path)
        meta_path = self.root/config["reference"]["metadata"]
        meta = json.loads(meta_path.read_text())
        meta["conditions"] = ["reference_timestamp_is_socketcan_receive_proxy"]
        meta_path.write_text(json.dumps(meta))
        config["thresholds"] = {"position_3d_rmse_m": 1}
        result = ae.evaluate(config, self.root/"conditional")["estimates"][0]
        self.assertEqual(result["acceptance"]["status"], "CONDITIONAL")
        self.assertEqual(set(result["conditions"]), {"calibration_field_measurement_epoch_unverified", "reference_timestamp_is_socketcan_receive_proxy"})

    def test_receive_time_basis_automatically_adds_condition(self):
        ref, _, _ = self.dataset("proxy", 50)
        meta_path = self.root/ref["metadata"]
        meta = json.loads(meta_path.read_text())
        meta["time_basis"] = "socketcan_receive_proxy"
        meta_path.write_text(json.dumps(meta))
        parsed = ae.read_metadata(meta_path, reference=True)
        self.assertTrue(any("measurement epoch" in c for c in parsed["conditions"]))

    def test_freeze_is_immutable_and_evaluation_cannot_fit(self):
        path, _, config = self.calibrate()
        with self.assertRaises(FileExistsError):
            ae.fit(config, "calibration", path)
        config["role"] = "evaluation"
        with self.assertRaisesRegex(ValueError, "role: calibration"):
            ae.fit(config, "calibration", self.root/"illegal.json")

    def test_same_map_must_use_same_frozen_artifact(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path)
        changed = json.loads(path.read_text())
        changed["comment"] = "different artifact even if numerically same"
        second = self.write_json("second.json", changed)
        config["estimates"].append({**config["estimates"][0], "name": "second", "alignment": second.name})
        with self.assertRaisesRegex(ValueError, "identical frozen alignment"):
            ae.evaluate(config, self.root/"bad")

    def test_before_after_uses_exact_common_population_and_keeps_outages(self):
        path, _, _ = self.calibrate()
        config = self.evaluation(path)
        original = self.root/config["estimates"][0]["tum"]
        shorter = self.root/"shorter.tum"
        # A missing middle output must not change the other stream's group population.
        rows = original.read_text().splitlines()
        shorter.write_text("\n".join(rows[:10]+rows[11:])+"\n")
        config["estimates"].append({**config["estimates"][0], "name": "optimized", "tum": shorter.name})
        config["common_sample_groups"] = [{"name": "before_after", "estimates": ["evaluation", "optimized"]}]
        report = ae.evaluate(config, self.root/"common")
        group = report["common_sample_groups"][0]
        self.assertEqual(group["paired_common_count"], 20)
        self.assertEqual(group["paired_rows_excluded_from_common"], {"evaluation": 1, "optimized": 0})
        self.assertGreater(report["estimates"][1]["coverage"]["estimate_unavailable_duration_sec"], 0)
        self.assertTrue((self.root/"common"/group["sample_csv"]).is_file())


if __name__ == "__main__":
    unittest.main()
