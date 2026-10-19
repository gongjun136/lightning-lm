"""A-only diagnostics must retain their fit bias, outages and non-acceptance boundary."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import accuracy_diagnostic as diagnostic
import accuracy_eval as ae


class MappingDiagnosticTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.interval = {"start_ns": 0, "end_ns": 12*ae.NS, "epoch": "Aclock"}
        base = {"schema_version": 1, "source_recording_id": "mapping_A", "body_frame_id": "rear_axle_flu",
                "time_basis": "socketcan_receive_proxy", "recorded_intervals": [self.interval],
                "conditions": ["SYNTHETIC TEST DATA", "HEIGHT_DATUM_UNCONFIRMED: altitude used conditionally."]}
        (self.root/"reference.json").write_text(json.dumps({**base, "world_frame_id": "enu:diagnostic-origin"}), encoding="utf-8")
        angle = np.radians(40)
        self.rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        self.translation = np.array([10., 20., 5.])
        reference, front, optimized = [], [], []
        for i in range(21):
            t = i*.5
            p = np.array([10*np.cos(i*.15), 10*np.sin(i*.15), i*.02])
            q = [0., 0., 0., 1.]
            valid = i >= 4 and i != 12
            reference.append([int(t*ae.NS), *(self.rotation @ p+self.translation), *ae.matrix_quat(self.rotation),
                              int(valid), "" if valid else "RTK_INVALID", "Aclock", "socketcan_receive_proxy"])
            optimized.append(" ".join([f"{t:.9f}", *map(str, p), *map(str, q)]))
            front.append(" ".join([f"{t:.9f}", *map(str, p+np.array([1., 0, 0])), *map(str, q)]))
        with (self.root/"reference.csv").open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["timestamp_ns", "x", "y", "z", "qx", "qy", "qz", "qw", "reference_valid", "reason", "epoch", "time_basis"])
            writer.writerows(reference)
        for name, lines, kind in (("frontend_kf", front, "frontend"), ("optimized_kf", optimized, "optimized_keyframes")):
            (self.root/(name+".tum")).write_text("\n".join(lines)+"\n", encoding="utf-8")
            (self.root/(name+".json")).write_text(json.dumps({**base, "world_frame_id": "map:Aasset",
                "trajectory_kind": kind, "clock_epoch": "Aclock"}), encoding="utf-8")
        self.raw = {"schema_version": 1, "role": "diagnostic_mapping_a",
            "reference": {"csv": "reference.csv", "metadata": "reference.json"},
            "estimates": [{"name": n, "tum": n+".tum", "metadata": n+".json"} for n in ("frontend_kf", "optimized_kf")],
            "analysis_intervals": [self.interval], "fit_estimate": "optimized_kf", "fit_common_sample_group": "keyframes",
            "common_sample_groups": [{"name": "keyframes", "estimates": ["frontend_kf", "optimized_kf"]}],
            "time": {"offset_status": "assumed", "allow_conditional": True, "estimate_to_reference_offset_ns": 0,
                     "reference_epoch": "Aclock", "evidence": "SYNTHETIC declared zero receive-proxy offset, not measurement time."},
            "pairing": {"max_reference_gap_sec": .6, "availability_hold_sec": .6},
            "fit": {"method": "position", "min_pairs": 3, "min_secondary_span_m": .1}, "thresholds": None,
            "conditions": ["SYNTHETIC TEST DATA"]}

    def config(self):
        path = self.root/"diagnostic.yaml"
        path.write_text(yaml.safe_dump(self.raw), encoding="utf-8")
        return diagnostic.load_config(path)

    def test_shared_fit_does_not_hide_frontend_bias_and_keeps_full_outages(self):
        report = diagnostic.diagnose(self.config(), self.root/"report")
        self.assertEqual(report["status"], diagnostic.STATUS)
        self.assertFalse(report["independent_accuracy_assessed"])
        self.assertFalse(report["localization_accuracy_assessed"])
        front, optimized = report["estimates"]
        self.assertAlmostEqual(front["metrics"]["position_3d_m"]["rmse"], 1)
        self.assertLess(optimized["metrics"]["position_3d_m"]["rmse"], 1e-12)
        self.assertEqual(front["diagnostic_alignment_sha256"], optimized["diagnostic_alignment_sha256"])
        self.assertEqual(report["common_sample_groups"][0]["paired_common_count"], 16)
        alignment = report["diagnostic_alignment"]["content"]
        self.assertEqual(alignment["fit_pair_count"], 16)
        self.assertEqual(alignment["transform"]["scale"], 1)
        np.testing.assert_allclose(alignment["transform"]["translation_m"], self.translation, atol=1e-12)
        for estimate in report["estimates"]:
            self.assertIsNone(estimate["acceptance"]["pass"])
            self.assertEqual(estimate["coverage"]["recorded_duration_sec"], 12)
            self.assertEqual(estimate["coverage"]["reference_valid_duration_sec"], 7)
            self.assertAlmostEqual(estimate["coverage"]["estimate_unavailable_duration_sec"], 1.4)
            self.assertEqual(estimate["sample_counts"]["unpaired"], 5)
        self.assertTrue(any("HEIGHT_DATUM_UNCONFIRMED" in c for c in report["conditions"]))
        self.assertTrue(any("measurement epoch" in c for c in report["conditions"]))
        self.assertEqual(ae.sha256(self.root/"report"/"diagnostic_alignment.json"), report["diagnostic_alignment"]["sha256"])
        self.assertTrue((self.root/"report"/"alignment_fit_samples.csv").is_file())

    def test_diagnostic_artifact_cannot_become_formal_frozen_alignment(self):
        config = self.config()
        diagnostic.diagnose(config, self.root/"report")
        _, ref, estimates = ae.load_inputs(config)
        with self.assertRaisesRegex(ValueError, "not a supported frozen alignment"):
            ae.load_alignment(self.root/"report"/"diagnostic_alignment.json", config, ref, estimates[0])

    def test_fit_and_primary_comparison_use_identical_common_mask(self):
        path = self.root/"frontend_kf.tum"
        rows = path.read_text().splitlines()
        path.write_text("\n".join(row for i, row in enumerate(rows) if i != 10)+"\n")
        report = diagnostic.diagnose(self.config(), self.root/"report")
        self.assertEqual(report["diagnostic_alignment"]["content"]["fit_pair_count"], 15)
        self.assertEqual(report["common_sample_groups"][0]["paired_common_count"], 15)
        self.assertEqual(report["estimates"][1]["sample_counts"]["paired"], 16)

    def test_no_valid_reference_does_not_become_zero_error(self):
        path = self.root/"reference.csv"
        with path.open(newline="") as stream:
            reader = csv.DictReader(stream)
            fields, rows = reader.fieldnames, list(reader)
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({**row, "reference_valid": 0, "reason": "NO_RTK"})
        with self.assertRaisesRegex(ValueError, "at least 3 matched poses"):
            diagnostic.diagnose(self.config(), self.root/"bad")
        self.assertFalse((self.root/"bad").exists())

    def test_rpe_is_gauge_invariant_and_never_crosses_invalid_interior(self):
        self.raw["relative_motion"] = {"deltas_sec": [1, 10], "max_estimate_gap_sec": .6}
        report = diagnostic.diagnose(self.config(), self.root/"report")
        for estimate in report["estimates"]:
            short, long = estimate["relative_motion"]
            self.assertEqual(short["start_samples_in_analysis"], 21)
            self.assertEqual(short["eligible_start_count"], 21)
            self.assertEqual(short["paired_count"], 12)
            self.assertLess(short["metrics"]["translation_3d_m"]["max"], 1e-12)
            self.assertLess(short["metrics"]["rotation_geodesic_deg"]["max"], 1e-12)
            self.assertFalse(short["global_alignment_applied"])
            self.assertIsNone(short["pass"])
            self.assertGreater(short["rejected_reasons"]["reference_segment_discontinuous"], 0)
            self.assertEqual(long["eligible_start_count"], 4)
            self.assertEqual(long["paired_count"], 0)
            self.assertEqual(long["status"], "NO_PAIRED_DATA")
            self.assertIsNone(long["metrics"]["translation_3d_m"]["rmse"])
        self.assertIn("固定时距相对位姿差", (self.root/"report"/"report.md").read_text(encoding="utf-8"))

    def test_rpe_sparse_keyframes_do_not_bridge_estimate_outages(self):
        config = self.config()
        config["relative_motion"] = {"deltas_sec": [1], "max_estimate_gap_sec": .25}
        reference, _, estimates = ae.load_inputs(config)
        result, _ = diagnostic.relative_motion(config, reference, estimates[0])
        self.assertEqual(result[0]["paired_count"], 0)
        self.assertIsNone(result[0]["metrics"]["translation_3d_m"]["rmse"])
        self.assertGreater(result[0]["rejected_reasons"]["estimate_segment_discontinuous"], 0)

    def test_rpe_standard_se3_error_uses_full_rotation_and_body_translation(self):
        # Distinct noncommuting start rotations make simple world displacement
        # subtraction incorrect and verify the standard inverse-relative chain.
        def rx(a):
            a = np.radians(a)
            return np.array([[1, 0, 0], [0, np.cos(a), -np.sin(a)], [0, np.sin(a), np.cos(a)]])
        def rz(a):
            a = np.radians(a)
            return np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        ref_r = [rz(30)@rx(20), rz(60)@rx(-10)]
        est_r = [rx(-15)@rz(40), rx(10)@rz(80)]
        ref_p = [np.array([1., 2, 3]), np.array([4., 6, 5])]
        est_p = [np.array([0., -1, 2]), np.array([3., 4, 8])]
        ref_rows = [ae.Pose(i*ae.NS, ref_p[i], ae.matrix_quat(ref_r[i]), "Aclock") for i in range(2)]
        est_rows = [ae.Pose(i*ae.NS, est_p[i], ae.matrix_quat(est_r[i]), "Aclock") for i in range(2)]
        interval = [{"start_ns": 0, "end_ns": 2*ae.NS, "epoch": "Aclock"}]
        config = {"analysis_intervals": interval, "relative_motion": {"deltas_sec": [1], "max_estimate_gap_sec": 1.1}}
        estimate = {"name": "full_se3", "rows": est_rows, "metadata": {"recorded_intervals": interval}, "offset": 0, "reference_epoch": "Aclock"}
        results, _ = diagnostic.relative_motion(config, ae.Reference(ref_rows, interval, 2*ae.NS), estimate)
        result = results[0]
        r_ref, r_est = ref_r[0].T@ref_r[1], est_r[0].T@est_r[1]
        t_ref, t_est = ref_r[0].T@(ref_p[1]-ref_p[0]), est_r[0].T@(est_p[1]-est_p[0])
        error_t, error_r = r_ref.T@(t_est-t_ref), r_ref.T@r_est
        self.assertEqual(result["paired_count"], 1)
        self.assertAlmostEqual(result["metrics"]["translation_3d_m"]["rmse"], np.linalg.norm(error_t))
        self.assertAlmostEqual(result["metrics"]["rotation_geodesic_deg"]["rmse"], ae.rotation_angle_deg(error_r))
        for key, expected in zip(("error_x_m", "error_y_m", "error_z_m"), error_t):
            self.assertAlmostEqual(result["metrics"][key]["bias"], expected)

    def test_plots_keep_diagnostic_warning_and_reject_acceptance_claim(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        import accuracy_plot_report as plot
        report = diagnostic.diagnose(self.config(), self.root/"report")
        source = self.root/"report"/"report.json"
        output = self.root/"plots"
        plot.render(source, output)
        self.assertIn("A-ONLY ALIGNED DIAGNOSTIC", (output/"index.html").read_text(encoding="utf-8"))
        self.assertEqual(json.loads((output/"plot_manifest.json").read_text())["report_type"],
                         "mapping_a_aligned_diagnostic")
        report["independent_accuracy_assessed"] = True
        source.write_text(json.dumps(report), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "explicitly exclude independent accuracy"):
            plot.render(source, self.root/"bad_plots")
        self.assertFalse((self.root/"bad_plots").exists())

    def test_rejects_localization_other_recording_and_threshold_claims(self):
        self.raw["thresholds"] = {"position_3d_rmse_m": .1}
        with self.assertRaisesRegex(ValueError, "do not accept accuracy thresholds"):
            self.config()
        self.raw["thresholds"] = None
        path = self.root/"frontend_kf.json"
        meta = json.loads(path.read_text())
        meta["trajectory_kind"] = "localization"
        path.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, "cannot claim localization"):
            diagnostic.diagnose(self.config(), self.root/"bad")
        meta["trajectory_kind"] = "frontend"
        meta["source_recording_id"] = "independent_C"
        path.write_text(json.dumps(meta))
        with self.assertRaisesRegex(ValueError, "same physical source_recording_id"):
            diagnostic.diagnose(self.config(), self.root/"bad")
        self.assertFalse((self.root/"bad").exists())


if __name__ == "__main__":
    unittest.main()
