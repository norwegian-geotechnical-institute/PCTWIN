"""Unit tests for step B0. Run from ps_scenario_rates:  python -m unittest discover tests"""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import copy
import numpy as np
from ptha_rates import branches, config


def v10(tree):
    """The b-value and coupling settings of version 1.0, for the tests of those rules."""
    t = copy.deepcopy(tree)
    t["levels"]["b_value"] = {"n": 10, "min": 0.7, "max": 1.2, "weights": "linear"}
    t["levels"]["coupling"].update(rule="offset", offset=0.2, floor=0.1, ceiling=1.0)
    t["options"]["budget_area"] = "mesh"
    return t


class TestBranches(unittest.TestCase):
    def setUp(self):
        self.cfg = config.load()
        self.tree = v10(self.cfg["tree"])

    def test_config_valid(self):
        errors, _ = config.validate(self.cfg)
        self.assertEqual(errors, [])

    def test_b_values(self):
        b = branches.b_values(self.tree)
        self.assertEqual(len(b), 10)
        self.assertAlmostEqual(b[0], 0.7); self.assertAlmostEqual(b[-1], 1.2)

    def test_b_weights_peak_and_sum(self):
        b = branches.b_values(self.tree)
        for b_hat in (0.67, 0.787, 0.86, 1.028, 1.3):
            w = branches.b_weights(b_hat, self.tree)
            self.assertAlmostEqual(w.sum(), 1.0)
            self.assertTrue((w > 0).all())
            self.assertEqual(np.argmax(w), np.argmin(np.abs(b - b_hat)))

    def test_b_weights_symmetric_at_centre(self):
        w = branches.b_weights(0.95, self.tree)
        np.testing.assert_allclose(w, w[::-1])

    def test_coupling_levels(self):
        c = branches.coupling_levels(0.2, self.tree)       # Sumba: floor applies
        self.assertAlmostEqual(c["low"], 0.1); self.assertAlmostEqual(c["high"], 0.4)
        c = branches.coupling_levels(1.0, self.tree)       # Arakan: ceiling applies
        self.assertAlmostEqual(c["low"], 0.8); self.assertAlmostEqual(c["high"], 1.0)

    def test_b_grid_v11(self):
        tree = self.cfg["tree"]; b = branches.b_values(tree); w = branches.b_weights(0.86, tree)
        np.testing.assert_allclose(b, [0.56, 0.71, 0.86, 1.01, 1.16])
        np.testing.assert_allclose(w, np.array([1, 2, 3, 2, 1]) / 9)
        self.assertAlmostEqual((b * w).sum(), 0.86)                     # weighted mean = catalogue b (decision 24)

    def test_coupling_factor_rule(self):
        tree = self.cfg["tree"]                                         # decision 25: x0.5 / x1 / x2 within 0.1 to 1.3
        self.assertEqual(branches.coupling_levels(0.8, tree), {"low": 0.4, "preferred": 0.8, "high": 1.3})
        self.assertEqual(branches.coupling_levels(0.176, tree), {"low": 0.1, "preferred": 0.176, "high": 0.352})

    def test_as_level_expands_the_branches(self):
        cfg = copy.deepcopy(self.cfg); cfg["tree"]["event_tree"]["as_level"] = ["scaling_law", "scaling_sample"]
        t = branches.branch_table(cfg, "sunda"); n0 = len(branches.branch_table(self.cfg, "sunda"))
        self.assertEqual(len(t), n0 * 2 * 10); self.assertAlmostEqual(t.weight.sum(), 1.0); self.assertTrue(t.branch_id.is_unique)

    def test_branch_table(self):
        t = branches.branch_table(self.cfg, "sunda")
        self.assertEqual(len(t), 240)
        self.assertEqual(int((t.rate_source == "geodynamic").sum()), 180)
        self.assertAlmostEqual(t.weight.sum(), 1.0)
        self.assertTrue(t.branch_id.is_unique)
        self.assertEqual(t.spatial_variant.nunique(), 8)


if __name__ == "__main__":
    unittest.main()
