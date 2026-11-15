"""Unit tests for step B0. Run from ps_scenario_rates:  python -m unittest discover tests"""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from ptha_rates import branches, config


class TestBranches(unittest.TestCase):
    def setUp(self):
        self.cfg = config.load()
        self.tree = self.cfg["tree"]

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

    def test_branch_table(self):
        t = branches.branch_table(self.cfg, "sunda")
        self.assertEqual(len(t), 480)
        self.assertEqual(int((t.rate_source == "geodynamic").sum()), 360)
        self.assertAlmostEqual(t.weight.sum(), 1.0)
        self.assertTrue(t.branch_id.is_unique)
        self.assertEqual(t.spatial_variant.nunique(), 8)


if __name__ == "__main__":
    unittest.main()
