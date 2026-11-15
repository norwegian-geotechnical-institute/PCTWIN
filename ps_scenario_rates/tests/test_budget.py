"""Unit tests for step B2 on a two-segment synthetic mesh."""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import budget, config


class TestBudget(unittest.TestCase):
    def setUp(self):
        self.tree = config.load()["tree"]
        self.zone = {"rigidity_pa": 3.0e10, "segments": [{"name": "A", "coupling": 0.5, "convergence": 40.0},
                                                         {"name": "B", "coupling": 0.2, "convergence": 60.0}]}
        self.tri = pd.DataFrame({"id": [1, 2, 3], "segment": [1, 1, 2], "area_km2": [100.0, 200.0, 50.0]})

    def test_loading(self):
        l = budget.loading_table(self.tri, self.zone, self.tree)
        np.testing.assert_allclose(l.loading_preferred_mm_yr, [20.0, 20.0, 12.0])
        np.testing.assert_allclose(l.loading_low_mm_yr, [12.0, 12.0, 6.0])        # B: floor 0.1 applies
        np.testing.assert_allclose(l.loading_high_mm_yr, [28.0, 28.0, 24.0])

    def test_budget(self):
        t = budget.budget_table(budget.loading_table(self.tri, self.zone, self.tree), self.zone).set_index("unit")
        a = 3.0e10 * 300e6 * 0.020        # rigidity x area (m2) x loading (m/yr)
        b = 3.0e10 * 50e6 * 0.012
        self.assertAlmostEqual(t.budget_preferred_nm_yr["A"] / a, 1.0)
        self.assertAlmostEqual(t.budget_preferred_nm_yr["B"] / b, 1.0)
        self.assertAlmostEqual(t.budget_preferred_nm_yr["whole_zone"] / (a + b), 1.0)
        self.assertAlmostEqual(t.coupling_preferred["whole_zone"], (a + b) / (3.0e10 * (300e6 * 0.040 + 50e6 * 0.060)))

    def test_doubling_coupling_doubles_budget(self):
        z2 = {"rigidity_pa": 3.0e10, "segments": [{"name": "A", "coupling": 0.25, "convergence": 40.0},
                                                  {"name": "B", "coupling": 0.1, "convergence": 60.0}]}
        t1 = budget.budget_table(budget.loading_table(self.tri, self.zone, self.tree), self.zone)
        t2 = budget.budget_table(budget.loading_table(self.tri, z2, self.tree), z2)
        np.testing.assert_allclose(t1.budget_preferred_nm_yr, 2 * t2.budget_preferred_nm_yr)


if __name__ == "__main__":
    unittest.main()
