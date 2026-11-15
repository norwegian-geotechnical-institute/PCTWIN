"""Unit tests for steps B5 to B7 on a small synthetic zone with two segments."""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import shares


def positions():
    # one cell (cell 0) with four positions: two inside segment 1, one inside segment 2, one crossing (60% / 40%)
    return pd.DataFrame({"position": [0, 1, 2, 3], "cell": [0, 0, 0, 0],
                         "seg_1": [1.0, 1.0, 0.0, 0.6], "seg_2": [0.0, 0.0, 1.0, 0.4]})


class TestShares(unittest.TestCase):
    def test_unconstrained_follow_loading(self):
        vals = {v: np.array([30.0, 10.0]) for v in shares.VARIANTS}
        U = shares.unconstrained_shares(positions(), vals, 2)[0]
        w = np.array([30.0, 30.0, 10.0, 0.6 * 30 + 0.4 * 10])
        np.testing.assert_allclose(U, w / w.sum())
        self.assertAlmostEqual(U.sum(), 1.0)

    def test_unconstrained_scale_invariant(self):
        a = shares.unconstrained_shares(positions(), {v: np.array([30.0, 10.0]) for v in shares.VARIANTS}, 2)
        b = shares.unconstrained_shares(positions(), {v: np.array([60.0, 20.0]) for v in shares.VARIANTS}, 2)
        np.testing.assert_allclose(a, b)        # doubling the loading everywhere leaves the shares unchanged

    def test_segmented(self):
        X = shares.segmented_shares(positions(), 2, 1.0)
        np.testing.assert_allclose(X[:, 0], np.array([1, 1, 0, 0.6]) / 2.6)
        np.testing.assert_allclose(X[:, 1], np.array([0, 0, 1, 0.4]) / 1.4)
        np.testing.assert_allclose(X.sum(axis=0), [1.0, 1.0])

    def test_segmented_rates_conserve_each_segment(self):
        X = shares.segmented_shares(positions(), 2, 1.0)
        r = np.array([0.02, 0.005])                      # bin rate of segment 1 and of segment 2
        rate = (X * r).sum(axis=1)
        self.assertAlmostEqual(rate.sum(), r.sum())
        self.assertAlmostEqual(rate[3], 0.02 * 0.6 / 2.6 + 0.005 * 0.4 / 1.4)

    def test_crossing_penalty(self):
        X1 = shares.segmented_shares(positions(), 2, 1.0); X2 = shares.segmented_shares(positions(), 2, 0.5)
        np.testing.assert_allclose(X2.sum(axis=0), [1.0, 1.0])
        self.assertLess(X2[3].sum(), X1[3].sum()); self.assertGreater(X2[0, 0], X1[0, 0])
        np.testing.assert_allclose(X2[:, 0], np.array([1, 1, 0, 0.3]) / 2.3)

    def test_position_table(self):
        scen = {"position": np.array([0, 0, 1]), "rig": np.array([0, 0, 0]), "scal": np.array([1, 1, 1]), "bin": np.array([0, 0, 0]),
                "tri": np.array([5, 5, 9]), "split": np.array([0.5, 0.5, 1.0]), "total_M0": np.array([1e19, 3e19, 2e19]),
                "seg_1": np.array([1.0, 0.5, 0.0]), "seg_2": np.array([0.0, 0.5, 1.0])}
        p = shares.position_table(scen, 2)
        np.testing.assert_allclose(p.seg_1, [0.75, 0.0]); np.testing.assert_allclose(p.seg_2, [0.25, 1.0])
        np.testing.assert_allclose(p.m0_nm, [2e19, 2e19]); self.assertEqual(p.n_scenarios.tolist(), [2, 1])


if __name__ == "__main__":
    unittest.main()
