"""Unit tests for step B1 on small synthetic cases."""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import inputs


class TestInputs(unittest.TestCase):
    def test_parse_names(self):
        p = inputs.parse_names(pd.Series(["PS_sun_00_6_8012_01_03_05747_002", "PS_mak_01_9_0260_02_00_00012_010"]))
        self.assertEqual(p.iloc[0].tolist(), [0, 68012, 1, 3, 5747, 2])
        self.assertEqual(p.iloc[1].tolist(), [1, 90260, 2, 0, 12, 10])

    def test_bin_table(self):
        b = inputs.bin_table(np.array([65000, 70000, 80000, 70000]))
        np.testing.assert_allclose(b.mw_centre, [6.5, 7.0, 8.0])
        np.testing.assert_allclose(b.mw_lower, [6.25, 6.75, 7.5])
        np.testing.assert_allclose(b.mw_upper, [6.75, 7.5, 8.5])
        np.testing.assert_allclose(b.m0_centre_nm, 10 ** (1.5 * b.mw_centre + 9.1))

    def test_counts_and_split(self):
        # one position with two scaling samples: sample 0 has three slip realisations, sample 1 has one;
        # a second position with a single scenario
        s = pd.DataFrame({"rig": 0, "scal": 1, "mw": [70000] * 4 + [80000], "tri": [5, 5, 5, 5, 9],
                          "sreal": [0, 0, 0, 1, 0], "slip": [1, 2, 3, 0, 0]})
        bins = inputs.bin_table(s.mw.to_numpy())
        s = inputs.add_bins_and_counts(s, bins)
        self.assertEqual(s.n_sreal.tolist(), [2, 2, 2, 2, 1])
        self.assertEqual(s.n_slip.tolist(), [3, 3, 3, 1, 1])
        np.testing.assert_allclose(s.split, [1 / 6, 1 / 6, 1 / 6, 1 / 2, 1.0])
        np.testing.assert_allclose(s.groupby("position").split.sum(), 1.0)

    def test_segment_map_fills_from_nearest(self):
        tri = pd.DataFrame({"id": [1, 2, 3], "lon": [90.0, 90.1, 95.0], "lat": [0.0, 0.0, 0.0]})
        s = pd.DataFrame({"bin": [0, 0, 1], "tri": [1, 3, 2], "seg_1": [1.0, 0.1, 0.5], "seg_2": [0.0, 0.9, 0.5]})
        seg, n_near = inputs.segment_map(s, tri, 2)
        self.assertEqual(seg.tolist(), [1, 1, 2])      # triangle 2 has no smallest scenario: takes triangle 1's segment
        self.assertEqual(n_near, 1)


if __name__ == "__main__":
    unittest.main()
