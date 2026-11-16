"""Unit tests for rupture footprints on a small synthetic strip of triangles."""
import os
import sys
import tempfile
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import footprints, rates, shares


def strip(n=10):
    # n triangles in a row along a meridian, 10 km2 each, about 11 km apart; the first half is segment 1
    return pd.DataFrame({"id": np.arange(1, n + 1), "lon": 100.0, "lat": 0.1 * np.arange(n), "depth_km": 10.0, "area_km2": 10.0,
                         "segment": np.r_[np.ones(n // 2, int), 2 * np.ones(n - n // 2, int)]})


class TestFootprints(unittest.TestCase):
    def setUp(self):
        self.tri = strip()
        # four positions in one cell: starting on triangles 1, 3, 5 (target 10 km2) and on 5 with target 30 km2
        self.pos = pd.DataFrame({"position": range(4), "cell": 0, "tri": [1, 3, 5, 5]})
        self.fp = footprints.estimate(self.tri, self.pos.tri.to_numpy(), np.array([10.0, 10.0, 10.0, 30.0]))

    def test_estimate_takes_nearest_triangles_up_to_the_area(self):
        n = np.diff(self.fp["indptr"])
        np.testing.assert_array_equal(n, [1, 1, 1, 3])
        np.testing.assert_array_equal(np.sort(self.fp["index"][3:]), [3, 4, 5])           # rows of triangles 4, 5, 6
        np.testing.assert_allclose(footprints.area_km2(self.fp, self.tri), [10, 10, 10, 30])
        np.testing.assert_allclose(np.add.reduceat(self.fp["frac"], self.fp["indptr"][:-1]), 1, rtol=1e-6)

    def test_area_is_limited_by_the_mesh(self):
        fp = footprints.estimate(self.tri, np.array([1]), np.array([1e6]))
        self.assertEqual(len(fp["index"]), len(self.tri))

    def test_segment_fractions(self):
        F = footprints.segment_fractions(self.fp, self.tri.segment.to_numpy(), 2)
        np.testing.assert_allclose(F, [[1, 0], [1, 0], [1, 0], [2 / 3, 1 / 3]], rtol=1e-6)

    def test_unconstrained_shares_follow_the_loading_under_the_footprint(self):
        load = np.where(self.tri.segment == 1, 30.0, 10.0); load[0] = 60.0                # triangle 1 is loaded twice as fast
        U = shares.unconstrained_shares_fp(self.pos, self.fp, {v: load for v in shares.VARIANTS})[0]
        w = np.array([60.0, 30.0, 30.0, (30 + 30 + 10) / 3])
        np.testing.assert_allclose(U, w / w.sum(), rtol=1e-6)

    def test_segmented_shares(self):
        load = np.where(self.tri.segment == 1, 30.0, 10.0)
        X, Ff = shares.segmented_shares_fp(self.pos, self.fp, {v: load for v in shares.VARIANTS}, self.tri.segment.to_numpy(), 2)
        self.assertEqual(X.shape, (4, 4, 2))
        np.testing.assert_allclose(X[0][:, 0], np.array([1, 1, 1, 2 / 3]) / (3 + 2 / 3), rtol=1e-6)   # constant loading cancels
        np.testing.assert_allclose(X[0][:, 1], [0, 0, 0, 1], atol=1e-12)
        np.testing.assert_allclose(X.sum(axis=1), 1)

    def test_position_rates_take_the_shares_of_the_variant(self):
        X = np.zeros((4, 2, 1)); X[0, :, 0] = [1, 0]; X[2, :, 0] = [0, 1]
        D = dict(R=dict(rate=np.ones((2, 4, 1, 1, 1, 1)), variants=shares.VARIANTS, forms=["f"], mmax=["m"]), pos=pd.DataFrame({"bin": [0, 0]}), pooled=np.ones(2), U=None, X=X)
        br = lambda src, c: pd.Series(dict(rate_source=src, coupling=c, segmentation="segmented", mfd_form="f", b_index=1, mmax="m"))
        np.testing.assert_allclose(rates.position_rates(D, br("seismic", "none")), [1, 0])
        np.testing.assert_allclose(rates.position_rates(D, br("geodynamic", "preferred")), [0, 1])

    def test_from_file_averages_the_geometries_of_a_position(self):
        # one position with two scaling samples (split 0.5 each): one covers triangle 1, the other triangles 1 and 2
        scen = dict(rig=np.array([0, 0]), mw=np.array([70000, 70000]), scal=np.array([1, 1]), sreal=np.array([0, 1]), tri=np.array([1, 1]),
                    position=np.array([0, 0]), split=np.array([0.5, 0.5]))
        with tempfile.TemporaryDirectory() as d:
            np.savez(os.path.join(d, "footprints_syn_00.npz"), rig=[0, 0], mw=[70000, 70000], scal=[1, 1], sreal=[0, 1], tri=[1, 1],
                     indptr=[0, 1, 3], triangle=[1, 1, 2])
            fp = footprints.from_file(d, "syn", scen, self.tri, 1)
            np.testing.assert_array_equal(fp["index"], [0, 1]); np.testing.assert_allclose(fp["frac"], [0.75, 0.25])
            scen["sreal"] = np.array([0, 2])
            self.assertRaises(ValueError, footprints.from_file, d, "syn", scen, self.tri, 1)


if __name__ == "__main__":
    unittest.main()
