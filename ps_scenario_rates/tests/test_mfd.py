"""Unit tests for step B4 against analytic results."""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from ptha_rates import mfd

LO = np.array([6.5, 7.0, 7.5, 8.0, 8.5, 9.0]); HI = np.array([7.0, 7.5, 8.0, 8.5, 9.0, 9.5])


class TestMfd(unittest.TestCase):
    def test_truncated_gr_bin_shapes(self):
        b, mmax = 0.9, 9.2
        g, mean_m0, _ = mfd.bin_shapes("truncated_gr", LO, HI, b, mmax, mmax)
        np.testing.assert_allclose(g[:5], 10 ** (-b * LO[:5]) - 10 ** (-b * HI[:5]))     # Mmax cancels below the top
        self.assertAlmostEqual(g[5], 10 ** (-b * 9.0) - 10 ** (-b * 9.2))                # top bin cut at Mmax
        # analytic mean moment of a bin: b/(1.5-b) * 10^9.1 * [10^((1.5-b)M)] / g
        an = b / (1.5 - b) * 10 ** 9.1 * (10 ** ((1.5 - b) * 7.5) - 10 ** ((1.5 - b) * 7.0)) / g[1]
        self.assertAlmostEqual(mean_m0[1] / an, 1.0, places=4)

    def test_bins_above_mmax_are_empty(self):
        r, *_ = mfd.bin_rates("truncated_gr", LO, HI, 1.0, 7.8, budget=1e20)
        self.assertTrue((r[3:] == 0).all() and (r[:3] > 0).all())

    def test_xi_analytic(self):
        b, mmax = 0.9, 9.2
        x = mfd.xi("truncated_gr", 6.5, b, mmax, mmax, 9.5)
        self.assertAlmostEqual(x, 1 - 10 ** (-(1.5 - b) * (mmax - 6.5)), places=4)

    def test_geodynamic_rates_release_the_budget(self):
        for form in ("truncated_gr", "tapered_pareto"):
            for tail in ("truncate_at_last_bin", "not_modelled"):
                r, A, x, mean_m0 = mfd.bin_rates(form, LO, HI, 0.85, 9.1, tail=tail, budget=3e20)
                self.assertAlmostEqual(np.nansum(r * mean_m0) / (x * 3e20), 1.0, places=9)
                self.assertTrue(0 < x <= 1)

    def test_seismic_rates_reproduce_the_observed_rate(self):
        for form in ("truncated_gr", "tapered_pareto"):
            r, A, *_ = mfd.bin_rates(form, LO, HI, 0.9, 9.3, rate_mc=120.0, mc=4.5)
            self.assertAlmostEqual(float(mfd.cumulative(form, 4.5, 0.9, 9.3, A, HI[-1])) / 120.0, 1.0)
            self.assertAlmostEqual(r.sum() / float(mfd.cumulative(form, 6.5, 0.9, 9.3, A, HI[-1])), 1.0)

    def test_budget_scales_rates_and_mmax_lowers_them(self):
        r1, *_ = mfd.bin_rates("truncated_gr", LO, HI, 0.9, 9.0, budget=1e20)
        r2, *_ = mfd.bin_rates("truncated_gr", LO, HI, 0.9, 9.0, budget=2e20)
        r3, *_ = mfd.bin_rates("truncated_gr", LO, HI, 0.9, 9.5, budget=1e20)
        np.testing.assert_allclose(r2, 2 * r1)
        self.assertTrue((r3[:5] < r1[:5]).all())
        # the scale is budget / moment of the whole distribution, which grows as 10^((1.5 - b) Mmax)
        self.assertAlmostEqual(r1[0] / r3[0], 10 ** ((1.5 - 0.9) * 0.5), places=3)

    def test_tapered_tends_to_gr_far_below_the_corner(self):
        g_t, *_ = mfd.bin_shapes("tapered_pareto", LO[:2], HI[:2], 0.9, 12.0, np.inf)
        g_g, *_ = mfd.bin_shapes("truncated_gr", LO[:2], HI[:2], 0.9, 12.0, 12.0)
        self.assertAlmostEqual((g_t[0] / g_t[1]) / (g_g[0] / g_g[1]), 1.0, places=3)

    def test_tapered_has_fewer_large_events_than_gr_with_the_same_small_ones(self):
        rt, At, *_ = mfd.bin_rates("tapered_pareto", LO, HI, 0.9, 9.0, rate_mc=100.0, mc=4.5)
        rg, Ag, *_ = mfd.bin_rates("truncated_gr", LO, HI, 0.9, 9.6, rate_mc=100.0, mc=4.5)
        self.assertLess(rt[-1], rg[-1])
        self.assertAlmostEqual(rt[0] / rg[0], 1.0, places=1)


if __name__ == "__main__":
    unittest.main()
