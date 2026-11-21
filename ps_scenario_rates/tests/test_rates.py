"""Unit tests for step B7: rigidity factor, scaling law and scaling sample as levels (decisions 27, 28)."""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import config, rates, shares


def D_two_positions():
    # two positions in one bin: rigidity 00 / scaling law 01, rigidity 01 / scaling law 02; unconstrained shares 1 each
    pos = pd.DataFrame({"bin": [0, 0], "rig": [0, 1], "scal": [1, 2]})
    R = dict(rate=np.ones((2, 4, 1, 1, 1, 1)), variants=shares.VARIANTS, forms=["f"], mmax=["m"])
    fac = {"preferred": np.array([[1.0, 1.25], [1.0, 1.25]])}            # unit x position
    return dict(R=R, pos=pos, w_rig=np.array([0.4, 0.6]), pooled=np.array([0.4 * 0.5, 0.6 * 0.5]), rig_factor=fac,
                U=np.ones((4, 2)), X=np.ones((2, 1)))


class TestRates(unittest.TestCase):
    br = staticmethod(lambda src, c, law="pooled": pd.Series(dict(rate_source=src, coupling=c, segmentation="unconstrained",
                                                                   mfd_form="f", b_index=1, mmax="m", scaling_law=law)))

    def test_rigidity_factor_on_geodynamic_only(self):
        D = D_two_positions()
        np.testing.assert_allclose(rates.position_rates(D, self.br("seismic", "none")), [0.2, 0.3])
        np.testing.assert_allclose(rates.position_rates(D, self.br("geodynamic", "preferred")), [0.2, 0.3 * 1.25])

    def test_scaling_law_as_level(self):
        D = D_two_positions()
        np.testing.assert_allclose(rates.position_rates(D, self.br("seismic", "none", "02")), [0.0, 0.6])

    def test_sample_split_as_level(self):
        # position 0: samples 0 (two slip realisations) and 1 (one); position 1: sample 1 only
        scen = dict(position=np.array([0, 0, 0, 1]), sreal=np.array([0, 0, 1, 1]), split=np.array([0.25, 0.25, 0.5, 1.0]))
        s0 = rates.sample_split(scen, pd.Series(dict(scaling_sample="00")))
        np.testing.assert_allclose(s0, [0.5, 0.5, 0.0, 1.0])            # position 1 has no sample 00: keeps its pooled split
        mean = (s0 + rates.sample_split(scen, pd.Series(dict(scaling_sample="01")))) / 2
        np.testing.assert_allclose(mean, scen["split"])                # equal weights over samples = pooled split

    def test_rigidity_in_use(self):
        tree = config.load()["tree"]; tree["options"]["rigidity"]["without_footprints"] = "skip"
        self.assertEqual(config.rigidity_in_use(tree, footprints_have_slip=False), {"00": 1.0})
        self.assertEqual(config.rigidity_in_use(tree, footprints_have_slip=True), {"00": 0.4, "01": 0.6})


if __name__ == "__main__":
    unittest.main()
