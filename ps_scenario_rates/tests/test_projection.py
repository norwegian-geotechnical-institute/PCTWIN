"""Unit tests for the provider's UTM areas (decision 29). Reference values from pyproj 3.8 (EPSG 32746 and 32646)."""
import os
import sys
import unittest
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from ptha_rates import projection


class TestProjection(unittest.TestCase):
    def test_against_pyproj(self):
        x, y = projection.utm_xy(110.0, -8.0, 46, south=True)
        self.assertAlmostEqual(float(x), 2400550.059705855, places=3); self.assertAlmostEqual(float(y), 9075569.283031084, places=3)
        x, y = projection.utm_xy(93.0, 10.0, 46)
        self.assertAlmostEqual(float(x), 500000.0, places=3); self.assertAlmostEqual(float(y), 1105412.4913010783, places=3)

    def test_area_stretch_grows_away_from_the_central_meridian(self):
        # small flat triangle (0.1 deg) at the central meridian and 20 deg east of it, on the equator
        tri = lambda lon: np.array([[[lon, 0.0, 0.0], [lon + 0.1, 0.0, 0.0], [lon, 0.1, 0.0]]])
        a0, a20 = projection.triangle_area_utm_km2(tri(93.0), 46)[0], projection.triangle_area_utm_km2(tri(113.0), 46)[0]
        self.assertAlmostEqual(a0, 61.50, places=1)                     # 0.5 x (11.1 km)^2 at the central meridian
        self.assertGreater(a20 / a0, 1.10)                               # transverse Mercator scale ~ 1 + x^2 / R^2


if __name__ == "__main__":
    unittest.main()
