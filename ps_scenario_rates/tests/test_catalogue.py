import unittest
import pandas as pd
from ptha_rates import catalogue


class TestCatalogueFilter(unittest.TestCase):
    def test_reasons(self):
        c = pd.DataFrame({"year": [2013, 2013, 2011, 1945, 1977], "month": [9, 9, 1, 11, 3], "day": [24, 28, 18, 27, 21],
                          "lat": [26.95, 27.18, 28.04, 24.98, 27.61], "lon": [65.50, 65.51, 63.85, 63.68, 56.39],
                          "mw": [7.8, 6.8, 7.4, 8.1, 6.7], "dist_km": [4.0, 4.9, 2.9, 1.2, 6.6]})
        flt = {"lon_min": 57.3, "max_dist_km": 25,
               "exclude": [{"name": "seq", "from": "2013-09-24", "to": "2013-12-31", "lat": 26.95, "lon": 65.5, "radius_deg": 1.0},
                           {"name": "slab", "from": "2011-01-18", "to": "2011-01-18", "lat": 28.04, "lon": 63.85, "radius_deg": 0.5, "mw_min": 7.0}]}
        why = catalogue.excluded(c, flt)
        self.assertEqual(list(why), ["seq", "seq", "slab", "", "west of 57.3 E"])

    def test_no_filter(self):
        c = pd.DataFrame({"year": [2000], "month": [1], "day": [1], "lat": [0.0], "lon": [0.0], "mw": [6.0], "dist_km": [100.0]})
        self.assertEqual(list(catalogue.excluded(c, None)), [""])


if __name__ == "__main__":
    unittest.main()
