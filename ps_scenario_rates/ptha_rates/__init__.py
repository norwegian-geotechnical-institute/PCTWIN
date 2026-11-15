"""Annual rates for PS scenarios on every branch of the logic tree.

The method, the steps and the decisions are described in docs/guide.md:
B0 configuration, B1 input tables, B2 loading and budget, B3 catalogue and
parameters, B4 bin rates, B5 shares, B6 segmented branches, B7 scenario rates,
B8 end-to-end checks.
"""
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # ps_scenario_rates
# Input files are read from ps_scenario_rates/data, or from the folder
# named by the environment variable PTHA_DATA (see data/README.md).
DATA = os.environ.get("PTHA_DATA", os.path.join(ROOT, "data"))
CONFIG = os.path.join(ROOT, "config")
PRODUCTS = os.path.join(ROOT, "products")
CHECKS = os.path.join(ROOT, "checks")
