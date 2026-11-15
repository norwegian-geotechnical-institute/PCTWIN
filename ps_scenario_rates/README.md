# PS scenario rates

Annual rates of the predominant-subduction (PS) tsunami scenarios on every branch of the logic tree.
This covers steps 2 to 5 of the workflow in the repository README (logic tree, MFD and bin rates,
distribution of the rates over the scenarios) and stops where the hazard aggregation starts.

Status: **Sunda, first run** (homogeneous rigidity model, 480 branches, 2,075,771 scenarios). Makran is not done yet.

**Read `docs/guide.md` first.** It explains the method, each step with its limitations, the decisions taken,
the results, and how to turn the rates into hazard curves.

## Use the rates (no input data needed)

```python
# from the folder ps_scenario_rates
import numpy as np, pandas as pd
from ptha_rates import rates

branches = pd.read_csv("products/B0_config/branch_table_sunda.csv")      # 480 branches with their weights
scen = np.load("products/B1_inputs/scenarios_sunda.npz")
names = rates.scenario_names("sun", scen)                                # scenario names, fixed order
r = rates.scenario_rates("sunda", branches.branch_id[0])                 # annual rate of every scenario on one branch
mean = np.load("products/B7_rates/scenario_rates_sunda.npz")["rate_mean"]  # weighted mean over all branches
```

`products/B7_rates/scenario_mean_rates_sunda.csv.gz` holds scenario name and mean rate as text.

## Run it

Needs Python 3 with numpy, pandas, matplotlib and pyyaml. Two of the three input files are in `data/`; the third is too large for the repository and has to be placed there (see `data/README.md`). Then:

```
python scripts/run_all.py
```

This runs the unit tests, steps B0 to B8, the sensitivity test, and refreshes `docs/figures/`. It takes about
three minutes for Sunda and stops at the first failed check. Each step can also be run alone
(`python scripts/run_B4_bin_rates.py`), in order.

To start again from nothing, delete everything the steps produced and rerun:

```
python scripts/clean.py
python scripts/run_all.py
```

`clean.py --dry-run` lists what would be deleted. It never touches `config/`, `data/`, the code or the documents, and it keeps
`products/B1_inputs` when the rupture fractions file is not there to rebuild it.

## Layout

| Folder | Content |
|---|---|
| `config/` | `logic_tree.yaml` (levels, weights, options) and `zones.yaml` (zone and segment values). Every modelling choice is set here. |
| `ptha_rates/` | The code, one module per task. |
| `scripts/` | `run_all.py`, `clean.py`, one run script per step (`run_B<n>_<name>.py`) and the sensitivity test (`sens_coupling.py`). |
| `tests/` | Unit tests on small synthetic cases with known answers. |
| `products/<step>/` | What each step produces and later steps read. Kept in the repository so the rates can be used without the input data. |
| `checks/<step>/` | Written by every run, not kept in the repository: `report.md` per step and `INDEX.md` listing every check with its result. |
| `data/` | Input files (the rupture fractions file has to be added; see its README). |
| `docs/` | `guide.md` and the figures. |

## Conventions

- A modelling choice is changed in `config/`, never in the code, and recorded in the decision table of `docs/guide.md`.
- A step is not built on until its checks pass. Checks marked INFO report a number for the record; PASS and FAIL are tests.
- Sensitivity tests change nothing in `config/` or in the products of the steps; they write to `products/sensitivity/`.
