"""Step B7. Scenario rates: the mean over all branches, and what is needed to rebuild any branch.

Run from ps_scenario_rates after steps B0 to B6:   python scripts/run_B7_rates.py [zone ...]     (default: sunda)
Writes products/B7_rates/ and checks/B7_rates/.

Products per zone:
  scenario_rates_<zone>.npz        one value per scenario, in the order of products/B1_inputs/scenarios_<zone>.npz:
                                   the version, the parts of the name, position, split, and the weighted mean annual rate over
                                   all branches and over four groups of branches
  scenario_mean_rates_<zone>.csv.gz  scenario name and mean annual rate, all scenarios
  position_rates_<zone>.csv        mean annual rate of every rupture position
  README.md                        how to rebuild the rate of any scenario on any branch
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, __version__, binrates, config, rates
from ptha_rates.report import CheckReport

README = """# Scenario rates ({zone})

Version {version}. Annual rates of the PS scenarios of rigidity model 00 (homogeneous) on the {n_branch} branches of the logic tree.
Scenarios of the other rigidity models have no rate in this run.

## Mean rates

`scenario_mean_rates_{zone}.csv.gz`: scenario name and annual rate, weighted mean over all branches.
`scenario_rates_{zone}.npz` holds the same with the parts of the name and the means over four groups of branches
(seismic, geodynamic, unconstrained, segmented).

## Any branch

A rate is not stored for every scenario on every branch. It is rebuilt from three small tables:

    rate(scenario, branch) = pooled weight x bin rate(branch) x share(position) x split(scenario)

- branch table: `../B0_config/branch_table_{zone}.csv` (branch id, alternative at each level, weight)
- bin rates: `../B4_bin_rates/bin_rates_{zone}.npz` (unit x variant x form x b x Mmax x bin)
- shares: `../B5_shares/shares_{zone}.npz` (unconstrained branches) and `../B6_segmented/shares_segmented_{zone}.npz`
- split and position of each scenario: `../B1_inputs/scenarios_{zone}.npz`

In Python, from the folder `ps_scenario_rates`:

    from ptha_rates import rates
    r = rates.scenario_rates("{zone}", "G-U-TG-b04-Marea-Cpref")    # one value per scenario

Scenario names in the same order: `rates.scenario_names("{code}", numpy.load("products/B1_inputs/scenarios_{zone}.npz"))`.
"""


def main(zones):
    cfg = config.load()
    out = os.path.join(PRODUCTS, "B7_rates"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B7_rates", "scenario rates")
    for zone in zones:
        z = cfg["zones"][zone]
        D = rates.load(zone, cfg); t = D["table"]; w = t.weight.to_numpy(); pos = D["pos"]
        scen = np.load(os.path.join(PRODUCTS, "B1_inputs", f"scenarios_{zone}.npz"))
        bins = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"bins_{zone}.csv"))
        groups = {"all": None, "seismic": (t.rate_source == "seismic").to_numpy(), "geodynamic": (t.rate_source == "geodynamic").to_numpy(),
                  "unconstrained": (t.segmentation == "unconstrained").to_numpy(), "segmented": (t.segmentation == "segmented").to_numpy()}
        PM = {g: rates.mean_position_rates(D, m) for g, m in groups.items()}
        Z = binrates.zone_bin_rates(D["R"], t)
        active = sum(wt for wt in cfg["tree"]["event_tree"]["rigidity_model"].values())
        expect = float((Z.sum(1) * w).sum()) * active
        rep.check(f"{zone}: mean rates add up to the mean bin rates", "Sum of the mean position rates must equal the weighted mean of the zone's total rate",
                  bool(np.isclose(PM["all"].sum(), expect, rtol=1e-9)), f"{PM['all'].sum():.6f} per year against {expect:.6f}")
        sr = PM["all"][scen["position"]] * scen["split"]
        rep.check(f"{zone}: scenario rates add up to the position rates", "The split divides a position's rate among its realisations without loss",
                  bool(np.isclose(sr.sum(), PM["all"].sum(), rtol=1e-9)), f"{len(sr):,} scenarios, total {sr.sum():.6f} per year")
        rep.check(f"{zone}: every scenario has a positive finite mean rate", "No scenario of the rigidity model in use is left without a rate",
                  bool(np.isfinite(sr).all() and (sr > 0).all()), f"smallest {sr.min():.2e}, largest {sr.max():.2e} per year")
        br = t.iloc[len(t) // 2]
        direct = rates.position_rates(D, br)[scen["position"]] * scen["split"]
        rep.check(f"{zone}: reader returns the rates of a single branch", "rates.scenario_rates must agree with the direct calculation",
                  bool(np.allclose(rates.scenario_rates(zone, br.branch_id), direct)), f"branch {br.branch_id}: total {direct.sum():.4f} per year")
        kb = pos.bin.to_numpy()
        byb = np.bincount(kb, weights=PM["all"], minlength=len(bins))
        rep.check(f"{zone}: mean rate by magnitude", "Annual rate of events at or above a magnitude, from the scenario rates", None,
                  "; ".join(f"Mw>={bins.mw_lower[k]:.2f}: " + (f"{byb[k:].sum():.2f} per yr" if byb[k:].sum() >= 1 else f"1 per {1 / byb[k:].sum():,.0f} yr") for k in (0, 5, 10, 16, 21, 26)))
        rep.check(f"{zone}: mean rate of all scenarios by group of branches", "Total annual rate of modelled events", None,
                  "; ".join(f"{g}: {v.sum():.2f} per yr" for g, v in PM.items()))
        i = int(np.argmax(sr)); nm = rates.scenario_names(z["code"], {k: scen[k][[i]] for k in ("rig", "mw", "scal", "sreal", "tri", "slip")})[0]
        rep.check(f"{zone}: scenario with the highest mean rate", "For orientation", None, f"{nm}: 1 per {1 / sr[i]:,.0f} yr")

        cols = {k: scen[k] for k in ("rig", "mw", "scal", "sreal", "tri", "slip", "bin", "position", "split")}
        for g, v in PM.items():
            cols["rate_mean" if g == "all" else f"rate_mean_{g}"] = v[scen["position"]] * scen["split"]
        np.savez_compressed(os.path.join(out, f"scenario_rates_{zone}.npz"), version=np.array(__version__), **cols)
        pd.DataFrame({"scenario": rates.scenario_names(z["code"], scen), "rate_per_yr": sr}).to_csv(
            os.path.join(out, f"scenario_mean_rates_{zone}.csv.gz"), index=False, float_format="%.6e", compression={"method": "gzip", "mtime": 0})   # mtime 0: the file is identical from run to run
        pp = pos[["position", "rig", "scal", "bin", "tri", "n_scenarios"]].copy(); pp["mw_centre"] = bins.mw_centre.to_numpy()[kb]
        for g, v in PM.items():
            pp["rate_mean" if g == "all" else f"rate_mean_{g}"] = v
        pp.to_csv(os.path.join(out, f"position_rates_{zone}.csv"), index=False, float_format="%.6e")
        with open(os.path.join(out, "README.md"), "w", encoding="utf-8") as f:
            f.write(README.format(zone=zone, code=z["code"], n_branch=len(t), version=__version__))
    rep.note("Mean rates are weighted over the branches with the weights of the branch table. No Bayesian update is applied (decision 13).")
    rep.note("Rates are given to scenarios of rigidity model 00 only (decision 14).")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda"]) else 1)
