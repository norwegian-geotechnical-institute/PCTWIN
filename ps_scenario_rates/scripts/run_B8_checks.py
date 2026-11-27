"""Step B8. End-to-end checks on the scenario rates.

Run from ps_scenario_rates after step B7:   python scripts/run_B8_checks.py [zone ...]     (default: sunda makran)
Writes checks/B8_end_to_end/ (report, tables and a figure).
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, binrates, config, plotting, rates
from ptha_rates.report import CheckReport
plt = plotting.plt


def wq(x, w, q):
    o = np.argsort(x); cw = np.cumsum(w[o]) / w.sum()
    return x[o][min(np.searchsorted(cw, q), len(x) - 1)]


def main(zones):
    cfg = config.load()
    rep = CheckReport("B8_end_to_end", "end-to-end checks on the scenario rates")
    for zone in zones:
        z = cfg["zones"][zone]; names = [s["name"] for s in z["segments"]]; n_seg = len(names)
        D = rates.load(zone, cfg); t = D["table"]; w = t.weight.to_numpy(); pos = D["pos"]; R = D["R"]
        bins = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"bins_{zone}.csv"))
        bud = pd.read_csv(os.path.join(PRODUCTS, "B2_budget", f"budget_{zone}.csv")).set_index("unit")
        tri = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"triangles_{zone}.csv")).set_index("id")
        cat = pd.read_csv(os.path.join(PRODUCTS, "B3_catalogue", f"catalogue_{zone}.csv"))
        P5 = pd.read_csv(os.path.join(PRODUCTS, "B5_shares", f"positions_{zone}.csv"))
        F = P5[[f"seg_{i + 1}" for i in range(n_seg)]].to_numpy()
        kb = pos.bin.to_numpy(); m0p = pos.m0_nm.to_numpy(); nb = len(bins)
        Z = binrates.zone_bin_rates(R, t)
        geo = (t.rate_source == "geodynamic").to_numpy(); unc = (t.segmentation == "unconstrained").to_numpy()

        err, mom_ratio, cen_ratio, seg_ratio = 0.0, np.full(len(t), np.nan), np.full(len(t), np.nan), np.full((len(t), n_seg), np.nan)
        for i, br in enumerate(t.itertuples()):
            pr = rates.position_rates(D, br)
            byb = np.bincount(kb, weights=pr, minlength=nb)
            ze = rates.expected_bin_rates(D, br)                               # bin rates x pooled rigidity factor (decision 27)
            ok = ze > 0
            err = max(err, np.abs(byb[ok] / ze[ok] - 1).max(), np.abs(byb[~ok]).max() if (~ok).any() else 0.0)
            i_f, ib, im = R["forms"].index(br.mfd_form), int(br.b_index) - 1, R["mmax"].index(br.mmax)
            us = [0] if br.segmentation == "unconstrained" else list(range(1, n_seg + 1))
            r_u = binrates.branch_bin_rates(R, br)[us]
            fu = (D["rig_factor_unit"][br.coupling] if br.rate_source == "geodynamic" else np.ones(len(R["units"])))[us][:, None]
            mfd_mom = np.nansum(fu * r_u * R["mean_m0"][us, i_f, ib, im])       # moment the MFD assigns to the modelled bins (x rigidity factor)
            mom_ratio[i] = (pr * m0p).sum() / mfd_mom
            cen_ratio[i] = (fu * r_u * bins.m0_centre_nm.to_numpy()).sum() / mfd_mom
            if br.rate_source == "geodynamic":
                lev = br.coupling
                fseg = D["rig_factor_unit"][lev][1:]                           # budget of each segment with the pooled rigidity models
                target = np.array([R["xi"][g + 1 if br.segmentation == "segmented" else 0, i_f, ib, im] * bud.loc[n, f"budget_{lev}_nm_yr"] * fseg[g] for g, n in enumerate(names)])
                seg_ratio[i] = (pr[:, None] * m0p[:, None] * F).sum(0) / target
        rep.check(f"{zone}: rate is conserved in every bin on every branch", "Position rates summed over a bin must equal the bin rate of the branch; nothing lost, nothing added",
                  bool(err < 1e-9), f"{len(t)} branches x {nb} bins; largest relative error {err:.1e}")
        rep.check(f"{zone}: moment released by the scenarios against the moment the MFD assigns", "Scenarios carry the moment of their bin centre (and total_M0 differs for some, issue H2), "
                  "while the MFD assigns the mean moment of each bin; the ratio shows the combined effect", None,
                  f"weighted mean {np.average(mom_ratio, weights=w):.3f}; 5% to 95% of branches {wq(mom_ratio, w, .05):.3f} to {wq(mom_ratio, w, .95):.3f}; "
                  f"from bin centres alone {np.average(cen_ratio, weights=w):.3f} ({wq(cen_ratio, w, .05):.3f} to {wq(cen_ratio, w, .95):.3f})")
        for name, m in (("unconstrained", geo & unc), ("segmented", geo & ~unc)):
            sr = seg_ratio[m]; ww = w[m]
            rep.check(f"{zone}: moment released in each segment against its budget, geodynamic {name}", "Spatial balance: over many cycles the slip released on a segment should match the slip deficit accumulating there (ratio 1)", None,
                      "; ".join(f"{n}: {np.average(sr[:, g], weights=ww):.2f} ({sr[:, g].min():.2f} to {sr[:, g].max():.2f})" for g, n in enumerate(names)))
        # events by segment: seismic unconstrained against the catalogue
        m = ~geo & unc
        pm = rates.mean_position_rates(D, m)
        k0 = int(np.searchsorted(bins.mw_lower.to_numpy(), z["seismic_anchor"]["magnitude"] - 1e-3))
        got = (pm[:, None] * F)[kb >= k0].sum(0); obs = np.array([((cat.segment == g + 1) & (cat.mw >= z["seismic_anchor"]["magnitude"]) & (cat.year >= z["seismic_anchor"]["start_year"])).sum() for g in range(n_seg)])
        rep.check(f"{zone}: events by segment, seismic unconstrained branches against the catalogue", "The seismic shares follow the observed events; modelled share of events above the anchor magnitude per segment", None,
                  "; ".join(f"{n}: {100 * a:.1f}% (observed {100 * b:.1f}%)" for n, a, b in zip(names, got / got.sum(), obs / obs.sum())))
        # figure: mean rate along the trench, with a map of how the arc is unrolled
        PMs = rates.mean_position_rates(D, ~geo); PMg = rates.mean_position_rates(D, geo)
        trench = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"trench_{zone}.csv"))
        plotting.rate_along_trench(os.path.join(rep.dir, f"rate_along_trench_{zone}.png"), f"{zone.capitalize()}: weighted mean rate along the trench", tri, trench,
                                   tri.along_trench_km.reindex(pos.tri.to_numpy()).to_numpy(), kb, bins, names,
                                   [(PMs, "--", "black", "seismic branches"), (PMg, "-", "black", "geodynamic branches")])
        pd.DataFrame({"branch_id": t.branch_id, "weight": w, "moment_scenarios_over_mfd": mom_ratio, "moment_bin_centres_over_mfd": cen_ratio,
                      **{f"segment_moment_over_budget_{n}": seg_ratio[:, g] for g, n in enumerate(names)}}).to_csv(os.path.join(rep.dir, f"moment_balance_{zone}.csv"), index=False, float_format="%.5g")
    rep.note("Modelled against observed numbers of large events are in the report of step B4 (they need only the bin rates).")
    rep.note("moment_balance_<zone>.csv holds the ratios for every branch; rate_along_trench_<zone>.png shows where the rate falls.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda", "makran"]) else 1)
