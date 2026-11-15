"""Write docs/results.md: the numbers and figures of the present run, read from products/ and checks/.

Run from ps_scenario_rates after the steps (run_all.py does it):   python scripts/make_results.py [zone ...]
Nothing in results.md is typed by hand; do not edit it. Change the configuration, run, and it is rewritten.
"""
import glob
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import CHECKS, PRODUCTS, ROOT, __version__, binrates, branches, config, rates


def table(header, rows):
    return ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)] + ["| " + " | ".join(str(x) for x in r) + " |" for r in rows] + [""]


def period(rate):
    return f"{rate:.2f} per year" if rate >= 1 else f"once per {1 / rate:,.0f} years"


def zone_section(zone, cfg):
    z = cfg["zones"][zone]; tree = cfg["tree"]; lv = tree["levels"]; names = [s["name"] for s in z["segments"]]; n_seg = len(names)
    P = lambda *a: os.path.join(PRODUCTS, *a)
    bins = pd.read_csv(P("B1_inputs", f"bins_{zone}.csv")); bud = pd.read_csv(P("B2_budget", f"budget_{zone}.csv")).set_index("unit")
    par = pd.read_csv(P("B3_catalogue", f"parameters_{zone}.csv")).set_index("unit"); big = pd.read_csv(P("B4_bin_rates", f"large_events_{zone}.csv"))
    D = rates.load(zone, cfg); t = D["table"]; w = t.weight.to_numpy(); pos = D["pos"]
    S = np.load(P("B5_shares", f"shares_{zone}.npz")); F = np.column_stack([S[f"seg_{i + 1}"] for i in range(n_seg)])
    geo = (t.rate_source == "geodynamic").to_numpy()
    L = [f"## {zone.capitalize()}", "",
         f"{int(bins.n_scenarios.sum()):,} scenarios in {len(pos):,} rupture positions; {len(bins)} magnitude bins, Mw {bins.mw_lower.iloc[0]:.2f} to {bins.mw_upper.iloc[-1]:.2f}; "
         f"{len(t)} branches ({int(geo.sum())} geodynamic, {int((~geo).sum())} seismic).", ""]
    idx = os.path.join(CHECKS, "INDEX.md")
    if os.path.exists(idx):
        txt = open(idx, encoding="utf-8").read()
        L += [f"Checks of this run: {txt.count('| PASS |')} passed, {txt.count('| FAIL |')} failed, {txt.count('| INFO |')} reported for information (`checks/INDEX.md`).", ""]

    L += ["### Configuration used", "", "Weights of the logic tree (`config/logic_tree.yaml`):", ""]
    bv, bw = branches.b_values(tree), branches.b_weights(z["whole_zone"]["b"], tree)
    rows = [[k.replace("_", " "), ", ".join(f"{a.replace('_', ' ')} {v:g}" for a, v in lv[k].items())] for k in ("rate_source", "segmentation", "mfd_form", "mmax")]
    rows.insert(3, ["b-value", ", ".join(f"{b:.2f}: {x:.3f}" for b, x in zip(bv, bw))])
    rows.append(["coupling (geodynamic)", ", ".join(f"{a} {v:.3f}" for a, v in lv["coupling"]["weights"].items())])
    L += table(["Level", "Alternatives and weights"], rows)
    L += ["Zone and segment values (`config/zones.yaml`; Mmax observed from the catalogue):", ""]
    L += table(["Unit", "b", "Coupling low / preferred / high", "Convergence (mm/yr)", "Mmax observed / area limit / mean", f"Events of Mw >= {par.anchor_mw.iloc[0]:g} since {int(par.anchor_start_year.iloc[0])}"],
               [[u.replace("_", " "), f"{par.b_catalogue[u]:.3f}", f"{bud.coupling_low[u]:.3g} / {bud.coupling_preferred[u]:.3g} / {bud.coupling_high[u]:.3g}", f"{bud.convergence_mm_yr[u]:.0f}",
                 f"{par.mmax_observed[u]:.1f} / {par.mmax_area_limit[u]:.1f} / {par.mmax_mean[u]:.2f}", int(par.n_ge_anchor[u])] for u in ["whole_zone"] + names])

    L += ["### Moment budget", ""]
    tot = bud.budget_preferred_nm_yr["whole_zone"]
    L += table(["Unit", "Low (N m/yr)", "Preferred", "High", "Share of the preferred budget"],
               [[u.replace("_", " "), f"{bud.budget_low_nm_yr[u]:.3g}", f"{bud.budget_preferred_nm_yr[u]:.3g}", f"{bud.budget_high_nm_yr[u]:.3g}", f"{100 * bud.budget_preferred_nm_yr[u] / tot:.0f}%"] for u in ["whole_zone"] + names])

    yrs = int(big.years.iloc[0]); thr = sorted(big.mw_threshold.unique())
    L += [f"### Large events in {yrs} years", "", "Weighted mean over the branches of each group, with the 5% to 95% range of branches in brackets. Reported only; no weight is changed.", ""]
    obs = big.drop_duplicates("mw_threshold").set_index("mw_threshold").observed
    rows = [["Observed"] + [int(obs[x]) for x in thr]]
    for g in big.group.unique():
        b = big[big.group == g].set_index("mw_threshold")
        rows.append([g[0].upper() + g[1:]] + [f"{b.modelled_mean[x]:.1f} ({b.modelled_p05[x]:.1f} to {b.modelled_p95[x]:.1f})" for x in thr])
    L += table([""] + [f"Mw >= {x:g}" for x in thr], rows)
    L += [f"![MFDs of all branches](figures/mfd_{zone}.png)", ""]

    Z = binrates.zone_bin_rates(D["R"], t); mean_z = (Z * w[:, None]).sum(0); lo = bins.mw_lower.to_numpy()
    L += ["### Mean rate of the zone", "", "Weighted mean over all branches."]
    L += [""] + table(["Magnitude", "Rate"], [[f"Mw >= {lo[k]:.2f}", period(mean_z[k:].sum())] for k in sorted({0} | {int(np.searchsorted(lo, x)) for x in (7.0, 7.5, 8.0, 8.5, 9.0, 9.25) if x < lo[-1]})])
    L += ["Total rate of modelled events: " + "; ".join(f"{lab} branches {(Z[m].sum(1) * w[m]).sum() / w[m].sum():.2f} per year" for lab, m in (("seismic", ~geo), ("geodynamic", geo))) + ".", ""]

    PMs, PMg = rates.mean_position_rates(D, ~geo), rates.mean_position_rates(D, geo)
    L += ["### Where the rate falls", "", "Events per year by segment, weighted mean. A rupture's rate is divided between segments by its segment fractions.", ""]
    L += table(["Segment", "Seismic branches", "Geodynamic branches", f"Observed (Mw >= {par.anchor_mw.iloc[0]:g})"],
               [[n, f"{(PMs * F[:, g]).sum():.3f}", f"{(PMg * F[:, g]).sum():.3f}", f"{par.rate_anchor_per_yr[n]:.3f}"] for g, n in enumerate(names)])
    L += [f"![Rate along the trench](figures/rate_along_trench_{zone}.png)", ""]

    mb = os.path.join(CHECKS, "B8_end_to_end", f"moment_balance_{zone}.csv")
    if os.path.exists(mb):
        m = pd.read_csv(mb); ww = m.weight.to_numpy(); g_ = m[f"segment_moment_over_budget_{names[0]}"].notna().to_numpy()
        unc = (t.segmentation == "unconstrained").to_numpy()
        L += ["### Consistency", "", f"Moment released by the scenarios against the moment the MFD assigns: {np.average(m.moment_scenarios_over_mfd, weights=ww):.3f} "
              f"(a scenario carries the moment of its bin centre: {np.average(m.moment_bin_centres_over_mfd, weights=ww):.3f} from that alone).", "",
              "Moment released in each segment against its budget, geodynamic branches (1 = balanced):", ""]
        L += table(["Segment", "Unconstrained", "Segmented"], [[n] + [f"{np.average(m[f'segment_moment_over_budget_{n}'][g_ & k], weights=ww[g_ & k]):.2f}" for k in (unc, ~unc)] for n in names])

    L += ["### Other figures", "", f"![Logic tree](figures/logic_tree_{zone}.png)", "", f"![Weights](figures/weights_panel_{zone}.png)", "", f"![Catalogue](figures/catalogue_fit_{zone}.png)", ""]
    return L


def sensitivity_section():
    L = []
    for rep in sorted(glob.glob(os.path.join(PRODUCTS, "sensitivity", "*", "report.md"))):
        tag = os.path.basename(os.path.dirname(rep)); body = open(rep, encoding="utf-8").read().strip().split("\n")
        L += ["### " + body[0].lstrip("# ")] + [("##" + x if x.startswith("## ") else x) for x in body[1:]] + [""]
        for f in ("rate_along_trench", "mfd", "mfd_segments"):
            if os.path.exists(os.path.join(os.path.dirname(rep), f + ".png")):
                L += [f"![{f.replace('_', ' ')}](figures/sensitivity_{tag}_{f}.png)", ""]
    return (["## Sensitivity tests", "", "In each test \"now\" or \"present\" is the configuration and \"test\" the changed input."
             " The lines about file names refer to `products/sensitivity/`.", ""] + L) if L else []


def main(zones):
    cfg = config.load()
    L = ["# Results", "", f"Version {__version__} of the rates. Written by `scripts/make_results.py` from the products of the run; not edited by hand.",
         "What the quantities mean, the limitations of each step and the decisions behind the inputs are in `guide.md`.", ""]
    for zone in zones:
        L += zone_section(zone, cfg)
    L += sensitivity_section()
    out = os.path.join(ROOT, "docs", "results.md")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L).rstrip("\n") + "\n")
    print("written", os.path.relpath(out, ROOT))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["sunda"]))
