"""Step B2. Loading on every triangle and the moment budget of each zone and segment.

Run from ps_scenario_rates after step B1:   python scripts/run_B2_budget.py [zone ...]     (default: sunda)
Writes products/B2_budget/ and checks/B2_budget/.

Products per zone:
  loading_<zone>.csv   per triangle: segment, area, convergence, coupling and loading at the three coupling levels
  budget_<zone>.csv    per segment and whole zone: area, loading and moment budget at the three coupling levels
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, budget, config
from ptha_rates.report import CheckReport


def main(zones):
    cfg = config.load(); tree = cfg["tree"]
    out = os.path.join(PRODUCTS, "B2_budget"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B2_budget", "loading and moment budget")
    for zone in zones:
        z = cfg["zones"][zone]
        tri = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"triangles_{zone}.csv"))
        load = budget.loading_table(tri, z, tree)
        t = budget.budget_table(load, z)
        seg, whole = t[t.unit != "whole_zone"], t[t.unit == "whole_zone"].iloc[0]
        c = tree["levels"]["coupling"]
        cc = load[[f"coupling_{l}" for l in budget.LEVELS]].to_numpy()
        rep.check(f"{zone}: coupling within floor and ceiling", "Decision 3: preferred -/+ 0.2, kept within 0.1 to 1",
                  bool((cc >= c["floor"] - 1e-12).all() and (cc <= c["ceiling"] + 1e-12).all()), f"{cc.min():.2f} to {cc.max():.2f}")
        rep.check(f"{zone}: low <= preferred <= high on every triangle", "The three levels must be ordered",
                  bool((cc[:, 0] <= cc[:, 1]).all() and (cc[:, 1] <= cc[:, 2]).all()), f"{len(load):,} triangles")
        rep.check(f"{zone}: whole-zone budget equals the sum of the segments", "The unconstrained and segmented branches must share one budget",
                  all(np.isclose(whole[f"budget_{l}_nm_yr"], seg[f"budget_{l}_nm_yr"].sum()) for l in budget.LEVELS),
                  "; ".join(f"{l}: {whole[f'budget_{l}_nm_yr']:.3e} N m/yr" for l in budget.LEVELS))
        direct = z["rigidity_pa"] * (load.area_km2 * 1e6 * load.loading_preferred_mm_yr * 1e-3).sum()
        rep.check(f"{zone}: budget recomputed directly from the triangles", "Independent sum over all triangles",
                  bool(np.isclose(direct, whole.budget_preferred_nm_yr)), f"{direct:.4e} N m/yr")
        rep.check(f"{zone}: area in the budget equals the mesh area", "No triangle left out", bool(np.isclose(whole.area_km2, tri.area_km2.sum())),
                  f"{whole.area_km2:,.0f} km2")
        for _, r in seg.iterrows():
            ref = next(s for s in z["segments"] if s["name"] == r.unit).get("ref_cs")
            ours = [r[f"loading_{l}_mm_yr"] for l in budget.LEVELS]
            if ref:
                inside = ref[0] <= ours[1] <= ref[2]
                txt = (f"ours {ours[0]:.1f} / {ours[1]:.1f} / {ours[2]:.1f}; reference {ref[0]:g} / {ref[1]:g} / {ref[2]:g}; "
                       f"preferred is {'inside' if inside else 'OUTSIDE'} the reference range ({ours[1] / ref[1]:.1f} x the reference preferred)")
            else:
                txt = f"ours {ours[0]:.1f} / {ours[1]:.1f} / {ours[2]:.1f}; no reference value"
            rep.check(f"{zone} {r.unit}: loading against Davies et al. (2018), mm/yr", "Coupled convergence, low / preferred / high", None, txt)
        rep.check(f"{zone}: share of the preferred budget by segment", "Shows which segments carry the geodynamic rates", None,
                  "; ".join(f"{r.unit} {100 * r.budget_preferred_nm_yr / whole.budget_preferred_nm_yr:.0f}%" for _, r in seg.iterrows()))
        m9 = 10 ** (1.5 * 9.0 + 9.1)
        rep.check(f"{zone}: size of the whole-zone budget", "For orientation: years needed to accumulate the moment of one Mw 9.0", None,
                  "; ".join(f"{l}: Mw {budget.mw_of(whole[f'budget_{l}_nm_yr']):.2f} per year, one Mw 9.0 per {m9 / whole[f'budget_{l}_nm_yr']:.0f} yr"
                            for l in budget.LEVELS))
        rep.check(f"{zone}: whole-zone effective coupling", "Budget divided by the budget with coupling 1", None,
                  " / ".join(f"{whole[f'coupling_{l}']:.2f}" for l in budget.LEVELS))
        dup = [r.unit for _, r in seg.iterrows() if len({round(r[f'coupling_{l}'], 6) for l in budget.LEVELS}) < 3]
        rep.check(f"{zone}: segments with identical coupling levels", "Their geodynamic branches coincide", None, ", ".join(dup) or "none")
        load.round(6).to_csv(os.path.join(out, f"loading_{zone}.csv"), index=False)
        t.to_csv(os.path.join(out, f"budget_{zone}.csv"), index=False, float_format="%.6g")
    rep.note(f"Rigidity {cfg['zones'][zones[0]]['rigidity_pa'] / 1e9:.0f} GPa throughout: the first run uses the homogeneous rigidity model only (decision 14).")
    rep.note("Coupling and convergence are constant within a segment (values of Table 1). Convergence is used as given (decision 6). No edge correction.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda"]) else 1)
