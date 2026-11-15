"""Step B5. Share of each bin's rate among rupture positions, unconstrained branches.

Run from ps_scenario_rates after steps B1 to B3:   python scripts/run_B5_shares.py [zone ...]     (default: sunda)
Writes products/B5_shares/ and checks/B5_shares/.

Products per zone:
  shares_<zone>.npz      the position table (one value per position) and share_unconstrained (variant x position);
                         variants: seismic, geodynamic_low, geodynamic_preferred, geodynamic_high
  positions_<zone>.csv   the position table as text
  segment_values_<zone>.csv   loading and seismicity density per segment used for the weights
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, config, shares
from ptha_rates.report import CheckReport


def main(zones):
    cfg = config.load()
    out = os.path.join(PRODUCTS, "B5_shares"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B5_shares", "shares among rupture positions, unconstrained branches")
    for zone in zones:
        z = cfg["zones"][zone]; n_seg = len(z["segments"]); names = [s["name"] for s in z["segments"]]
        scen = np.load(os.path.join(PRODUCTS, "B1_inputs", f"scenarios_{zone}.npz"))
        bud = pd.read_csv(os.path.join(PRODUCTS, "B2_budget", f"budget_{zone}.csv")).set_index("unit").loc[names]
        par = pd.read_csv(os.path.join(PRODUCTS, "B3_catalogue", f"parameters_{zone}.csv")).set_index("unit").loc[names]
        bins = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"bins_{zone}.csv"))
        pos = shares.position_table(scen, n_seg)
        vals = {"seismic": (par.n_ge_anchor / bud.area_km2).to_numpy(),
                "geodynamic_low": bud.loading_low_mm_yr.to_numpy(), "geodynamic_preferred": bud.loading_preferred_mm_yr.to_numpy(),
                "geodynamic_high": bud.loading_high_mm_yr.to_numpy()}
        U = shares.unconstrained_shares(pos, vals, n_seg)
        cell = pos.cell.to_numpy(); ncell = cell.max() + 1
        ev = cfg["tree"]["event_tree"]
        n_expected = len(bins) * sum(w > 0 for w in ev["scaling_law"].values()) * sum(w > 0 for w in ev["rigidity_model"].values())
        rep.check(f"{zone}: positions", "One row per rigidity model + scaling law + magnitude bin + initiation triangle", len(pos) == scen["position"].max() + 1,
                  f"{len(pos):,} positions from {len(scen['position']):,} scenarios")
        rep.check(f"{zone}: every bin has positions for every scaling law", "Otherwise part of a bin's rate would have no scenario to go to",
                  ncell == n_expected, f"{ncell} cells of {n_expected}; positions per cell {np.bincount(cell).min():,} to {np.bincount(cell).max():,}")
        F = pos[[f"seg_{i + 1}" for i in range(n_seg)]].to_numpy()
        rep.check(f"{zone}: expected segment fractions of a position sum to 1", "They are averages of the scenario fractions", bool(np.allclose(F.sum(1), 1, atol=1e-3)),
                  f"{F.sum(1).min():.5f} to {F.sum(1).max():.5f}")
        tot = np.array([np.bincount(cell, weights=u, minlength=ncell) for u in U])
        rep.check(f"{zone}: shares sum to 1 in every cell", "The whole of a bin's rate is given to its positions", bool(np.allclose(tot, 1)), f"largest deviation {np.abs(tot - 1).max():.1e}")
        rep.check(f"{zone}: shares are not negative", "A share is a fraction", bool((U >= 0).all()), f"smallest {U.min():.2e}, largest {U.max():.2e}")
        # uniform loading must give equal shares within a cell
        Uu = shares.unconstrained_shares(pos, {v: np.ones(n_seg) for v in shares.VARIANTS}, n_seg)[0]
        rep.check(f"{zone}: equal loading gives equal shares", "With the same loading everywhere a position's share must be 1 / (positions in its cell)",
                  bool(np.allclose(Uu, 1 / np.bincount(cell)[cell], rtol=1e-3)),
                  f"largest relative deviation {np.abs(Uu * np.bincount(cell)[cell] - 1).max():.1e} (segment fractions sum to 1 within 4e-5)")
        # where does the rate go? share of a bin's rate by dominant segment, compared with the target
        dom = F.argmax(1)
        for v, target, label in (("geodynamic_preferred", (bud.area_km2 * bud.loading_preferred_mm_yr), "area x loading"), ("seismic", par.n_ge_anchor, "events above the anchor")):
            u = U[shares.VARIANTS.index(v)]
            k0 = (pos.bin == 0).to_numpy()
            got = np.array([(u[k0] * F[k0, g]).sum() for g in range(n_seg)]); got = got / got.sum()
            tg = (target / target.sum()).to_numpy()
            rep.check(f"{zone}: where the rate of the smallest bin goes, {v}", f"Share of the bin's rate falling in each segment, against the segment's share of {label}", None,
                      "; ".join(f"{n}: {100 * a:.1f}% (target {100 * b:.1f}%)" for n, a, b in zip(names, got, tg)))
        rep.check(f"{zone}: positions per segment in the smallest bin", "Positions are initiation triangles; an uneven density of positions shifts rate between segments", None,
                  "; ".join(f"{n}: {int(((dom == g) & (pos.bin == 0)).sum()):,} positions, {bud.area_km2.iloc[g] / max(((dom == g) & (pos.bin == 0)).sum(), 1):.0f} km2 each" for g, n in enumerate(names)))
        np.savez_compressed(os.path.join(out, f"shares_{zone}.npz"), share_unconstrained=U, variants=np.array(shares.VARIANTS),
                            **{c: pos[c].to_numpy() for c in pos.columns})
        pos.round(6).to_csv(os.path.join(out, f"positions_{zone}.csv"), index=False)
        pd.DataFrame(vals, index=names).rename_axis("segment").to_csv(os.path.join(out, f"segment_values_{zone}.csv"), float_format="%.6g")
    rep.note("First-run approximation (decision 9): loading and seismicity density under a rupture are taken from its segment fractions. "
             "Inside a segment every position of a cell therefore gets the same weight.")
    rep.note("Seismicity density of a segment = events at or above the anchor magnitude since the anchor year, per km2 of the segment.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda"]) else 1)
