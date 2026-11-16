"""Step B6. Shares for the segmented branches, with boundary-crossing ruptures.

Run from ps_scenario_rates after step B5:   python scripts/run_B6_segmented.py [zone ...]     (default: sunda)
Writes products/B6_segmented/ and checks/B6_segmented/.

Product per zone:
  shares_segmented_<zone>.npz   share_segmented (position x segment): share of segment g's bin rate, inside the
                                position's cell, that goes to the position
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, config, footprints, shares
from ptha_rates.report import CheckReport


def main(zones):
    cfg = config.load(); pen = cfg["tree"]["options"]["crossing_penalty"]
    out = os.path.join(PRODUCTS, "B6_segmented"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B6_segmented", "shares for the segmented branches")
    for zone in zones:
        z = cfg["zones"][zone]; n_seg = len(z["segments"]); names = [s["name"] for s in z["segments"]]
        pos = pd.read_csv(os.path.join(PRODUCTS, "B5_shares", f"positions_{zone}.csv"))
        bins = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"bins_{zone}.csv"))
        mode = cfg["tree"]["options"]["footprints"]
        if mode == "segment_fraction_approximation":
            F = pos[[f"seg_{i + 1}" for i in range(n_seg)]].to_numpy()
            make = lambda p: shares.segmented_shares(pos, n_seg, p)
            X = make(pen); layers = [X]
        else:                                               # with rupture footprints: one set of shares per rate variant
            tri = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"triangles_{zone}.csv"))
            load = pd.read_csv(os.path.join(PRODUCTS, "B2_budget", f"loading_{zone}.csv"))
            par = pd.read_csv(os.path.join(PRODUCTS, "B3_catalogue", f"parameters_{zone}.csv")); bud = pd.read_csv(os.path.join(PRODUCTS, "B2_budget", f"budget_{zone}.csv"))
            fp = footprints.load(shares.footprint_path(zone)); tv = shares.triangle_values(load, par, bud, names)
            make = lambda p: shares.segmented_shares_fp(pos, fp, tv, tri.segment.to_numpy(), n_seg, p)[0]
            X, F = shares.segmented_shares_fp(pos, fp, tv, tri.segment.to_numpy(), n_seg, pen); layers = list(X)
        cell = pos.cell.to_numpy(); ncell = cell.max() + 1; X0 = layers[0]
        tots = [np.column_stack([np.bincount(cell, weights=L[:, g], minlength=ncell) for g in range(n_seg)]) for L in layers]; tot = tots[0]
        rep.check(f"{zone}: shares on a segment sum to 1 or 0 in every cell", "1 where some rupture of the cell touches the segment, 0 where none does",
                  bool(all(np.all(np.isclose(t_, 1) | np.isclose(t_, 0)) for t_ in tots)), f"{int(np.isclose(tot, 1).sum())} cell-segment pairs with ruptures, {int(np.isclose(tot, 0).sum())} without")
        cells = pos.drop_duplicates("cell").set_index("cell").sort_index()
        empty = [(names[g], bins.mw_centre[cells.bin[c]], int(cells.scal[c])) for c, g in zip(*np.where(np.isclose(tot, 0)))]
        rep.check(f"{zone}: cells in which a segment has no rupture", "The segment's rate in that bin and scaling law has no scenario and is lost (reported in step B8)", None,
                  "; ".join(f"{n} Mw {m:g} law {s:02d}" for n, m, s in empty[:12]) + (f"; and {len(empty) - 12} more" if len(empty) > 12 else "") or "none")
        single = (F > 0).sum(1) == 1
        g1 = F.argmax(1)
        other = np.max([L[single].sum(axis=1) - L[single, g1[single]] for L in layers], axis=0)
        rep.check(f"{zone}: positions inside one segment take rate from that segment only", "Their share on every other segment must be zero",
                  bool((other == 0).all()), f"{int(single.sum()):,} positions inside one segment, {int((~single).sum()):,} crossing a boundary")
        Xp = make(0.5); Xp = Xp if Xp.ndim == 2 else Xp[0]
        tp = np.column_stack([np.bincount(cell, weights=Xp[:, g], minlength=ncell) for g in range(n_seg)])
        rep.check(f"{zone}: a crossing penalty keeps every segment's rate", "Shares are normalised after the penalty, so a penalty moves rate from crossing to non-crossing ruptures without losing any",
                  bool(np.all(np.isclose(tp, 1) | np.isclose(tp, 0)) and (Xp[~single].sum() < X0[~single].sum())), f"tested with penalty 0.5; penalty in use: {pen}")
        cross_share = np.array([[(X0[(pos.bin == k).to_numpy() & ~single, g]).sum() / max(X0[(pos.bin == k).to_numpy(), g].sum(), 1e-300) for g in range(n_seg)] for k in (0, 10, 16, 30)])
        rep.check(f"{zone}: part of a segment's rate taken by boundary-crossing ruptures", "Rises with magnitude; at the largest magnitudes most ruptures cross a boundary", None,
                  "; ".join(f"Mw {bins.mw_centre[k]:g}: " + ", ".join(f"{n} {100 * c:.0f}%" for n, c in zip(names, row)) for k, row in zip((0, 10, 16, 30), cross_share)))
        np.savez_compressed(os.path.join(out, f"shares_segmented_{zone}.npz"), share_segmented=X, segments=np.array(names), position=pos.position.to_numpy())
    if cfg["tree"]["options"]["footprints"] == "segment_fraction_approximation":
        rep.note("In the first run the loading and the seismicity density are constant inside a segment, so they cancel from these shares: "
                 "one set of segmented shares serves the seismic branches and all three coupling levels.")
    else:
        rep.note(f"footprints = {cfg['tree']['options']['footprints']}: share_segmented is variant x position x segment (variants: " + ", ".join(shares.VARIANTS) + ").")
    rep.note(f"crossing_penalty = {pen} (1 = no penalty).")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda"]) else 1)
