"""Steps B5 and B6. How the rate of a magnitude bin is shared among rupture positions.

A position is rigidity model + scaling law + magnitude bin + initiation triangle (decision 19).
Shares are normalised inside each cell = rigidity model + scaling law + magnitude bin, because the
bin rate is first divided between rigidity models and scaling laws by their pooled weights.

First run (decision 9): rupture footprints are not available. The loading (or seismicity density)
under a rupture is taken as the sum over segments of the segment value x the rupture's segment
fraction. Rupture area and released slip are the same for every position of a cell at this level
of approximation, so they cancel from the shares.

Unconstrained branches (step B5): one share per position for each of four variants
    geodynamic low / preferred / high   weight = sum_g fraction_g x loading_g
    seismic                             weight = sum_g fraction_g x (events above the anchor per km2 in g)
Segmented branches (step B6): one share per position and segment
    weight on segment g = fraction_g x penalty, with penalty = crossing_penalty for positions that
    touch more than one segment, else 1. The segment's own loading or density is the same for all
    its positions and cancels, so in the first run this share is the same for all four variants.
"""
import os
import numpy as np
import pandas as pd
from . import DATA, PRODUCTS, footprints

VARIANTS = ["seismic", "geodynamic_low", "geodynamic_preferred", "geodynamic_high"]


def position_table(scen, n_seg):
    """One row per position: its cell, initiation triangle, expected segment fractions (mean over
    its realisations, weighted by their split) and expected moment."""
    seg = [f"seg_{i + 1}" for i in range(n_seg)]
    d = pd.DataFrame({k: scen[k] for k in ["position", "rig", "scal", "bin", "tri"]})
    first = d.drop_duplicates("position").set_index("position").sort_index()
    w = scen["split"]
    for c in seg:
        first[c] = np.bincount(scen["position"], weights=w * scen[c], minlength=len(first))
    first["m0_nm"] = np.bincount(scen["position"], weights=w * scen["total_M0"], minlength=len(first))
    first["n_scenarios"] = np.bincount(scen["position"], minlength=len(first))
    first["cell"] = first.groupby(["rig", "scal", "bin"], sort=False).ngroup()
    return first.reset_index()


def _normalise(weight, cell):
    tot = np.bincount(cell, weights=weight)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(tot[cell] > 0, weight / tot[cell], 0.0)


def unconstrained_shares(pos, seg_values, n_seg):
    """seg_values: dict variant -> array of one value per segment (loading, or seismicity density).
    Returns array variant x position."""
    F = pos[[f"seg_{i + 1}" for i in range(n_seg)]].to_numpy()
    return np.array([_normalise(F @ np.asarray(seg_values[v], float), pos.cell.to_numpy()) for v in VARIANTS])


def segmented_shares(pos, n_seg, crossing_penalty=1.0):
    """Returns array position x segment: the share of segment g's bin rate (inside the position's
    cell) that goes to the position."""
    F = pos[[f"seg_{i + 1}" for i in range(n_seg)]].to_numpy()
    pen = np.where((F > 0).sum(axis=1) > 1, crossing_penalty, 1.0)
    return np.column_stack([_normalise(F[:, g] * pen, pos.cell.to_numpy()) for g in range(n_seg)])


# ---------------------------------------------------------------- with rupture footprints
# Options `estimated` and `from_file`. The weight of a position is the area-weighted mean, over the triangles of its
# footprint, of the loading (geodynamic) or of the seismicity density (seismic) on those triangles. Loading that varies
# from triangle to triangle is therefore followed inside a segment. On the segmented branches the same sum is taken
# over the triangles of each segment separately, so these shares differ between the variants.

def triangle_values(load, par, bud, names):
    """dict variant -> one value per triangle, in the order of the loading table of step B2.
    Seismic: events above the anchor per km2 of the triangle's segment (as in the approximation)."""
    seg = load.segment.to_numpy() - 1
    dens = (par.set_index("unit").loc[names, "n_ge_anchor"].to_numpy() / bud.set_index("unit").loc[names, "area_km2"].to_numpy())[seg]
    return {"seismic": dens, **{f"geodynamic_{lv}": load[f"loading_{lv}_mm_yr"].to_numpy() for lv in ("low", "preferred", "high")}}


def get_footprints(cfg, zone, tri, pos, bins, scen=None):
    """Footprint table for the option in the configuration (None for the approximation)."""
    opt = cfg["tree"]["options"]; mode = opt["footprints"]
    if mode == "segment_fraction_approximation":
        return None
    if mode == "estimated":
        mw = bins.mw_centre.to_numpy()[pos.bin.to_numpy()]; area = np.zeros(len(pos))
        for k, coef in opt["footprint_area"].items():
            m = (pos.scal == int(k)).to_numpy(); area[m] = footprints.law_area(mw[m], coef)
        return footprints.estimate(tri, pos.tri.to_numpy(), area)
    return footprints.from_file(os.path.join(DATA, "footprints"), cfg["zones"][zone]["code"], scen, tri, len(pos))


def footprint_path(zone):
    return os.path.join(PRODUCTS, "B5_shares", f"footprints_{zone}.npz")


def unconstrained_shares_fp(pos, fp, tri_values):
    """As unconstrained_shares, with the weight taken under the footprint. Returns variant x position."""
    return np.array([_normalise(footprints.row_sums(fp, tri_values[v]), pos.cell.to_numpy()) for v in VARIANTS])


def segmented_shares_fp(pos, fp, tri_values, tri_segment, n_seg, crossing_penalty=1.0):
    """As segmented_shares, with the weight on segment g = sum over the footprint's triangles in g of
    area fraction x value. Returns (variant x position x segment, footprint segment fractions)."""
    seg = np.asarray(tri_segment); cell = pos.cell.to_numpy()
    Ff = np.column_stack([footprints.row_sums(fp, (seg == g + 1).astype(float)) for g in range(n_seg)])
    pen = np.where((Ff > 0).sum(axis=1) > 1, crossing_penalty, 1.0)
    X = np.array([np.column_stack([_normalise(footprints.row_sums(fp, np.where(seg == g + 1, tri_values[v], 0.0)) * pen, cell) for g in range(n_seg)]) for v in VARIANTS])
    return X, Ff
