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
import numpy as np
import pandas as pd

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
