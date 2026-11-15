"""Step B0. Logic tree levels, their weights and the branch table."""
import itertools
import numpy as np
import pandas as pd

SHORT = {"seismic": "S", "geodynamic": "G", "unconstrained": "U", "segmented": "X",
         "truncated_gr": "TG", "tapered_pareto": "TP",
         "observed": "Mobs", "area_limit": "Marea", "mean": "Mmean",
         "low": "Clow", "preferred": "Cpref", "high": "Chigh", "none": "C-"}


def b_values(tree):
    b = tree["levels"]["b_value"]
    return np.linspace(b["min"], b["max"], b["n"])


def b_weights(b_hat, tree):
    """Weights of the b-values: highest at the catalogue estimate b_hat and falling
    linearly with distance from it. L is the distance from b_hat to the farther end
    of the range plus one grid step, so the end values keep a small weight."""
    b = b_values(tree)
    if tree["levels"]["b_value"]["weights"] != "linear":
        raise ValueError("only 'linear' b-value weights are implemented")
    step = b[1] - b[0]
    L = max(abs(b_hat - b[0]), abs(b[-1] - b_hat)) + step
    w = 1 - np.abs(b - b_hat) / L
    return w / w.sum()


def coupling_levels(preferred, tree):
    """Low, preferred and high coupling: preferred -/+ offset, kept within floor and ceiling."""
    c = tree["levels"]["coupling"]
    clip = lambda x: float(min(max(x, c["floor"]), c["ceiling"]))
    return {"low": clip(preferred - c["offset"]), "preferred": clip(preferred),
            "high": clip(preferred + c["offset"])}


def branch_table(cfg, zone):
    """One row per branch: alternative at each level, weight and spatial variant.

    The b-value weight of a branch is taken from the whole-zone catalogue estimate,
    for unconstrained and segmented branches alike: a branch uses one b-value in
    every segment, so it can carry only one b-value weight.
    """
    tree, z = cfg["tree"], cfg["zones"][zone]
    lv = tree["levels"]
    bv, bw = b_values(tree), b_weights(z["whole_zone"]["b"], tree)
    rows = []
    for src, seg, form, ib, mm in itertools.product(lv["rate_source"], lv["segmentation"],
                                                     lv["mfd_form"], range(len(bv)), lv["mmax"]):
        coup = lv["coupling"]["weights"] if src == "geodynamic" else {"none": 1.0}
        for cl, wc in coup.items():
            w = lv["rate_source"][src] * lv["segmentation"][seg] * lv["mfd_form"][form] * bw[ib] * lv["mmax"][mm] * wc
            rows.append(dict(
                branch_id=f"{SHORT[src]}-{SHORT[seg]}-{SHORT[form]}-b{ib + 1:02d}-{SHORT[mm]}-{SHORT[cl]}",
                rate_source=src, segmentation=seg, mfd_form=form, b_index=ib + 1, b_value=round(float(bv[ib]), 6),
                mmax=mm, coupling=cl, weight=w, spatial_variant=f"{src}|{seg}|{cl}"))
    return pd.DataFrame(rows)
