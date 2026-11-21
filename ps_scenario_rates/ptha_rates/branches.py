"""Step B0. Logic tree levels, their weights and the branch table."""
import itertools
import numpy as np
import pandas as pd

SHORT = {"seismic": "S", "geodynamic": "G", "unconstrained": "U", "segmented": "X",
         "truncated_gr": "TG", "tapered_pareto": "TP",
         "observed": "Mobs", "area_limit": "Marea", "mean": "Mmean",
         "low": "Clow", "preferred": "Cpref", "high": "Chigh", "none": "C-"}


def b_values(tree):
    """The b-values of the logic tree: the listed values, or n values from min to max."""
    b = tree["levels"]["b_value"]
    if "values" in b:
        return np.asarray(b["values"], float)
    return np.linspace(b["min"], b["max"], b["n"])


def b_weights(b_hat, tree):
    """Weights of the b-values. A list in the configuration is used as given (normalised to sum to 1;
    decision 24: centred grid with symmetric weights, so the weighted mean is the catalogue b).
    'linear' (version 1.0): highest at the catalogue estimate b_hat and falling linearly with distance
    from it; L is the distance from b_hat to the farther end of the range plus one grid step, so the
    end values keep a small weight."""
    b = b_values(tree)
    rule = tree["levels"]["b_value"]["weights"]
    if isinstance(rule, (list, tuple)):
        w = np.asarray(rule, float)
        if len(w) != len(b):
            raise ValueError("b_value: one weight per value is needed")
        return w / w.sum()
    if rule != "linear":
        raise ValueError("b_value weights must be a list or 'linear'")
    step = b[1] - b[0]
    L = max(abs(b_hat - b[0]), abs(b[-1] - b_hat)) + step
    w = 1 - np.abs(b - b_hat) / L
    return w / w.sum()


def coupling_levels(preferred, tree):
    """Low, preferred and high coupling, kept within floor and ceiling.
    rule factor (decision 25): preferred x the factor of each level; rule offset (version 1.0): preferred -/+ offset."""
    c = tree["levels"]["coupling"]
    clip = lambda x: float(min(max(x, c["floor"]), c["ceiling"]))
    if c.get("rule", "offset") == "factor":
        return {lev: clip(preferred * f) for lev, f in c["factors"].items()}
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
    t = pd.DataFrame(rows)
    return expand_event_tree_levels(t, tree)


def expand_event_tree_levels(t, tree):
    """Event-tree items made logic-tree levels (event_tree: as_level, decision 28). Each branch is repeated for every
    scaling law (weights of event_tree: scaling_law) and/or every scaling sample (equal weights, samples 00 to
    n_scaling_samples - 1). Without as_level the columns hold 'pooled'."""
    ev = tree["event_tree"]; lev = ev.get("as_level", [])
    t["scaling_law"] = "pooled"; t["scaling_sample"] = "pooled"
    if "scaling_law" in lev:
        laws = {k: w for k, w in ev["scaling_law"].items() if w > 0}
        t = pd.concat([t.assign(scaling_law=k, weight=t.weight * w, branch_id=t.branch_id + f"-L{k}") for k, w in laws.items()], ignore_index=True)
    if "scaling_sample" in lev:
        n = int(ev.get("n_scaling_samples", 10))
        t = pd.concat([t.assign(scaling_sample=f"{k:02d}", weight=t.weight / n, branch_id=t.branch_id + f"-S{k:02d}") for k in range(n)], ignore_index=True)
    return t
