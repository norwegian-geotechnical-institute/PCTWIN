"""Step B7. Annual rate of every rupture position and scenario on any branch.

    position rate = pooled weight (rigidity model x scaling law)
                    x  bin rate of the branch
                    x  share of the position                     (unconstrained branches)
                  = pooled weight x sum over segments g of
                    (bin rate of segment g  x  share of the position on g)   (segmented branches)
    scenario rate = position rate x split (1 / (scaling samples x slip realisations that exist))
"""
import os
import numpy as np
import pandas as pd
from . import PRODUCTS, binrates, config, shares


def load(zone, cfg=None):
    """Everything needed to compute rates for a zone, read from the products of steps B0 to B6."""
    cfg = cfg or config.load()
    P = PRODUCTS
    R = dict(np.load(os.path.join(P, "B4_bin_rates", f"bin_rates_{zone}.npz"), allow_pickle=False))
    R = {k: (v.tolist() if k in ("units", "variants", "forms", "mmax") else v) for k, v in R.items()}
    S = np.load(os.path.join(P, "B5_shares", f"shares_{zone}.npz"), allow_pickle=False)
    ev = cfg["tree"]["event_tree"]
    pos = pd.DataFrame({k: S[k] for k in ("position", "rig", "scal", "bin", "tri", "cell", "m0_nm", "n_scenarios")})
    pooled = np.array([ev["rigidity_model"][f"{r:02d}"] * ev["scaling_law"][f"{s:02d}"] for r, s in zip(pos.rig, pos.scal)])
    return dict(zone=zone, cfg=cfg, R=R, pos=pos, pooled=pooled, U=S["share_unconstrained"],
                X=np.load(os.path.join(P, "B6_segmented", f"shares_segmented_{zone}.npz"))["share_segmented"],
                table=pd.read_csv(os.path.join(P, "B0_config", f"branch_table_{zone}.csv")))


def position_rates(D, branch):
    """Annual rate of every position on one branch (a row of the branch table)."""
    r = binrates.branch_bin_rates(D["R"], branch)            # unit x bin
    k = D["pos"].bin.to_numpy()
    if branch.segmentation == "unconstrained":
        v = shares.VARIANTS.index(binrates.variant_of(branch.rate_source, branch.coupling))
        return D["pooled"] * r[0][k] * D["U"][v]
    return D["pooled"] * (r[1:][:, k].T * D["X"]).sum(axis=1)


def mean_position_rates(D, mask=None):
    """Weighted mean rate of every position over the branches selected by mask (all by default)."""
    t = D["table"]; out = np.zeros(len(D["pos"])); wsum = 0.0
    for i, br in enumerate(t.itertuples()):
        if mask is None or mask[i]:
            out += br.weight * position_rates(D, br); wsum += br.weight
    return out / wsum


def scenario_names(zone_code, scen):
    """Scenario names PS_<zone>_<rig>_<Mw int>_<Mw dec>_<scal>_<sreal>_<tri>_<slip> from the columns of the scenario table."""
    return [f"PS_{zone_code}_{r:02d}_{m // 10000}_{m % 10000:04d}_{s:02d}_{q:02d}_{t:05d}_{l:03d}"
            for r, m, s, q, t, l in zip(scen["rig"], scen["mw"], scen["scal"], scen["sreal"], scen["tri"], scen["slip"])]


def scenario_rates(zone, branch_id):
    """Annual rate of every scenario of the zone on one branch, in the order of the scenario table
    (products/B1_inputs/scenarios_<zone>.npz). For use by the hazard calculation."""
    D = load(zone)
    br = D["table"][D["table"].branch_id == branch_id].iloc[0]
    scen = np.load(os.path.join(PRODUCTS, "B1_inputs", f"scenarios_{zone}.npz"))
    return position_rates(D, br)[scen["position"]] * scen["split"]
