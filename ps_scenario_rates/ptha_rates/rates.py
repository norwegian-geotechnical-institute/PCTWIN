"""Step B7. Annual rate of every rupture position and scenario on any branch.

    position rate = pooled weight (rigidity model x scaling law; decision 4)
                    x  rigidity factor of the position's model (geodynamic branches; budget with that model's
                       mu(depth) / budget with the homogeneous rigidity; decision 27)
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
    ev = cfg["tree"]["event_tree"]; rig_w = config.rigidity_in_use(cfg["tree"])
    pos = pd.DataFrame({k: S[k] for k in ("position", "rig", "scal", "bin", "tri", "cell", "m0_nm", "n_scenarios")})
    w_rig = np.array([rig_w.get(f"{r:02d}", 0.0) for r in pos.rig])
    w_scal = np.array([ev["scaling_law"][f"{s:02d}"] for s in pos.scal])
    bud = pd.read_csv(os.path.join(P, "B2_budget", f"budget_{zone}.csv")).set_index("unit")
    units = R["units"]
    fac = {}                                                     # level -> unit x position
    for lev in ("low", "preferred", "high"):
        fac[lev] = np.ones((len(units), len(pos)))
        for r in np.unique(pos.rig):
            col = f"rigidity_factor_{r:02d}_{lev}"
            if col in bud:
                fac[lev][:, (pos.rig == r).to_numpy()] = bud.loc[units, col].to_numpy()[:, None]
    fac_unit = {lev: sum(w * (bud.loc[units, f"rigidity_factor_{k}_{lev}"].to_numpy() if f"rigidity_factor_{k}_{lev}" in bud else np.ones(len(units)))
                         for k, w in rig_w.items()) for lev in ("low", "preferred", "high")}
    return dict(zone=zone, cfg=cfg, R=R, pos=pos, w_rig=w_rig, w_scal=w_scal, pooled=w_rig * w_scal, rig_factor=fac, rig_factor_unit=fac_unit, U=S["share_unconstrained"],
                X=np.load(os.path.join(P, "B6_segmented", f"shares_segmented_{zone}.npz"))["share_segmented"],
                table=pd.read_csv(os.path.join(P, "B0_config", f"branch_table_{zone}.csv")))


def expected_bin_rates(D, branch):
    """Total rate per bin that the positions of a branch must carry: the bin rates of the branch, scaled on geodynamic
    branches by the pooled rigidity factor of each unit (sum over rigidity models of weight x factor, decision 27)."""
    r = binrates.branch_bin_rates(D["R"], branch)                    # unit x bin
    f = D["rig_factor_unit"][branch.coupling] if branch.rate_source == "geodynamic" else np.ones(len(r))
    return r[0] * f[0] if branch.segmentation == "unconstrained" else (r[1:] * f[1:, None]).sum(axis=0)


def pooled_weights(D, branch):
    """Pooled weight of every position on one branch: rigidity model x scaling law. When the scaling law is a
    logic-tree level (event_tree: as_level, decision 28) the branch keeps only the positions of its own law."""
    law = getattr(branch, "scaling_law", "pooled")
    if law in ("pooled", None) or (isinstance(law, float) and np.isnan(law)):
        return D["pooled"]
    return D["w_rig"] * (D["pos"].scal.to_numpy() == int(law))


def position_rates(D, branch):
    """Annual rate of every position on one branch (a row of the branch table)."""
    r = binrates.branch_bin_rates(D["R"], branch)            # unit x bin
    k = D["pos"].bin.to_numpy()
    w = pooled_weights(D, branch)
    geo = branch.rate_source == "geodynamic"
    F = D["rig_factor"][branch.coupling] if geo else None   # unit x position
    if branch.segmentation == "unconstrained":
        v = shares.VARIANTS.index(binrates.variant_of(branch.rate_source, branch.coupling))
        return w * r[0][k] * D["U"][v] * (F[0] if geo else 1.0)
    X = D["X"] if D["X"].ndim == 2 else D["X"][shares.VARIANTS.index(binrates.variant_of(branch.rate_source, branch.coupling))]   # per variant with footprints
    seg = r[1:][:, k].T * X                                  # position x segment
    if geo:
        seg = seg * F[1:].T
    return w * seg.sum(axis=1)


def sample_split(scen, branch):
    """Split of each scenario on a branch. Pooled samples (default): 1 / (samples x slip realisations at the
    position). With scaling_sample as a level (decision 28) the branch's sample takes the whole position rate,
    shared over its slip realisations; a position without that sample keeps the pooled split."""
    k = getattr(branch, "scaling_sample", "pooled")
    if k in ("pooled", None) or (isinstance(k, float) and np.isnan(k)):
        return scen["split"]
    mine = scen["sreal"] == int(k)
    n_mine = np.bincount(scen["position"], weights=mine, minlength=int(scen["position"].max()) + 1)[scen["position"]]
    return np.where(n_mine > 0, mine / np.maximum(n_mine, 1), scen["split"])


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
    return position_rates(D, br)[scen["position"]] * sample_split(scen, br)
