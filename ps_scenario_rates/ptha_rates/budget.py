"""Step B2. Loading (slip-deficit rate) on every triangle and the moment budget."""
import numpy as np
import pandas as pd
from . import branches

LEVELS = ["low", "preferred", "high"]


def loading_table(tri, zone_cfg, tree):
    """Loading = coupling x convergence (mm/yr) on each triangle, for the three coupling levels.
    First run: coupling and convergence are constant within a segment. The convergence rate is
    used as given, with no projection onto the interface (decision 6)."""
    if tree["options"]["convergence_projection"] != "none":
        raise NotImplementedError("only convergence_projection: none is implemented")
    out = tri[["id", "segment", "area_km2"]].copy()
    conv = np.array([s["convergence"] for s in zone_cfg["segments"]])[tri.segment.to_numpy() - 1]
    out["convergence_mm_yr"] = conv
    for lev in LEVELS:
        c = np.array([branches.coupling_levels(s["coupling"], tree)[lev] for s in zone_cfg["segments"]])[tri.segment.to_numpy() - 1]
        out[f"coupling_{lev}"] = c
        out[f"loading_{lev}_mm_yr"] = c * conv
    return out


def budget_table(load, zone_cfg):
    """Moment budget (N m/yr) = rigidity x sum of area x loading, per segment and for the whole
    zone (the sum of its segments), for each coupling level. Also the budget with coupling 1."""
    mu = zone_cfg["rigidity_pa"]
    area_m2 = load.area_km2.to_numpy() * 1e6
    rows = []
    for k, s in enumerate(zone_cfg["segments"], start=1):
        m = (load.segment == k).to_numpy()
        r = dict(unit=s["name"], area_km2=load.area_km2[m].sum(), convergence_mm_yr=s["convergence"])
        for lev in LEVELS:
            r[f"coupling_{lev}"] = load.loc[m, f"coupling_{lev}"].iloc[0]
            r[f"loading_{lev}_mm_yr"] = load.loc[m, f"loading_{lev}_mm_yr"].iloc[0]
            r[f"budget_{lev}_nm_yr"] = mu * (area_m2[m] * load.loc[m, f"loading_{lev}_mm_yr"].to_numpy() * 1e-3).sum()
        r["budget_uncoupled_nm_yr"] = mu * (area_m2[m] * load.convergence_mm_yr.to_numpy()[m] * 1e-3).sum()
        rows.append(r)
    t = pd.DataFrame(rows)
    whole = dict(unit="whole_zone", area_km2=t.area_km2.sum(),
                 convergence_mm_yr=(t.convergence_mm_yr * t.area_km2).sum() / t.area_km2.sum())
    for lev in LEVELS:
        whole[f"budget_{lev}_nm_yr"] = t[f"budget_{lev}_nm_yr"].sum()
        whole[f"loading_{lev}_mm_yr"] = (t[f"loading_{lev}_mm_yr"] * t.area_km2).sum() / t.area_km2.sum()
        whole[f"coupling_{lev}"] = whole[f"budget_{lev}_nm_yr"] / t.budget_uncoupled_nm_yr.sum()
    whole["budget_uncoupled_nm_yr"] = t.budget_uncoupled_nm_yr.sum()
    return pd.concat([pd.DataFrame([whole]), t], ignore_index=True)


def mw_of(m0):
    return (np.log10(m0) - 9.1) / 1.5
