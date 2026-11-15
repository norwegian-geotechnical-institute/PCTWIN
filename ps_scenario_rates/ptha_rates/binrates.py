"""Step B4. Bin rates of every unit (whole zone and segments) for every combination of
rate variant, MFD form, b-value and Mmax alternative.

Array axes:  unit x variant x form x b x mmax x bin
  variant: seismic, geodynamic_low, geodynamic_preferred, geodynamic_high
"""
import numpy as np
import pandas as pd
from . import branches, mfd

VARIANTS = ["seismic", "geodynamic_low", "geodynamic_preferred", "geodynamic_high"]
MMAX_COL = {"observed": "mmax_observed", "area_limit": "mmax_area_limit", "mean": "mmax_mean"}


def compute(cfg, par, budget, bins):
    """Returns a dict with the rate array, the scale A, xi, the mean moment per bin and the axis labels."""
    tree = cfg["tree"]; lv = tree["levels"]
    forms, mmaxs, bv = list(lv["mfd_form"]), list(lv["mmax"]), branches.b_values(tree)
    tail = tree["options"]["tapered_tail"]
    lo, hi = bins.mw_lower.to_numpy(), bins.mw_upper.to_numpy()
    units = par.unit.tolist()
    bud = budget.set_index("unit")
    shp = (len(units), len(VARIANTS), len(forms), len(bv), len(mmaxs))
    rate = np.zeros(shp + (len(lo),)); A = np.zeros(shp); XI = np.zeros(shp[0:1] + shp[2:])
    M0 = np.full(shp[0:1] + shp[2:] + (len(lo),), np.nan)
    for iu, u in par.reset_index(drop=True).iterrows():
        for i_f, form in enumerate(forms):
            for ib, b in enumerate(bv):
                for im, mm in enumerate(mmaxs):
                    mmax = u[MMAX_COL[mm]]
                    r, a, x, m0k = mfd.bin_rates(form, lo, hi, b, mmax, tail, rate_mc=u.rate_anchor_per_yr, mc=u.anchor_mw)
                    rate[iu, 0, i_f, ib, im], A[iu, 0, i_f, ib, im] = r, a
                    XI[iu, i_f, ib, im], M0[iu, i_f, ib, im] = x, m0k
                    for iv, lev in enumerate(("low", "preferred", "high"), start=1):
                        r, a, _, _ = mfd.bin_rates(form, lo, hi, b, mmax, tail, budget=bud.loc[u.unit, f"budget_{lev}_nm_yr"])
                        rate[iu, iv, i_f, ib, im], A[iu, iv, i_f, ib, im] = r, a
    return dict(rate=rate, scale=A, xi=XI, mean_m0=M0, units=units, variants=VARIANTS, forms=forms,
                b_values=bv, mmax=mmaxs, tail=tail)


def variant_of(rate_source, coupling):
    return "seismic" if rate_source == "seismic" else f"geodynamic_{coupling}"


def branch_bin_rates(R, branch):
    """Bin rates of one branch (a row of the branch table): array unit x bin. For an unconstrained
    branch use the row of 'whole_zone'; for a segmented branch the rows of the segments."""
    iv = R["variants"].index(variant_of(branch.rate_source, branch.coupling))
    return R["rate"][:, iv, R["forms"].index(branch.mfd_form), int(branch.b_index) - 1, R["mmax"].index(branch.mmax)]


def zone_bin_rates(R, table):
    """Total bin rate of the zone on every branch (branch x bin): the whole-zone unit for
    unconstrained branches, the sum of the segments for segmented ones."""
    out = np.zeros((len(table), R["rate"].shape[-1]))
    for i, br in enumerate(table.itertuples()):
        r = branch_bin_rates(R, br)
        out[i] = r[0] if br.segmentation == "unconstrained" else r[1:].sum(axis=0)
    return out


def long_table(R, bins):
    idx = pd.MultiIndex.from_product([R["units"], R["variants"], R["forms"], np.round(R["b_values"], 6), R["mmax"], bins.mw_centre],
                                     names=["unit", "variant", "mfd_form", "b_value", "mmax", "mw_centre"])
    return pd.DataFrame({"rate_per_yr": R["rate"].ravel()}, index=idx).reset_index()
