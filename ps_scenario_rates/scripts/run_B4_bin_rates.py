"""Step B4. Annual rate of every magnitude bin on every branch.

Run from ps_scenario_rates after steps B0 to B3:   python scripts/run_B4_bin_rates.py [zone ...]     (default: sunda makran)
Writes products/B4_bin_rates/ and checks/B4_bin_rates/.

Products per zone:
  bin_rates_<zone>.npz        rate array (unit x variant x form x b x mmax x bin) with its axis labels, scale and xi
  bin_rates_<zone>.csv        the same rates as a long table
  zone_bin_rates_<zone>.csv   total rate of the zone in each bin, one row per branch
  mfd_<zone>.png              cumulative rates of all branches with the observed large events
  large_events_<zone>.csv     modelled against observed numbers of large events
"""
import math
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, binrates, catalogue, config, inputs, mfd, plotting
from ptha_rates.report import CheckReport
plt = plotting.plt


def wquantile(x, w, q):
    o = np.argsort(x); cw = np.cumsum(w[o]) / w.sum()
    return x[o][np.searchsorted(cw, q)]


def poisson_sf(k, mu):
    """P(N >= k) for a Poisson mean mu."""
    return 1.0 - sum(math.exp(-mu) * mu ** i / math.factorial(i) for i in range(k))


def zone_cumulative(R, table, par, m, m_top):
    """Annual rate of events >= m for the whole zone on every branch (branch x len(m))."""
    out = np.zeros((len(table), len(m)))
    for i, br in enumerate(table.itertuples()):
        iv = R["variants"].index(binrates.variant_of(br.rate_source, br.coupling))
        i_f, ib, im = R["forms"].index(br.mfd_form), int(br.b_index) - 1, R["mmax"].index(br.mmax)
        us = [0] if br.segmentation == "unconstrained" else range(1, len(R["units"]))
        for iu in us:
            mmax = par.iloc[iu][binrates.MMAX_COL[br.mmax]]
            out[i] += mfd.cumulative(br.mfd_form, m, R["b_values"][ib], mmax, R["scale"][iu, iv, i_f, ib, im], m_top, R["tail"])
    return out


def main(zones):
    cfg = config.load()
    out = os.path.join(PRODUCTS, "B4_bin_rates"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B4_bin_rates", "bin rates for every branch")
    for zone in zones:
        z = cfg["zones"][zone]
        bins = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"bins_{zone}.csv"))
        par = pd.read_csv(os.path.join(PRODUCTS, "B3_catalogue", f"parameters_{zone}.csv"))
        bud = pd.read_csv(os.path.join(PRODUCTS, "B2_budget", f"budget_{zone}.csv"))
        table = pd.read_csv(os.path.join(PRODUCTS, "B0_config", f"branch_table_{zone}.csv"))
        cat = pd.read_csv(os.path.join(PRODUCTS, "B3_catalogue", f"catalogue_{zone}.csv"))
        R = binrates.compute(cfg, par, bud, bins)
        rate, m_top = R["rate"], bins.mw_upper.iloc[-1]
        w = table.weight.to_numpy()

        rep.check(f"{zone}: rates are finite and not negative", "Basic sanity of every bin rate",
                  bool(np.isfinite(rate).all() and (rate >= 0).all()), f"{rate.size:,} values; {int((rate == 0).sum()):,} are zero (bins above Mmax in truncated G-R)")
        # geodynamic: sum of rate x mean moment = xi x budget
        err = 0.0
        for iu, unit in enumerate(R["units"]):
            for iv, lev in enumerate(("low", "preferred", "high"), start=1):
                b_ = bud.set_index("unit").loc[unit, f"budget_{lev}_nm_yr"]
                rel = np.nansum(rate[iu, iv] * R["mean_m0"][iu], axis=-1) / (R["xi"][iu] * b_) - 1
                err = max(err, np.abs(rel).max())
        rep.check(f"{zone}: geodynamic rates release the moment budget", "Sum of bin rate x mean moment must equal xi x budget for every unit, form, b, Mmax and coupling level",
                  bool(err < 1e-9), f"largest relative error {err:.1e}")
        # seismic: cumulative rate at Mc equals the observed rate
        err = 0.0
        for iu, u in par.iterrows():
            for i_f, form in enumerate(R["forms"]):
                for ib, b in enumerate(R["b_values"]):
                    for im, mm in enumerate(R["mmax"]):
                        mmax_ = u[binrates.MMAX_COL[mm]]
                        if form == "truncated_gr" and mmax_ <= u.anchor_mw:
                            continue
                        n = float(mfd.cumulative(form, u.anchor_mw, b, mmax_, R["scale"][iu, 0, i_f, ib, im], m_top, R["tail"]))
                        err = max(err, abs(n / u.rate_anchor_per_yr - 1))
        rep.check(f"{zone}: seismic rates reproduce the observed rate at the anchor magnitude", "The seismic a-value is anchored on the catalogue (decision 20)", bool(err < 1e-9), f"largest relative error {err:.1e}")
        ipref = R["variants"].index("geodynamic_preferred")
        rep.check(f"{zone}: higher coupling gives higher rates", "Rates are proportional to the budget",
                  bool((rate[:, 1] <= rate[:, 2] + 1e-30).all() and (rate[:, 2] <= rate[:, 3] + 1e-30).all()), "low <= preferred <= high in every bin")
        rep.check(f"{zone}: xi", "Fraction of the moment released above the smallest modelled magnitude; must lie between 0 and 1",
                  bool((R["xi"] > 0).all() and (R["xi"] <= 1).all()), f"{R['xi'].min():.3f} to {R['xi'].max():.3f}")
        Z = binrates.zone_bin_rates(R, table)
        rep.check(f"{zone}: one set of zone bin rates per branch", "Every branch of the branch table must find its rates", bool(Z.shape == (len(table), len(bins)) and (Z.sum(1) > 0).all()),
                  f"{Z.shape[0]} branches x {Z.shape[1]} bins")

        # ---- information
        seis = (table.rate_source == "seismic").to_numpy(); unc = (table.segmentation == "unconstrained").to_numpy()
        groups = {"seismic, unconstrained": seis & unc, "seismic, segmented": seis & ~unc,
                  "geodynamic, unconstrained": ~seis & unc, "geodynamic, segmented": ~seis & ~unc, "all branches": np.ones(len(table), bool)}
        thr = np.array(catalogue.THRESHOLDS); years = cat.year.max() - cat.year.min() + 1
        C = zone_cumulative(R, table, par, thr, m_top)
        obs = [int((cat.mw >= t).sum()) for t in thr]
        rows = []
        for name, m in groups.items():
            for j, t in enumerate(thr):
                mean = float((C[m, j] * w[m]).sum() / w[m].sum())
                rows.append(dict(group=name, mw_threshold=t, observed=obs[j], years=years, modelled_mean=mean * years,
                                 modelled_p05=wquantile(C[m, j], w[m], 0.05) * years, modelled_p95=wquantile(C[m, j], w[m], 0.95) * years,
                                 prob_at_least_observed=poisson_sf(obs[j], mean * years)))
        L = pd.DataFrame(rows); L.to_csv(os.path.join(out, f"large_events_{zone}.csv"), index=False, float_format="%.4g")
        for name in groups:
            g = L[L.group == name]
            rep.check(f"{zone}: large events in {years} years, {name}", "Modelled (weighted mean, 5% to 95% of branches) against observed. Reported only; no weights are changed (decision 13)", None,
                      "; ".join(f"Mw>={r.mw_threshold:g}: observed {r.observed}, modelled {r.modelled_mean:.1f} ({r.modelled_p05:.1f} to {r.modelled_p95:.1f})" for r in g.itertuples()))
        # implied coupling of the seismic branches
        unc_b = bud.set_index("unit").budget_uncoupled_nm_yr
        txt = []
        for name, m in (("unconstrained", seis & unc), ("segmented", seis & ~unc)):
            mom = np.array([np.nansum(binrates.branch_bin_rates(R, br)[[0] if br.segmentation == "unconstrained" else slice(1, None)]
                                      * R["mean_m0"][[0] if br.segmentation == "unconstrained" else slice(1, None), R["forms"].index(br.mfd_form), int(br.b_index) - 1, R["mmax"].index(br.mmax)]
                                      / R["xi"][[0] if br.segmentation == "unconstrained" else slice(1, None), R["forms"].index(br.mfd_form), int(br.b_index) - 1, R["mmax"].index(br.mmax)][:, None])
                            for br in table[m].itertuples()])
            c = mom / unc_b["whole_zone"]; ww = w[m]
            txt.append(f"{name}: weighted mean {float((c * ww).sum() / ww.sum()):.2f}, 5% to 95% {wquantile(c, ww, 0.05):.2f} to {wquantile(c, ww, 0.95):.2f}")
        rep.check(f"{zone}: coupling implied by the seismic branches", "Moment rate of a seismic branch divided by the budget with coupling 1; compare with the whole-zone coupling of the geodynamic branches (" + " / ".join(f"{bud.set_index('unit').loc['whole_zone', f'coupling_{lv}']:.2f}" for lv in ("low", "preferred", "high")) + ")", None, "; ".join(txt))
        mean_z = (Z * w[:, None]).sum(0)
        rep.check(f"{zone}: weighted mean rate of the zone", "Return period of events at or above a magnitude, all branches", None,
                  "; ".join(f"Mw>={bins.mw_lower[k]:.2f}: " + (f"{mean_z[k:].sum():.1f} per yr" if mean_z[k:].sum() >= 1 else f"1 per {1 / mean_z[k:].sum():,.0f} yr")
                            for k in inputs.summary_bins(bins)))
        tp = np.array([R["forms"].index(f) for f in table.mfd_form]) == R["forms"].index("tapered_pareto")
        k9 = int(min(np.searchsorted(bins.mw_lower.to_numpy(), 9.0), len(bins) - 1))
        rep.check(f"{zone}: tapered Pareto against truncated G-R above Mw {bins.mw_lower[k9]:.2f}", "Effect of the MFD form on the largest events (weighted mean rate of the zone)", None,
                  f"tapered {float((Z[tp, k9:].sum(1) * w[tp]).sum() / w[tp].sum()):.2e}/yr, truncated {float((Z[~tp, k9:].sum(1) * w[~tp]).sum() / w[~tp].sum()):.2e}/yr")

        # ---- products
        np.savez_compressed(os.path.join(out, f"bin_rates_{zone}.npz"), rate=rate, scale=R["scale"], xi=R["xi"], mean_m0=R["mean_m0"],
                            units=np.array(R["units"]), variants=np.array(R["variants"]), forms=np.array(R["forms"]),
                            b_values=R["b_values"], mmax=np.array(R["mmax"]), mw_centre=bins.mw_centre.to_numpy(),
                            mw_lower=bins.mw_lower.to_numpy(), mw_upper=bins.mw_upper.to_numpy())
        binrates.long_table(R, bins).to_csv(os.path.join(out, f"bin_rates_{zone}.csv"), index=False, float_format="%.6g")
        zt = pd.DataFrame(Z, columns=[f"mw_{c:.4f}" for c in bins.mw_centre]); zt.insert(0, "weight", w); zt.insert(0, "branch_id", table.branch_id)
        zt.to_csv(os.path.join(out, f"zone_bin_rates_{zone}.csv"), index=False, float_format="%.6g")

        # ---- figure
        mm = np.linspace(bins.mw_lower.iloc[0], m_top, 120)
        CC = zone_cumulative(R, table, par, mm, m_top)
        fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), sharey=True)
        for ax, (src, m) in zip(axes, (("Seismic branches", seis), ("Geodynamic branches", ~seis))):
            for i in np.where(m)[0]:
                ax.semilogy(mm, np.where(CC[i] > 0, CC[i], np.nan), color="0.75", lw=0.3, ls="-" if unc[i] else "--")
            ax.semilogy(mm, (CC[m] * w[m, None]).sum(0) / w[m].sum(), color="black", lw=1.6, label="weighted mean")
            anc = z["seismic_anchor"]; ok = (cat.mw >= anc["magnitude"]) & (cat.year >= anc["start_year"])
            big = np.sort(cat.mw[ok].to_numpy())[::-1]
            ax.semilogy(big, np.arange(1, len(big) + 1) / (z["catalogue_end_year"] - anc["start_year"] + 1), "o", ms=3, mfc="white", mec="black",
                        label=f"observed, {anc['start_year']}-{z['catalogue_end_year']}")
            ax.set_title(f"{zone.capitalize()}: {src} ({int(m.sum())})"); ax.set_xlabel("Mw")
            ax.set_ylim(1e-6, 30); ax.legend(frameon=False)
        axes[0].set_ylabel("events per year, magnitude >= Mw, whole zone")
        fig.text(0.01, 0.01, "Grey: one line per branch (solid: unconstrained; dashed: segmented, summed over segments).\n"
                 f"Observed: events at or above the anchor magnitude (Mw {z['seismic_anchor']['magnitude']:g}) since {z['seismic_anchor']['start_year']}, taken as complete.", fontsize=plotting.NOTE, va="bottom")
        fig.tight_layout(rect=(0, 0.09, 1, 1)); fig.savefig(os.path.join(out, f"mfd_{zone}.png")); plt.close(fig)
    rep.note("Scenario data are not used in this step. Rates are for the whole zone or segment; step B5 shares them among scenarios.")
    rep.note("Tapered Pareto: corner magnitude = the Mmax alternative; the shape is cut at the upper edge of the last bin "
             f"(tapered_tail = {cfg['tree']['options']['tapered_tail']}), so that the budget is balanced over the bins that have scenarios.")
    rep.note("The b-value of a branch is applied to every unit; each unit uses its own Mmax alternatives, catalogue rate and budget.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda", "makran"]) else 1)
