"""Sensitivity test: bin rates with a different preferred coupling in some segments.

Run from ps_scenario_rates after step B4:
    python scripts/sens_coupling.py sunda Java=0.176 Arakan=0.45

Nothing in config/ or in the products of steps B0 to B4 is changed. The coupling is replaced in
memory, steps B2 and B4 are recomputed, and the result is compared with the present products.
Writes products/sensitivity/coupling_<tag>/: report.md, large_events.csv, budget.csv, and (after step B7)
rate_by_segment.csv and rate_along_trench.png, the rate along the trench with and without the change;
mfd.png (whole zone) and mfd_segments.png (the segments that were changed).
"""
import copy
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, binrates, branches, budget, catalogue, config, footprints, mfd, plotting, rates, shares
import run_B4_bin_rates as b4
plt = plotting.plt


def unit_cumulative(R, table, par, iu, m, m_top):
    """Annual rate of events >= m in one unit (segment) on every branch (branch x len(m))."""
    out = np.zeros((len(table), len(m)))
    for i, br in enumerate(table.itertuples()):
        iv = R["variants"].index(binrates.variant_of(br.rate_source, br.coupling))
        i_f, ib, im = R["forms"].index(br.mfd_form), int(br.b_index) - 1, R["mmax"].index(br.mmax)
        out[i] = mfd.cumulative(br.mfd_form, m, R["b_values"][ib], par.iloc[iu][binrates.MMAX_COL[br.mmax]], R["scale"][iu, iv, i_f, ib, im], m_top, R["tail"])
    return out


def mfd_figures(out, zone, what, Rs, table, par, cat, bins, names, changed):
    w = table.weight.to_numpy(); m_top = bins.mw_upper.iloc[-1]; mm = np.linspace(bins.mw_lower.iloc[0], m_top, 120)
    geo = (table.rate_source == "geodynamic").to_numpy(); unc = (table.segmentation == "unconstrained").to_numpy()
    years = cat.year.max() - cat.year.min() + 1
    mean = lambda C, m: (C[m] * w[m, None]).sum(0) / w[m].sum()
    obs = lambda c: (np.sort(c.mw[c.mw >= mm[0]].to_numpy())[::-1], np.arange(1, int((c.mw >= mm[0]).sum()) + 1) / years)
    CC = {k: b4.zone_cumulative(Rs[k], table, par, mm, m_top) for k in ("present", "test")}
    # whole zone: every geodynamic branch, as configured and in the test
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.6), sharey=True)
    for ax, k, ttl in zip(axes, ("present", "test"), ("coupling as configured", what)):
        for i in np.where(geo)[0]:
            ax.semilogy(mm, np.where(CC[k][i] > 0, CC[k][i], np.nan), color="0.75", lw=0.3, ls="-" if unc[i] else "--")
        ax.semilogy(mm, mean(CC[k], geo), color="black", lw=1.6, label="geodynamic, weighted mean")
        if k == "test":
            ax.semilogy(mm, mean(CC["present"], geo), color="black", lw=1.0, ls=":", label="geodynamic as configured, weighted mean")
        ax.semilogy(mm, mean(CC["present"], ~geo), color="black", lw=1.2, ls="--", label="seismic, weighted mean")
        ax.semilogy(*obs(cat), "o", ms=3, mfc="white", mec="black", label=f"observed, {cat.year.min()}-{cat.year.max()}")
        ax.set_title(f"{zone.capitalize()}: geodynamic branches, {ttl}"); ax.set_xlabel("Mw"); ax.set_ylim(1e-6, 30); ax.legend(frameon=False, loc="lower left")
    axes[0].set_ylabel("events per year, magnitude >= Mw, whole zone")
    fig.text(0.01, 0.01, "Grey: one line per geodynamic branch (solid: unconstrained; dashed: segmented, summed over segments).\n"
             "The seismic branches do not depend on coupling. Observed rates assume the catalogue is complete above Mw 6.35 for its whole length.", fontsize=plotting.NOTE, va="bottom")
    fig.tight_layout(rect=(0, 0.09, 1, 1)); fig.savefig(os.path.join(out, "mfd.png")); plt.close(fig)
    # the segments that were changed: segmented branches only (the only ones with an MFD per segment)
    fig, axes = plt.subplots(1, len(changed), figsize=(5.5 * len(changed), 5.6), sharey=True, squeeze=False)
    for ax, n in zip(axes[0], changed):
        iu = names.index(n) + 1; seg = ~unc
        Cp, Ct = (unit_cumulative(Rs[k], table, par, iu, mm, m_top) for k in ("present", "test"))
        for i in np.where(geo & seg)[0]:
            ax.semilogy(mm, np.where(Ct[i] > 0, Ct[i], np.nan), color="0.75", lw=0.3)
        ax.semilogy(mm, mean(Ct, geo & seg), color="black", lw=1.6, label="geodynamic in the test, weighted mean")
        ax.semilogy(mm, mean(Cp, geo & seg), color="black", lw=1.0, ls=":", label="geodynamic as configured, weighted mean")
        ax.semilogy(mm, mean(Cp, ~geo & seg), color="black", lw=1.2, ls="--", label="seismic, weighted mean")
        c = cat[cat.segment == iu]
        ax.semilogy(*obs(c), "o", ms=3, mfc="white", mec="black", label=f"observed in {n}, {cat.year.min()}-{cat.year.max()}")
        ax.set_title(f"{n}: segmented branches, {what}" if len(changed) == 1 else f"{n}: segmented branches"); ax.set_xlabel("Mw"); ax.set_ylim(1e-6, 30); ax.legend(frameon=False, loc="lower left")
    axes[0][0].set_ylabel("events per year, magnitude >= Mw, in the segment")
    fig.text(0.01, 0.01, f"Grey: one line per geodynamic segmented branch in the test ({what}).", fontsize=plotting.NOTE, va="bottom")
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(os.path.join(out, "mfd_segments.png")); plt.close(fig)


def summary(R, table, par, cat, m_top, label):
    w = table.weight.to_numpy(); thr = np.array(catalogue.THRESHOLDS); years = cat.year.max() - cat.year.min() + 1
    C = b4.zone_cumulative(R, table, par, np.r_[R["mw_lower0"], thr], m_top)
    geo = (table.rate_source == "geodynamic").to_numpy(); unc = (table.segmentation == "unconstrained").to_numpy()
    rows = []
    for name, m in (("geodynamic, unconstrained", geo & unc), ("geodynamic, segmented", geo & ~unc), ("seismic", ~geo), ("all branches", np.ones(len(w), bool))):
        r = dict(case=label, group=name, rate_first_bin_per_yr=float((C[m, 0] * w[m]).sum() / w[m].sum()))
        for j, t in enumerate(thr, start=1):
            r[f"n_ge_{t:g}"] = float((C[m, j] * w[m]).sum() / w[m].sum()) * years
        rows.append(r)
    return pd.DataFrame(rows)


def main(zone, overrides):
    cfg = config.load(); P = PRODUCTS
    tri = pd.read_csv(os.path.join(P, "B1_inputs", f"triangles_{zone}.csv")); bins = pd.read_csv(os.path.join(P, "B1_inputs", f"bins_{zone}.csv"))
    par = pd.read_csv(os.path.join(P, "B3_catalogue", f"parameters_{zone}.csv")); cat = pd.read_csv(os.path.join(P, "B3_catalogue", f"catalogue_{zone}.csv"))
    table = pd.read_csv(os.path.join(P, "B0_config", f"branch_table_{zone}.csv")); m_top = bins.mw_upper.iloc[-1]
    cfg2 = copy.deepcopy(cfg)
    for s in cfg2["zones"][zone]["segments"]:
        if s["name"] in overrides:
            s["coupling"] = overrides[s["name"]]
    out_rows, buds, Rs, loads = [], {}, {}, {}
    for label, c in (("present", cfg), ("test", cfg2)):
        z = c["zones"][zone]
        loads[label] = budget.loading_table(tri, z, c["tree"]); bud = budget.budget_table(loads[label], z)
        R = binrates.compute(c, par, bud, bins); R["mw_lower0"] = bins.mw_lower.iloc[0]
        out_rows.append(summary(R, table, par, cat, m_top, label)); buds[label] = bud.set_index("unit"); Rs[label] = R
    S = pd.concat(out_rows, ignore_index=True)
    tag = "_".join(f"{k}{v:g}" for k, v in overrides.items()).replace(".", "p")
    out = os.path.join(P, "sensitivity", f"coupling_{tag}"); os.makedirs(out, exist_ok=True)
    S.to_csv(os.path.join(out, "large_events.csv"), index=False, float_format="%.4g")
    B = pd.DataFrame({"budget_preferred_present": buds["present"].budget_preferred_nm_yr, "budget_preferred_test": buds["test"].budget_preferred_nm_yr,
                      "coupling_present": buds["present"].coupling_preferred, "coupling_test": buds["test"].coupling_preferred})
    B["ratio"] = B.budget_preferred_test / B.budget_preferred_present
    B.to_csv(os.path.join(out, "budget.csv"), float_format="%.4g")
    thr = catalogue.THRESHOLDS; obs = [int((cat.mw >= t).sum()) for t in thr]
    obs0 = float((cat.mw >= bins.mw_lower.iloc[0]).sum()) / (cat.year.max() - cat.year.min() + 1)
    L = [f"# Sensitivity test: coupling {', '.join(f'{k} = {v:g}' for k, v in overrides.items())} ({zone})", "",
         "Nothing in the configuration or in the products of steps B0 to B4 was changed. Low and high coupling follow the coupling rule "
         "of the logic tree (decision 25).", "", "## Preferred coupling and budget", "",
         "| Unit | Coupling now | Coupling in test | Budget now (N m/yr) | Budget in test | Ratio |", "|---|---|---|---|---|---|"]
    L += [f"| {u} | {r.coupling_present:.2f} | {r.coupling_test:.2f} | {r.budget_preferred_present:.3e} | {r.budget_preferred_test:.3e} | {r.ratio:.2f} |" for u, r in B.iterrows()]
    L += ["", f"## Events in {cat.year.max() - cat.year.min() + 1} years, weighted mean over branches", "",
          f"| Group | Case | Rate at Mw >= {bins.mw_lower.iloc[0]:.2f} (per yr) | " + " | ".join(f"Mw >= {t:g}" for t in thr) + " |", "|---|---|---|" + "---|" * len(thr),
          f"| Observed | | {obs0:.2f} | " + " | ".join(str(o) for o in obs) + " |"]
    for g in S.group.unique():
        for _, r in S[S.group == g].iterrows():
            L.append(f"| {g} | {r.case} | {r.rate_first_bin_per_yr:.2f} | " + " | ".join(f"{r[f'n_ge_{t:g}']:.1f}" for t in thr) + " |")
    # where the rate falls along the trench: bin rates and unconstrained shares both change with the coupling
    z = cfg["zones"][zone]; names = [s["name"] for s in z["segments"]]; n_seg = len(names); what = ", ".join(f"{k} = {v:g}" for k, v in overrides.items())
    D = rates.load(zone, cfg); pos = shares.position_table(np.load(os.path.join(P, "B1_inputs", f"scenarios_{zone}.npz")), n_seg)
    assert (pos.tri.to_numpy() == D["pos"].tri.to_numpy()).all()

    if cfg["tree"]["options"]["footprints"] == "segment_fraction_approximation":
        def geo_shares(bud):
            b = bud.loc[names]
            vals = {"seismic": np.ones(n_seg), **{f"geodynamic_{lv}": b[f"loading_{lv}_mm_yr"].to_numpy() for lv in ("low", "preferred", "high")}}
            return shares.unconstrained_shares(pos, vals, n_seg)[1:]
        assert np.allclose(geo_shares(buds["present"]), D["U"][1:]), "present shares are not reproduced"
        U2 = D["U"].copy(); U2[1:] = geo_shares(buds["test"]); D2 = dict(D, R=Rs["test"], U=U2)
    else:                                                   # with rupture footprints the segmented shares change with the loading as well
        fp = footprints.load(shares.footprint_path(zone))
        tv = {k: shares.triangle_values(loads[k], par, buds[k].reset_index(), names) for k in ("present", "test")}
        ratio = shares.mu_ratio(cfg, zone, tri, np.unique(pos.rig))
        assert np.allclose(shares.unconstrained_shares_fp(pos, fp, tv["present"], ratio), D["U"]), "present shares are not reproduced"
        D2 = dict(D, R=Rs["test"], U=shares.unconstrained_shares_fp(pos, fp, tv["test"], ratio),
                  X=shares.segmented_shares_fp(pos, fp, tv["test"], tri.segment.to_numpy(), n_seg, cfg["tree"]["options"]["crossing_penalty"], ratio)[0])
    geo = (table.rate_source == "geodynamic").to_numpy()
    PMs, PMg, PMg2 = rates.mean_position_rates(D, ~geo), rates.mean_position_rates(D, geo), rates.mean_position_rates(D2, geo)
    tri_i = tri.set_index("id"); kb = D["pos"].bin.to_numpy(); F = pos[[f"seg_{i + 1}" for i in range(n_seg)]].to_numpy()
    plotting.rate_along_trench(os.path.join(out, "rate_along_trench.png"), f"{zone.capitalize()}: weighted mean rate along the trench, coupling test", tri_i,
                               pd.read_csv(os.path.join(P, "B1_inputs", f"trench_{zone}.csv")), tri_i.along_trench_km.reindex(D["pos"].tri.to_numpy()).to_numpy(), kb, bins, names,
                               [(PMs, "--", "black", "seismic branches"), (PMg, "-", "0.6", "geodynamic, coupling as configured"), (PMg2, "-", "black", f"geodynamic, {what}")],
                               note=f"Test: preferred coupling {what}.")
    years = cat.year.max() - cat.year.min() + 1; big = (bins.mw_centre >= 7.5).to_numpy()[kb]
    G = pd.DataFrame({"segment": names})
    for lab, pm in (("seismic", PMs), ("geodynamic_present", PMg), ("geodynamic_test", PMg2)):
        G[f"{lab}_per_yr"] = (pm[:, None] * F).sum(0); G[f"{lab}_n_ge_7.5"] = (pm[big, None] * F[big]).sum(0) * years
    G.to_csv(os.path.join(out, "rate_by_segment.csv"), index=False, float_format="%.4g")
    L += ["", f"## Rate by segment, weighted mean over branches (all modelled magnitudes per year; Mw >= 7.5 in {years} years)", "",
          "| Segment | Seismic | Geodynamic now | Geodynamic in test | Seismic, Mw >= 7.5 | Geodynamic now, Mw >= 7.5 | Geodynamic in test, Mw >= 7.5 |", "|---|---|---|---|---|---|---|"]
    L += [f"| {r.segment} | {r.seismic_per_yr:.3f} | {r.geodynamic_present_per_yr:.3f} | {r.geodynamic_test_per_yr:.3f} | {r['seismic_n_ge_7.5']:.1f} | {r['geodynamic_present_n_ge_7.5']:.1f} | {r['geodynamic_test_n_ge_7.5']:.1f} |" for _, r in G.iterrows()]
    L += ["", "A rupture's rate is divided between segments by its segment fractions. rate_along_trench.png shows the same along the trench."]
    mfd_figures(out, zone, what, Rs, table, par, cat, bins, names, [n for n in names if n in overrides])
    L += ["", "mfd.png: cumulative rates of the whole zone on every geodynamic branch, as configured and in the test. mfd_segments.png: the same for each segment that was changed (segmented branches)."]
    with open(os.path.join(out, "report.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main(sys.argv[1], {k: float(v) for k, v in (a.split("=") for a in sys.argv[2:])})
