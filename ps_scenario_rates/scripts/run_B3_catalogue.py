"""Step B3. Catalogue by segment, largest observed events and the observed rate that anchors the seismic branches.

Run from ps_scenario_rates after step B1:   python scripts/run_B3_catalogue.py [zone ...]     (default: sunda makran)
Writes products/B3_catalogue/ and checks/B3_catalogue/.

Products per zone:
  parameters_<zone>.csv   per unit (whole zone and segments): catalogue b, anchor magnitude and observed rate
                          above it, the three Mmax alternatives, counts of large events, Mc and its rate
  catalogue_<zone>.csv    the catalogue with the segment of each event
  b_values_<zone>.csv     the b-values of the logic tree and their weights
  catalogue_fit_<zone>.png  observed cumulative rates with the Gutenberg-Richter lines of every b-value
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, branches, catalogue, config, mfd, plotting
from ptha_rates.report import CheckReport
plt = plotting.plt


def figure(cat, par, cfg, zone, path):
    tree = cfg["tree"]; bv = branches.b_values(tree); z = cfg["zones"][zone]
    nrow = int(np.ceil(len(par) / 3))
    fig, axes = plt.subplots(nrow, 3, figsize=(11, 3.7 * nrow + 0.7), sharex=True, squeeze=False)
    for ax, (_, u) in zip(axes.ravel(), par.iterrows()):
        e = cat if u.segment == 0 else cat[cat.segment == u.segment]
        e1 = e[(e.year >= u.mc_start_year) & (e.mw >= u.mc)]
        m = np.sort(e1.mw.to_numpy())[::-1]
        ax.semilogy(m, np.arange(1, len(m) + 1) / u.years_complete, color="0.45", lw=0.9, ls="--")
        e2 = e[(e.year >= u.anchor_start_year) & (e.mw >= u.anchor_mw)]
        m = np.sort(e2.mw.to_numpy())[::-1]
        ax.semilogy(m, np.arange(1, len(m) + 1) / u.anchor_years, color="black", lw=1.4)
        mm = np.linspace(u.anchor_mw, u.mmax_mean, 200)
        for b in bv:
            A = u.rate_anchor_per_yr / float(mfd.shape("truncated_gr", u.anchor_mw, b, u.mmax_mean, u.mmax_mean))
            ax.semilogy(mm, A * mfd.shape("truncated_gr", mm, b, u.mmax_mean, u.mmax_mean), color="0.6", lw=0.5)
        ax.axvline(u.anchor_mw, color="black", lw=0.6, ls=":")
        name = u.unit.replace("_", " ")
        ax.set_title(f"{name[0].upper() + name[1:]}: {u.n_ge_anchor} events")
        ax.set_ylim(1e-4, 1e3)
    for ax in axes.ravel()[len(par):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("Mw")
        ax.xaxis.set_tick_params(labelbottom=True)
    for ax in axes[:, 0]:
        ax.set_ylabel("events per year, magnitude >= Mw")
    a = z["seismic_anchor"]
    fig.text(0.01, 0.01, f"Black: observed rate of events above the anchor magnitude (Mw {a['magnitude']:g}, since {a['start_year']}). Grey: truncated G-R for the b-values of the logic tree, anchored there.\n"
             "Dashed: observed rate above Mc in its complete period (mc_start_year in zones.yaml); not used for the rates.", fontsize=plotting.NOTE, va="bottom")
    fig.tight_layout(rect=(0, 0.16 / nrow, 1, 1)); fig.savefig(path); plt.close(fig)


def main(zones):
    cfg = config.load(); tree = cfg["tree"]
    out = os.path.join(PRODUCTS, "B3_catalogue"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B3_catalogue", "catalogue rates and parameter sets")
    for zone in zones:
        z = cfg["zones"][zone]
        tri = pd.read_csv(os.path.join(PRODUCTS, "B1_inputs", f"triangles_{zone}.csv"))
        cat = catalogue.load_catalogue(z, tri)
        par = catalogue.parameter_table(cat, z, cfg["tree"])
        seg = par[par.segment > 0]; whole = par[par.segment == 0].iloc[0]
        allc = catalogue.load_catalogue(z, tri, keep_excluded=True); x = allc[allc.excluded != ""]
        big = x[x.mw >= 6.0].sort_values("mw", ascending=False)
        rep.check(f"{zone}: events removed by the catalogue filter", "zones.yaml catalogue_filter: interface events only, as far as location allows", None,
                  ("; ".join(f"{k}: {v:,}" for k, v in x.excluded.value_counts().items()) + ". Removed events of Mw >= 6: " +
                   ", ".join(f"{int(e.year)} Mw {e.mw:.1f} ({e.excluded})" for _, e in big.iterrows())) if len(x) else "no filter")
        rep.check(f"{zone}: events in the catalogue", "Whole catalogue is read (after the filter)", None,
                  f"{len(cat):,} events, {cat.year.min()} to {cat.year.max()}, Mw {cat.mw.min():.2f} to {cat.mw.max():.2f}")
        rep.check(f"{zone}: distance from each event to the mesh", "Events are assigned to the segment of the nearest triangle; no depth filter is applied",
                  bool(cat.dist_km.max() < 50), f"median {cat.dist_km.median():.1f} km, largest {cat.dist_km.max():.1f} km")
        rep.check(f"{zone}: every event has a segment", "Segment counts must add up to the catalogue", int(seg.n_events.sum()) == len(cat),
                  "; ".join(f"{u.unit}: {u.n_events:,}" for _, u in seg.iterrows()))
        few = par[par.n_ge_anchor < 10]
        rep.check(f"{zone}: units anchored on fewer than 10 events", "The anchor rate of these units is uncertain (one standard error is rate / sqrt(n))", None,
                  "; ".join(f"{u.unit}: {u.n_ge_anchor} events, rate uncertain by about {100 / np.sqrt(u.n_ge_anchor):.0f}%" for _, u in few.iterrows()) or "none")
        rep.check(f"{zone}: events above the anchor magnitude", "The seismic a-value is anchored on this count (decision 20); every unit needs at least one event",
                  bool((par.n_ge_anchor >= 1).all()),
                  "; ".join(f"{u.unit}: {u.n_ge_anchor} events of Mw >= {u.anchor_mw:g} in {u.anchor_years} yr = {u.rate_anchor_per_yr:.3f}/yr" for _, u in par.iterrows()))
        a = z["seismic_anchor"]; e_ = cat[cat.mw >= a["magnitude"]]
        r1 = (e_.year < 1964).sum() / (1964 - cat.year.min()); r2 = ((e_.year >= 1964) & (e_.year < 2000)).sum() / 36; r3 = (e_.year >= 2000).sum() / (z["catalogue_end_year"] - 1999)
        rep.check(f"{zone}: rate above the anchor magnitude by period", "Tests whether the catalogue is complete above the anchor from its start, and shows the 2004-2012 cluster", None,
                  f"before 1964: {r1:.2f}/yr; 1964-1999: {r2:.2f}/yr; from 2000: {r3:.2f}/yr")
        rep.check(f"{zone}: rate above Mc in its complete period", "For the record; not used for the rates",
                  None,
                  "; ".join(f"{u.unit}: Mc {u.mc} from {int(u.mc_start_year)}, {u.n_ge_mc:,} events in {u.years_complete} yr = {u.rate_mc_per_yr:.1f}/yr" for _, u in par.iterrows()))
        o = tree["options"].get("mmax_observed", {"n_sigma": 0, "sigma": 0})
        rep.check(f"{zone}: Mmax alternatives", f"Observed = largest event + {o['n_sigma']:g} x {o['sigma']:g} (decision 23), at most the area limit",
                  bool((par.mmax_observed <= par.mmax_area_limit).all() and (par.mmax_observed >= par.largest_event_mw).all()),
                  "; ".join(f"{u.unit}: largest event {u.largest_event_mw:.2f}, observed {u.mmax_observed:.2f}, area limit {u.mmax_area_limit:.2f}, mean {u.mmax_mean:.2f}" for _, u in par.iterrows()))
        rep.check(f"{zone}: whole-zone largest event is the largest of the segments", "Consistency of the two segmentation alternatives",
                  bool(np.isclose(whole.mmax_observed, seg.mmax_observed.max())), f"{whole.mmax_observed:.2f}")
        rep.check(f"{zone}: large events since {cat.year.min()}", "Used later to compare modelled with observed numbers", None,
                  "; ".join(f"Mw>={t:g}: {int(whole[f'n_ge_{t:g}'])}" for t in catalogue.THRESHOLDS))
        bv = branches.b_values(tree); bw = branches.b_weights(z["whole_zone"]["b"], tree)
        rep.check(f"{zone}: b-values of the logic tree", "Values and weights of the b-value level (decisions 2 and 24); weighted mean",
                  bool(len(bv) == len(bw) and np.all(np.diff(bv) > 0)),
                  f"{bv[0]:.3f} to {bv[-1]:.3f} in steps of {bv[1] - bv[0]:.4f}; weighted mean {(bv * bw).sum():.3f}")
        # how well does each unit's own catalogue b reproduce the observed number of M >= 6 events in the complete period?
        txt = []
        for _, u in par.iterrows():
            e = cat if u.segment == 0 else cat[cat.segment == u.segment]
            obs6 = int(((e.mw >= 6.0) & (e.year >= u.mc_start_year)).sum())
            A = u.rate_mc_per_yr / float(mfd.shape("truncated_gr", u.mc, u.b_catalogue, u.mmax_mean, u.mmax_mean))
            txt.append(f"{u.unit}: observed {obs6}, G-R with catalogue b {A * float(mfd.shape('truncated_gr', 6.0, u.b_catalogue, u.mmax_mean, u.mmax_mean)) * u.years_complete:.0f}")
        rep.check(f"{zone}: events of Mw >= 6 since Mc became complete, observed against an anchor at Mc", "Why the rates are not anchored at Mc: extrapolating the rate at Mc with the catalogue b gives too many larger events", None, "; ".join(txt))
        par.to_csv(os.path.join(out, f"parameters_{zone}.csv"), index=False, float_format="%.6g")
        cat.round(4).to_csv(os.path.join(out, f"catalogue_{zone}.csv"), index=False)
        pd.DataFrame({"b_index": np.arange(1, len(bv) + 1), "b_value": bv.round(6), "weight": bw}).to_csv(os.path.join(out, f"b_values_{zone}.csv"), index=False)
        figure(cat, par, cfg, zone, os.path.join(out, f"catalogue_fit_{zone}.png"))
    rep.note("The catalogue is not declustered, as in the method document. Aftershocks of 2004 to 2012 are inside the complete period.")
    rep.note("Seismic branches are anchored on the observed rate of events at or above the anchor magnitude over the whole catalogue (decision 20), "
             "not on the rate at Mc: the period in which the catalogue is complete at Mc (from 2000) holds the 2004-2012 cluster.")
    rep.note("Start years of the complete period at Mc are not in Table 1. They were read from the number of events above Mc per five years and are kept for the record only.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda", "makran"]) else 1)
