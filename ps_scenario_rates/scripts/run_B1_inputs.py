"""Step B1. Freeze the inputs into checked tables.

Run from ps_scenario_rates:   python scripts/run_B1_inputs.py [zone ...]     (default: sunda)
Reads Data (never written to). Writes products/B1_inputs/ and checks/B1_inputs/.

Products per zone:
  triangles_<zone>.csv        id, centroid, depth, area, dip, segment, distance along the trench
  trench_<zone>.csv           the trench line (lon, lat, distance in km) used for plots along the arc
  bins_<zone>.csv             bin centre, edges, moment, counts of scenarios, positions and geometries
  scenarios_<zone>.npz        one entry per column, one value per scenario (rigidity models in use)
  realisations_<zone>.csv     numbers of scaling samples and slip realisations, by magnitude bin
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import pandas as pd
from ptha_rates import PRODUCTS, config, geometry, inputs
from ptha_rates.report import CheckReport

V2_ROWS = {"sunda": 7953344, "makran": 1018327}        # rows in the v2 fractions files that were checked


def main(zones):
    cfg = config.load()
    out = os.path.join(PRODUCTS, "B1_inputs"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B1_inputs", "frozen input tables")
    models = [int(k) for k in config.rigidity_in_use(cfg["tree"])]
    skipped = config.skipped_rigidity(cfg["tree"])
    for zone in zones:
        z = cfg["zones"][zone]; n_seg = len(z["segments"]); seg_cols = [f"seg_{i + 1}" for i in range(n_seg)]
        tri = inputs.triangle_table(z)
        s, n_all = inputs.read_scenarios(z, models)
        bins = inputs.bin_table(s.mw.to_numpy(), cfg["tree"]["options"]["bin_edges"])
        s = inputs.add_bins_and_counts(s, bins, cfg["tree"]["event_tree"]["scaling_sample_00_counts_as_sample"])
        tri["segment"], n_near = inputs.segment_map(s, tri, n_seg)
        tri["segment_name"] = [z["segments"][k - 1]["name"] for k in tri.segment]
        line = geometry.trench_line(tri, z["arc_centre"])
        tri["along_trench_km"] = geometry.along_trench(tri.lon, tri.lat, line)

        # ---- checks
        rep.check(f"{zone}: rows in the v2 fractions file", "The file read must be the v2 file that was checked",
                  n_all == V2_ROWS[zone], f"{n_all:,} read, {V2_ROWS[zone]:,} expected")
        rep.check(f"{zone}: scenarios kept", f"Rigidity model(s) in use {models} (decision 27)" + (f"; left out for want of footprints with slip: {skipped}" if skipped else ""), None,
                  f"{len(s):,} of {n_all:,}")
        key = s.rig.astype("int64") * 10 ** 17 + s.mw.astype("int64") * 10 ** 12 + s.scal.astype("int64") * 10 ** 10 \
            + s.sreal.astype("int64") * 10 ** 8 + s.tri.astype("int64") * 10 ** 3 + s.slip
        rep.check(f"{zone}: scenario names are unique", "Rates are keyed by scenario name", key.is_unique, f"{key.nunique():,} names")
        fs = s[seg_cols].sum(axis=1)
        rep.check(f"{zone}: segment fractions sum to 1", "Fractions divide a rupture's moment among segments",
                  bool((np.abs(fs - 1) < 1e-3).all()), f"min {fs.min():.5f}, max {fs.max():.5f}")
        rep.check(f"{zone}: total_M0 is positive", "Zero-moment rows were removed in v2", bool((s.total_M0 > 0).all()),
                  f"min {s.total_M0.min():.3e} N m")
        rep.check(f"{zone}: initiation triangles exist in the mesh", "Scenario position is its initiation triangle",
                  bool(s.tri.isin(tri.id).all()), f"triangle ids {s.tri.min()} to {s.tri.max()}, mesh has {len(tri):,}")
        rep.check(f"{zone}: every scenario is in one bin", "Rates are assigned by magnitude bin",
                  bool((bins.mw_centre.to_numpy()[s.bin] * 1e4).round().astype(int).tolist() == s.mw.tolist()), f"{len(bins)} bins")
        e = bins
        rep.check(f"{zone}: bin edges are continuous and contain their centres", "Bins must tile the magnitude range",
                  bool(np.allclose(e.mw_upper[:-1], e.mw_lower[1:]) and ((e.mw_lower < e.mw_centre) & (e.mw_centre < e.mw_upper)).all()),
                  f"Mw {e.mw_lower.iloc[0]:.4f} to {e.mw_upper.iloc[-1]:.4f}; widths {(e.mw_upper - e.mw_lower).min():.3f} to {(e.mw_upper - e.mw_lower).max():.3f}")
        tot = s.groupby("position").split.sum()
        rep.check(f"{zone}: split sums to 1 for every position", "No rate is lost or duplicated when a position's rate is divided among its realisations",
                  bool(np.allclose(tot, 1)), f"{len(tot):,} positions; sum of split = {s.split.sum():,.3f}")
        rep.check(f"{zone}: triangles without a starting scenario", "Their segment is taken from the nearest triangle (approximate map)",
                  None, f"{n_near:,} of {len(tri):,}")
        rep.check(f"{zone}: every triangle has a segment", "Loading is assigned by segment in step B2",
                  bool(tri.segment.between(1, n_seg).all()), ", ".join(f"{k}: {v:,}" for k, v in tri.segment_name.value_counts(sort=False).items()))
        g = tri.groupby("segment_name", sort=False).agg(triangles=("id", "size"), area_km2=("area_km2", "sum"),
                                                        lat_min=("lat", "min"), lat_max=("lat", "max"))
        order = [x["name"] for x in z["segments"]]
        rep.check(f"{zone}: area by segment (km2)", "Used for the moment budget; to compare with the provider's values", None,
                  "; ".join(f"{k}: {g.area_km2[k]:,.0f}" for k in order) + f"; total {tri.area_km2.sum():,.0f}")
        at = tri.groupby("segment").along_trench_km.median()
        rep.check(f"{zone}: distance along the trench", "Used for plots along the arc; segments must follow one another along it",
                  bool(at.is_monotonic_increasing), f"trench length {line.dist_km.iloc[-1]:,.0f} km; " + "; ".join(
                      f"{z['segments'][k - 1]['name']}: {tri.along_trench_km[tri.segment == k].quantile(.02):,.0f} to {tri.along_trench_km[tri.segment == k].quantile(.98):,.0f} km" for k in range(1, n_seg + 1)))
        lat_c = tri.groupby("segment").lat.mean()
        rep.check(f"{zone}: segments are ordered along the arc", "seg_1..seg_n must follow the order in zones.yaml",
                  bool(lat_c.is_monotonic_decreasing or lat_c.is_monotonic_increasing),
                  "; ".join(f"{k}: lat {g.lat_min[k]:.1f} to {g.lat_max[k]:.1f}" for k in order))
        dmw = (np.log10(s.total_M0) - 9.1) / 1.5 - s.mw / 1e4
        rep.check(f"{zone}: total_M0 against the named magnitude", "Known issue H2: some scenarios differ from the bin magnitude", None,
                  f"within 0.01: {100 * (dmw.abs() <= 0.01).mean():.1f}%; range {dmw.min():+.3f} to {dmw.max():+.3f}")

        # ---- products
        bins["n_scenarios"] = s.groupby("bin").size().reindex(bins.bin).to_numpy()
        bins["n_positions"] = s.groupby("bin").position.nunique().reindex(bins.bin).to_numpy()
        bins["n_geometries"] = s.groupby("bin").geometry.nunique().reindex(bins.bin).to_numpy()
        bins["pct_crossing"] = (100 * ((s[seg_cols] > 0).sum(axis=1) > 1).groupby(s.bin).mean()).reindex(bins.bin).round(2).to_numpy()
        bins.round(6).to_csv(os.path.join(out, f"bins_{zone}.csv"), index=False)
        tri.round(6).to_csv(os.path.join(out, f"triangles_{zone}.csv"), index=False)
        line.round(4).to_csv(os.path.join(out, f"trench_{zone}.csv"), index=False)
        p1 = s.drop_duplicates("position"); g1 = s.drop_duplicates("geometry")
        real = pd.DataFrame({"mw_centre": bins.mw_centre,
                             "positions": p1.groupby("bin").size().reindex(bins.bin).to_numpy(),
                             "sreal_per_position_mean": p1.groupby("bin").n_sreal.mean().reindex(bins.bin).round(2).to_numpy(),
                             "positions_with_1_sreal": (p1.n_sreal == 1).groupby(p1.bin).sum().reindex(bins.bin).to_numpy(),
                             "positions_with_10_sreal": (p1.n_sreal == 10).groupby(p1.bin).sum().reindex(bins.bin).to_numpy(),
                             "geometries": g1.groupby("bin").size().reindex(bins.bin).to_numpy(),
                             "geometries_with_1_slip": (g1.n_slip == 1).groupby(g1.bin).sum().reindex(bins.bin).to_numpy(),
                             "geometries_with_10_slip": (g1.n_slip == 10).groupby(g1.bin).sum().reindex(bins.bin).to_numpy(),
                             "slip_per_geometry_max": g1.groupby("bin").n_slip.max().reindex(bins.bin).to_numpy()})
        real.to_csv(os.path.join(out, f"realisations_{zone}.csv"), index=False)
        np.savez_compressed(os.path.join(out, f"scenarios_{zone}.npz"), **{c: s[c].to_numpy() for c in s.columns})
        rep.check(f"{zone}: realisations per position and geometry", "Uneven numbers of realisations (issue H3); the split uses the counts that exist", None,
                  f"positions {len(p1):,}: {100 * (p1.n_sreal == 1).mean():.0f}% have 1 scaling sample, {100 * (p1.n_sreal == 10).mean():.0f}% have 10; "
                  f"geometries {len(g1):,}: {100 * (g1.n_slip == 1).mean():.0f}% have 1 slip realisation, {100 * (g1.n_slip == 10).mean():.0f}% have 10")
    rule = cfg["tree"]["options"]["bin_edges"]
    rep.note("Bin edges from config/magnitude_intervals.csv, the intervals of the scenario provider (decision 22)." if rule == "provider"
             else "Bin edges are midpoints between bin centres (decision 5).")
    rep.note("The triangle-to-segment map is derived from the scenarios, not supplied; it can be wrong by one segment next to a boundary.")
    rep.note("scenarios_<zone>.npz is read with numpy.load; columns: rig, mw (x 1e4), scal, sreal, tri, slip, slip_ratio, total_M0, "
             "seg_1.., bin, position, geometry, n_slip, n_sreal, split.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda"]) else 1)
