"""Step B3. Catalogue: events by segment, largest observed events and the observed rate that
anchors the seismic branches."""
import os
import numpy as np
import pandas as pd
from . import DATA

KM_PER_DEG = 111.19
THRESHOLDS = [7.5, 8.0, 8.5, 9.0]


def load_catalogue(zone_cfg, tri, keep_excluded=False):
    """Catalogue with the segment of each event (segment of the nearest triangle) and the
    distance to that triangle in km. Events removed by the zone's catalogue_filter are dropped, or kept
    with the reason in column 'excluded' if keep_excluded."""
    c = pd.read_csv(os.path.join(DATA, zone_cfg["catalogue"]))
    c = c.rename(columns={"YEAR": "year", "MONTH": "month", "DAY": "day", "LAT": "lat", "LON": "lon", "DEPTH": "depth_km",
                          "Unified_Mw": "mw", "Mw": "mw"})
    c = c[(c.mw > 0) & (c.year <= zone_cfg["catalogue_end_year"])].reset_index(drop=True)   # rows without a magnitude (Makran: Mw blank or 0)
    tlon, tlat = tri.lon.to_numpy()[None], tri.lat.to_numpy()[None]
    j, d = np.empty(len(c), int), np.empty(len(c))
    for i in range(0, len(c), 2000):
        la, lo = c.lat.to_numpy()[i:i + 2000, None], c.lon.to_numpy()[i:i + 2000, None]
        dd = KM_PER_DEG * np.hypot((lo - tlon) * np.cos(np.deg2rad(la)), la - tlat)
        j[i:i + 2000], d[i:i + 2000] = dd.argmin(1), dd.min(1)
    c["segment"] = tri.segment.to_numpy()[j]
    c["dist_km"] = d
    c["excluded"] = excluded(c, zone_cfg.get("catalogue_filter"))
    cols = ["year", "lat", "lon", "depth_km", "mw", "segment", "dist_km"]
    return c[cols + ["excluded"]] if keep_excluded else c[c.excluded == ""][cols].reset_index(drop=True)


def excluded(c, flt):
    """Events removed by the zone's catalogue_filter (zones.yaml), with the reason. The filter keeps interface
    events only, as far as location allows: no focal mechanisms are in the catalogue."""
    why = pd.Series("", index=c.index)
    if not flt:
        return why
    if flt.get("lon_min") is not None:
        why[c.lon < flt["lon_min"]] = f"west of {flt['lon_min']} E"
    if flt.get("max_dist_km") is not None:
        why[(why == "") & (c.dist_km > flt["max_dist_km"])] = f"more than {flt['max_dist_km']} km from the mesh"
    for e in flt.get("exclude", []):
        t = c.year * 10000 + c.month * 100 + c.day
        m = (t >= int(e["from"].replace("-", ""))) & (t <= int(e["to"].replace("-", ""))) & \
            (np.hypot(c.lat - e["lat"], c.lon - e["lon"]) <= e["radius_deg"]) & (c.mw >= e.get("mw_min", 0))
        why[(why == "") & m] = e["name"]
    return why


def units(zone_cfg):
    """Whole zone first, then the segments. Each: name, segment number (0 = whole zone), parameters."""
    out = [dict(unit="whole_zone", segment=0, **zone_cfg["whole_zone"])]
    for k, s in enumerate(zone_cfg["segments"], start=1):
        out.append(dict(unit=s["name"], segment=k, **{f: s[f] for f in ("b", "mc", "mc_start_year", "mmax_area", "mmax_observed")}))
    return out


def observed_mmax(largest, area_limit, tree=None):
    """The 'observed' Mmax alternative: the largest event plus n_sigma x sigma (decision 23; version 1.0 had
    n_sigma 0), kept at or below the area limit."""
    o = (tree or {}).get("options", {}).get("mmax_observed", {"n_sigma": 0, "sigma": 0})
    return float(min(largest + o["n_sigma"] * o["sigma"], area_limit))


def parameter_table(cat, zone_cfg, tree=None):
    """Per unit: catalogue b; the observed annual rate at the anchor magnitude, which sets the a-value
    of the seismic branches (decision 20); the three Mmax alternatives (largest observed, area limit,
    mean of the two); and, for the record, Mc with the rate above it in its complete period."""
    end = zone_cfg["catalogue_end_year"]
    anc = zone_cfg["seismic_anchor"]
    anc_years = end - anc["start_year"] + 1
    rows = []
    for u in units(zone_cfg):
        e = cat if u["segment"] == 0 else cat[cat.segment == u["segment"]]
        years = end - u["mc_start_year"] + 1
        n = int(((e.mw >= u["mc"]) & (e.year >= u["mc_start_year"])).sum())
        largest = u["mmax_observed"] if u["mmax_observed"] is not None else float(e.mw.max())
        obs = observed_mmax(largest, u["mmax_area"], tree)
        n_anc = int(((e.mw >= anc["magnitude"]) & (e.year >= anc["start_year"])).sum())
        r = dict(unit=u["unit"], segment=u["segment"], b_catalogue=u["b"],
                 anchor_mw=anc["magnitude"], anchor_start_year=anc["start_year"], anchor_years=anc_years,
                 n_ge_anchor=n_anc, rate_anchor_per_yr=n_anc / anc_years,
                 mc=u["mc"], mc_start_year=u["mc_start_year"],
                 years_complete=years, n_ge_mc=n, rate_mc_per_yr=n / years, n_events=len(e),
                 largest_event_mw=largest, mmax_observed=obs, mmax_area_limit=u["mmax_area"], mmax_mean=(obs + u["mmax_area"]) / 2)
        for t in THRESHOLDS:
            r[f"n_ge_{t:g}"] = int((e.mw >= t).sum())
        rows.append(r)
    return pd.DataFrame(rows)
