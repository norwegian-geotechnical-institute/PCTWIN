"""Step B1. Freeze the inputs into checked tables: triangles, magnitude bins, scenarios."""
import os
import zipfile
import numpy as np
import pandas as pd
from . import DATA

KM_PER_DEG = 111.19


def fractions_handle(zone_cfg):
    """Open the zone's fractions csv inside Rupture_fractions_v2.zip in the data folder."""
    z = zipfile.ZipFile(os.path.join(DATA, "Rupture_fractions_v2.zip"))
    name = [n for n in z.namelist() if os.path.basename(n) == zone_cfg["fractions_csv"]][0]
    return z.open(name)


def triangle_table(zone_cfg):
    """id, centroid (lon, lat, depth km), area (km2) and dip (deg) of every triangle.
    Area and dip come from the three vertices on a local flat-earth projection."""
    a = np.loadtxt(os.path.join(DATA, zone_cfg["triangles"]))
    v = a[:, 1:].reshape(-1, 3, 3)
    lat0 = np.deg2rad(v[:, :, 1].mean(axis=1))[:, None]
    P = np.stack([v[:, :, 0] * KM_PER_DEG * np.cos(lat0), v[:, :, 1] * KM_PER_DEG, np.abs(v[:, :, 2])], axis=2)
    n = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
    norm = np.linalg.norm(n, axis=1)
    return pd.DataFrame({"id": a[:, 0].astype(int), "lon": v[:, :, 0].mean(1), "lat": v[:, :, 1].mean(1),
                         "depth_km": np.abs(v[:, :, 2]).mean(1), "area_km2": 0.5 * norm,
                         "dip_deg": np.degrees(np.arccos(np.abs(n[:, 2]) / norm))})


def parse_names(fn):
    """Split scenario names PS_<zone>_<rig>_<Mw int>_<Mw dec>_<scal>_<sreal>_<tri>_<slip>."""
    return pd.DataFrame({"rig": fn.str.slice(7, 9).astype("int8"),
                         "mw": (fn.str.slice(10, 11) + fn.str.slice(12, 16)).astype("int32"),   # Mw x 1e4
                         "scal": fn.str.slice(17, 19).astype("int8"),
                         "sreal": fn.str.slice(20, 22).astype("int8"),
                         "tri": fn.str.slice(23, 28).astype("int32"),
                         "slip": fn.str.slice(29, 32).astype("int16")}, index=fn.index)


def read_scenarios(zone_cfg, rigidity_models, chunk=1_000_000):
    """Scenarios of the given rigidity models (list of ints) from the v2 fractions file.
    Returns the table and the number of rows read in total."""
    parts, n_all = [], 0
    with fractions_handle(zone_cfg) as fh:
        for c in pd.read_csv(fh, chunksize=chunk):
            n_all += len(c)
            c = c[c.filename.str.slice(7, 9).astype("int8").isin(rigidity_models)]
            if len(c) == 0:
                continue
            fn = c.pop("filename")
            parts.append(pd.concat([parse_names(fn), c], axis=1))
    s = pd.concat(parts, ignore_index=True)
    seg = [k for k in s.columns if k.startswith("seg_")]
    s[seg] = s[seg].astype("float32")
    return s, n_all


def bin_table(mw_int, rule="midpoints"):
    """Magnitude bins from the bin centres found in the data. Edges are midpoints between
    centres; the first and last bins are made symmetric about their centres."""
    if rule != "midpoints":
        raise ValueError("only 'midpoints' bin edges are implemented")
    c = np.sort(np.unique(mw_int)) / 1e4
    mid = (c[:-1] + c[1:]) / 2
    lo = np.r_[c[0] - (mid[0] - c[0]), mid]
    hi = np.r_[mid, c[-1] + (c[-1] - mid[-1])]
    return pd.DataFrame({"bin": np.arange(len(c)), "mw_centre": c, "mw_lower": lo, "mw_upper": hi,
                         "m0_centre_nm": 10 ** (1.5 * c + 9.1)})


def add_bins_and_counts(s, bins, sample00_counts=True):
    """Add the bin index and the counts needed to split a position's rate over its realisations.

    position  = rigidity model + scaling law + magnitude bin + initiation triangle
    geometry  = position + scaling sample
    n_sreal   = scaling samples that exist for the position
    n_slip    = slip realisations that exist for the geometry
    split     = 1 / (n_sreal x n_slip); sums to 1 over the scenarios of a position
    """
    if not sample00_counts:
        raise NotImplementedError("scaling sample 00 is counted as one sample (decision 10)")
    s["bin"] = np.searchsorted((bins.mw_centre.to_numpy() * 1e4).round().astype("int32"), s.mw.to_numpy()).astype("int8")
    pos = ["rig", "scal", "bin", "tri"]
    s["position"] = s.groupby(pos, sort=False).ngroup().astype("int32")
    s["geometry"] = s.groupby(["position", "sreal"], sort=False).ngroup().astype("int32")
    s["n_slip"] = s.groupby("geometry").slip.transform("size").astype("int16")
    s["n_sreal"] = s.groupby("position").sreal.transform("nunique").astype("int16")
    s["split"] = 1.0 / (s.n_sreal.astype("float64") * s.n_slip)
    return s


def segment_map(s, tri, n_seg):
    """Segment of each triangle. No triangle-to-segment table has been supplied, so each
    triangle takes the dominant segment of the smallest scenarios that start on it;
    triangles on which no such scenario starts take the segment of the nearest triangle
    that has one. Returns (segment per triangle, number assigned from the nearest)."""
    seg_cols = [f"seg_{i + 1}" for i in range(n_seg)]
    small = s[s.bin == s.bin.min()]
    dom = small.groupby("tri")[seg_cols].mean()
    seg = pd.Series(0, index=tri.id.to_numpy(), dtype="int8")
    seg.loc[dom.index] = (dom.to_numpy().argmax(axis=1) + 1).astype("int8")
    have = seg.to_numpy() > 0
    lon, lat = tri.lon.to_numpy(), tri.lat.to_numpy()
    out = seg.to_numpy().copy()
    for i in np.where(~have)[0]:
        j = np.where(have)[0][np.argmin(((lon[have] - lon[i]) * np.cos(np.deg2rad(lat[i]))) ** 2 + (lat[have] - lat[i]) ** 2)]
        out[i] = out[j]
    return out, int((~have).sum())
