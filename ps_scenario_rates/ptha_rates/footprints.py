"""Rupture footprints: the triangles that each rupture position covers.

A footprint table is a sparse matrix position x triangle held as three arrays:
    indptr   start of each position's entries (length = positions + 1)
    index    row of the triangle in the triangle table (0-based)
    frac     fraction of the rupture's area (or, from slip files with slip, of its moment) on that triangle; sums to 1 for every position
Two sources (option `footprints` in config/logic_tree.yaml):
    estimated   grown around the initiation triangle until the area of the scaling law is reached
    from_file   read from the files written by tools/reduce_slip_files.py from the slip files
With `segment_fraction_approximation` (the default) no footprint table is used.
"""
import glob
import os
import numpy as np

R_EARTH = 6371.0


def xyz(tri):
    """Cartesian coordinates (km) of the triangle centroids; straight-line distances between
    neighbours are then independent of latitude."""
    lon, lat, r = np.deg2rad(tri.lon.to_numpy()), np.deg2rad(tri.lat.to_numpy()), R_EARTH - tri.depth_km.to_numpy()
    return np.column_stack([r * np.cos(lat) * np.cos(lon), r * np.cos(lat) * np.sin(lon), r * np.sin(lat)])


def law_area(mw, coef):
    """Median rupture area (km2) of a scaling law: log10(A) = a + b Mw."""
    return 10.0 ** (coef[0] + coef[1] * np.asarray(mw, float))


def estimate(tri, pos_tri, pos_area):
    """Footprint of every position: the triangles nearest to its initiation triangle, taken in order of
    distance until their summed area reaches the target area (at least one triangle; at most the mesh).
    tri: triangle table; pos_tri: initiation triangle id per position; pos_area: target area (km2) per position."""
    P, area = xyz(tri), tri.area_km2.to_numpy()
    row_of = np.full(int(tri.id.max()) + 1, -1); row_of[tri.id.to_numpy()] = np.arange(len(tri))
    start = row_of[np.asarray(pos_tri)]
    if (start < 0).any():
        raise ValueError("an initiation triangle is not in the mesh")
    n_pos = len(start); count = np.zeros(n_pos, np.int64); chunks = [None] * n_pos
    order_pos = np.argsort(start, kind="stable"); bounds = np.flatnonzero(np.r_[True, np.diff(start[order_pos]) > 0, True])
    for a, b in zip(bounds[:-1], bounds[1:]):
        members = order_pos[a:b]; i = start[members[0]]
        d = ((P - P[i]) ** 2).sum(axis=1)
        need = pos_area[members].max()
        k = min(len(area), max(8, int(2.5 * need / area.mean())))          # enough neighbours, then sort only those
        while True:
            near = np.argpartition(d, k - 1)[:k] if k < len(area) else np.arange(len(area))
            near = near[np.argsort(d[near], kind="stable")]; cum = np.cumsum(area[near])
            if cum[-1] >= need or k >= len(area):
                break
            k = min(len(area), 2 * k)
        n = np.minimum(np.searchsorted(cum, pos_area[members]) + 1, len(near))
        for p, m in zip(members, n):
            chunks[p] = (near[:m], area[near[:m]] / cum[m - 1]); count[p] = m
    indptr = np.r_[0, np.cumsum(count)]
    return dict(indptr=indptr, index=np.concatenate([c[0] for c in chunks]).astype(np.int16 if len(tri) < 32768 else np.int32),
                frac=np.concatenate([c[1] for c in chunks]).astype(np.float32))


def row_sums(fp, tri_values):
    """For every position: sum over its triangles of frac x value (the area-weighted mean of a triangle
    quantity under the rupture). tri_values: one value per triangle, or triangles x k."""
    v = np.asarray(tri_values, float)[fp["index"]]
    w = fp["frac"].astype(float) if v.ndim == 1 else fp["frac"].astype(float)[:, None]
    return np.add.reduceat(w * v, fp["indptr"][:-1], axis=0)


def segment_fractions(fp, tri_segment, n_seg):
    """Fraction of every rupture's area in each segment (position x segment)."""
    onehot = (np.asarray(tri_segment)[:, None] == np.arange(1, n_seg + 1)[None, :]).astype(float)
    return row_sums(fp, onehot)


def area_km2(fp, tri):
    """Footprint area of every position."""
    a = tri.area_km2.to_numpy()[fp["index"]]
    return np.add.reduceat(a, fp["indptr"][:-1])


def from_file(folder, zone_code, scen, tri, n_pos, mu_of=None):
    """Footprints of every position from the files written by tools/reduce_slip_files.py (one or more
    footprints_<zone>_*.npz in `folder`). Each file lists rupture geometries by scenario name with their
    triangles. A position's footprint is the mean over its scenarios, weighted by their split.
    scen: the scenario table of step B1 (rig, mw, scal, sreal, tri, position, split).
    When the files hold slip (tools/reduce_slip_files.py --slip all) and mu_of is given, frac is the moment fraction of
    each triangle, mu x area x mean slip (decision 27); otherwise the area fraction. mu_of(rig, rows) returns the rigidity
    of the triangles in the given rows of tri for rigidity model rig (int)."""
    files = sorted(glob.glob(os.path.join(folder, f"footprints_{zone_code}_*.npz")))
    if not files:
        raise FileNotFoundError(f"no footprints_{zone_code}_*.npz in {folder}; see data/README.md")
    key = lambda rig, mw, scal, sreal, t: ((np.asarray(rig, np.int64) * 100000 + mw) * 100 + scal) * 100 * 100000 + np.asarray(sreal, np.int64) * 100000 + t
    k_scen = key(scen["rig"], scen["mw"], scen["scal"], scen["sreal"], scen["tri"])
    o = np.argsort(k_scen, kind="stable"); ks = k_scen[o]
    w_geo_pos = np.add.reduceat(scen["split"][o], np.flatnonzero(np.r_[True, np.diff(ks) > 0]))     # weight of each geometry in its position
    k_u = ks[np.r_[True, np.diff(ks) > 0]]; pos_u = scen["position"][o][np.r_[True, np.diff(ks) > 0]]
    row_of = np.full(int(tri.id.max()) + 1, -1); row_of[tri.id.to_numpy()] = np.arange(len(tri)); area = tri.area_km2.to_numpy()
    P, T, W = [], [], []; seen = np.zeros(len(k_u), bool)
    for f in files:
        d = np.load(f)
        k_f = key(d["rig"], d["mw"], d["scal"], d["sreal"], d["tri"])
        j = np.searchsorted(k_u, k_f); ok = (j < len(k_u)) & (k_u[np.minimum(j, len(k_u) - 1)] == k_f)
        ip = d["indptr"]; has_slip = "slip" in d.files and mu_of is not None
        for g in np.flatnonzero(ok):
            if seen[j[g]]:
                continue
            seen[j[g]] = True
            t = row_of[d["triangle"][ip[g]:ip[g + 1]]]; a = area[t]
            if has_slip:
                a = a * d["slip"][ip[g]:ip[g + 1]] * mu_of(int(d["rig"][g]), t)
            P.append(np.full(len(t), pos_u[j[g]])); T.append(t); W.append(w_geo_pos[j[g]] * a / a.sum())
    if not seen.all():
        raise ValueError(f"{int((~seen).sum()):,} of {len(seen):,} rupture geometries have no footprint in {folder}")
    P, T, W = np.concatenate(P), np.concatenate(T), np.concatenate(W)
    code = P.astype(np.int64) * len(tri) + T                         # merge the geometries of a position
    u, inv = np.unique(code, return_inverse=True); frac = np.bincount(inv, weights=W)
    p_u, t_u = u // len(tri), u % len(tri)
    return dict(indptr=np.r_[0, np.cumsum(np.bincount(p_u, minlength=n_pos))], index=t_u.astype(np.int16 if len(tri) < 32768 else np.int32),
                frac=frac.astype(np.float32))


def save(path, fp):
    np.savez_compressed(path, **fp)


def load(path):
    d = np.load(path)
    return {k: d[k] for k in ("indptr", "index", "frac")}
