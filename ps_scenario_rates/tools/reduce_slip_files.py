"""Reduce slip files to rupture footprints: which triangles each rupture covers and, by default, the mean slip on each.

To be run where the slip files are. Needs only Python 3 and numpy; it does not need the rest of this repository.

    python reduce_slip_files.py SLIP_FOLDER MESH_FILE ZONE OUTPUT_FOLDER [--rig 00] [--label NAME]

    SLIP_FOLDER    folder searched, with its subfolders, for slip files  PS_<zone>_<rig>_<Mw>_<scal>_<sreal>_<tri>_<slip>.dat
                   (one header line, then one row per ruptured triangle: lon, lat, depth of the three vertices, rake, slip)
    MESH_FILE      the mesh as text: triangle id and lon, lat, depth of its three vertices (e.g. Sunda_triangles.txt)
    ZONE           zone code used in the file names (sun, mak)
    OUTPUT_FOLDER  where footprints_<zone>_<label>.npz is written
    --rig          only this rigidity model (two digits); default: all
    --label        name of the output file; default: the rigidity model, or "all"
    --slip         all (default): read every slip realisation and keep the mean slip on each triangle of a rupture geometry
                   (needed for rigidity models other than 00, decision 27); none: triangles only, first realisation read

A rupture geometry is rigidity model, magnitude, scaling law, scaling sample and initiation triangle; its slip realisations
(the last number of the name) differ in slip, not in the area they cover. With --slip all the triangles of a geometry are
the union over its realisations and the slip on a triangle is the mean over the realisations (0 where a realisation has none).
It can be run several times on parts of the scenario set (for example one magnitude folder at a time, with a
different --label each time); all files footprints_<zone>_*.npz are read together afterwards.

The output goes into  ps_scenario_rates/data/footprints/  and is used with  footprints: from_file  in
config/logic_tree.yaml. Arrays: rig, mw (Mw x 10000), scal, sreal, tri (one entry per rupture geometry);
indptr and triangle (ids of the triangles of each geometry: triangle[indptr[i]:indptr[i+1]]); with --slip all also
slip (m, aligned with triangle) and n_real (realisations read per geometry).

Check printed at the end (models 00 and 01): the moment of each geometry, sum of mu x area x slip, against the moment of
its named Mw. Areas are computed as the scenario generator does (vertices on one UTM zone, --utm, default 46 for sun and 41
for mak) and mu is that of the rigidity model (00: 30 GPa; 01: the depth-dependent function with Fact_mu 0.5). The ratio
should be 1 within about 1%; anything else means the files or the mesh are not the ones expected.
"""
import os
import sys
import numpy as np

TOL = 1e-4          # degrees / km: agreement required between a slip-file triangle and a mesh triangle
UTM_ZONE = {"sun": 46, "mak": 41}
PREM_MU = np.array([[0, 26.5], [3, 26.5], [15, 43.9], [24.4, 67.8], [40, 67.6], [60, 67.3], [80, 64.6], [115, 64.2]])


def utm_xy(lon, lat, zone, k0=0.9996):
    """WGS84 -> UTM easting, northing (m); as ptha_rates/projection.py (north and south differ only by a constant)."""
    a, f = 6378137.0, 1 / 298.257223563; n = f / (2 - f); A = a / (1 + n) * (1 + n ** 2 / 4 + n ** 4 / 64)
    al = [n / 2 - 2 * n ** 2 / 3 + 5 * n ** 3 / 16 + 41 * n ** 4 / 180, 13 * n ** 2 / 48 - 3 * n ** 3 / 5 + 557 * n ** 4 / 1440,
          61 * n ** 3 / 240 - 103 * n ** 4 / 140, 49561 * n ** 4 / 161280]
    lam = np.deg2rad(lon - (6 * zone - 183)); phi = np.deg2rad(lat); s = 2 * np.sqrt(n) / (1 + n)
    t = np.sinh(np.arctanh(np.sin(phi)) - s * np.arctanh(s * np.sin(phi)))
    xi = np.arctan2(t, np.cos(lam)); eta = np.arctanh(np.sin(lam) / np.sqrt(1 + t ** 2))
    x = eta + sum(c * np.cos(2 * j * xi) * np.sinh(2 * j * eta) for j, c in enumerate(al, 1))
    y = xi + sum(c * np.sin(2 * j * xi) * np.cosh(2 * j * eta) for j, c in enumerate(al, 1))
    return k0 * A * x, k0 * A * y


def mu_pa(rig, z):
    """Rigidity (Pa) at depths z (km): 00 homogeneous 30 GPa; 01 as ptha_rates/rigidity.py; None for other models."""
    if rig == 0:
        return np.full(len(z), 3e10)
    if rig != 1:
        return None
    d, m = PREM_MU[:, 0], np.log10(PREM_MU[:, 1]); j = np.clip(np.searchsorted(d, z, side="right") - 1, 0, len(d) - 2)
    p = np.where((z >= d[0]) & (z <= d[-1]), 10 ** (m[j] + (m[j + 1] - m[j]) / (d[j + 1] - d[j]) * (z - d[j])), 0.0)
    mu = 10 ** (0.5631 + 0.0437 * z)
    return np.where(mu >= p, np.where(z > 10, p, mu), mu + 0.5 * (p - mu)) * 1e9


def parse(name):
    """PS_sun_00_6_8012_01_03_05747_002 -> (rig, mw x 1e4, scal, sreal, tri)."""
    p = name.split("_")
    return int(p[2]), int(p[3]) * 10000 + int(p[4]), int(p[5]), int(p[6]), int(p[7])


def main(argv):
    args, opt, i = [], {}, 0
    while i < len(argv):
        if argv[i].startswith("--") and i + 1 < len(argv):
            opt[argv[i][2:]] = argv[i + 1]; i += 2
        else:
            args.append(argv[i]); i += 1
    if len(args) != 4:
        print(__doc__); return 1
    folder, mesh_file, zone, out = args
    rig_only = int(opt["rig"]) if "rig" in opt else None
    label = opt.get("label", f"{rig_only:02d}" if rig_only is not None else "all")
    want_slip = opt.get("slip", "all") == "all"
    utm = int(opt.get("utm", UTM_ZONE.get(zone, 46)))
    mesh = np.loadtxt(mesh_file); ids = mesh[:, 0].astype(np.int64); mc = mesh[:, 1:10].reshape(-1, 3, 3).mean(axis=1)
    order = np.argsort(mc[:, 0]); mcs = mc[order]

    def match(a9):
        q = a9.reshape(-1, 3, 3).mean(axis=1); lo = np.searchsorted(mcs[:, 0], q[:, 0] - TOL); hit_row = np.full(len(q), -1)
        for k in range(8):
            c = np.clip(lo + k, 0, len(mcs) - 1); hit = (np.abs(mcs[c] - q) < TOL).all(axis=1) & (hit_row < 0); hit_row[hit] = order[c[hit]]
        return hit_row

    files = []
    for root, _, names in os.walk(folder):
        files += [os.path.join(root, n) for n in names if n.startswith(f"PS_{zone}_") and n.endswith(".dat")]
    files.sort()
    groups = {}
    for f in files:
        k = parse(os.path.basename(f)[:-4])
        if rig_only is not None and k[0] != rig_only:
            continue
        if k not in groups:
            groups[k] = [f]
        elif want_slip:
            groups[k].append(f)
    keys, tris, slips, nreal, ratio, n_unmatched = [], [], [], [], [], 0
    for i, (k, fs) in enumerate(sorted(groups.items())):
        ids_all, slip_all, m0 = [], [], 0.0
        mu_ok = mu_pa(k[0], np.zeros(1)) is not None
        for f in fs:
            a = np.loadtxt(f, skiprows=1, ndmin=2)
            j = match(a[:, :9]); n_unmatched += int((j < 0).sum()); ok = j >= 0
            ids_all.append(ids[j[ok]])
            if want_slip:
                s_ = a[ok, 10]; slip_all.append(s_)
                if mu_ok:
                    v = a[ok, :9].reshape(-1, 3, 3); x, y = utm_xy(v[:, :, 0], v[:, :, 1], utm)
                    P = np.stack([x, y, -np.abs(v[:, :, 2]) * 1e3], axis=2)
                    area = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1)
                    m0 += (mu_pa(k[0], np.abs(v[:, :, 2]).mean(1)) * area * s_).sum()
        u, inv = np.unique(np.concatenate(ids_all), return_inverse=True)
        keys.append(k); tris.append(u.astype(np.int32)); nreal.append(len(fs))
        if want_slip:
            slips.append((np.bincount(inv, weights=np.concatenate(slip_all), minlength=len(u)) / len(fs)).astype(np.float32))
            if mu_ok:
                ratio.append(m0 / len(fs) / 10 ** (1.5 * k[1] / 1e4 + 9.1))
        if (i + 1) % 5000 == 0:
            print(f"{i + 1:,} of {len(groups):,} rupture geometries read", flush=True)
    if not keys:
        print("no slip files found"); return 1
    K = np.array(keys, dtype=np.int64)
    os.makedirs(out, exist_ok=True); path = os.path.join(out, f"footprints_{zone}_{label}.npz")
    extra = dict(slip=np.concatenate(slips), n_real=np.array(nreal, np.int16)) if want_slip else {}
    np.savez_compressed(path, rig=K[:, 0].astype(np.int8), mw=K[:, 1].astype(np.int32), scal=K[:, 2].astype(np.int8), sreal=K[:, 3].astype(np.int8),
                        tri=K[:, 4].astype(np.int32), indptr=np.r_[0, np.cumsum([len(t) for t in tris])].astype(np.int64), triangle=np.concatenate(tris), **extra)
    print(f"{len(keys):,} rupture geometries from {sum(len(f) for f in groups.values()):,} files written to {path}")
    if want_slip and ratio:
        q = np.quantile(ratio, [0.05, 0.5, 0.95]); ok = 0.98 < q[0] and q[2] < 1.02
        print(f"check: moment (mu x UTM-{utm} area x slip) / moment of the named Mw, {len(ratio):,} geometries: "
              f"5% {q[0]:.3f}, median {q[1]:.3f}, 95% {q[2]:.3f}: " + ("OK" if ok else "NOT within 2% of 1; check the files and the mesh"))
    if n_unmatched:
        print(f"WARNING: {n_unmatched:,} slip-file triangles were not found in the mesh (wrong mesh file?)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
