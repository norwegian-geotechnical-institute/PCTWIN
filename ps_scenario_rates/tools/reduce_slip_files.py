"""Reduce slip files to rupture footprints: which triangles each rupture covers.

To be run where the slip files are. Needs only Python 3 and numpy; it does not need the rest of this repository.

    python reduce_slip_files.py SLIP_FOLDER MESH_FILE ZONE OUTPUT_FOLDER [--rig 00] [--label NAME]

    SLIP_FOLDER    folder searched, with its subfolders, for slip files  PS_<zone>_<rig>_<Mw>_<scal>_<sreal>_<tri>_<slip>.dat
                   (one header line, then one row per ruptured triangle: lon, lat, depth of the three vertices, rake, slip)
    MESH_FILE      the mesh as text: triangle id and lon, lat, depth of its three vertices (e.g. Sunda_triangles.txt)
    ZONE           zone code used in the file names (sun, mak)
    OUTPUT_FOLDER  where footprints_<zone>_<label>.npz is written
    --rig          only this rigidity model (two digits); default: all
    --label        name of the output file; default: the rigidity model, or "all"

The slip realisations of a rupture (the last number of the name) share one footprint, so only the first file of each
rupture geometry (rigidity model, magnitude, scaling law, scaling sample, initiation triangle) is read.
It can be run several times on parts of the scenario set (for example one magnitude folder at a time, with a
different --label each time); all files footprints_<zone>_*.npz are read together afterwards.

The output goes into  ps_scenario_rates/data/footprints/  and is used with  footprints: from_file  in
config/logic_tree.yaml. Arrays: rig, mw (Mw x 10000), scal, sreal, tri (one entry per rupture geometry);
indptr and triangle (ids of the triangles of each geometry: triangle[indptr[i]:indptr[i+1]]).
"""
import os
import sys
import numpy as np

TOL = 1e-4          # degrees / km: agreement required between a slip-file triangle and a mesh triangle


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
    seen, keys, tris, n_unmatched = set(), [], [], 0
    for i, f in enumerate(files):
        k = parse(os.path.basename(f)[:-4])
        if (rig_only is not None and k[0] != rig_only) or k in seen:
            continue
        seen.add(k)
        a = np.loadtxt(f, skiprows=1, usecols=range(9), ndmin=2)
        j = match(a); n_unmatched += int((j < 0).sum())
        keys.append(k); tris.append(ids[j[j >= 0]].astype(np.int32))
        if len(keys) % 5000 == 0:
            print(f"{len(keys):,} ruptures read ({i + 1:,} of {len(files):,} files)", flush=True)
    if not keys:
        print("no slip files found"); return 1
    K = np.array(keys, dtype=np.int64)
    os.makedirs(out, exist_ok=True); path = os.path.join(out, f"footprints_{zone}_{label}.npz")
    np.savez_compressed(path, rig=K[:, 0].astype(np.int8), mw=K[:, 1].astype(np.int32), scal=K[:, 2].astype(np.int8), sreal=K[:, 3].astype(np.int8),
                        tri=K[:, 4].astype(np.int32), indptr=np.r_[0, np.cumsum([len(t) for t in tris])].astype(np.int64), triangle=np.concatenate(tris))
    print(f"{len(keys):,} rupture geometries from {len(files):,} files written to {path}")
    if n_unmatched:
        print(f"WARNING: {n_unmatched:,} slip-file triangles were not found in the mesh (wrong mesh file?)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
