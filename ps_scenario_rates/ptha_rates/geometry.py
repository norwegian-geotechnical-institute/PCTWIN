"""Distance along the trench, used to plot rates along the arc.

The trench line is traced through the shallowest triangles of the mesh. They are put in order
along the arc by their bearing from a point on the concave side of the arc (arc_centre in
zones.yaml; it is used for ordering only). Every triangle is then given the distance along that
line of the nearest point of the line, so a triangle down-dip of the trench has the same
coordinate as the trench point it faces.
"""
import numpy as np
import pandas as pd

KM_PER_DEG = 111.19


def trench_line(tri, arc_centre, shallow_km=10.0, n=400, smooth=7):
    """Polyline along the trench: DataFrame lon, lat, dist_km (0 at the first end met going
    anticlockwise from the centre's point of view, i.e. the north-western end for Sunda)."""
    s = tri[tri.depth_km <= shallow_km]
    ang = np.unwrap(np.arctan2(s.lat.to_numpy() - arc_centre[1], (s.lon.to_numpy() - arc_centre[0]) * np.cos(np.deg2rad(arc_centre[1]))))
    ang = np.mod(ang, 2 * np.pi)
    edges = np.linspace(ang.min(), ang.max() + 1e-9, n + 1)
    k = np.digitize(ang, edges) - 1
    g = pd.DataFrame({"k": k, "lon": s.lon.to_numpy(), "lat": s.lat.to_numpy(), "d": s.depth_km.to_numpy()})
    # in each bearing bin, the trench is at the shallowest triangles
    g = g[g.d <= g.groupby("k").d.transform("min") + 2.0].groupby("k")[["lon", "lat"]].mean().sort_index()
    line = g.rolling(smooth, center=True, min_periods=1).mean().reset_index(drop=True)
    dx = np.diff(line.lon) * KM_PER_DEG * np.cos(np.deg2rad((line.lat.to_numpy()[1:] + line.lat.to_numpy()[:-1]) / 2))
    dy = np.diff(line.lat) * KM_PER_DEG
    line["dist_km"] = np.r_[0.0, np.cumsum(np.hypot(dx, dy))]
    return line


def along_trench(lon, lat, line, step_km=5.0):
    """Distance along the trench line of the point of the line nearest to each (lon, lat)."""
    d = np.arange(0.0, line.dist_km.iloc[-1] + step_km, step_km)
    llon, llat = np.interp(d, line.dist_km, line.lon), np.interp(d, line.dist_km, line.lat)
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    out = np.empty(len(lon))
    for i in range(0, len(lon), 4000):
        la, lo = lat[i:i + 4000, None], lon[i:i + 4000, None]
        dd = np.hypot((lo - llon[None]) * np.cos(np.deg2rad(la)), la - llat[None])
        out[i:i + 4000] = d[dd.argmin(axis=1)]
    return out
