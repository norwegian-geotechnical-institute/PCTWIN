"""Triangle areas as the scenario provider computes them (ANTI-FASc; issue H1, decision 29).

The provider projects the vertices of every triangle to one UTM zone for the whole zone (Sunda: zone 46, EPSG 32646
north of the equator and 32746 south of it; Makran: zone 41) and takes the area of the triangle in (easting, northing,
-depth). Far from the central meridian (93 E for zone 46) the projection stretches distances by up to about 10% (east
of 115 E), so these areas are larger than the true ones. The slip of each scenario was set from its named moment with
these areas, so its true moment is smaller than named by the same factor.

utm_xy is the forward transverse Mercator projection on WGS84 (Krueger series to fourth order); it agrees with pyproj
(EPSG 32646 / 32746) to a few micrometres over the Sunda arc.
"""
import numpy as np

A_WGS, F_WGS = 6378137.0, 1 / 298.257223563


def utm_xy(lon, lat, zone, south=False, k0=0.9996):
    """WGS84 longitude, latitude (degrees) -> UTM easting, northing (m) in the given zone."""
    lon, lat = np.asarray(lon, float), np.asarray(lat, float)
    n = F_WGS / (2 - F_WGS); A = A_WGS / (1 + n) * (1 + n ** 2 / 4 + n ** 4 / 64)
    al = [n / 2 - 2 * n ** 2 / 3 + 5 * n ** 3 / 16 + 41 * n ** 4 / 180, 13 * n ** 2 / 48 - 3 * n ** 3 / 5 + 557 * n ** 4 / 1440,
          61 * n ** 3 / 240 - 103 * n ** 4 / 140, 49561 * n ** 4 / 161280]
    lam = np.deg2rad(lon - (6 * zone - 183)); phi = np.deg2rad(lat)
    s = 2 * np.sqrt(n) / (1 + n)
    t = np.sinh(np.arctanh(np.sin(phi)) - s * np.arctanh(s * np.sin(phi)))
    xi = np.arctan2(t, np.cos(lam)); eta = np.arctanh(np.sin(lam) / np.sqrt(1 + t ** 2))
    x = eta + sum(a * np.cos(2 * j * xi) * np.sinh(2 * j * eta) for j, a in enumerate(al, 1))
    y = xi + sum(a * np.sin(2 * j * xi) * np.cosh(2 * j * eta) for j, a in enumerate(al, 1))
    return 500000 + k0 * A * x, (10000000 if south else 0) + k0 * A * y


def triangle_area_utm_km2(v, zone):
    """Area (km2) of triangles given as an array (n, 3, 3) of (lon, lat, depth km) per vertex, in the provider's way:
    vertices to UTM (x, y) in metres and -depth in metres, area from the cross product. North or south (false northing)
    does not change the area."""
    x, y = utm_xy(v[:, :, 0], v[:, :, 1], zone)
    P = np.stack([x, y, -np.abs(v[:, :, 2]) * 1e3], axis=2)
    return 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1) / 1e6
