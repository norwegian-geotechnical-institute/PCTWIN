"""Rigidity of the interface for each rigidity model of the scenarios (decision 27).

  00  homogeneous: the zone's rigidity_pa (30 GPa)
  01  theoretical: mu = 10^(0.5631 + 0.0437 z) GPa (z in km). Where this is below the PREM value it is raised
      towards PREM by fact_mu x (PREM - mu); below 10 km depth it is capped at PREM. PREM is read between its
      tabulated depths by interpolation linear in log10(mu). Same as assign_rigidity() of the scenario provider.
"""
import numpy as np

PREM_MU = np.array([[0, 26.5], [3, 26.5], [15, 43.9], [24.4, 67.8], [40, 67.6], [60, 67.3], [80, 64.6], [115, 64.2]])   # km, GPa
A, B = 0.5631, 0.0437
KINDS = ["homogeneous", "theoretical"]


def prem(depth_km):
    """PREM rigidity (GPa), log-linear between tabulated depths; 0 outside the table (as in the provider's code)."""
    z = np.asarray(depth_km, float); out = np.zeros(z.shape)
    d, m = PREM_MU[:, 0], np.log10(PREM_MU[:, 1])
    j = np.clip(np.searchsorted(d, z, side="right") - 1, 0, len(d) - 2)
    inside = (z >= d[0]) & (z <= d[-1])
    v = 10 ** (m[j] + (m[j + 1] - m[j]) / (d[j + 1] - d[j]) * (z - d[j]))
    out[inside] = v[inside]
    return out


def theoretical(depth_km, fact_mu=0.5):
    """Model 01 rigidity (GPa) at the given depths (km)."""
    z = np.asarray(depth_km, float)
    mu = 10 ** (A + B * z); p = prem(z)
    return np.where(mu >= p, np.where(z > 10, p, mu), mu + fact_mu * (p - mu))


def mu_pa(model, depth_km, zone_cfg, tree):
    """Rigidity (Pa) of one rigidity model ("00", "01") at the given depths."""
    r = tree["options"]["rigidity"]; kind = r["models"][model]
    if kind == "homogeneous":
        return np.full(np.shape(depth_km), float(zone_cfg["rigidity_pa"]))
    if kind == "theoretical":
        return theoretical(depth_km, r.get("fact_mu", 0.5)) * 1e9
    raise ValueError(f"rigidity model {model}: kind must be one of {KINDS}")
