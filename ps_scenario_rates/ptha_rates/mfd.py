"""Step B4. Magnitude-frequency distributions and the rate of each magnitude bin.

Two forms, both written as a cumulative shape S(M), proportional to the number of events of
magnitude >= M:

  truncated Gutenberg-Richter   S(M) = 10^(-bM) - 10^(-b Mmax)                    for M < Mmax, else 0
  tapered Pareto                S(M) = (M0/M0_ref)^(-beta) exp((M0_ref - M0)/M0_c)   beta = 2b/3

with M0 = 10^(1.5 M + 9.1) N m and the corner moment M0_c that of the branch's Mmax (decision 7).
With tapered_tail = truncate_at_last_bin the tapered shape is cut at the upper edge of the last
bin (S(M) - S(M_top)), because there are no scenarios above it.

Bin rate r_k = A x g_k, with g_k = S(lower edge) - S(upper edge).
  seismic branches     A = observed rate at the anchor magnitude / S(anchor magnitude)
  geodynamic branches  A = xi x budget / sum_k(g_k x mean moment of bin k)
xi is the fraction of the moment of the whole distribution that is released above the smallest
modelled magnitude (decision 8).
"""
import numpy as np

M_REF = 6.0        # reference magnitude of the tapered Pareto shape; cancels in every rate
M_LOW = 0.0        # the "whole distribution" is integrated from here; moment below is negligible
DM = 0.002         # magnitude step of the numerical moment integrals


def m0(mw):
    return 10 ** (1.5 * np.asarray(mw, dtype=float) + 9.1)


def _raw(form, m, b, mmax):
    m = np.asarray(m, dtype=float)
    if form == "truncated_gr":
        return 10 ** (-b * m)
    if form == "tapered_pareto":
        return (m0(m) / m0(M_REF)) ** (-2 * b / 3) * np.exp((m0(M_REF) - m0(m)) / m0(mmax))
    raise ValueError(form)


def cut_magnitude(form, mmax, m_top, tail="truncate_at_last_bin"):
    """Magnitude above which the shape is zero."""
    if form == "truncated_gr":
        return mmax
    return m_top if tail == "truncate_at_last_bin" else np.inf


def shape(form, m, b, mmax, m_cut):
    """Cumulative shape S(M): relative number of events of magnitude >= M."""
    m = np.asarray(m, dtype=float)
    tail = 0.0 if np.isinf(m_cut) else _raw(form, m_cut, b, mmax)
    return np.where(m < m_cut, _raw(form, np.minimum(m, m_cut), b, mmax) - tail, 0.0)


def moment_integral(form, m1, m2, b, mmax, m_cut):
    """Sum of moment x number of events between m1 and m2 (relative units)."""
    m2 = min(m2, m_cut) if not np.isinf(m_cut) else m2
    if m2 <= m1:
        return 0.0
    g = np.linspace(m1, m2, max(int(np.ceil((m2 - m1) / DM)), 2) + 1)
    s = shape(form, g, b, mmax, m_cut)
    return float(((s[:-1] - s[1:]) * m0((g[:-1] + g[1:]) / 2)).sum())


def bin_shapes(form, lower, upper, b, mmax, m_cut):
    """g_k and the mean moment of an event in each bin (N m; nan where g_k = 0)."""
    g = shape(form, lower, b, mmax, m_cut) - shape(form, upper, b, mmax, m_cut)
    mom = np.array([moment_integral(form, lo, hi, b, mmax, m_cut) for lo, hi in zip(lower, upper)])
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_m0 = np.where(g > 0, mom / g, np.nan)
    return g, mean_m0, mom


def xi(form, m_min, b, mmax, m_cut, m_top, tail="truncate_at_last_bin"):
    """Fraction of the moment of the whole distribution released in the modelled range,
    from m_min to the cut (or to m_top when the tail is not modelled)."""
    whole_top = m_cut if not np.isinf(m_cut) else mmax + 3.0
    upper = m_cut if not np.isinf(m_cut) else m_top
    return moment_integral(form, m_min, upper, b, mmax, m_cut) / moment_integral(form, M_LOW, whole_top, b, mmax, m_cut)


def bin_rates(form, lower, upper, b, mmax, tail="truncate_at_last_bin", budget=None, rate_mc=None, mc=None):
    """Annual rate of each bin. Give either budget (N m/yr; geodynamic) or, for the seismic branches,
    rate_mc and mc: the observed annual rate of events at or above an anchor magnitude, and that magnitude.
    Returns (rates, scale A, xi, mean moment per bin)."""
    lower, upper = np.asarray(lower, float), np.asarray(upper, float)
    m_cut = cut_magnitude(form, mmax, upper[-1], tail)
    g, mean_m0, mom = bin_shapes(form, lower, upper, b, mmax, m_cut)
    x = xi(form, lower[0], b, mmax, m_cut, upper[-1], tail)
    if budget is not None:
        A = x * budget / mom.sum()
    else:
        A = rate_mc / float(shape(form, mc, b, mmax, m_cut))
    return A * g, A, x, mean_m0


def cumulative(form, m, b, mmax, A, m_top, tail="truncate_at_last_bin"):
    """Annual rate of events of magnitude >= m for a scale A."""
    return A * shape(form, m, b, mmax, cut_magnitude(form, mmax, m_top, tail))
