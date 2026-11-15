"""Figures: the logic tree and its weights (step B0) and the rate along the trench (step B8, sensitivity tests)."""
import os
import tempfile
os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp())
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# One style for every figure: text sized to be read when the figure is placed at page width.
matplotlib.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "xtick.labelsize": 10,
                            "ytick.labelsize": 10, "legend.fontsize": 10, "figure.titlesize": 12, "savefig.dpi": 200})
NOTE = 10          # size of notes written under a figure
from matplotlib.patches import FancyBboxPatch
from . import branches

LABEL = {"rate_source": "Rate source", "segmentation": "Segmentation", "mfd_form": "MFD form",
         "b_value": "b-value", "mmax": "Mmax", "coupling": "Coupling",
         "seismic": "Seismic", "geodynamic": "Geodynamic", "unconstrained": "Unconstrained",
         "segmented": "Segmented", "truncated_gr": "Truncated G-R", "tapered_pareto": "Tapered Pareto",
         "observed": "Observed", "area_limit": "Area limit", "mean": "Mean of the two",
         "low": "Low", "preferred": "Preferred", "high": "High"}


def logic_tree(cfg, zone, table, path):
    """One column per level, one box per alternative with its weight. As in Figure 1 of
    the method document, only the first alternative of each level is expanded."""
    tree, z = cfg["tree"], cfg["zones"][zone]
    lv = tree["levels"]
    bv, bw = branches.b_values(tree), branches.b_weights(z["whole_zone"]["b"], tree)
    off = lv["coupling"]["offset"]
    cols = [("rate_source", [(LABEL[k], w) for k, w in lv["rate_source"].items()]),
            ("segmentation", [(LABEL[k], w) for k, w in lv["segmentation"].items()]),
            ("mfd_form", [(LABEL[k], w) for k, w in lv["mfd_form"].items()]),
            ("b_value", [(f"b = {b:.3f}", w) for b, w in zip(bv, bw)]),
            ("mmax", [(LABEL[k], w) for k, w in lv["mmax"].items()]),
            ("coupling", [(f"{LABEL[k]} ({'%+.1f' % d})" if d else LABEL[k], w)
                          for (k, w), d in zip(lv["coupling"]["weights"].items(), (-off, 0, off))])]
    nmax = max(len(a) for _, a in cols)
    bw_, bh, dx, dy = 2.9, 0.66, 3.3, 0.9
    fig, ax = plt.subplots(figsize=(16, 0.62 * nmax + 2.0))
    ax.set_xlim(-0.2, dx * len(cols) - 0.4); ax.set_ylim(-dy * nmax - 0.9, 1.0); ax.axis("off")
    prev = None
    for j, (key, alts) in enumerate(cols):
        x = j * dx
        ax.text(x + bw_ / 2, 0.62, LABEL[key], ha="center", va="center", fontsize=13, fontweight="bold")
        if key == "coupling":
            ax.text(x + bw_ / 2, 0.28, "geodynamic branches only", ha="center", va="center", fontsize=10, style="italic")
        top = -(nmax - len(alts)) * dy / 2
        ys = [top - i * dy for i in range(len(alts))]
        for (name, w), y in zip(alts, ys):
            ax.add_patch(FancyBboxPatch((x, y - bh), bw_, bh, boxstyle="round,pad=0,rounding_size=0.06",
                                        fc="white", ec="black", lw=0.8))
            ax.text(x + 0.1, y - bh / 2, name, ha="left", va="center", fontsize=10.5)
            ax.text(x + bw_ - 0.1, y - bh / 2, f"{w:.3f}", ha="right", va="center", fontsize=10.5, fontweight="bold")
            if prev is not None:
                ax.plot([prev[0], x], [prev[1], y - bh / 2], color="black", lw=0.4)
        prev = (x + bw_, ys[0] - bh / 2)
    n_g = int((table.rate_source == "geodynamic").sum()); n_s = int((table.rate_source == "seismic").sum())
    ax.text(-0.1, -dy * nmax - 0.55,
            f"{zone.capitalize()}: {len(table)} branches ({n_g} geodynamic, {n_s} seismic). Numbers in bold are the weights of each alternative "
            f"within its level.\nb-value weights use the whole-zone catalogue estimate b = {z['whole_zone']['b']}. "
            f"The weight of a branch is the product of its level weights.", fontsize=NOTE, va="top")
    fig.savefig(path, dpi=200, bbox_inches="tight"); plt.close(fig)


def weights_panel(cfg, zone, path):
    """b-value weights for the whole zone (used for branch weights) and, for reference,
    what they would be with each segment's own estimate; and the Mmax weights."""
    tree, z = cfg["tree"], cfg["zones"][zone]
    bv = branches.b_values(tree)
    units = [("Whole zone (used)", z["whole_zone"]["b"])] + [(s["name"], s["b"]) for s in z["segments"]]
    fig, axes = plt.subplots(2, 4, figsize=(11, 6.2), sharey=True); axes = axes.ravel()
    for ax, (name, bh) in zip(axes, units):
        ax.bar(bv, branches.b_weights(bh, tree), width=0.04, fc="0.75" if "used" not in name else "0.3", ec="black", lw=0.5)
        ax.axvline(bh, color="black", lw=0.8, ls="--")
        ax.set_title(f"{name}\ncatalogue b = {bh}"); ax.set_xlabel("b-value")
    axes[0].set_ylabel("weight"); axes[4].set_ylabel("weight")
    for ax in axes[len(units):-1]:
        ax.axis("off")
    mm = tree["levels"]["mmax"]
    axes[-1].bar(range(len(mm)), list(mm.values()), fc="0.3", ec="black", lw=0.5)
    axes[-1].set_xticks(range(len(mm))); axes[-1].set_xticklabels([LABEL[k].replace(" ", "\n", 1) for k in mm], fontsize=9.5)
    axes[-1].set_title("Mmax weights\n(all zones and segments)")
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def rate_along_trench(path, title, tri, trench, s_pos, kb, bins, names, curves, width=100.0, note=""):
    """Rate per `width` km of trench for three magnitude classes, with a map of how the arc is unrolled.
    tri: triangle table indexed by id (lon, lat, segment, along_trench_km); trench: trench line (lon, lat, dist_km);
    s_pos, kb: distance along the trench and magnitude bin of every position;
    curves: list of (rate per position, line style, colour, label)."""
    n_seg = len(names)
    L = trench.dist_km.iloc[-1]; edges = np.arange(0.0, np.ceil(L / width) * width + width, width)
    bounds = [(tri.along_trench_km[tri.segment == g].quantile(.98) + tri.along_trench_km[tri.segment == g + 1].quantile(.02)) / 2 for g in range(1, n_seg)]
    mids = np.diff(np.r_[0.0, bounds, L]) / 2 + np.r_[0.0, bounds]
    classes = [("Mw 6.35 to 7.5", (bins.mw_centre < 7.5).to_numpy()), ("Mw 7.5 to 8.5", ((bins.mw_centre >= 7.5) & (bins.mw_centre < 8.5)).to_numpy()), ("Mw 8.5 and above", (bins.mw_centre >= 8.5).to_numpy())]
    fig = plt.figure(figsize=(11, 7.6))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.75, 1], hspace=0.12, wspace=0.16)
    axes = [fig.add_subplot(gs[i, 0]) for i in range(3)]
    for ax, (lab, kk) in zip(axes, classes):
        sel = kk[kb]
        for pm_, ls, col, nm in curves:
            h = np.histogram(s_pos[sel], bins=edges, weights=pm_[sel])[0]
            ax.step(edges, np.r_[h, h[-1]], where="post", color=col, ls=ls, lw=1.1, label=nm)
        for b_ in bounds:
            ax.axvline(b_, color="0.6", lw=0.6)
        ax.set_xlim(0, L); ax.set_ylim(bottom=0); ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        ax.text(0.012, 0.9, lab, transform=ax.transAxes, fontsize=11, fontweight="bold", va="top")
        if ax is not axes[-1]:
            ax.set_xticklabels([])
    for m_, n_ in zip(mids, names):
        axes[0].text(m_, 1.075, n_, transform=axes[0].get_xaxis_transform(), ha="center", va="bottom", fontsize=10)
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False, loc="upper center", bbox_to_anchor=(0.775, 0.89))   # above the map
    axes[1].set_ylabel(f"events per year per {width:g} km of trench"); axes[2].set_xlabel("distance along the trench (km)")
    am = fig.add_subplot(gs[:, 1])
    am.scatter(tri.lon, tri.lat, s=0.4, color="0.8", rasterized=True)
    am.plot(trench.lon, trench.lat, color="black", lw=1.6)
    for dk in np.arange(0, L + 1, 1000):
        x_, y_ = np.interp(dk, trench.dist_km, trench.lon), np.interp(dk, trench.dist_km, trench.lat)
        am.plot(x_, y_, "o", ms=5, mfc="white", mec="black")
        am.annotate(f"{dk:,.0f}", (x_, y_), xytext=(-7, -3), textcoords="offset points", ha="right", va="center", fontsize=9.5)
    for b_ in bounds:
        x_, y_ = np.interp(b_, trench.dist_km, trench.lon), np.interp(b_, trench.dist_km, trench.lat)
        am.plot(x_, y_, "|", ms=14, mew=1.2, color="0.35")
    for m_, n_ in zip(mids, names):
        j = (tri.along_trench_km - m_).abs() < 150
        flat = tri.lat[j].mean() < -7.5                    # where the arc runs east-west, put the name above it
        am.text(tri.lon[j].mean() + (0.0 if flat else 0.8), tri.lat[j].mean() + (1.6 if flat else 0.0), n_, fontsize=9.5,
                ha="center" if flat else "left", va="bottom" if flat else "center", style="italic")
    am.set_aspect(1 / np.cos(np.deg2rad(tri.lat.mean()))); am.set_xlabel("longitude"); am.set_ylabel("latitude"); am.yaxis.tick_right(); am.yaxis.set_label_position("right")
    am.set_xlim(tri.lon.min() - 6.5, tri.lon.max() + 1)
    fig.suptitle(title, y=0.97)
    fig.text(0.01, 0.01, "Right: the trench line (black) with distance along it in km; grey, the mesh. A rupture is placed at the trench point nearest to its initiation triangle.\n"
             "Left: grey lines are approximate segment boundaries." + (" " + note if note else ""), fontsize=NOTE, va="bottom")
    fig.subplots_adjust(left=0.07, right=0.94, top=0.9, bottom=0.14); fig.savefig(path); plt.close(fig)
