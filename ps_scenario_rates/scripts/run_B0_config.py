"""Step B0. Validate the configuration, write the branch table and draw the logic tree.

Run from ps_scenario_rates:   python scripts/run_B0_config.py [zone ...]     (default: sunda makran)
Writes products/B0_config/ and checks/B0_config/.
"""
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from ptha_rates import PRODUCTS, branches, config, plotting
from ptha_rates.report import CheckReport


def main(zones):
    cfg = config.load()
    out = os.path.join(PRODUCTS, "B0_config"); os.makedirs(out, exist_ok=True)
    rep = CheckReport("B0_config", "configuration and branch table")
    errors, missing = config.validate(cfg)
    rep.check("Configuration is valid", "Weights of every level must sum to 1 and no level may be absent",
              not errors, "; ".join(errors) or "no errors")
    tree = cfg["tree"]; lv = tree["levels"]
    for zone in zones:
        z = cfg["zones"][zone]
        miss = [f"{w}: {f}" for zz, w, f in missing if zz == zone]
        rep.check(f"{zone}: values still to be filled", "Null entries in zones.yaml; each is needed by a later step",
                  None, "; ".join(miss) or "none")
        t = branches.branch_table(cfg, zone)
        t.to_csv(os.path.join(out, f"branch_table_{zone}.csv"), index=False)
        nb = len(branches.b_values(tree)); ev = tree["event_tree"]; al = ev.get("as_level", [])
        n_ev = (len(config.active(ev["scaling_law"])) if "scaling_law" in al else 1) * (int(ev.get("n_scaling_samples", 10)) if "scaling_sample" in al else 1)
        n_seis = len(lv["segmentation"]) * len(lv["mfd_form"]) * nb * len(lv["mmax"]) * n_ev
        n_geo = n_seis * len(lv["coupling"]["weights"])
        rep.check(f"{zone}: number of branches", "Must equal the product of the numbers of alternatives",
                  (t.rate_source == "geodynamic").sum() == n_geo and (t.rate_source == "seismic").sum() == n_seis,
                  f"{n_geo} geodynamic + {n_seis} seismic = {len(t)}")
        rep.check(f"{zone}: branch weights sum to 1", "Weights are probabilities over the whole tree",
                  abs(t.weight.sum() - 1) < 1e-9, f"sum = {t.weight.sum():.12f}")
        rep.check(f"{zone}: branch ids are unique", "Each branch is looked up by its id",
                  t.branch_id.is_unique, f"{t.branch_id.nunique()} ids")
        g = t.groupby("rate_source").weight.sum()
        rep.check(f"{zone}: weight of each rate source", "Seismic branches have no coupling level and must still carry the level weight",
                  all(abs(g[k] - w) < 1e-9 for k, w in lv["rate_source"].items()), ", ".join(f"{k} {v:.4f}" for k, v in g.items()))
        rep.check(f"{zone}: spatial variants", "Share sets needed in steps B5 and B6", t.spatial_variant.nunique() == 8,
                  f"{t.spatial_variant.nunique()} variants")
        bw = branches.b_weights(z["whole_zone"]["b"], tree); bv = branches.b_values(tree)
        rep.check(f"{zone}: b-value weights", "Highest at the value nearest the catalogue estimate, all positive",
                  bool(np.argmax(bw) == np.argmin(np.abs(bv - z["whole_zone"]["b"])) and (bw > 0).all()),
                  "; ".join(f"{b:.3f}: {w:.3f}" for b, w in zip(bv, bw)))
        cl = {s["name"]: branches.coupling_levels(s["coupling"], tree) for s in z["segments"] if s["coupling"] is not None}
        rep.check(f"{zone}: coupling levels", "Low / preferred / high from the coupling rule, within floor and ceiling (decision 25)", None,
                  "; ".join(f"{k}: {v['low']:.2f}/{v['preferred']:.2f}/{v['high']:.2f}" for k, v in cl.items()) or "no coupling values")
        dup = [k for k, v in cl.items() if len({v['low'], v['preferred'], v['high']}) < 3]
        rep.check(f"{zone}: segments with fewer than three distinct coupling levels",
                  "Identical levels give identical branches for that segment", None, ", ".join(dup) or "none")
        plotting.logic_tree(cfg, zone, t, os.path.join(out, f"logic_tree_{zone}.png"))
        plotting.weights_panel(cfg, zone, os.path.join(out, f"weights_panel_{zone}.png"))
    rep.note("The b-value weight of a branch uses the whole-zone catalogue estimate, also in the segmented branches: "
             "a branch applies one b-value to every segment, so it can carry only one b-value weight. "
             "The weights panel shows, for reference, what the weights would be with each segment's own estimate.")
    rep.note("Products: branch_table_<zone>.csv, logic_tree_<zone>.png, weights_panel_<zone>.png in products/B0_config.")
    return rep.write()


if __name__ == "__main__":
    sys.exit(0 if main(sys.argv[1:] or ["sunda", "makran"]) else 1)
