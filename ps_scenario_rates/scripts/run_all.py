"""Run the unit tests, every step in order and the sensitivity tests; refresh docs/figures and docs/results.md.

Run from ps_scenario_rates:   python scripts/run_all.py [zone ...]     (default: sunda makran)
Stops at the first step that fails a check. Step B1 is skipped when the rupture fractions file is not in the data
folder and its tables are already in products/B1_inputs. For a fresh start, run scripts/clean.py first.
"""
import glob
import os
import shutil
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STEPS = ["B0_config", "B1_inputs", "B2_budget", "B3_catalogue", "B4_bin_rates", "B5_shares", "B6_segmented", "B7_rates", "B8_checks"]
# Sensitivity tests kept with the run. Sunda: the coupling of Table 1 for Java and Arakan (decision 21).
SENSITIVITY = {"sunda": [["Java=0.76", "Arakan=1.0"]]}


def run(args):
    print(">", " ".join(args), flush=True)
    return subprocess.run([sys.executable] + args, cwd=ROOT).returncode


def main(zones):
    if run(["-m", "unittest", "discover", "tests"]):
        return 1
    sys.path.insert(0, ROOT)
    from ptha_rates import DATA, PRODUCTS
    have_fractions = os.path.exists(os.path.join(DATA, "Rupture_fractions_v2.zip"))
    for s in STEPS:
        if s == "B1_inputs" and not have_fractions:
            missing = [z for z in zones if not os.path.exists(os.path.join(PRODUCTS, "B1_inputs", f"scenarios_{z}.npz"))]
            if missing:
                print(f"Step B1 needs Rupture_fractions_v2.zip in {DATA} (see data/README.md); no stored tables for: {', '.join(missing)}"); return 1
            print("Step B1 skipped: Rupture_fractions_v2.zip is not in the data folder; using the tables in products/B1_inputs.", flush=True)
            continue
        if run([os.path.join("scripts", f"run_{s}.py")] + zones):
            print(f"Step {s} failed a check; see checks/{s}/report.md"); return 1
    for z in zones:
        for case in SENSITIVITY.get(z, []):
            if run([os.path.join("scripts", "sens_coupling.py"), z] + case):
                return 1
    fig = os.path.join(ROOT, "docs", "figures"); os.makedirs(fig, exist_ok=True)
    found = glob.glob(os.path.join(ROOT, "products", "B*", "*.png")) + glob.glob(os.path.join(ROOT, "checks", "B*", "*.png"))
    for f in found:
        shutil.copyfile(f, os.path.join(fig, os.path.basename(f)))
    for f in glob.glob(os.path.join(ROOT, "products", "sensitivity", "*", "*.png")):
        shutil.copyfile(f, os.path.join(fig, f"sensitivity_{os.path.basename(os.path.dirname(f))}_{os.path.basename(f)}"))
    if run([os.path.join("scripts", "make_results.py")] + zones):
        return 1
    print(f"Done. Check results: checks/INDEX.md. Numbers and figures of this run: docs/results.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["sunda", "makran"]))
