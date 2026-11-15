"""Delete everything the steps produce, so that the chain can be run again from a fresh start.

Run from ps_scenario_rates:
    python scripts/clean.py              delete the outputs
    python scripts/clean.py --dry-run    only list what would be deleted
    python scripts/clean.py --all        also delete products/B1_inputs when the rupture fractions file is missing

Deleted: products/ (tables and figures of every step and of the sensitivity tests), checks/ (check reports),
docs/figures/ and the __pycache__ folders. Never touched: config/, data/, the code and the documents.

Step B1 needs the rupture fractions file, which is not in the repository (data/README.md). If that file is not
in the data folder, products/B1_inputs is kept, because it could not be made again; run_all.py then starts from
the tables that are there. After cleaning:   python scripts/run_all.py
"""
import os
import shutil
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ptha_rates import CHECKS, DATA, PRODUCTS, ROOT

FRACTIONS = os.path.join(DATA, "Rupture_fractions_v2.zip")


def targets(keep_b1):
    out = []
    if os.path.isdir(PRODUCTS):
        for name in sorted(os.listdir(PRODUCTS)):
            if not (keep_b1 and name == "B1_inputs"):
                out.append(os.path.join(PRODUCTS, name))
    out += [CHECKS, os.path.join(ROOT, "docs", "figures")]
    for folder, dirs, _ in os.walk(ROOT):
        if ".git" in dirs:
            dirs.remove(".git")
        out += [os.path.join(folder, d) for d in dirs if d == "__pycache__"]
    return [t for t in out if os.path.exists(t)]


def main(args):
    unknown = [a for a in args if a not in ("--dry-run", "--all")]
    if unknown:
        print(__doc__); return 1
    dry = "--dry-run" in args
    keep_b1 = not os.path.exists(FRACTIONS) and "--all" not in args
    todo = targets(keep_b1)
    for t in todo:
        print(("would delete " if dry else "deleting ") + os.path.relpath(t, ROOT))
        if not dry:
            shutil.rmtree(t) if os.path.isdir(t) else os.remove(t)
    if not todo:
        print("Nothing to delete.")
    if keep_b1 and os.path.isdir(os.path.join(PRODUCTS, "B1_inputs")):
        print(f"Kept products/B1_inputs: {os.path.basename(FRACTIONS)} is not in {DATA}, so step B1 could not be run again. "
              "Use --all to delete it as well.")
    if not dry:
        print("Clean. To run everything again:   python scripts/run_all.py")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
