"""Check reports. Each step writes checks/<step>/report.md and status.json;
checks/INDEX.md lists every check of every step."""
import datetime
import glob
import json
import os
from . import CHECKS


class CheckReport:
    def __init__(self, step, title):
        self.step, self.title = step, title
        self.dir = os.path.join(CHECKS, step)
        os.makedirs(self.dir, exist_ok=True)
        self.rows, self.notes = [], []

    def check(self, name, why, passed, detail=""):
        """Record one check. passed is True, False, or None for information only."""
        status = "INFO" if passed is None else "PASS" if passed else "FAIL"
        self.rows.append(dict(name=name, why=why, status=status, detail=str(detail)))
        print(f"[{status}] {name}: {detail}", flush=True)
        return passed

    def note(self, text):
        self.notes.append(text)

    @property
    def ok(self):
        return not any(r["status"] == "FAIL" for r in self.rows)

    def write(self):
        when = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        n_fail = sum(r["status"] == "FAIL" for r in self.rows)
        L = [f"# Checks for step {self.step}: {self.title}", "",
             f"Run {when}. Result: **{'all checks passed' if self.ok else f'{n_fail} check(s) FAILED'}**.", "",
             "| Check | Why | Result | Detail |", "|---|---|---|---|"]
        L += [f"| {r['name']} | {r['why']} | {r['status']} | {r['detail']} |" for r in self.rows]
        if self.notes:
            L += ["", "## Notes", ""] + [f"- {n}" for n in self.notes]
        with open(os.path.join(self.dir, "report.md"), "w", encoding="utf-8") as f:
            f.write("\n".join(L) + "\n")
        with open(os.path.join(self.dir, "status.json"), "w", encoding="utf-8") as f:
            json.dump(dict(step=self.step, title=self.title, run=when, ok=self.ok, checks=self.rows), f, indent=1)
        write_index()
        return self.ok


def write_index():
    L = ["# Index of checks", "",
         "One row per check. Reports are in the folder of each step. Rebuilt every time a step is run.", "",
         "| Step | Check | Why | Result | Last run |", "|---|---|---|---|---|"]
    for p in sorted(glob.glob(os.path.join(CHECKS, "*", "status.json"))):
        with open(p, encoding="utf-8") as f:
            s = json.load(f)
        L += [f"| {s['step']} | {r['name']} | {r['why']} | {r['status']} | {s['run']} |" for r in s["checks"]]
    with open(os.path.join(CHECKS, "INDEX.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
