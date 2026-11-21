"""Step B0. Load and validate the configuration files."""
import os
import yaml
from . import CONFIG

LEVELS = ["rate_source", "segmentation", "mfd_form", "b_value", "mmax", "coupling"]
TOL = 1e-9
FOOTPRINTS = ["segment_fraction_approximation", "estimated", "from_file"]


def load(config_dir=None):
    """Return {'tree': ..., 'zones': ...} read from logic_tree.yaml and zones.yaml."""
    d = config_dir or CONFIG
    with open(os.path.join(d, "logic_tree.yaml"), encoding="utf-8") as f:
        tree = yaml.safe_load(f)
    with open(os.path.join(d, "zones.yaml"), encoding="utf-8") as f:
        zones = yaml.safe_load(f)
    return {"tree": tree, "zones": zones}


def _sum_to_one(name, weights, errors):
    s = sum(weights.values())
    if abs(s - 1) > TOL:
        errors.append(f"weights of {name} sum to {s:.6f}, not 1")
    if any(w < 0 for w in weights.values()):
        errors.append(f"negative weight in {name}")


def validate(cfg):
    """Return (errors, missing). errors stop the run; missing lists values that are
    still null, each as (zone, where, field)."""
    errors, missing = [], []
    lv = cfg["tree"]["levels"]
    for k in LEVELS:
        if k not in lv:
            errors.append(f"level {k} is not defined")
    for k in ("rate_source", "segmentation", "mfd_form", "mmax"):
        _sum_to_one(k, lv[k], errors)
    _sum_to_one("coupling", lv["coupling"]["weights"], errors)
    b = lv["b_value"]
    if "values" in b:
        if len(b["values"]) < 2 or (isinstance(b["weights"], list) and len(b["weights"]) != len(b["values"])):
            errors.append("b_value needs at least 2 values and one weight per value")
    elif not (b["n"] >= 2 and b["min"] < b["max"]):
        errors.append("b_value needs n >= 2 and min < max")
    c = lv["coupling"]
    if not (0 <= c["floor"] < c["ceiling"]):
        errors.append("coupling floor and ceiling must satisfy 0 <= floor < ceiling")
    if c.get("rule", "offset") == "factor" and set(c["factors"]) != set(c["weights"]):
        errors.append("coupling: factors and weights must name the same levels")
    for k in ("rigidity_model", "scaling_law"):
        _sum_to_one(k, cfg["tree"]["event_tree"][k], errors)
    r = cfg["tree"]["options"].get("rigidity", {})
    for k in active(cfg["tree"]["event_tree"]["rigidity_model"]):
        if k not in r.get("models", {"00": "homogeneous"}):
            errors.append(f"rigidity model {k} has weight but no entry in options: rigidity: models")
    if r.get("without_footprints", "skip") not in ("skip", "approximate"):
        errors.append("options: rigidity: without_footprints must be skip or approximate")
    for k in cfg["tree"]["event_tree"].get("as_level", []):
        if k not in ("scaling_law", "scaling_sample"):
            errors.append(f"event_tree: as_level: {k} is not one of scaling_law, scaling_sample")
    fp = cfg["tree"]["options"].get("footprints")
    if fp not in FOOTPRINTS:
        errors.append(f"options: footprints must be one of {', '.join(FOOTPRINTS)}")
    if fp == "estimated" and any(k not in cfg["tree"]["options"].get("footprint_area", {}) for k in active(cfg["tree"]["event_tree"]["scaling_law"])):
        errors.append("options: footprint_area needs an entry for every scaling law in use")
    for zone, z in cfg["zones"].items():
        for f in ("triangles", "catalogue"):
            if z.get(f) is None:
                missing.append((zone, "zone", f))
        for f, v in z["whole_zone"].items():
            if v is None:
                missing.append((zone, "whole_zone", f))
        for s in z["segments"]:
            for f, v in s.items():
                if v is None and f != "ref_cs":
                    missing.append((zone, s["name"], f))
            if s.get("coupling") is not None and not (0 <= s["coupling"] <= cfg["tree"]["levels"]["coupling"]["ceiling"]):
                errors.append(f"{zone} {s['name']}: coupling outside 0 to the ceiling")
    return errors, missing


def rigidity_in_use(tree, footprints_have_slip=None):
    """Weights of the rigidity models that are used in this run (decision 27). A model listed under
    options: rigidity: needs_footprints is left out (and the other weights renormalised) unless
    footprints: from_file is set or without_footprints is 'approximate'."""
    w = {k: v for k, v in tree["event_tree"]["rigidity_model"].items() if v > 0}
    r = tree["options"].get("rigidity", {})
    has_fp = tree["options"]["footprints"] == "from_file" if footprints_have_slip is None else footprints_have_slip
    if r.get("without_footprints", "skip") == "skip" and not has_fp:
        w = {k: v for k, v in w.items() if k not in r.get("needs_footprints", [])}
    if not w:
        raise ValueError("no rigidity model left to use")
    s = sum(w.values())
    return {k: v / s for k, v in w.items()}


def skipped_rigidity(tree):
    """Rigidity models with weight that this run leaves out."""
    return [k for k, v in tree["event_tree"]["rigidity_model"].items() if v > 0 and k not in rigidity_in_use(tree)]


def active(weights):
    """Keys of a weight table that have non-zero weight."""
    return [k for k, w in weights.items() if w > 0]
