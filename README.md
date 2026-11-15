# PTHA integration

Code for the probabilistic tsunami hazard assessment (PTHA) in PCTWIN WP1.

The repository brings the pieces developed across the work package (source models, magnitude-frequency distributions, rupture realisations, tsunami propagation) into one PTHA workflow. The first version covers predominant subduction (PS) sources for the Sunda and Makran subduction zones. Background seismicity (BS) will be added later.

## Where things stand

| Step of the workflow | Status | Where |
|---|---|---|
| 1. Source model and scenario set | Done by the scenario provider (meshes, scenarios, rupture fractions) | inputs, see `ps_scenario_rates/data/README.md` |
| 2. Epistemic model / logic tree | Sunda: first run done | `ps_scenario_rates/` |
| 3. MFD and annual rates per magnitude bin | Sunda: first run done | `ps_scenario_rates/` |
| 4. Rates of the aleatory realisations (position, scaling sample, slip) | Sunda: first run done | `ps_scenario_rates/` |
| 5. Scenario catalogue with annual rates on every branch | Sunda: first run done | `ps_scenario_rates/products/` |
| 6. Unit-source combination / initial conditions | not in this repository yet | |
| 7. Tsunami propagation (Tsunami-HySEA) | not in this repository yet | |
| 8. Intensities at points of interest and inundation | not in this repository yet | |
| 9. Hazard aggregation | to do; how to use the rates is in `ps_scenario_rates/docs/guide.md`, section 2 | |

Makran has not been run yet; what it needs is listed in the guide.

## PS scenario rates

`ps_scenario_rates/` gives every PS scenario an annual rate on each of the 480 branches of the logic tree. It is self-contained: its own configuration, code, unit tests, result tables and figures.

- Start with `ps_scenario_rates/README.md` (how to use the rates, how to run).
- The method, a table of every step with its limitations and planned improvements, the decisions taken and the open points are in `ps_scenario_rates/docs/guide.md`.

- The numbers and figures of the current version are in `ps_scenario_rates/docs/results.md`; every modelling choice is set in `ps_scenario_rates/config/`.
- Current version of the rates: 1.0 (tag `v1.0`).

## Repository structure

```
ps_scenario_rates/
    config/       logic tree and zone parameters
    ptha_rates/   the code
    scripts/      one run script per step, run_all.py, sensitivity test
    tests/        unit tests
    products/     result tables and figures of every step
    data/         input files
    docs/         guide.md, results.md and figures
```

Folders for the later steps (initial conditions, propagation, hazard aggregation) are to be added as that work comes in.
