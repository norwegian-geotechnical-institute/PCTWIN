# PS scenario rates: guide

Version 1.0 (15 November 2026). Sunda, first run is done. Makran not yet done.

This guide explains the method, the steps, the decisions and the open points. The numbers and figures of the run are in
[`results.md`](results.md), which the run writes itself. The inputs are in `config/`.

## 1. Introduction

This folder computes the annual rate of every scenario, on every branch of a logic tree that carries the uncertainty in how often earthquakes of each magnitude occur and where. This run (version 1.0) follows the discussion on 11th and 12th november 2026 meeting within WP1. The choices are simplified (quite a lot, I would say) , below is the list of main features that are in this first run:

| | In this run |
|---|---|
| Zone | Sunda megathrust, five segments: Arakan, Aceh-Andaman, Sumatra, Java, Sumba |
| Scenarios | The 2,075,771 scenarios of rigidity model 00 (homogeneous, 30 GPa). The other rigidity models have no rate yet. (I don't have rigidity models!!) |
| Magnitudes | 31 bins, Mw 6.35 to 9.52 (This won't change, unless dipankar has used some other bin edges) |
| Branches | 480 (360 geodynamic, 120 seismic) |
| Weights | Several are equal weights. No Bayesian update yet (must needed later). |

Note that, I have laid out the code with triangle-level detail in mind. The moment budget is already computed triangle by triangle, so depth-dependent rigidity or coupling that varies in space can be added to the budget with little effort. The placement of the rates within a segment can also be made triangle-level through rupture footprints, as an optional setting (end of section 5); version 1.0 does not use it.

## 2. Computing hazard from the rates

For computing hazard, the tables in `products/` are enough; the input data are not needed.

**Mean hazard curve.** Use the mean rate of each scenario over all branches:

```python
import numpy as np
from ptha_rates import rates

scen  = np.load("products/B1_inputs/scenarios_sunda.npz")
names = rates.scenario_names("sun", scen)                       # PS_sun_00_..., one per scenario
mean  = np.load("products/B7_rates/scenario_rates_sunda.npz")["rate_mean"]   # events per year, same order

# h[i] = tsunami height of scenario names[i] at one point of interest (from the simulations)
thresholds = np.logspace(-2, 1.3, 50)
annual_rate = np.array([mean[h >= x].sum() for x in thresholds])     # rate of exceedance
prob_50yr   = 1 - np.exp(-annual_rate * 50)                          # Poisson probability in 50 years
```

**Curve of one branch, and percentiles.** Replace `mean` by the rates of a branch and repeat for all 480:

```python
import pandas as pd
branches = pd.read_csv("products/B0_config/branch_table_sunda.csv")   # branch_id, alternative at each level, weight
D = rates.load("sunda")                                               # read the tables once
curves = []
for br in branches.itertuples():
    r = rates.position_rates(D, br)[scen["position"]] * scen["split"]  # rate of every scenario on this branch
    curves.append([r[h >= x].sum() for x in thresholds])
curves = np.array(curves)                                             # branch x threshold
# weighted mean over branches equals the curve from the mean rates; weighted percentiles give the epistemic spread
```

`rates.scenario_rates("sunda", branch_id)` returns the same for a single branch.

Things to know before using them:

- **Order.** Rates are in the order of `products/B1_inputs/scenarios_sunda.npz`. Match your simulation results to
  `names`; do not rely on the order of any other file. `products/B7_rates/scenario_mean_rates_sunda.csv.gz`
  has name and mean rate side by side.
- **Coverage.** Only scenarios of rigidity model 00 have a rate. A simulated scenario of another rigidity model must
  be left out for now, not given a zero.
- **Units.** Events per year. The rates of all scenarios add up to the total rate of Mw 6.35 and larger.
- **Version.** `scenario_rates_sunda.npz` carries the version of the rates (`version`); quote it with any hazard result.
- **Weights.** `weight` in the branch table sums to 1 over the 480 branches. The split between rigidity models
  and scaling laws is already inside the rate (they are pooled, not branches).
- **Means of groups.** `scenario_rates_sunda.npz` also holds the mean over the seismic, geodynamic, unconstrained and
  segmented branches (`rate_mean_seismic` and so on), for a quick view of what drives the hazard.
- **Provisional.** Section 5 below lists what each step approximates and section 9 the open points on the inputs. Rates
  will change when those are settled; rerun and replace, do not patch!.

The formula behind a rate, should you want to rebuild it from the tables yourself:

    rate(scenario, branch) = pooled weight x bin rate(branch) x share(position) x split(scenario)

| Table | File |
|---|---|
| Branches and weights | `products/B0_config/branch_table_sunda.csv` |
| Scenarios: parts of the name, bin, position, split | `products/B1_inputs/scenarios_sunda.npz` |
| Bin rates (unit x variant x form x b x Mmax x bin) | `products/B4_bin_rates/bin_rates_sunda.npz` |
| Shares, unconstrained branches (variant x position) | `products/B5_shares/shares_sunda.npz` |
| Shares, segmented branches (position x segment) | `products/B6_segmented/shares_segmented_sunda.npz` |

## 3. Method in brief
The method is detailed in this doc file: https://docs.google.com/document/d/1egjXVHOM4ZHP2JqCGmNSFqsZXr3nyb-Q/edit?rtpof=true&tab=t.0 However, the key ideas are included here as well.

Two terms you will meet below and all over the code:

- **Position**: where and how a rupture is placed: rigidity model + scaling law + magnitude bin + initiation triangle.
- **Split**: a position usually has several scenarios (scaling samples and slip realisations); its rate is divided
  equally among those that exist.

A **unit** is the whole zone or one segment, and a **bin** is one of the magnitudes at which scenarios exist.

Bin rates:

- Geodynamic branches: the MFD shape (truncated Gutenberg-Richter or tapered Pareto, with the b-value and Mmax of
  the branch) is scaled so that the moment it releases above the smallest modelled magnitude equals that part of
  the budget.
- Seismic branches: the same shape is scaled so that the rate above Mw 6.35 equals the rate observed since 1906.

Where the rate goes:

- **Unconstrained** branches use one MFD for the whole zone. A bin's rate is shared among positions in proportion
  to the loading (geodynamic) or the density of observed events (seismic) under the rupture.
- **Segmented** branches use one MFD per segment. A rupture takes rate from each segment it touches in proportion
  to the fraction of its moment in that segment.

Each scenario then receives its position's rate times its split.

## 4. The logic tree

![Logic tree](figures/logic_tree_sunda.png)

Six levels: rate source, segmentation, MFD form, b-value, Mmax and, on the geodynamic branches, coupling. The weight of a
branch is the product of the weights of its alternatives. Rigidity model and scaling law (Strasser et al. 2010;
Thingbaijam et al. 2017) are pooled inside each bin and are not branches.

The alternatives and weights are set in `config/logic_tree.yaml`, the values of each zone and segment (b, coupling,
convergence, Mmax) in `config/zones.yaml`. Both are listed as used, with the figure of the weights, in `results.md`.

## 5. The steps

One run script per step (`scripts/run_B<n>_<name>.py`); each writes `products/<step>/` and a check report.

| Step | What it does | Main products | Limitation in this run | Improvement |
|---|---|---|---|---|
| B0 Configuration | Reads and validates the two configuration files; builds the 480 branches with their weights | `branch_table`, logic tree figure, weights panel | Weights of rate source and coupling levels are equal by default | Agreed weights; Bayesian update against the record |
| B1 Input tables | Reads the scenario table and the mesh; assigns bins, positions, splits; gives every triangle a segment and a distance along the trench | `scenarios.npz`, `bins`, `triangles`, `trench` | Bin edges are midpoints between bin centres. The triangle-to-segment map is derived from the smallest scenarios (approximate at boundaries). Realisations per rupture are uneven (1 to 10). | Bin edges and segment table are needed from Dipankar to correct this |
| B2 Loading and budget | Coupling x convergence per triangle; moment budget per segment and for the zone at three coupling levels | `budget` | Coupling and convergence are constant inside a segment. Convergence is used as given (not projected onto the interface). One rigidity (30 GPa). | Values per triangle; a budget per rigidity model |
| B3 Catalogue | Assigns catalogue events to segments; counts events above the anchor magnitude; largest observed event per unit | `parameters`, catalogue fit figure | Events go to the segment of the nearest triangle with no depth or mechanism filter. Arakan (3 events) and Sumba (8) have uncertain anchor rates. b and Mc are taken from the Madan's estimates, not refitted. | Interface-only catalogue (thrust-mechanism events only); uncertainty on the anchor rate as a branch? |
| B4 Bin rates | Annual rate of every bin for every unit, rate variant, MFD form, b-value and Mmax | `bin_rates.npz`, MFD figure, modelled against observed large events | The tapered Pareto is cut at the top of the last bin (no scenarios above Mw 9.52). In segments with no large event on record, the "largest observed" Mmax forces the whole budget into moderate events. | Review the Mmax alternatives for quiet segments; decide the tail treatment |
| B5 Shares, unconstrained | Share of each bin's rate for every position, for the seismic variant and the three coupling levels | `shares.npz` | Rupture footprints are not used by default: loading under a rupture is taken from its segment fractions, so inside a segment all positions of a bin get the same weight. Can be run with Dipankar's slip models (end of this section). | Footprint-based shares (needs a reduction of the slip files on Dipankar's side) |
| B6 Shares, segmented | Share of each segment's bin rate for every position, including ruptures that cross a boundary | `shares_segmented.npz` | No penalty on boundary-crossing ruptures (`crossing_penalty: 1`). Same shares for all rate variants, for the reason given under B5. | Penalty from evidence on segment boundaries, if wanted |
| B7 Scenario rates | Mean rate of every scenario over all branches and over four groups; reader for any single branch | `scenario_rates.npz`, `scenario_mean_rates.csv.gz` | Rates per branch are rebuilt on request, not stored (480 x 2 million values) | |
| B8 End-to-end checks | Rate conserved in every bin on every branch; moment released against the budget, in total and per segment; events per segment against the catalogue; rate along the trench | check report, rate-along-trench figure | Reports, does not correct | |

The unit tests (`python -m unittest discover tests`) check the building blocks on small cases with known answers.

### Rupture footprints (optional setting)

A footprint is the list of triangles a rupture covers. With footprints, the weight of a rupture in steps B5 and B6 is the mean
loading (or seismicity density) over its own triangles, so loading that varies inside a segment is followed. The setting is
`footprints` in `config/logic_tree.yaml`:

| Setting | Footprint comes from | Status |
|---|---|---|
| `segment_fraction_approximation` | not used; segment fractions of the fractions file | default; version 1.0 |
| `estimated` | the triangles nearest to the initiation triangle, up to the median area of the scaling law | runs now, needs no new input |
| `from_file` | the slip files, reduced by `tools/reduce_slip_files.py` and placed in `data/footprints/` | runs once those files exist |

Checked against the slip sample available (105 ruptures, Mw 9.03, scaling sample 00): the area in the slip files equals the
scaling-law area within 0.1%, and the estimated footprint shares 96% of its triangles with the real one (median). Other
magnitudes and the other scaling samples, whose areas differ from the median, could not be checked. With the present
loading, constant inside each segment, `estimated` changes the scenario rates very little (0.2% of the rate moves, 2.5%
for Mw 8.5 and above), so it matters only once loading per triangle is supplied.

`reduce_slip_files.py` needs only Python and numpy and is run where the slip files are; its header explains the arguments.

## 6. Results

See [`results.md`](results.md): the configuration used, the moment budget, modelled against observed large events, the mean
rate of the zone, where the rate falls along the trench, the consistency checks and the figures. It is rewritten by every
run (`scripts/make_results.py`), so it always matches the products.

How to read it: the seismic branches match the record at moderate magnitudes by construction. The geodynamic branches are
an independent estimate; where they give more moderate events than the record, the likely cause is the "largest observed"
and "mean" Mmax alternatives in segments whose record holds no large event (step B4 in section 5).

## 7. Sensitivity to coupling

The coupling of Java and Arakan in Table 1 of the source characterisation is 0.76 and 1.00. With those values the
geodynamic branches gave several times more moderate events in the two segments than the catalogue holds, so the
configuration uses 0.176 and 0.45, for which coupling x convergence equals the preferred coupled convergence of
Davies et al. (2018) (12 and 8.1 mm/yr). The Table 1 values are kept as a sensitivity test:

```
python scripts/sens_coupling.py sunda Java=0.76 Arakan=1.0
```

The test recomputes budget, bin rates and shares in memory and writes to `products/sensitivity/`; nothing in the
configuration or the main products changes. Any segment and value can be tested the same way. `run_all.py` runs the
tests listed at its top, and their tables and figures appear at the end of `results.md`.

## 8. Decisions

Taken by Himanshu; open to review by the group. The numbers are those used in the code and the configuration.

| # | Question | Decision | Status |
|---|---|---|---|
| 1 | Weights of seismic against geodynamic | Equal, both kept | Open |
| 2 | b-values | Ten values over 0.7 to 1.2 for every unit; weight falls linearly with distance from the catalogue estimate | Decided |
| 3 | Low and high coupling | Preferred -/+ 0.2, floor 0.1, ceiling 1; equal weights | Decided; weights open |
| 4 | Rigidity model and scaling law | Pooled inside the bin with equal weights; not branches | Decided |
| 5 | Bin edges | Midpoints between bin centres | Stand-in |
| 6 | Convergence | Used as given; not projected onto the interface | Decided |
| 7 | Tapered Pareto corner | Equal to the Mmax alternative of the branch | Decided |
| 8 | Smallest magnitude for the budget | Lower edge of the first bin | Decided |
| 9 | Rupture footprints | Segment-fraction approximation in the first run | Decided |
| 10 | Scaling sample 00 | Counted as one sample among the others | Decided |
| 11 | Zero-slip triangles in a footprint | Not needed for the first run | Open |
| 12 | Output for the hazard calculation | Bin rates and shares for all branches, with a reader | Decided |
| 13 | Weight updating | Product of the level weights; no Bayesian update; modelled against observed reported only | Decided |
| 14 | Rigidity models | Homogeneous model only in the first run | Decided |
| 15 | b-value weight on segmented branches | From the whole-zone estimate | Decided |
| 16 | Sunda whole-zone area-limit Mmax | 9.5 | Decided |
| 17 | Segment Mmax | Area limit from the source characterisation; observed from the catalogue; third value their mean | Decided |
| 18 | First and last bins | Symmetric about their centres | Decided |
| 19 | Position | Rigidity model + scaling law + magnitude bin + initiation triangle | Decided |
| 20 | Seismic anchor | Observed rate of Mw >= 6.35 since 1906, per unit. The rate at the completeness magnitude since 2000 is not used: that period holds the 2004 to 2012 cluster. | Decided |
| 21 | Coupling of Java and Arakan | 0.176 and 0.45 (section 7); Table 1 values kept as a sensitivity test | Decided 15 Nov 2026 |

Choices made in the code and not yet confirmed: the tapered Pareto is cut at the top of the last bin with the
budget balanced over the modelled bins (`tapered_tail`); catalogue events are assigned to the nearest segment
without a filter; the seismic shares use events above the anchor per km2 of segment.

## 9. Some questions and potential consideration on the inputs

Here are few things I found while checking the scenario tables that might have some consequences. Flagged here in case we decide to do something about these: 

| # | Point | Effect |
|---|---|---|
| H1 | In the one slip sample checked (Mw 9.03), slip x area x 30 GPa is about 12% below the named moment for ruptures centred south of 8°S (Java, Sumba) | Tsunami sources there may carry less moment than the bin they get a rate for |
| H2 | 6% of Sunda scenarios have a total moment 0.01 to 0.07 magnitude units above the named magnitude (mostly slip realisation 000, Mw 7.93 and below) | No scenario changes bin; about 0.2% in the moment balance |
| H3 | Realisations are uneven: most rupture geometries have one slip realisation, others ten; scaling samples per triangle vary from 1 to 10 | The split uses the counts that exist |




## 10. Next steps

- **Makran.** The code takes any zone listed in `config/zones.yaml`. Makran needs its mesh as a text file, a
  catalogue, and the coupling and Mmax values that are still empty there. Then `python scripts/run_all.py makran`.
- **Depth-dependent rigidity models** (01 to 03): set their weights under `event_tree: rigidity_model` once a
  budget per model is coded.
- **Real rupture footprints.** The option is in the code (end of section 5). What is missing is the footprints themselves,
  reduced from Dipankar's slip files with `tools/reduce_slip_files.py`.
- **Weights from data** (Bayesian update against the record of large events).
- **Variable coupling and convergence rate.** Read them per triangle in step B2, not per segment. The footprints only
  start to matter once this is in.

## 11. Versions

The version is the number in the file `VERSION`; it is written into the scenario rates and into `results.md`, and each
version is a tag in the repository. Raise it whenever the rates change.

| Version | Date | What |
|---|---|---|
| 1.0 | 15 November 2026 | Sunda, first run: homogeneous rigidity model, 480 branches, decisions 1 to 21 |

To change an input: edit `config/`, run `python scripts/run_all.py`, add a row to the decisions in section 8 and to this
table, raise `VERSION`, commit and tag.
