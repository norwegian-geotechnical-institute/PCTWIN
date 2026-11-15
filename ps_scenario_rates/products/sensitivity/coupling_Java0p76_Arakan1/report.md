# Sensitivity test: coupling Java = 0.76, Arakan = 1 (sunda)

Nothing in the configuration or in the products of steps B0 to B4 was changed. Low and high coupling follow the usual rule (preferred -/+ 0.2, within 0.1 to 1).

## Preferred coupling and budget

| Unit | Coupling now | Coupling in test | Budget now (N m/yr) | Budget in test | Ratio |
|---|---|---|---|---|---|
| whole_zone | 0.42 | 0.65 | 6.579e+20 | 1.010e+21 | 1.54 |
| Arakan | 0.45 | 1.00 | 6.097e+19 | 1.355e+20 | 2.22 |
| Aceh-Andaman | 0.70 | 0.70 | 1.552e+20 | 1.552e+20 | 1.00 |
| Sumatra | 0.80 | 0.80 | 2.836e+20 | 2.836e+20 | 1.00 |
| Java | 0.18 | 0.76 | 8.379e+19 | 3.618e+20 | 4.32 |
| Sumba | 0.20 | 0.20 | 7.430e+19 | 7.430e+19 | 1.00 |

## Events in 120 years, weighted mean over branches

| Group | Case | Rate at Mw >= 6.35 (per yr) | Mw >= 7.5 | Mw >= 8 | Mw >= 8.5 | Mw >= 9 |
|---|---|---|---|---|---|---|
| Observed | | 1.06 | 18 | 5 | 3 | 1 |
| geodynamic, unconstrained | present | 1.72 | 16.3 | 5.5 | 1.8 | 0.4 |
| geodynamic, unconstrained | test | 2.53 | 24.0 | 8.1 | 2.6 | 0.6 |
| geodynamic, segmented | present | 2.78 | 24.9 | 7.5 | 1.8 | 0.3 |
| geodynamic, segmented | test | 4.76 | 41.1 | 11.5 | 2.5 | 0.4 |
| seismic | present | 1.06 | 12.1 | 4.3 | 1.4 | 0.3 |
| seismic | test | 1.06 | 12.1 | 4.3 | 1.4 | 0.3 |
| all branches | present | 1.66 | 16.3 | 5.4 | 1.6 | 0.4 |
| all branches | test | 2.35 | 22.3 | 7.1 | 2.0 | 0.4 |

## Rate by segment, weighted mean over branches (all modelled magnitudes per year; Mw >= 7.5 in 120 years)

| Segment | Seismic | Geodynamic now | Geodynamic in test | Seismic, Mw >= 7.5 | Geodynamic now, Mw >= 7.5 | Geodynamic in test, Mw >= 7.5 |
|---|---|---|---|---|---|---|
| Arakan | 0.026 | 0.291 | 0.601 | 0.4 | 2.7 | 5.4 |
| Aceh-Andaman | 0.177 | 0.389 | 0.390 | 2.4 | 4.3 | 4.4 |
| Sumatra | 0.594 | 0.833 | 0.838 | 7.6 | 9.0 | 9.2 |
| Java | 0.195 | 0.436 | 1.505 | 2.6 | 4.6 | 15.6 |
| Sumba | 0.068 | 0.306 | 0.312 | 0.9 | 3.5 | 3.7 |

A rupture's rate is divided between segments by its segment fractions. rate_along_trench.png shows the same along the trench.

mfd.png: cumulative rates of the whole zone on every geodynamic branch, as configured and in the test. mfd_segments.png: the same for each segment that was changed (segmented branches).
