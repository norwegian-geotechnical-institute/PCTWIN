# Sensitivity test: coupling Java = 0.76, Arakan = 1 (sunda)

Nothing in the configuration or in the products of steps B0 to B4 was changed. Low and high coupling follow the coupling rule of the logic tree (decision 25).

## Preferred coupling and budget

| Unit | Coupling now | Coupling in test | Budget now (N m/yr) | Budget in test | Ratio |
|---|---|---|---|---|---|
| whole_zone | 0.41 | 0.63 | 6.816e+20 | 1.055e+21 | 1.55 |
| Arakan | 0.45 | 1.00 | 6.076e+19 | 1.350e+20 | 2.22 |
| Aceh-Andaman | 0.70 | 0.70 | 1.546e+20 | 1.546e+20 | 1.00 |
| Sumatra | 0.80 | 0.80 | 2.874e+20 | 2.874e+20 | 1.00 |
| Java | 0.18 | 0.76 | 9.015e+19 | 3.893e+20 | 4.32 |
| Sumba | 0.20 | 0.20 | 8.875e+19 | 8.875e+19 | 1.00 |

## Events in 120 years, weighted mean over branches

| Group | Case | Rate at Mw >= 6.34 (per yr) | Mw >= 7.5 | Mw >= 8 | Mw >= 8.5 | Mw >= 9 |
|---|---|---|---|---|---|---|
| Observed | | 1.06 | 18 | 5 | 3 | 1 |
| geodynamic, unconstrained | present | 1.65 | 15.3 | 5.3 | 1.8 | 0.5 |
| geodynamic, unconstrained | test | 2.47 | 23.0 | 7.9 | 2.7 | 0.7 |
| geodynamic, segmented | present | 2.53 | 23.6 | 7.4 | 2.0 | 0.4 |
| geodynamic, segmented | test | 4.38 | 40.6 | 11.9 | 2.9 | 0.5 |
| seismic | present | 1.08 | 14.0 | 5.5 | 2.0 | 0.5 |
| seismic | test | 1.08 | 14.0 | 5.5 | 2.0 | 0.5 |
| all branches | present | 1.74 | 17.8 | 6.1 | 1.9 | 0.5 |
| all branches | test | 2.60 | 25.8 | 8.4 | 2.5 | 0.6 |

## Rate by segment, weighted mean over branches (all modelled magnitudes per year; Mw >= 7.5 in 120 years)

| Segment | Seismic | Geodynamic now | Geodynamic in test | Seismic, Mw >= 7.5 | Geodynamic now, Mw >= 7.5 | Geodynamic in test, Mw >= 7.5 |
|---|---|---|---|---|---|---|
| Arakan | 0.027 | 0.303 | 0.567 | 0.4 | 3.1 | 5.8 |
| Aceh-Andaman | 0.183 | 0.384 | 0.386 | 2.8 | 4.2 | 4.3 |
| Sumatra | 0.612 | 0.786 | 0.795 | 8.9 | 8.5 | 8.8 |
| Java | 0.195 | 0.394 | 1.549 | 2.9 | 4.4 | 16.9 |
| Sumba | 0.066 | 0.309 | 0.316 | 1.0 | 3.5 | 3.8 |

A rupture's rate is divided between segments by its segment fractions. rate_along_trench.png shows the same along the trench.

mfd.png: cumulative rates of the whole zone on every geodynamic branch, as configured and in the test. mfd_segments.png: the same for each segment that was changed (segmented branches).
