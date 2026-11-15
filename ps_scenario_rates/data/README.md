# Input files

Read from this folder, or from the folder named by the environment variable `PTHA_DATA`. The inputs are only read, never written to.

| File | In the repository | Content | Read by |
|---|---|---|---|
| `Sunda_triangles.txt` | yes | The Sunda mesh: one row per triangle | step B1 |
| `Seismicity Distribution_Sunda.csv` | yes | Earthquake catalogue with unified Mw, 1906 to 2025 (columns YEAR, MONTH, DAY, HOUR, MIN, SEC, LAT, LON, DEPTH, Unified_Mw, Sigma) | step B3 |
| `Rupture_fractions_v2.zip` | no (135 MB, too large for one file on GitHub) | Second version of the rupture fractions tables (`Sunda_rupture_fractions.csv`, `Makran_rupture_fractions.csv`): one row per scenario with its name, total moment and the fraction of its moment in each segment | step B1 |

To rerun the whole chain, place `Rupture_fractions_v2.zip` here. Without it, steps B2 to B8 and the sensitivity tests still
run, because the scenario table that step B1 makes from it is kept in `products/B1_inputs/`.

File names are set per zone in `config/zones.yaml`. The rates can be used without any of these files (see `docs/guide.md`, section 2).
