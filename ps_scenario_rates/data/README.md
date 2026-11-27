# Input files

Read from this folder, or from the folder named by the environment variable `PTHA_DATA`. The inputs are only read, never written to.

| File | In the repository | Content | Read by |
|---|---|---|---|
| `Sunda_triangles.txt` | yes | The Sunda mesh: one row per triangle | step B1 |
| `Seismicity Distribution_Sunda.csv` | yes | Earthquake catalogue with unified Mw, 1906 to 2025 (columns YEAR, MONTH, DAY, HOUR, MIN, SEC, LAT, LON, DEPTH, Unified_Mw, Sigma) | step B3 |
| `Makran_triangles.txt` | yes | The Makran mesh: one row per triangle (Dipankar's text file) | step B1 |
| `Seismicity Distribution_Makran.csv` | yes | Makran earthquake catalogue, 1907 to 2026 (columns YEAR ... DEPTH, Mw, SE); rows without Mw are dropped and the filter of `zones.yaml` applied | step B3 |
| `Rupture_fractions_v2.zip` | no (135 MB, too large for one file on GitHub) | Second version of the rupture fractions tables (`Sunda_rupture_fractions.csv`, `Makran_rupture_fractions.csv`): one row per scenario with its name, total moment and the fraction of its moment in each segment | step B1 |
| `footprints/footprints_<zone>_*.npz` | no; optional | Rupture footprints and mean slip reduced from the slip files with `tools/reduce_slip_files.py`. Only read with `footprints: from_file` in `config/logic_tree.yaml`; needed for rigidity model 01 (`tools/README.md`) | steps B5, B6 |

To rerun the whole chain, place `Rupture_fractions_v2.zip` here. Without it, steps B2 to B8 and the sensitivity tests still
run, because the scenario table that step B1 makes from it is kept in `products/B1_inputs/`.

The magnitude intervals of the scenario provider are a modelling choice and are kept in `config/magnitude_intervals.csv`.

File names are set per zone in `config/zones.yaml`. The rates can be used without any of these files (see `docs/guide.md`, section 2).
