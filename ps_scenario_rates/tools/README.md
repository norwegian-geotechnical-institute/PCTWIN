# Running version 1.1 with rigidity model 01 (on the server that holds the slip files)

Rigidity model 01 is placed by moment under each rupture, so it needs the slip of every scenario (decision 27).
Without it the chain leaves 01 out and gives model 00 all the weight.

1. **Reduce the slip files** (Python 3 and numpy only). One run per rigidity model, or per magnitude folder to run several at once:

       python tools/reduce_slip_files.py PS/sun/00 data/Sunda_triangles.txt sun data/footprints --rig 00 --label 00
       python tools/reduce_slip_files.py PS/sun/01 data/Sunda_triangles.txt sun data/footprints --rig 01 --label 01

   Each run reads every slip realisation and writes `data/footprints/footprints_sun_<label>.npz` (triangles and mean slip of
   each rupture geometry). It prints a moment check for model 00: about 1 in the west and lower in the east (issue H1).
   Text files read at roughly 400 per second per process, so split the work by magnitude folder (`--label 00_9_0260` and so on)
   if it is slow; all `footprints_sun_*.npz` files are read together.

2. **Switch the placement to the footprints** in `config/logic_tree.yaml`:  `footprints: from_file`.

3. **Run the chain** (with `Rupture_fractions_v2.zip` in `data/`):

       python scripts/run_all.py sunda

   Check `checks/B1_inputs` (rigidity models in use: 00 and 01), `checks/B5_shares` (footprint files hold slip) and
   `checks/B7_rates` (rigidity weights 0.4 / 0.6). The rates are in `products/B7_rates/`.

Send back `products/` and `checks/` (or the whole folder without `data/`).
