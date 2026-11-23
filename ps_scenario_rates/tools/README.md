# Running version 1.1 with rigidity model 01 (on the server that holds the slip files)

Rigidity model 01 is placed by moment under each rupture, so it needs the slip of every scenario (decision 27).
Without it the chain leaves 01 out and gives model 00 all the weight. Run everything from the `ps_scenario_rates` folder.

1. **Reduce the slip files** (Python 3 and numpy only), one run per rigidity model:

       python tools/reduce_slip_files.py PS/sun/00 data/Sunda_triangles.txt sun data/footprints --rig 00 --label 00
       python tools/reduce_slip_files.py PS/sun/01 data/Sunda_triangles.txt sun data/footprints --rig 01 --label 01

   `PS/sun/00` is the folder that holds the slip files of model 00 (searched with its subfolders); change it to where they are.
   Each run writes `data/footprints/footprints_sun_<label>.npz` (triangles and mean slip of each rupture geometry) and ends
   with a check of the moment of each rupture against its named Mw, computed with the same areas and rigidity as the scenario
   generator. It should say OK.

   The two models hold about 4 million slip files. One process reads about 400 per second, so about 3 hours in all. To go
   faster, run one magnitude folder per process, with its own label, for example:

       python tools/reduce_slip_files.py PS/sun/00/9_0260 data/Sunda_triangles.txt sun data/footprints --rig 00 --label 00_9_0260

   All `footprints_sun_*.npz` files in `data/footprints/` are read together afterwards.

2. **Switch the placement to the footprints** in `config/logic_tree.yaml`:  `footprints: from_file`.

3. **Run the chain**, with `Rupture_fractions_v2.zip` in `data/` (needs several GB of memory; more than 4 GB):

       python scripts/run_all.py sunda

4. **Check** that `checks/INDEX.md` shows no FAIL, that `checks/B1_inputs` lists rigidity models 00 and 01 in use, and that
   `checks/B5_shares` confirms the footprint files hold slip.

If all is fine, the scenario rates of both models are in `products/B7_rates/`; section 2 of `docs/guide.md` shows how to
compute hazard curves from them.
