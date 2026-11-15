# Scenario rates (sunda)

Annual rates of the PS scenarios of rigidity model 00 (homogeneous) on the 480 branches of the logic tree.
Scenarios of the other rigidity models have no rate in this run.

## Mean rates

`scenario_mean_rates_sunda.csv.gz`: scenario name and annual rate, weighted mean over all branches.
`scenario_rates_sunda.npz` holds the same with the parts of the name and the means over four groups of branches
(seismic, geodynamic, unconstrained, segmented).

## Any branch

A rate is not stored for every scenario on every branch. It is rebuilt from three small tables:

    rate(scenario, branch) = pooled weight x bin rate(branch) x share(position) x split(scenario)

- branch table: `../B0_config/branch_table_sunda.csv` (branch id, alternative at each level, weight)
- bin rates: `../B4_bin_rates/bin_rates_sunda.npz` (unit x variant x form x b x Mmax x bin)
- shares: `../B5_shares/shares_sunda.npz` (unconstrained branches) and `../B6_segmented/shares_segmented_sunda.npz`
- split and position of each scenario: `../B1_inputs/scenarios_sunda.npz`

In Python, from the folder `ps_scenario_rates`:

    from ptha_rates import rates
    r = rates.scenario_rates("sunda", "G-U-TG-b04-Marea-Cpref")    # one value per scenario

Scenario names in the same order: `rates.scenario_names("sun", numpy.load("products/B1_inputs/scenarios_sunda.npz"))`.
