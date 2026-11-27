# Scenario rates

Zones of this run: sunda (code sun), makran (code mak). Replace <zone> and <code> below by one of them.

Version 1.2. Annual rates of the PS scenarios on the 240 / 240 branches of the logic tree, for the rigidity models
listed in checks/B7_rates (decision 27). Scenarios of the other rigidity models have no rate in this run.

## Mean rates

`scenario_mean_rates_<zone>.csv.gz`: scenario name and annual rate, weighted mean over all branches.
`scenario_rates_<zone>.npz` holds the same with the parts of the name and the means over four groups of branches
(seismic, geodynamic, unconstrained, segmented).

## Any branch

A rate is not stored for every scenario on every branch. It is rebuilt from three small tables:

    rate(scenario, branch) = pooled weight x bin rate(branch) x share(position) x split(scenario)

- branch table: `../B0_config/branch_table_<zone>.csv` (branch id, alternative at each level, weight)
- bin rates: `../B4_bin_rates/bin_rates_<zone>.npz` (unit x variant x form x b x Mmax x bin)
- shares: `../B5_shares/shares_<zone>.npz` (unconstrained branches) and `../B6_segmented/shares_segmented_<zone>.npz`
- split and position of each scenario: `../B1_inputs/scenarios_<zone>.npz`

In Python, from the folder `ps_scenario_rates`:

    from ptha_rates import rates
    r = rates.scenario_rates("<zone>", "G-U-TG-b04-Marea-Cpref")    # one value per scenario

Scenario names in the same order: `rates.scenario_names("<code>", numpy.load("products/B1_inputs/scenarios_<zone>.npz"))`.
