# PTHA integration

Code for the probabilistic tsunami hazard assessment (PTHA) in PCTWIN WP1.

This repository brings together the pieces being developed across the work package (source models, magnitude-frequency distributions, rupture realisations, tsunami propagation) into one PTHA workflow. The first version covers predominant subduction (PS) sources only, for the Makran and Sunda subduction zones. Background seismicity (BS) will be added later.

The implementation is still in an early stage, so expect things to move around.

## Workflow

The PTHA is built up in the following steps:

1. Source model
2. Epistemic model / logic tree
3. MFD and annual rates per magnitude bin
4. Aleatory event generation
5. Scenario catalogue
6. Unit-source combination / initial conditions
7. Tsunami propagation (Tsunami-HySEA)
8. Intensities at points of interest (POIs) and inundation
9. Hazard aggregation

## Modelling choices

Epistemic uncertainty is handled through logic-tree branches, while aleatory variability is sampled during event generation. The main components currently considered are:

- MFD type: truncated and tapered Pareto
- b-value, including constraints from geodynamic modelling
- Mmax: global priors with Bayesian fitting
- Non-uniform magnitude bins (Makran Mw 6.5-9.02, Sunda Mw 6.8-9.5)
- Segmentation: combined and split
- Scaling relations: Strasser et al. (2010) and Thingbaijam et al. (2017), with scaling-law uncertainty
- Rigidity: depth-variable and homogeneous profiles
- Stochastic slip realisations
- Focal mechanism and depth

Scenario rates are computed by distributing the annual rate of each magnitude bin over the aleatory realisations (slip, scaling, etc.) within each logic-tree branch.
