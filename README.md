# LichenReceptorModels

Example code to fit Bayesian multivariate receptor models presented in "Multivariate Receptor Modeling with Widely Dispersed Lichens as Bioindicators of Air Quality" by Heiner, Grimm, Smith, Leavitt, Christensen, Carling, and St. Clair.

## Folders

**data**: Contains `rhizoplaca96.csv` with raw lichen elemental data and `hardZeroId2.csv` with the Id2 structured zero pattern for profile identifiability.

**helperScripts**: R scripts that support the analysis including prior specification of pollutant-source profiles, simulation from the profile priors, elemental detection limits, and identifiability tests.

**plots**: Empty folder to collect plots generated by analysis.

**postsim**: Empty folder to collect MCMC output.

**progress**: Empty folder to collect progress logs output by the `.sh` scripts.

**stan_code**: Contains STAN model files that are called by `1_run_models.R` and `3_simulation_study.R`.

## Analysis scripts
R scripts that comprise the primary analysis workflow are named with ordered numeric prefixes. Most should be run interactively. MCMC settings in the `stan()` function calls are currently set to a low number of iterations. `2_analyze_results.R` has commented options whose use is dictated by the specific model output (i.e., base or sparse). The `_summaries.R` files collect results from multiple runs and summarize them across different model settings.

The simulation study can be run with the `run_simulations.sh` shell script.
