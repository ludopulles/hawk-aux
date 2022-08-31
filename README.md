# hawk-aux

# HAWK signature scheme
This repository contains:

* a reference to a C implementation of HAWK,
* multiple sage files used in the cryptanalysis and used for generating all the tables in the C implementation,
* data files that were used to generate all the figures in the paper.

[DPPvW21] **NTRU Fatigue: How Stretched is Overstretched ?**
by _Léo Ducas, Eamonn W. Postlethwaite, Ludo Pulles, Wessel van Woerden_

# Contributers

* Léo Ducas,
* Eamonn W. Postlethwaite,
* Ludo Pulles,
* Wessel van Woerden,

# Requirements

* C/C++ compiler, preferable `gcc`,
* [SageMath](https://www.sagemath.org/).

# Description of files
Short description of the files.

In code:

* BKZ_simulator.sage contains a the code used to simualate BKZ reduction, adapted from the [leaky-LWE-estimator](https://github.com/lducas/leaky-LWE-Estimator)
* _attack_Zn.sage contains the functions that generate and apply lattice reduction to forms, _attack_Zn.py is created via sage --preparse _attack_Zn.sage and renaming. This preparsing must take place if _attack_Zn.sage is altered
* attack_Zn.sage contains the high level functions for running experiments in parallel
* bkz_strat.json contains the strategies used by fpylll for BKZ reduction
* convert_sigma_to_stddev.py converts a Gaussian parameter sigma to standard deviation, as these vary below smoothing
* coprimality.sage contains the functions that estimate whether a sampled f, g can be completed using TowerSolve
* exp_varying_n.sage runs experiments where the dimension n varies
* exp_varying_sigma.sage runs experiments where the Gaussian parameter sigma varies
* find_params.sage contains the functions used to find the three parameter sets in the paper, and to compare security in [FALCON](https://falcon-sign.info/)'s security methodology
* generate_C_tables.sage contains the functions that determine tables for sampling from discrete Gaussians with sufficiently small Rényi divergence
* hawk.sage contains a slow but simple implementation of hawk for didactic purposes
* load_strategies.sage loads the BKZ strategies
* pred_varying_n.sage runs the BKZ_simualtor on instances in dimension n
* proba_utils.py contains functions that generate and examine discrete distributions used in find_params.sage

In data:

* XXX_varying_sigma_80trials.txt contain the average block sizes required to solve XXX dimensional instances with the given Gaussian parameters sigma
* XXX_varying_stddev_80trials.txt as above but with sigma converted to standard deviation via conver_sigma_to_stddev.sage
* estimates_tower_solve.txt estimates the number of completable samples and the proportion of these that are completed by TowerSolve for various powers of 2 using coprimality.sage
* estimates_tower_solve_hawk.txt as above but further stipulating HAWK KGen conditions
* experimental_tower_solve_1.5.txt data on how many sampled (f, g) with sigma = 1.5 are completable, and how many of these are completable by TowerSolve (to compare to estimates_tower_solve.txt)
* experimental_tower_solve_1.5_hawk.txt as above but further stipulating HAWK KGen conditions
* pred_varying_n.txt data from pred_varying_n.sage
* varying_n_40trials.txt experimental data on the block size required to solve instances for varying dimensions
* varying_structured_n_40trials.txt as above but with structured lattices for dimensions we can experiment on
