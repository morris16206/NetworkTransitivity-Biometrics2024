README.md
================

## Overview

This folder contains the simulation code used in the manuscript titled
*Estimating Longitudinal Change in Network Transitivity: Application to
Viral Genetic Linkage Networks*.

The code is used to evaluate the performance of the proposed method
under cross-sectional and longitudinal network settings, as well as
varying levels of transitivity.

## Requirements

- R (\>= 4.3.1)
- Required R packages:
  - `EnvStats`
  - `igraph`
  - `RcppAlgos`
  - `data.table`
  - `plyr`
  - `tidyverse`
  - `locfit`
  - `qqplotr`
  - `gtable`
  - `gridExtra`
  - `Rcpp`
  - `RcppEigen`

You can install them via:

``` r
install.packages(c(
  "EnvStats", "igraph", "RcppAlgos", "data.table", "plyr", "tidyverse", "locfit",
  "qqplotr", "gtable", "gridExtra",
  "Rcpp", "RcppEigen"
), repos = "https://cloud.r-project.org")
```

## How to Run

All required functions are provided in `Simulation_function.R`. To
generate simulation datasets and perform the proposed UGEE method, run
the corresponding `.R` file for each scenario:

- **Cross-sectional networks**: `Simulation_cs_cluster.R`
- **Longitudinal networks**: `Simulation_l_cluster.R`
- **Networks with varying transitivity levels** (cross-sectional):
  `Simulation_diff_trans_cluster.R`

Each of the three `.R` files listed above begins with the following line
to load the necessary functions. Be sure to run it before proceeding:

``` r
source("Simulation_function.R")
```

**Note**: The simulation may take a long time (from several hours up to
days or even weeks, depending on the sample size and the number of time
points). The provided simulation code is designed for scenarios where
computational resources are manageable. For applications involving very
large sample sizes or long time series, further optimization may be
required.

The scripts will generate and save simulation results in `.RData` format
in the working directory. The results are provided as Monte Carlo
replicates. Use the `Simulation.R` file to summarize the Monte Carlo
results.

## Reproducibility

To ensure reproducibility, a seed is set by default (`set.seed(...)`).
You can change the randomness by modifying any `set.seed()` call found
in the `.R` files.

## Contact

For questions about the code or simulations, please contact:

**Tsung-Chin Wu**  
Email: <tswu@health.ucsd.edu>
