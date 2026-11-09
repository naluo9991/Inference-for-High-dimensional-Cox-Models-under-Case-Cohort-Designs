# Inference for High-dimensional Cox Models under Case-Cohort Designs

R code accompanying the manuscript **Inference for High-dimensional Cox Models under Case-Cohort Designs**.

## Overview

This project studies estimation and inference for sparse high-dimensional Cox proportional hazards models under case-cohort sampling. A fixed-size subcohort is selected by simple random sampling without replacement, and all cases outside the subcohort are also included.

The repository contains simulation and real-data analysis scripts for Lasso estimation and coefficient-wise hypothesis testing, accounting for the additional variability introduced by subcohort sampling. The real-data application uses gene expression and overall survival data from the METABRIC study.

## Methods

The scripts implement three tests of a zero regression coefficient:

| Code | Method |
| --- | --- |
| `S` | Decorrelated score test |
| `W` | Decorrelated Wald test |
| `L` | Decorrelated partial likelihood ratio test (PLRT) |

For each tested coefficient, the remaining coefficients are treated as high-dimensional nuisance parameters. See the manuscript for the theoretical results and assumptions.

## Repository contents

| Path | Description |
| --- | --- |
| `simulation/rej_n050_s2.R` | Simulation script for evaluating rejection rates across a grid of target coefficient values. |
| `real/real.R` | METABRIC analysis for the full cohort and three subcohort sizes. |
| `real/data/YMETABRIC.RData` | Supplied survival-data file. |
| `real/result/` | Saved coefficient estimates, test p-values, and subcohort membership indicators. |

Both R scripts contain their own inference functions and experiment code.

## Requirements and execution

The scripts use R and load the following packages:

`glmnet`, `survival`, `MASS`, `mvtnorm`, `glasso`, `lpSolve`, `iterators`, `foreach`, `Rmpi`, `doMPI`, `doParallel`, and `nleqslv`.

The `parallel` package is included with R. The real-data script optionally uses `matrixStats` to accelerate cumulative sums.

Both scripts initialize an MPI backend using `startMPIcluster()` and `registerDoMPI()`. Run them as batch jobs in an MPI-capable environment using the launcher or scheduler configured for that environment. The final `mpi.quit()` call terminates the R session.

Input and output paths are relative to the R working directory. Set the working directory and check file paths before execution. New output files are written to that directory; they are not automatically placed in `real/result/`.

## Simulation study

The supplied script uses the following default configuration:

| Parameter | Meaning | Default |
| --- | --- | --- |
| `n` | Full-cohort size | `200` |
| `d` | Number of covariates | `500` |
| `n0` | Subcohort size | `50` |
| `p0` | Subcohort sampling fraction, `n0 / n` | `0.25` |
| `s` | Number of nonzero nuisance coefficients | `2` |
| `rho` | Correlation parameter in the Toeplitz covariance matrix | `0.4` |
| `iters` | Repetitions per method and target coefficient value | `1000` |
| `rej_data_list` | Values of the target coefficient | `seq(0, 0.55, by = 0.01)` |
| `method1_list` | Tests to evaluate | `c("S", "W", "L")` |

The target is the first regression coefficient. The `s` nonzero nuisance coefficients are selected from coordinates `2:d` and set to one. Thus, the total number of nonzero coefficients is `s` under the null and `s + 1` under the alternatives.

Each simulation CSV contains the p-values from the repeated experiments, with one repetition per column. The proportion below `0.05` gives the empirical type-I error rate when the target coefficient is zero, and empirical power otherwise.

For example, `rej_d500_n0050_rho40_Dirac_S_s2_00.csv` corresponds to the score test with `d = 500`, `n0 = 50`, `rho = 0.4`, `s = 2`, and a zero target coefficient.

Other settings can be investigated by changing the parameters and signal generation in the script. For a short exploratory run, reduce `iters` and the number of values in `rej_data_list`; the MPI requirement still applies.

## METABRIC application

The manuscript analyzes 1,979 patients and 20,603 gene expression variables, using overall survival as the outcome. The original data were obtained from cBioPortal. The manuscript describes patient matching, exclusions, and gene-wise median imputation of missing expression measurements.

The current script expects two input files in the R working directory:

| File | Object expected by the script |
| --- | --- |
| `XMETABRIC.RData` | `X_dlbcl`: numeric matrix with patients in rows and genes in columns. |
| `YMETABRIC.RData` | `Y_dlbcl`: outcomes in the same patient order; column 1 is event status (`1` = death, `0` = censored), and column 2 is observed time in months. |

The object names above are those used in the current script. The analysis concerns METABRIC.

**Data availability:** The processed gene expression matrix ('XMETABRIC.RData') is not included in this repository because of its large file size. The survival data are provided in 'real/data/YMETABRIC.RData'.

The script analyzes the full cohort and subcohorts of sizes `495`, `990`, and `1484`, adding all cases outside each subcohort. It performs one draw per subcohort size per run. The repeated-sampling comparisons and summary table in the manuscript require repeated runs and separate post-processing.

### Outputs

| Filename pattern | Contents |
| --- | --- |
| `beta_hat_<n0>.csv` | Lasso coefficient estimates in input gene order. |
| `n0<n0>_S.csv` | Unadjusted score-test p-values, one gene per column. |
| `n0<n0>_W.csv` | Unadjusted Wald-test p-values, one gene per column. |
| `n0<n0>_L.csv` | Unadjusted PLRT p-values, one gene per column. |
| `n0_<n0>.rds` | Logical subcohort membership vector in full-cohort patient order. |

For example, `n0495_S.csv` contains score-test p-values for `n0 = 495`. The p-values are not BH-adjusted. In the manuscript, BH adjustment is applied separately within each method, and the reported significant gene set is the intersection of the three sets with adjusted p-values below `0.05`.

## Reproducibility notes

- Subcohort selection and cross-validation involve randomness. Reusing a saved membership vector fixes the sampled patients; reproducing the complete analysis also requires controlling cross-validation folds and random-number generation, including parallel workers where relevant.
- Use separate output directories for repeated runs, since the current filenames do not include a repetition identifier.
- The manuscript assumes untied failure times. The optimized cumulative-risk-set calculations in `real/real.R` do not include explicit tie handling; tied observations require a consistent, documented treatment.
