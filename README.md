# Early–Late Window Simulation Code

This repository contains R code for the simulation analyses used in our manuscript on Early–Late (E–L) treatment initiation comparisons under clone–censor–weighting (CCW).

## Overview

This repository provides code for the primary simulation analyses examining how alternative E–L initiation-window specifications affect estimated risk ratios under CCW.

The current version includes:

- simulation data generation for Sim0–Sim3
- five E–L window families
- primary analyses based on weighted empirical 12-month risks
- weighted risk ratio estimation
- support-overlap assessment using the baseline covariate L0
- effective sample size (ESS) calculations
- code for primary-analysis figures

## Current contents

- `EvsL_primary_analysis_public.R`  
  Main R script for the primary simulation analyses.

## Notes

- In this repository, the overlap coefficient (OVL) is defined as the overlap of the empirical distributions of the baseline covariate L0 between the Early and Late groups.
- The primary contrast is summarized as a weighted risk ratio based on weighted empirical 12-month risks.
- The current public version includes the primary analyses only.

## Planned additions

The following components may be added later:

- code for the Monte Carlo bias / RMSE analyses
- code for the secondary analyses excluding pre-treatment events
- additional reproducibility notes and figure-specific scripts, if needed

## Software requirements

The code was written in R and uses common packages for data manipulation, modeling, and plotting.

## Reproducibility

Users should run the script in R after installing the required packages. The script is organized to generate the simulated datasets, fit the weighting models, estimate weighted risks and risk ratios, and produce the primary analysis outputs.

## Contact

For questions regarding the code, please contact the repository owner.
