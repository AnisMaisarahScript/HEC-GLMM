# HEC Land-Cover Transition GLMM Workflow

This repository contains the R workflow used to analyse associations between
land-cover transitions and human–elephant conflict (HEC) intensity in Johor,
Peninsular Malaysia.

## Analysis workflow

The script implements:

1. Temporal alignment of two study periods:
   - T1: land-cover transitions (2010–2015) with HEC (2013–2017)
   - T2: land-cover transitions (2015–2021) with HEC (2018–2023)
2. Screening of land-cover transition predictors using variance inflation factors (VIF).
3. Generalized linear mixed-effects modelling (GLMM) of KDE-derived HEC intensity.
4. Estimation of marginal and conditional R².
5. Identification of significant land-cover transitions (p < 0.05).
6. Calculation of cumulative HEC impact (ΔHEC = β × n).
7. Translation of cumulative impact into predicted HEC presence-cell counts.
8. Verification using R², p-value, RMSE, and normalized RMSE (NRMSE).

## Main script

`HEC_LandCoverTransition_GLMM_Workflow.R`

The script contains detailed instructions describing the required input files,
data structure, parameters, analysis workflow, and generated outputs.

## Software

The analysis was conducted in R using the following main packages:

- glmmTMB
- performance
- car
- dplyr
- tidyr
- readr
- broom
- ggplot2

## Data availability

The analytical script is provided to support reproducibility of the statistical
workflow. Input data are not included in this repository where access is subject
to data-provider restrictions.

## Citation

If using this workflow, please cite the associated publication and this repository.

## Author

Anis Maisarah Fakhrulanuar
