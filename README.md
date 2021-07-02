# BMIM
Code to run simulations and analysis for "Bayesian Multiple Index Models for Environmental Mixtures"

## Functions
### bsmim2
Rcpp package to fit a BMIM in R
- bsmim2.R: contains main function bsmim2() for fitting BMIMs
- bsmim_predict2.R: contains functions for predicting response surface
- plot_univar_hnew2.R: contains functions for plotting exposure response curves (componentwise and indexwise)

## Simulations
Code to replicate simulation scenarios 
- NHANES_simA2.R: single index scenario
- NHANES_simB.R: 3-index scenario


as well as supplementary simulations:
- NHANES_simA1.R: linear single index scenario
- NHANES_simC.R: mis-specified index groupings
- NHANES_simD.R: generating under full BKMR

and Rmarkdown files to report results from main & supplementary simulations:
- report_NHANES_sims.Rmd
- report_NHANES_supp_sims.Rmd

## Analysis
Code to replicate NHANES data analysis.
- NHANES_analysis.R: code to run NHANES analysis
- NHANES_analysis_CV.R: code to run 5-fold cross-validation on NHANES data (on cluster looping iter_no 1-5)
- report_NHANES_analysis.Rmd: code to report results of NHANES analysis
