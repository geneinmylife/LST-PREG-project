

# Sex-Stratified Genetic Analyses Mapping the Influences of Sedentary Behaviors and Physical Activity on Female Reproductive Health

## Overview
This repository contains genome-wide association study (GWAS) and Mendelian Randomization (MR) code for a study investigating the causal relationships between sedentary behaviours, physical activity and risks of pregnancy-related disorders. 

![docs](vignettes/figures/Graphical_abstract.png)  

***Key findings*** of our study include:

- 18 novel sex-stratified variants were identified in our sex-stratified GWAS. 
- LST exhibits higher heritability in females compared to males.
- Female-specific genetic instruments reveal that higher LST increases risks of menorrhagia, endometriosis, PCOS and ectopic pregnancy.
- LST exceeding 4 hours/day elevates menorrhagia and endometriosis risk.

##

## Installation
The following working environments are required: 
- R (≥4.0)
- Required R packages:
  ```R
  install.packages("remotes")
  remotes::install_github("MRCIEU/TwoSampleMR")  # Package for core MR analysis
  install.packages("dplyr")                      # Package for data manipulation
   #other packages are listed in each code.
  ```

## Statistical analysis scripts

The following code were used in this study:

- ***gwas_code.sh*** script to conduct a sex-stratified GWAS on sedentary behaviours and physical activity using REGENIE.
- ***GCTA_cojo_code.sh*** script to identify the instrument variables of exposures.
- ***Post-GWAS analysis.R*** To calculate heritability and cross-trait genetic correlations.
- ***MR analysis.R*** To identify causal effect between female LST and reproductive outcomes, including univariate MR, Multivariable MR and nonlinear MR.
